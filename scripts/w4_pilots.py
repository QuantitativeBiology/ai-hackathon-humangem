"""W4 / M3 pilots (label-free). Subcommands:

  repair12  F0 (R1, RNA, M0) on the 12 lines of the compute-review probe (the first 12 default rows of the
            expression file): WT growth before/after the fixed repair rule, re-opened reactions, task LPs.
  lossrule  loss rule + order-independent drop on the 30 panel lines; 6 input orders per line must agree.
  cost      per-line cost: F0/F2/F3/F4 x R1 x {real, permuted} RNA (+ genomic re-sweep), one process per
            family (peak RSS per family), Gurobi; CORDA (F3b) separately with --families F3b.
  glpk      GLPK vs Gurobi per family: (a) full pipeline on GLPK; (b) the Gurobi-built final model swept on GLPK.

Outputs in results/w4/ (training lines by ModelID; non-training lines of the probe are masked).
Never opens a label file. Usage: .venv/bin/python scripts/w4_pilots.py <cmd> [--lines N] [--families ...]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import resource
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from nextgem import context as cx
from nextgem import context_inputs as ci
from nextgem import context_run as cr
from nextgem import paths
from nextgem.config import SEEDS

OUT = paths.results_dir() / "w4"
FASTCC = OUT / "fastcc_human2_M0.json"
PROBE12 = ["ACH-001113", "ACH-001289", "ACH-001339", "ACH-001979", "ACH-002438", "ACH-003165",
           "ACH-001538", "ACH-000242", "ACH-000708", "ACH-003477", "ACH-002989", "ACH-000327"]


# ----------------------------------------------------------------------------- inputs (parent process)
class Inputs:
    def __init__(self, lines, need_prot=False):
        from nextgem.model import load_model  # noqa: F401  (genes table only needs the SBML genes)
        import cobra

        self.split = pd.read_csv(paths.splits_dir() / "lineage_split.tsv", sep="\t")
        self.panel = pd.read_csv(OUT / "stage_a_panel.tsv", sep="\t")
        m = cobra.io.read_sbml_model(str(paths.data_root() / "Human-GEM-main/model/Human-GEM.xml"))
        self.genes = ci.gene_table(m)
        del m
        st = json.loads((OUT / "threshold_stats.json").read_text())
        gm = pd.read_parquet(OUT / "threshold_gene_means.parquet").set_index("gene")
        self.rna_stats = ci.ThresholdStats("rna", st["rna"]["p25"], st["rna"]["p75"], gm.rna_mean.dropna(),
                                           st["rna"]["floor_value"], st["rna"]["n_lines"])
        self.prot_stats = ci.ThresholdStats("protein", st["protein"]["p25"], st["protein"]["p75"],
                                            gm.protein_mean.dropna(), st["protein"]["floor_value"], st["protein"]["n_lines"])
        allines = sorted(set(lines) | set(self.panel.ModelID))
        self.rna = ci.read_expression(self.genes, allines)
        self.rna_perm = ci.permute_across_lines(self.rna.loc[self.panel.ModelID], SEEDS[0])
        self.cn = ci.read_cn(self.genes, allines)
        self.dmg = ci.read_damaging(self.genes, allines)
        self.prot = ci.read_procan(self.genes, allines) if need_prot else None

    def scores(self, line, rule="R1", modality="rna", variant="real"):
        rna = self.rna_perm if variant == "perm" else self.rna
        return ci.signed_scores(line, rule, modality, self.genes, rna, self.rna_stats, self.prot, self.prot_stats)

    def lost(self, line, **kw):
        return ci.lost_genes(line, self.dmg, self.cn, self.rna, **kw)


def mask(split: pd.DataFrame, line: str, k: int) -> tuple[str, str]:
    tier = split.set_index("ModelID").tier.get(line, "unscreened")
    return (line if tier in ("train", "unscreened") else f"masked_{tier}_{k}"), tier


def group(init_kwargs: dict, fn_name: str, payload) -> dict:
    """Run in a fresh worker process: init the parent, then call a pilot function on the payload."""
    t0, c0 = time.time(), time.process_time()
    cr.init_worker(**init_kwargs)
    res = globals()[fn_name](payload)
    return {"result": res, "load_wall_s": cr._W["load_wall_s"], "load_cpu_s": cr._W["load_cpu_s"],
            "wall_s": time.time() - t0, "cpu_s": time.process_time() - c0, "peak_rss_gb": cr.peak_rss_gb()}


def run_groups(groups, max_workers):
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx) as ex:
        futs = [ex.submit(group, *g) for g in groups]
        return [f.result() for f in futs]


# ----------------------------------------------------------------------------- repair12
def _repair12_worker(items):
    parent = cr._W["parent"]
    rows = []
    for label, tier, s in items:
        lm = cx.extract(parent, "F0", s.reindex(parent.genes))
        pre = cx._growth_with_closed(parent, lm.closed)
        c0 = time.process_time()
        cx.finish(parent, lm)
        rows.append({"line": label, "tier": tier, "n_off_genes": len(lm.off), "n_closed": len(lm.closed),
                     "wt_before_repair": pre, "growth_lp": lm.info["repair_growth_lp"],
                     "n_reopened": lm.info["n_reopened"], "reopen_frac_of_closed": lm.info["n_reopened"] / len(lm.closed),
                     "wt_after_repair": lm.info["wt_growth"], "wt_ratio_generic": lm.info["wt_ratio_generic"],
                     "failed": lm.info["failed"], "repair_cpu_s": time.process_time() - c0,
                     "task_lps": ";".join(f"{k}={v}" for k, v in lm.info["repair_task_lps"].items())})
    return rows


def cmd_repair12(a):
    hdr = pd.read_csv(paths.depmap_dir() / ci.EXPR_FILE, usecols=["ModelID", "IsDefaultEntryForModel"], nrows=40)
    first12 = list(hdr[hdr.IsDefaultEntryForModel == "Yes"].ModelID[:12])
    assert first12 == PROBE12, first12  # the compute-review probe's lines (probe2.py: first 12 default rows)
    I = Inputs(PROBE12)
    items = []
    for k, line in enumerate(PROBE12):
        label, tier = mask(I.split, line, k)
        items.append((label, tier, I.scores(line)))
    res = run_groups([({"fastcc_cache": None}, "_repair12_worker", items)], 1)[0]
    df = pd.DataFrame(res["result"])
    df.to_csv(OUT / "repair_probe12.tsv", sep="\t", index=False)
    print(df.drop(columns="task_lps").to_string())
    print({k: v for k, v in res.items() if k != "result"})


# ----------------------------------------------------------------------------- loss rule
def _loss_worker(items):
    parent = cr._W["parent"]
    rows = []
    for line, lost, expr, flags in items:
        rng = np.random.default_rng(SEEDS[0])
        orders = [lost, lost[::-1]] + [list(rng.permutation(lost)) for _ in range(4)]
        c0 = time.process_time()
        results = [cx.accept_losses(parent, o, expr) for o in orders]
        cpu = (time.process_time() - c0) / len(orders)
        key = lambda r: (tuple(r["contradicted"]), tuple(r["dropped_joint"]), tuple(r["accepted"]))  # noqa: E731
        same = all(key(r) == key(results[0]) for r in results)
        r0 = results[0]
        rows.append({"line": line, "n_lost_model_genes": len(r0["lost"]), "n_contradicted": len(r0["contradicted"]),
                     "n_dropped_joint": len(r0["dropped_joint"]), "n_accepted": len(r0["accepted"]),
                     "order_independent": same, "n_orders": len(orders), "cpu_s_per_run": cpu,
                     "joint_growth_ratio": r0["joint_growth"] / parent.wt_generic,
                     "contradicted": ";".join(r0["contradicted"]), "dropped_joint": ";".join(r0["dropped_joint"]),
                     "cn_missing": flags["cn_missing"], "_accepted": r0["accepted"]})
    return rows


def cmd_lossrule(a):
    I = Inputs([])
    items = []
    for line in I.panel.ModelID:
        lost, flags = I.lost(line)
        lost = [g for g in lost if g in set(I.genes.index)]
        items.append((line, lost, I.rna.loc[line], flags))
    res = run_groups([({"fastcc_cache": None}, "_loss_worker", items)], 1)[0]
    df = pd.DataFrame(res["result"])
    (OUT / "tmp" / "accepted_losses_panel.json").write_text(json.dumps(dict(zip(df.line, df.pop("_accepted")))))
    sym = I.genes.symbol.to_dict()
    for c in ("contradicted", "dropped_joint"):
        df[c] = df[c].map(lambda s: ";".join(sym.get(g, g) for g in s.split(";") if g))
    df.to_csv(OUT / "loss_rule_panel.tsv", sep="\t", index=False)
    print(df.to_string())
    print("all order-independent:", bool(df.order_independent.all()), {k: v for k, v in res.items() if k != "result"})


# ----------------------------------------------------------------------------- cost / glpk
def _jobs_worker(jobs):
    out = []
    for j in jobs:
        r = cr.run_job(j)
        out.append({"rows": r["rows"], "ko": {k: v for k, v in r["ko"].items()}})
    return out


def _sweep_fixed_worker(items):
    """GLPK sweep of a Gurobi-built final model (identical reaction set)."""
    parent = cr._W["parent"]
    out = []
    for family, line, final_closed, off, lost, eflux, lin, cap_max in items:
        lm = cx.LineModel(line=line, family=family, closed=set(final_closed), off=set(off), lost=set(lost),
                          eflux=dict(eflux), lin=dict(lin), cap_max=cap_max)
        lm.lost_closed = cx.loss_closure(parent, lm.lost)
        c0 = time.process_time()
        ko = cx.ko_sweep(parent, lm)
        out.append({"family": family, "line": line, "ko": ko[["gene", "ratio_r", "status"]],
                    "sweep_cpu_s": time.process_time() - c0})
    return out


def make_jobs(I, lines, families, rule="R1", variants=("real",), genomic=("none", "loss"), tasks=True, accepted=None):
    jobs = {f: [] for f in families}
    for f in families:
        for line in lines:
            for v in variants:
                g = genomic if v == "real" else ("none",)
                acc = tuple((accepted or {}).get(line, ()))
                jobs[f].append(cr.Job(build_id=f"{line}|{f}|{rule}|rna|{v}", line=line, family=f, rule=rule,
                                      modality="rna", variant=v, scores=I.scores(line, rule, "rna", v),
                                      genomic=g, accepted=acc, tasks=tasks))
    return jobs


def accepted_losses(I, lines):
    p = OUT / "tmp" / "accepted_losses_panel.json"
    if not p.exists():
        raise SystemExit("run `lossrule` first")
    d = json.loads(p.read_text())
    return {line: d[line] for line in lines}


def cmd_cost(a):
    I = Inputs([])
    lines = list(I.panel.ModelID[: a.lines]) if not a.line_ids else a.line_ids
    acc = accepted_losses(I, lines)
    families = a.families
    variants = ("real", "perm") if a.perm else ("real",)
    jobs = make_jobs(I, lines, families, variants=variants, accepted=acc, tasks=not a.no_tasks, genomic=tuple(a.genomic))
    if a.split_lines:  # one process per (family, line): isolates per-line RSS and parallelises long builds
        groups = [({"fastcc_cache": str(FASTCC), "solver": a.solver}, "_jobs_worker", [j]) for f in families for j in jobs[f]]
    else:
        groups = [({"fastcc_cache": str(FASTCC), "solver": a.solver}, "_jobs_worker", jobs[f]) for f in families]
    t0 = time.time()
    res = run_groups(groups, min(a.workers, len(groups)))
    rows, kos = [], []
    for r in res:
        for jr in r["result"]:
            for row in jr["rows"]:
                rows.append({**row, "worker_peak_rss_gb": r["peak_rss_gb"], "worker_load_cpu_s": r["load_cpu_s"]})
            for gvar, ko in jr["ko"].items():
                ko = ko.copy()
                ko["build_id"] = f"{jr['rows'][0]['build_id'].rsplit('|', 1)[0]}|{gvar}"
                kos.append(ko)
    df = pd.DataFrame(rows)
    tag = a.tag or f"{'_'.join(families)}_{a.solver}"
    df.to_csv(OUT / f"pilot_builds_{tag}.tsv", sep="\t", index=False)
    if kos:
        pd.concat(kos).to_parquet(OUT / "tmp" / f"pilot_ko_{tag}.parquet", index=False)
    cols = ["line", "family", "variant", "genomic", "failed", "n_closed_extract", "n_reopened", "repair_frac_reactions",
            "n_open", "wt_ratio_generic", "sweep_n_solved", "sweep_n_lethal", "tasks_passed", "extract_cpu_s",
            "variant_cpu_s", "job_cpu_s", "peak_rss_gb"]
    print(df[[c for c in cols if c in df]].to_string())
    print(f"wall {time.time() - t0:.0f} s")


def cmd_glpk(a):
    """(a) full pipeline on GLPK for the given families/line; (b) Gurobi final models swept on GLPK."""
    I = Inputs([])
    line = a.line_ids[0] if a.line_ids else I.panel.ModelID.iloc[0]
    ref = pd.concat([pd.read_csv(OUT / f"pilot_builds_{t}.tsv", sep="\t") for t in a.ref_tag])
    ref_ko = pd.concat([pd.read_parquet(OUT / "tmp" / f"pilot_ko_{t}.parquet") for t in a.ref_tag])
    acc = accepted_losses(I, [line])
    families = a.families
    jobs = make_jobs(I, [line], families, genomic=("none",), tasks=False, accepted=acc)
    # (b) payload: Gurobi final reaction sets
    import cobra
    m = cobra.io.read_sbml_model(str(paths.data_root() / "Human-GEM-main/model/Human-GEM.xml"))
    rxn_ids = [r.id for r in m.reactions]
    items = []
    for f in families:
        row = ref[(ref.line == line) & (ref.family == f) & (ref.variant == "real") & (ref.genomic == "none")].iloc[0]
        bits = np.unpackbits(np.frombuffer(bytes.fromhex(row.open_mask), dtype=np.uint8))[: len(rxn_ids)].astype(bool)
        final_closed = [r for r, b in zip(rxn_ids, bits) if not b]
        s = I.scores(line)
        off = list(s.index[(s < 0).to_numpy()])
        eflux, lin, capmax = ({}, {}, 1.0)
        items.append((f, line, final_closed, off, [], eflux, lin, capmax))
    del m
    groups = [({"fastcc_cache": str(FASTCC), "solver": "glpk"}, "_jobs_worker", jobs[f]) for f in families]
    sweep_items = [it for it in items if it[0] != "F4"]  # F4 has no closed set: (a) covers it
    groups.append(({"fastcc_cache": str(FASTCC), "solver": "glpk"}, "_sweep_fixed_worker", sweep_items))
    res = run_groups(groups, min(a.workers, len(groups)))
    rows, comp = [], []
    for f, r in zip(families, res[:-1]):
        jr = r["result"][0]
        row = jr["rows"][0]
        rows.append({**row, "worker_peak_rss_gb": r["peak_rss_gb"]})
        g_row = ref[(ref.line == line) & (ref.family == f) & (ref.variant == "real") & (ref.genomic == "none")].iloc[0]
        gk = ref_ko[ref_ko.build_id == f"{line}|{f}|R1|rna|real|none"].set_index("gene").ratio_r
        lk = jr["ko"]["none"].set_index("gene").ratio_r if "none" in jr["ko"] else pd.Series(dtype=float)
        a_open = np.unpackbits(np.frombuffer(bytes.fromhex(g_row.open_mask), dtype=np.uint8)).astype(bool)
        b_open = np.unpackbits(np.frombuffer(bytes.fromhex(row["open_mask"]), dtype=np.uint8)).astype(bool) if "open_mask" in row else a_open
        jac = (a_open & b_open).sum() / max(1, (a_open | b_open).sum())
        common = gk.index.intersection(lk.index)
        comp.append({"family": f, "line": line, "mode": "full_pipeline", "reaction_set_jaccard": jac,
                     "n_genes": len(common), "n_ratio_r_differ": int((gk[common] - lk[common]).abs().gt(1e-6).sum()),
                     "max_abs_diff": float((gk[common] - lk[common]).abs().max()) if len(common) else np.nan,
                     "glpk_job_cpu_s": row.get("job_cpu_s"), "gurobi_job_cpu_s": g_row.job_cpu_s,
                     "glpk_sweep_cpu_s": row.get("sweep_cpu_s"), "gurobi_sweep_cpu_s": g_row.sweep_cpu_s})
    for it in res[-1]["result"]:
        f = it["family"]
        gk = ref_ko[ref_ko.build_id == f"{line}|{f}|R1|rna|real|none"].set_index("gene").ratio_r
        lk = it["ko"].set_index("gene").ratio_r
        common = gk.index.intersection(lk.index)
        g_row = ref[(ref.line == line) & (ref.family == f) & (ref.variant == "real") & (ref.genomic == "none")].iloc[0]
        comp.append({"family": f, "line": line, "mode": "same_model_sweep", "reaction_set_jaccard": 1.0,
                     "n_genes": len(common), "n_ratio_r_differ": int((gk[common] - lk[common]).abs().gt(1e-6).sum()),
                     "max_abs_diff": float((gk[common] - lk[common]).abs().max()),
                     "glpk_sweep_cpu_s": it["sweep_cpu_s"], "gurobi_sweep_cpu_s": g_row.sweep_cpu_s})
    cdf = pd.DataFrame(comp)
    cdf.to_csv(OUT / "glpk_vs_gurobi.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(OUT / f"pilot_builds_glpk_{line}.tsv", sep="\t", index=False)
    print(cdf.to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["repair12", "lossrule", "cost", "glpk"])
    ap.add_argument("--lines", type=int, default=3)
    ap.add_argument("--line-ids", nargs="*")
    ap.add_argument("--families", nargs="*", default=["F0", "F2", "F3", "F4"])
    ap.add_argument("--solver", default="gurobi")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--perm", action="store_true")
    ap.add_argument("--split-lines", action="store_true")
    ap.add_argument("--genomic", nargs="*", default=["none", "loss"])
    ap.add_argument("--no-tasks", action="store_true")
    ap.add_argument("--tag")
    ap.add_argument("--ref-tag", nargs="*", default=["F0_F2_F3_F4_gurobi"])
    a = ap.parse_args()
    (OUT / "tmp").mkdir(parents=True, exist_ok=True)
    t0, c0 = time.time(), time.process_time()
    {"repair12": cmd_repair12, "lossrule": cmd_lossrule, "cost": cmd_cost, "glpk": cmd_glpk}[a.cmd](a)
    print(f"[{a.cmd}] parent wall {time.time() - t0:.0f} s, parent CPU {time.process_time() - c0:.0f} s, "
          f"parent peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB")


if __name__ == "__main__":
    main()
