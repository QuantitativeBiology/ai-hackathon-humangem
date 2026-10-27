"""Stage-A analysis, label-free part (DESIGN.md §9 W4, §13; amendments A4, A7, A10). Reads no label file.

Inputs: results/w4/stage_a/ (builds.tsv, ko_*.parquet, a6prime_lines.tsv; Stage A run at c280413), the split
file (lineage), Human-GEM (reaction order, gene symbols, task list), DepMap expression (panel lines, for the
expression reference of the variance partition).

Outputs (results/w4/stage_a_analysis/):
  build_gates.tsv           per build x genomic variant: label-blind gates and build statistics
  config_summary.tsv        per configuration (row, family, rule, modality, variant, genomic): statistics and gates
  a7_stability.tsv          A7: reaction-set and lethal-set Jaccard across the 3 seeds and vs the unperturbed build
  a6prime_comparison.tsv    A6' (A7 amendment): CORDA vs F0/F2/F3/F4 on the same lines and thresholds (descriptive)
  permuted_control.tsv      real vs permuted-expression builds (A1), per family x threshold
  positive_controls.tsv     TYR/DCT reactions in Skin vs other lines; urea-cycle reactions open
  variance_partition.tsv    marginal R^2 (method, threshold, lineage, line) per level x build set, permutation p
  h4_tests.json             H4.1, H4.2, H4.3 verdicts with bootstrap CIs and permuted-control nulls
  MANIFEST.tsv              provenance

Usage: .venv/bin/python scripts/w4_stage_a_labelfree.py [--perms 999] [--boot 2000] [--allow-dirty]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from nextgem import paths
from nextgem import w4_analysis as wa
from nextgem.config import SEEDS

W4 = paths.results_dir() / "w4"
STAGE = W4 / "stage_a"
OUT = W4 / "stage_a_analysis"
LP = ["F0", "F2", "F3", "F4"]
CFG_KEYS = ["row", "family", "rule", "modality", "variant", "genomic"]
TYR_DCT = ["TYR", "DCT"]
UREA = ["CPS1", "OTC", "ARG1", "ASS1", "ASL"]


def load():
    import cobra
    from raven_toolbox.tasks.tasklist import parse_task_list

    from nextgem.context_inputs import gene_table
    from nextgem.model import GEMS

    cobra.Configuration().solver = "glpk"  # structure only; no Gurobi session
    m = cobra.io.read_sbml_model(str(GEMS["human2"].sbml))
    rxn_ids = [r.id for r in m.reactions]
    gt = gene_table(m)
    sym_rxns = {}
    for s in TYR_DCT + UREA:
        gs = list(gt.index[gt.symbol == s])
        sym_rxns[s] = sorted({r.id for g in gs for r in m.genes.get_by_id(g).reactions})
    tasks = parse_task_list(str(GEMS["human2"].essential_tasks))
    split = pd.read_csv(paths.splits_dir() / "lineage_split.tsv", sep="\t").set_index("ModelID")
    b = pd.read_csv(STAGE / "builds.tsv", sep="\t", low_memory=False)
    b["lineage"] = b.line.map(split.lineage)
    b["build_cpu_s"] = b.extract_cpu_s.fillna(0) + b.variant_cpu_s.fillna(0)
    return dict(model=m, rxn_ids=rxn_ids, gt=gt, sym_rxns=sym_rxns, tasks=tasks, split=split, builds=b)


def ko_lethal(rows: list[str], gene_order) -> dict:
    """build_id -> boolean lethal vector (ratio_r < 0.01) over ``gene_order``."""
    out = {}
    for r in rows:
        k = pd.read_parquet(STAGE / f"ko_{r}.parquet", columns=["build_id", "gene", "ratio_r"])
        k["build_id"] = k.build_id.astype(str)
        for bid, g in k.groupby("build_id"):
            out[bid] = g.set_index("gene").ratio_r.reindex(gene_order).fillna(1.0).to_numpy() < 0.01
    return out


def part1(C, a):
    b = wa.build_gates(C["builds"], wa.viability_task_indices(C["tasks"]))
    cols = ["build_id", "line", "lineage", *CFG_KEYS, "failed", "gate_growth", "gate_viability", "gate_repair",
            "gate_all", "viability_passed", "tasks_passed", "wt_ratio_generic", "n_reopened", "repair_frac_reactions",
            "n_open", "sweep_n_lethal", "unreachable_core", "build_cpu_s", "job_wall_s", "status_tag", "p1_eligible"]
    b[[c for c in cols if c in b]].to_csv(OUT / "build_gates.tsv", sep="\t", index=False)
    g = b.groupby(CFG_KEYS, dropna=False)
    summ = g.agg(n_builds=("build_id", "size"), n_failed=("failed", "sum"),
                 wt_ratio_median=("wt_ratio_generic", "median"), wt_ratio_min=("wt_ratio_generic", "min"),
                 viability_pass_rate=("gate_viability", "mean"), tasks_passed_mean=("tasks_passed", "mean"),
                 reopened_median=("n_reopened", "median"), reopened_max=("n_reopened", "max"),
                 repair_frac_max=("repair_frac_reactions", "max"), n_open_median=("n_open", "median"),
                 lethal_median=("sweep_n_lethal", "median"), unreachable_core_median=("unreachable_core", "median"),
                 build_cpu_min_mean=("build_cpu_s", lambda x: x.mean() / 60)).reset_index()
    gates = wa.config_gates(b, CFG_KEYS)
    summ = summ.merge(gates, on=CFG_KEYS + [], how="left", suffixes=("", "_g"))
    summ["status_tag"] = np.where(summ.family == "F3b", "descriptive, not eligible for P1 (amendment A7)",
                                  np.where(summ.family == "F5", "label-calibrated (F5)", ""))
    summ.to_csv(OUT / "config_summary.tsv", sep="\t", index=False)
    return b, summ


def masks_for(b: pd.DataFrame, n_rxn: int) -> np.ndarray:
    return np.vstack([wa.unpack_mask(h, n_rxn) for h in b.open_mask])


def jac(x: np.ndarray, y: np.ndarray) -> float:
    u = (x | y).sum()
    return float((x & y).sum() / u) if u else 1.0


def part1_extra(C, b, a):
    n_rxn = len(C["rxn_ids"])
    genes = sorted(g.id for g in C["model"].genes)
    leth = ko_lethal(["A1", "A6p", "A7"], genes)
    # ---- A7 stability: F2 (GIMME cost noise) and F3 (FASTCORE expression noise), 3 seeds, vs the A1 R1 real build
    rows = []
    a7 = b[(b.row == "A7")]
    for (fam, line), grp in a7.groupby(["family", "line"]):
        M = masks_for(grp, n_rxn)
        L = [leth[x] for x in grp.build_id]
        base = b[(b.row == "A1") & (b.family == fam) & (b.line == line) & (b.rule == "R1") & (b.variant == "real")
                 & (b.genomic == "none")].iloc[0]
        mb = wa.unpack_mask(base.open_mask, n_rxn)
        lb = leth[base.build_id]
        pr = [(i, j) for i in range(len(grp)) for j in range(i + 1, len(grp))]
        rows.append({"family": fam, "line": line, "lineage": grp.lineage.iloc[0], "n_seeds": len(grp),
                     "rxn_jaccard_between_seeds": np.mean([jac(M[i], M[j]) for i, j in pr]),
                     "rxn_jaccard_vs_unperturbed": np.mean([jac(M[i], mb) for i in range(len(grp))]),
                     "lethal_jaccard_between_seeds": np.mean([jac(L[i], L[j]) for i, j in pr]),
                     "lethal_jaccard_vs_unperturbed": np.mean([jac(L[i], lb) for i in range(len(grp))]),
                     "n_open_sd": float(grp.n_open.std()), "n_lethal_sd": float(grp.sweep_n_lethal.std())})
    pd.DataFrame(rows).to_csv(OUT / "a7_stability.tsv", sep="\t", index=False)
    # ---- A6' comparison (label-free part), written into the analysis directory
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from w4_summarise import a6prime_comparison

    a6prime_comparison(STAGE, OUT / "a6prime_comparison.tsv")
    # ---- permuted-expression control (A1): real vs permuted builds, per family x threshold
    rows = []
    a1 = b[(b.row == "A1") & (b.genomic == "none")]
    for (fam, rule), grp in a1.groupby(["family", "rule"]):
        for line, gl in grp.groupby("line"):
            r_ = gl[gl.variant == "real"].iloc[0]
            p_ = gl[gl.variant == "perm"].iloc[0]
            rows.append({"family": fam, "rule": rule, "line": line, "lineage": r_.lineage,
                         "rxn_jaccard_real_perm": jac(wa.unpack_mask(r_.open_mask, n_rxn), wa.unpack_mask(p_.open_mask, n_rxn)),
                         "lethal_jaccard_real_perm": jac(leth[r_.build_id], leth[p_.build_id]),
                         "n_open_real": r_.n_open, "n_open_perm": p_.n_open,
                         "n_lethal_real": r_.sweep_n_lethal, "n_lethal_perm": p_.sweep_n_lethal,
                         "wt_real": r_.wt_ratio_generic, "wt_perm": p_.wt_ratio_generic,
                         "gate_all_real": r_.gate_all, "gate_all_perm": p_.gate_all})
    pc = pd.DataFrame(rows)
    pc.to_csv(OUT / "permuted_control.tsv", sep="\t", index=False)
    # ---- positive controls: TYR/DCT reactions open (Skin vs other lines); urea-cycle reactions open
    idx = {r: i for i, r in enumerate(C["rxn_ids"])}
    rows = []
    real = b[(b.row.isin(["A1", "A2"])) & (b.variant == "real") & (b.genomic == "none")]
    for (fam, rule), grp in real.groupby(["family", "rule"]):
        M = masks_for(grp, n_rxn)
        skin = (grp.lineage == "Skin").to_numpy()
        for sym, rx in C["sym_rxns"].items():
            cols = [idx[r] for r in rx]
            frac = M[:, cols].mean(1) if cols else np.zeros(len(grp))
            rows.append({"family": fam, "rule": rule, "gene": sym, "n_reactions": len(cols),
                         "frac_open_skin": float(frac[skin].mean()) if skin.any() else np.nan,
                         "frac_open_other": float(frac[~skin].mean()),
                         "lines_any_open": int((frac > 0).sum()), "n_lines": len(grp)})
    pd.DataFrame(rows).to_csv(OUT / "positive_controls.tsv", sep="\t", index=False)
    return leth, genes


def part2(C, b, leth, genes, a):
    """Method-vs-tissue variance partition (§9 W4): F0, F2, F3, F4 x {R1, R3} x {real, permuted}, 30 lines."""
    n_rxn = len(C["rxn_ids"])
    a1 = b[(b.row == "A1") & (b.genomic == "none")].sort_values(["variant", "family", "rule", "line"]).reset_index(drop=True)
    levels = {
        "reactions": lambda d: masks_for(d, n_rxn),
        "tasks57": lambda d: np.vstack([np.array([c == "1" for c in s]) for s in d.task_profile]),
        "ko_binary": lambda d: np.vstack([leth[x] for x in d.build_id]),
    }
    rows, h4 = [], {"levels": {}}
    for lev, fn in levels.items():
        H = {}
        for var in ("real", "perm"):
            d = a1[a1.variant == var].reset_index(drop=True)
            D2 = wa.jaccard_distance_sq(fn(d))
            meta = d.rename(columns={"family": "method", "rule": "threshold"})[["method", "threshold", "line", "lineage"]]
            res = {}
            for f, kw in (("method", {"block": "line"}), ("threshold", {"block": "line"}),
                          ("lineage", {"unit": "line"}), ("line", {})):
                pm = wa.permanova(D2, meta, f, perms=a.perms, seed=SEEDS[0], **kw)
                res[f] = pm
                rows.append({"level": lev, "build_set": f"{var} (all methods)", **{k: pm[k] for k in ("factor", "r2", "p_perm", "null_q95")}})
            H[var] = {"D2": D2, "meta": meta, "res": res}
            # within method: lineage R^2 (H4.2)
            for mth in LP:
                sel = (meta.method == mth).to_numpy()
                sub, sm = D2[np.ix_(sel, sel)], meta[sel].reset_index(drop=True)
                pm = wa.permanova(sub, sm, "lineage", perms=a.perms, seed=SEEDS[0], unit="line")
                rows.append({"level": lev, "build_set": f"{var} {mth}", "factor": "lineage", "r2": pm["r2"],
                             "p_perm": pm["p_perm"], "null_q95": pm["null_q95"]})
                H[var][f"lin_{mth}"] = pm
        # H4.1: R2(method) > R2(lineage), line bootstrap of the difference (real builds)
        bt = wa.bootstrap_r2_diff(H["real"]["D2"], H["real"]["meta"], "method", "lineage", B=a.boot, seed=SEEDS[0])
        h41 = {"r2_method": H["real"]["res"]["method"]["r2"], "r2_lineage": H["real"]["res"]["lineage"]["r2"],
               "diff": bt["diff"], "lo": bt["lo"], "hi": bt["hi"], "pass": bool(bt["lo"] > 0),
               # resampled duplicate lines are identical points inside one lineage, which inflates R2(lineage) in
               # every replicate; the basic (pivot) interval corrects that location bias. Both are reported.
               "basic_lo": 2 * bt["diff"] - bt["hi"], "basic_hi": 2 * bt["diff"] - bt["lo"]}
        # H4.2: within each method, R2(lineage | real) > 95th percentile of R2(lineage | permuted builds),
        # the permuted-build null = R2(lineage) on the permuted-expression builds under 999 lineage relabellings
        h42 = {}
        for mth in LP:
            real_r2 = H["real"][f"lin_{mth}"]["r2"]
            perm = H["perm"][f"lin_{mth}"]
            thr = max(perm["null_q95"], perm["r2"])
            h42[mth] = {"r2_lineage_real": real_r2, "r2_lineage_perm": perm["r2"], "perm_null_q95": perm["null_q95"],
                        "threshold": thr, "pass": bool(real_r2 > thr),
                        "label": "tissue-informed" if real_r2 > thr else "method-dominated"}
        h4["levels"][lev] = {"H4.1": h41, "H4.2": h42}
    h4["H4.3"] = {"H4.1": h4["levels"]["ko_binary"]["H4.1"]["pass"],
                  "H4.2_methods_passing": [m for m, v in h4["levels"]["ko_binary"]["H4.2"].items() if v["pass"]]}
    # ---- expression reference: same partition on the input expression (30 panel lines)
    from nextgem.context_inputs import permute_across_lines, read_expression

    lines = sorted(a1.line.unique())
    X = read_expression(C["gt"], lines)
    lin = C["split"].lineage.reindex(X.index)
    for nm, M_ in (("expression_real", X), ("expression_perm", permute_across_lines(X, SEEDS[0]))):
        Z = M_.to_numpy()
        D2 = ((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1)
        pm = wa.permanova(D2, pd.DataFrame({"lineage": lin.to_numpy()}), "lineage", perms=a.perms, seed=SEEDS[0])
        rows.append({"level": "input", "build_set": nm, "factor": "lineage", "r2": pm["r2"], "p_perm": pm["p_perm"],
                     "null_q95": pm["null_q95"]})
        h4[nm + "_r2_lineage"] = pm["r2"]
    # ---- CORDA sensitivity (A7): F0, F2, F3, F4, F3b x {R1, R3} real on the 10 A6' lines (1 line per lineage, so
    # lineage = line here: only R2(method) vs R2(line) is identifiable)
    a6 = pd.read_csv(STAGE / "a6prime_lines.tsv", sep="\t")
    sub = b[(b.line.isin(a6.ModelID)) & (b.genomic == "none") & (b.variant == "real")
            & (((b.row == "A1") & b.family.isin(LP)) | (b.row == "A6p"))].reset_index(drop=True)
    for lev, fn in levels.items():
        for with_corda in (False, True):
            d = sub if with_corda else sub[sub.family != "F3b"].reset_index(drop=True)
            D2 = wa.jaccard_distance_sq(fn(d))
            meta = d.rename(columns={"family": "method", "rule": "threshold"})[["method", "threshold", "line", "lineage"]]
            for f, kw in (("method", {"block": "line"}), ("line", {})):
                pm = wa.permanova(D2, meta, f, perms=a.perms, seed=SEEDS[0], **kw)
                rows.append({"level": lev, "build_set": f"A6' lines, {'with' if with_corda else 'without'} CORDA (sensitivity)",
                             "factor": f, "r2": pm["r2"], "p_perm": pm["p_perm"], "null_q95": pm["null_q95"]})
    vp = pd.DataFrame(rows)
    vp.to_csv(OUT / "variance_partition.tsv", sep="\t", index=False)
    (OUT / "h4_tests.json").write_text(json.dumps(h4, indent=1, default=float))
    return vp, h4


def main(argv=None) -> int:
    from nextgem.run import provenance

    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=999)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--allow-dirty", action="store_true", help="uncommitted code (recorded as dirty=True)")
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = [STAGE / f for f in ("builds.tsv", "ko_A1.parquet", "ko_A6p.parquet", "ko_A7.parquet", "a6prime_lines.tsv",
                                   "MANIFEST.tsv")] + [paths.splits_dir() / "lineage_split.tsv"]
    with provenance(OUT / "MANIFEST.tsv", inputs=inputs, seeds=[SEEDS[0]], solver="none (analysis)", workers=1,
                    allow_dirty=a.allow_dirty, author="W4 engineer agent") as prov:
        t0 = time.time()
        C = load()
        b, summ = part1(C, a)
        leth, genes = part1_extra(C, b, a)
        vp, h4 = part2(C, b, leth, genes, a)
        for f in ("build_gates.tsv", "config_summary.tsv", "a7_stability.tsv", "a6prime_comparison.tsv",
                  "permuted_control.tsv", "positive_controls.tsv", "variance_partition.tsv", "h4_tests.json"):
            prov.output(OUT / f)
        prov.note(script_part="label-free", perms=a.perms, boot=a.boot, wall_s=round(time.time() - t0))
    print(json.dumps({k: v for k, v in h4.items() if k != "levels"}, indent=1, default=float))
    for lev, v in h4["levels"].items():
        print(lev, "H4.1", {k: round(x, 4) if isinstance(x, float) else x for k, x in v["H4.1"].items()},
              "H4.2", {m: r["label"] for m, r in v["H4.2"].items()})
    return 0


if __name__ == "__main__":
    sys.exit(main())
