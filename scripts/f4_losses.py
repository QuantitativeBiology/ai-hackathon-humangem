"""F4a step 1 (label-free): the pre-registered genomic loss rule and the order-independent drop on the P1 parent, per
training line (A15 post hoc; RESEARCH_TREE D19; DESIGN §W4 genomic constraint).

Lines: the W1/M6/F1 primary training set (training tier, non-rare, with expression; 686 lines, 12 lineages).
Per line:
  1. lost model genes: ``nextgem.context_inputs.lost_genes`` (rule r2-11) on DepMap damaging mutations, WGS relative CN
     and log2(TPM+1) (default rows; Human1 genes by Entrez, ``context_inputs.gene_table``); no CN -> mutations only
     (flag ``cn_missing``);
  2. ``nextgem.context.accept_losses`` on the P1 parent (Human1 + W3 B-loo, the line's M2(f=1) bounds; GLPK, A10 LP
     limit): each loss alone (a loss that alone kills growth, < 1% of the parent's WT, is model-contradicted and
     dropped), the rest jointly, then the expression-ordered joint drop;
  3. ``context.loss_closure`` of the accepted set = the reactions F4a closes; the line's F4a sweep job is
     (signature, accepted set), mode 'full' if >= 1 reaction is closed, else 'partial' (f4_common.sweep_mode).
One process per M2(f=1) signature (one parent load each), <= 12 GLPK workers.
Check: each signature parent's WT growth equals the pinned P1 sweep's WT growth (results/m6/sweeps/sweeps.tsv or W3).

No label is read. Outputs (results/f4/losses/): line_losses.tsv, gene_counts.tsv, summary.json; MANIFEST row.

Usage: .venv/bin/python scripts/f4_losses.py --workers 12 [--allow-dirty] [--lines N --out results/f4/tmp/...]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import multiprocessing
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f4_common as fc  # noqa: E402

fc.blas_single_thread()


def _sig_job(sig: str, bounds: dict, items: list) -> dict:
    """Worker: one P1 parent on ``sig``; accept_losses for every (line, lost, expr) of that signature."""
    import signal

    from nextgem import context as cx

    def _on_prof(signum, frame):
        raise TimeoutError("job CPU limit")

    signal.signal(signal.SIGPROF, _on_prof)
    signal.setitimer(signal.ITIMER_PROF, float(fc.JOB_CPU_LIMIT_S))
    t0, c0 = time.time(), time.process_time()
    try:
        parent = fc.p1_parent(bounds, sig)
        load_cpu = time.process_time() - c0
        rows = []
        single = {}
        for line, lost, expr in items:
            c1 = time.process_time()
            r = cx.accept_losses(parent, lost, pd.Series(expr, dtype="float64"))
            closed = cx.loss_closure(parent, set(r["accepted"]))
            for g in r["lost"]:
                if g not in single:
                    single[g] = len(cx.loss_closure(parent, {g}))
            rows.append({"ModelID": line, "w2_sweep_id": sig, "lost": r["lost"], "contradicted": r["contradicted"],
                         "dropped_joint": r["dropped_joint"], "accepted": r["accepted"],
                         "joint_growth_ratio": r["joint_growth"] / parent.wt_generic, "closed_rxns": sorted(closed),
                         "cpu_s": time.process_time() - c1})
        return {"sig": sig, "rows": rows, "single_closure": single, "wt": parent.wt_generic, "solver": parent.solver,
                "load_cpu_s": load_cpu, "cpu_s": time.process_time() - c0, "wall_s": time.time() - t0}
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0.0)


def p1_wt_by_sig(sigs) -> dict:
    """Pinned P1 WT growth per signature: M6 new sweeps (sweeps.tsv) or the W3 B-loo sweep (classic-RPMI signature)."""
    m6 = pd.read_csv(fc.PROJECT / "results" / "m6" / "sweeps" / "sweeps.tsv", sep="\t")
    m6 = m6[(m6.parent == fc.PARENT) & (m6.variant == fc.BLOO) & (m6.status == "ok")]
    out = m6.groupby("w2_sweep_id").wt_growth.first().to_dict()
    w3dir = fc.PROJECT / "results" / "w3"
    ss = json.loads((w3dir / "structure" / "structure_summary.json").read_text())["w2_sweep_ids"][fc.PARENT]
    w3s = pd.read_csv(w3dir / "sweeps" / f"sweeps_{fc.PARENT}.tsv", sep="\t")
    hit = w3s[(w3s["medium"] == "RPMI_M2f1") & (w3s.variant == fc.BLOO)]
    if len(hit) == 1:
        out.setdefault(ss["RPMI_M2f1"], float(hit.wt_growth.iloc[0]))
    return {s: out.get(s, np.nan) for s in sigs}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=fc.MAX_WORKERS)
    ap.add_argument("--lines", type=int, help="pilot: first N lines only (needs --out under tmp/)")
    ap.add_argument("--out", help="output directory (default results/f4/losses)")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    if not 1 <= a.workers <= fc.MAX_WORKERS:
        raise SystemExit(f"--workers must be in 1..{fc.MAX_WORKERS}")
    out = fc.LOSSDIR if not a.out else (Path(a.out) if a.out.startswith("/") else fc.PROJECT / a.out)
    if a.lines and not fc.in_tmp(out):
        raise SystemExit("--lines writes only under a tmp/ directory")
    from nextgem import context_inputs as ci
    from nextgem import paths
    from nextgem.model import GEMS

    dm = paths.depmap_dir()
    inputs = [dm / ci.EXPR_FILE, dm / ci.CN_FILE, dm / ci.DAMAGING_FILE, GEMS[fc.PARENT].sbml,
              GEMS[fc.PARENT].genes_tsv, paths.splits_dir() / "lineage_split.tsv",
              fc.PROJECT / "results" / "m6" / "scores" / "line_assignment.tsv",
              fc.PROJECT / "results" / "w2" / "controls" / "line_sweeps.tsv",
              fc.PROJECT / "results" / "w3" / "structure" / "variants_human1.json",
              fc.PROJECT / "results" / "m6" / "sweeps" / "sweeps.tsv"]
    manifest = fc.MANIFEST if not fc.in_tmp(out) else out / "MANIFEST.tsv"
    with fc.provenance(inputs=inputs, solver="glpk", workers=a.workers, allow_dirty=a.allow_dirty or bool(a.lines),
                       manifest=manifest) as prov:
        outs, summ = run(a, out)
        prov.output(*outs)
        prov.note(step="f4_losses (F4a label-free loss rule, D19)",
                  **{k: summ[k] for k in ("n_lines", "n_lines_lost_ge1", "n_lines_accepted_ge1", "n_lines_closing_ge1",
                                          "n_sweep_jobs", "job_cpu_s")})
    return 0


def run(a, out: Path):
    import cobra

    from nextgem import context_inputs as ci
    from nextgem import labels, w1
    from nextgem.model import GEMS

    t_start = time.time()
    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    E = fc.scored_training_lines(split, lf)
    if a.lines:
        E = E[: a.lines]
    sig = fc.line_signatures().reindex(E)
    if sig.isna().any():
        raise RuntimeError(f"{int(sig.isna().sum())} lines without an M2(f=1) signature")
    cobra.Configuration().solver = "glpk"
    m = cobra.io.read_sbml_model(str(GEMS[fc.PARENT].sbml))
    genes = ci.gene_table(m, GEMS[fc.PARENT].genes_tsv)
    del m
    rna = ci.read_expression(genes, E)
    cn = ci.read_cn(genes, E)
    dmg = ci.read_damaging(genes, E)
    print(f"{len(E)} lines; RNA {rna.shape}, CN {cn.shape}, damaging {dmg.shape}; {time.time() - t_start:.0f} s",
          flush=True)
    lost, flags = {}, {}
    for l in E:
        lost[l], flags[l] = ci.lost_genes(l, dmg, cn, rna, rule=fc.LOSS_RULE)
    bounds, info = fc.signature_bounds(sig)
    groups = {}
    for l in E:
        expr = rna.loc[l].dropna().to_dict() if l in rna.index else {}
        groups.setdefault(sig[l], []).append((l, lost[l], expr))
    res = []
    ctx = multiprocessing.get_context("spawn")
    with cf.ProcessPoolExecutor(max_workers=min(a.workers, len(groups)), mp_context=ctx) as pool:
        futs = {pool.submit(_sig_job, s, bounds[s], items): s for s, items in groups.items()}
        for f in cf.as_completed(futs):
            r = f.result()
            res.append(r)
            print(f"  {r['sig']}: {len(r['rows'])} lines, WT {r['wt']:.6g}, {r['cpu_s']:.0f} s CPU, "
                  f"{time.time() - t_start:.0f} s", flush=True)
    # ---- tables
    rows, single = [], {}
    for r in res:
        single.update(r["single_closure"])
        for x in r["rows"]:
            rows.append(x)
    L = pd.DataFrame(rows).set_index("ModelID").reindex(E)
    L.insert(0, "lineage", split.lineage.reindex(E))
    L["cn_missing"] = [flags[l]["cn_missing"] for l in E]
    L["mut_missing"] = [flags[l]["mut_missing"] for l in E]
    for c in ("lost", "contradicted", "dropped_joint", "accepted", "closed_rxns"):
        L[f"n_{c}"] = L[c].map(len)
    L["job_key"] = [fc.job_key(L.at[l, "w2_sweep_id"], L.at[l, "accepted"]) if L.at[l, "n_accepted"] else ""
                    for l in E]
    L["sweep_mode"] = [fc.sweep_mode(n) if acc else "P1" for n, acc in zip(L.n_closed_rxns, L.n_accepted)]
    for c in ("lost", "contradicted", "dropped_joint", "accepted", "closed_rxns"):
        L[c] = L[c].map(fc.join_genes)
    cols = ["lineage", "w2_sweep_id", "cn_missing", "mut_missing", "n_lost", "n_contradicted", "n_dropped_joint",
            "n_accepted", "n_closed_rxns", "joint_growth_ratio", "sweep_mode", "job_key", "lost", "contradicted",
            "dropped_joint", "accepted", "closed_rxns", "cpu_s"]
    L = L[cols].rename_axis("ModelID").reset_index()
    # gene table
    cnt = {}
    for kind in ("lost", "accepted", "contradicted", "dropped_joint"):
        s = pd.Series([g for v in L[kind] for g in fc.split_genes(v)], dtype=str).value_counts()
        cnt[f"n_lines_{kind}"] = s
    G = pd.DataFrame(cnt).fillna(0).astype(int)
    G.index.name = "gene"
    G.insert(0, "symbol", genes.symbol.reindex(G.index))
    G.insert(1, "entrez", genes.entrez.reindex(G.index))
    G["n_rxns_closed_alone"] = [single.get(g, np.nan) for g in G.index]
    G = G.sort_values(["n_lines_lost", "symbol"], ascending=[False, True]).reset_index()
    # ---- checks
    pin = p1_wt_by_sig(groups)
    wt = {r["sig"]: r["wt"] for r in res}
    chk = []
    for s in sorted(groups):
        ref = pin[s]
        ok = bool(np.isfinite(ref) and abs(wt[s] - ref) <= 1e-5 * max(1.0, abs(ref)))
        chk.append({"w2_sweep_id": s, "medium_keys": info[s]["medium_keys"], "n_lines": len(groups[s]),
                    "parent_wt": wt[s], "pinned_p1_wt": ref, "pass": ok if np.isfinite(ref) else None})
    chk = pd.DataFrame(chk)
    if (chk["pass"] == False).any():  # noqa: E712  (None = no pinned M6 row: the W3 classic-RPMI sweep, checked below)
        raise RuntimeError("P1 parent WT differs from the pinned P1 sweep:\n" + chk.to_string())
    jobs = L[L.n_accepted > 0].drop_duplicates("job_key")
    top = G.head(25)
    summ = {
        "post_hoc": "A15 post hoc / exploratory", "rule": fc.LOSS_RULE, "n_lines": len(E),
        "n_lineages": int(L.lineage.nunique()), "n_signatures": len(groups),
        "n_lines_cn_missing": int(L.cn_missing.sum()), "n_lines_mut_missing": int(L.mut_missing.sum()),
        "n_lines_lost_ge1": int((L.n_lost > 0).sum()), "n_lines_accepted_ge1": int((L.n_accepted > 0).sum()),
        "n_lines_closing_ge1": int((L.n_closed_rxns > 0).sum()),
        "n_lines_contradicted_ge1": int((L.n_contradicted > 0).sum()),
        "n_lines_dropped_joint_ge1": int((L.n_dropped_joint > 0).sum()),
        "n_lost_per_line": {"median": float(L.n_lost.median()), "mean": float(L.n_lost.mean()), "max": int(L.n_lost.max())},
        "n_closed_per_line": {"median": float(L.n_closed_rxns.median()), "max": int(L.n_closed_rxns.max())},
        "n_distinct_lost_genes": int((G.n_lines_lost > 0).sum()),
        "n_sweep_jobs": int(len(jobs)), "n_sweep_jobs_full": int((jobs.sweep_mode == "full").sum()),
        "n_sweep_jobs_partial": int((jobs.sweep_mode == "partial").sum()),
        "top_lost_genes": [{"gene": r.gene, "symbol": r.symbol, "n_lines_lost": int(r.n_lines_lost),
                            "n_lines_accepted": int(r.n_lines_accepted),
                            "n_rxns_closed_alone": None if pd.isna(r.n_rxns_closed_alone) else int(r.n_rxns_closed_alone)}
                           for r in top.itertuples()],
        "contradicted_genes": {r.symbol or r.gene: int(r.n_lines_contradicted) for r in G.itertuples()
                               if r.n_lines_contradicted > 0},
        "parent_wt_check": {"n_checked": int(chk["pass"].notna().sum()),
                            "n_pass": int((chk["pass"] == True).sum()),  # noqa: E712
                            "unchecked": list(chk[chk["pass"].isna()].w2_sweep_id)},
        "job_cpu_s": round(float(sum(r["cpu_s"] for r in res)), 1), "wall_s": round(time.time() - t_start, 1),
    }
    outs = [fc.atomic_write_tsv(L, out / "line_losses.tsv"), fc.atomic_write_tsv(G, out / "gene_counts.tsv"),
            fc.atomic_write_tsv(chk, out / "parent_check.tsv"), fc.atomic_write_json(summ, out / "summary.json")]
    print(json.dumps({k: v for k, v in summ.items() if k != "top_lost_genes"}, indent=1, default=float))
    print(top.to_string(index=False))
    return outs, summ


if __name__ == "__main__":
    raise SystemExit(main())
