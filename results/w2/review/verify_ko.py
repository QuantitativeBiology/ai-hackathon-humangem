"""Independent check of the W2 KO sweep engine (review M3-W2). Label-free.

1. M0 sweep row of results/w2/ko/ratio_r_{parent}.parquet == results/baseline/{parent}_ko.parquet ratio_r, bitwise.
2. Sweep definitions rebuilt with W2's own sweep_context/build_sweeps: every parquet sweep id is reproduced.
3. Random gene x sweep pairs re-solved WITHOUT any shortcut (plain Gene.knock_out + FBA) on GLPK (a different
   solver from W2's Gurobi), stratified over shortcut classes; |ratio_r(review) - ratio_r(W2)| reported.
4. --full: one complete no-shortcut sweep per parent (all genes) on Gurobi, 1 session via nextgem.gurobi_pool.

Run: .venv/bin/python results/w2/review/verify_ko.py [--pairs 24] [--full]
Writes results/w2/review/ko_verify_pairs.tsv, ko_verify_full.tsv, ko_verify_summary.tsv.
"""
from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"  # before any model is loaded

import numpy as np
import pandas as pd

from nextgem.knockouts import RATIO_DECIMALS, WT_GROWTH_TOL
from nextgem.model import BIOMASS_REACTION

ROOT = Path(__file__).resolve().parents[3]
W2 = ROOT / "results" / "w2"
OUT = Path(__file__).resolve().parent


def rr(x: float) -> float:
    return float(np.clip(np.round(x, RATIO_DECIMALS), 0.0, 1.0) + 0.0)


def solve(m) -> tuple[float, str]:
    from optlang.interface import INFEASIBLE, OPTIMAL
    m.slim_optimize(error_value=math.nan)
    st = m.solver.status
    if st == OPTIMAL:
        return m.solver.objective.value, st
    return (0.0 if st == INFEASIBLE else math.nan), st


def closures(model) -> dict[str, frozenset]:
    out = {}
    for g in model.genes:
        before = {r.id: r.bounds for r in g.reactions}
        with model:
            g.knock_out()
            out[g.id] = frozenset(r.id for r in g.reactions if r.bounds == (0, 0) and before[r.id] != (0, 0))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=3, help="genes per sampled sweep")
    ap.add_argument("--sweeps", type=int, default=8, help="sweeps sampled per parent")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--seed", type=int, default=20261027)
    a = ap.parse_args()
    summary = []

    # 1. M0 bitwise
    labels = pd.read_csv(W2 / "ko" / "sweep_labels.tsv", sep="\t")
    for p in ("human1", "human2"):
        R = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
        sid = labels[(labels.parent == p) & (labels.variant == "M0")].sweep_id.unique()
        assert len(sid) == 1
        base = pd.read_parquet(ROOT / "results" / "baseline" / f"{p}_ko.parquet").set_index("gene").ratio_r
        w = R.loc[sid[0]]
        same_genes = set(w.index) == set(base.index)
        x, y = w.reindex(base.index).to_numpy(), base.to_numpy()
        bitwise = same_genes and np.array_equal(x, y, equal_nan=True) and x.tobytes() == y.astype(np.float64).tobytes()
        summary.append(dict(check=f"M0_bitwise_{p}", value=bitwise,
                            detail=f"{len(base)} genes; sweep {sid[0]}; same gene set {same_genes}"))

    # 2. rebuild sweeps
    from nextgem import media
    from nextgem.w2_engine import build_sweeps, sweep_context
    lines, sheets, doc, phys, sens, pcs, exchanges, hams = sweep_context()
    sweeps = build_sweeps(exchanges, hams, pcs, doc, phys, sens)
    by_id = {s.sweep_id: s for s in sweeps}
    for p in ("human1", "human2"):
        R = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
        missing = set(R.index) - set(by_id)
        extra = {k for k in by_id if k.startswith(p)} - set(R.index)
        summary.append(dict(check=f"sweeps_rebuilt_{p}", value=not missing and not extra,
                            detail=f"{len(R)} in parquet; missing {len(missing)}; extra {len(extra)}"))

    # 3. random pairs, GLPK, no shortcuts
    from nextgem.w2_runtime import load_parent
    rng = np.random.default_rng(a.seed)
    rows = []
    for p in ("human1", "human2"):
        m = load_parent(p, solver="glpk")
        m.objective = BIOMASS_REACTION
        clos = closures(m)
        R = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
        sw_ids = rng.choice(sorted(R.index), size=a.sweeps, replace=False)
        for sid in sw_ids:
            bounds = by_id[sid].bounds
            row = R.loc[sid]
            with m:
                m.medium = bounds
                wt, wst = solve(m)
            # strata: (i) ratio 1 with non-empty closure (mostly the zero-flux shortcut), (ii) ratio < 1,
            # (iii) a gene sharing its closure set with another gene (shortcut c), (iv) uniform
            g_all = list(row.index)
            s1 = [g for g in g_all if row[g] == 1.0 and clos[g]]
            s2 = [g for g in g_all if row[g] < 1.0]
            setsize = pd.Series([clos[g] for g in g_all]).value_counts()
            s3 = [g for g in s2 if clos[g] and setsize[clos[g]] > 1]
            picks = []
            for pool, tag in ((s1, "ratio1_closure"), (s2, "ratio_lt1"), (s3, "shared_closure"), (g_all, "uniform")):
                pool = [g for g in pool if g not in [q for q, _ in picks]]
                if pool:
                    k = a.pairs if tag in ("ratio1_closure", "ratio_lt1") else 1
                    for g in rng.choice(pool, size=min(k, len(pool)), replace=False):
                        picks.append((g, tag))
            for g, tag in picks:
                t0 = time.perf_counter()
                with m:
                    m.medium = bounds
                    m.genes.get_by_id(g).knock_out()
                    kg, kst = solve(m)
                ratio = kg / wt if wt > WT_GROWTH_TOL else math.nan
                rows.append(dict(parent=p, sweep_id=sid, gene=g, stratum=tag, n_closed=len(clos[g]),
                                 wt_glpk=wt, ko_glpk=kg, ko_status=kst, ratio_r_review=rr(ratio),
                                 ratio_r_w2=float(row[g]), abs_diff=abs(rr(ratio) - float(row[g])),
                                 secs=round(time.perf_counter() - t0, 2)))
        del m
    pairs = pd.DataFrame(rows)
    pairs.to_csv(OUT / "ko_verify_pairs.tsv", sep="\t", index=False)
    summary.append(dict(check="pairs_glpk_no_shortcut", value=bool((pairs.abs_diff <= 1e-5).all()),
                        detail=f"{len(pairs)} pairs over {pairs.sweep_id.nunique()} sweeps; max |diff| "
                               f"{pairs.abs_diff.max():.2e}; strata {pairs.stratum.value_counts().to_dict()}"))

    # 4. full no-shortcut sweep per parent on Gurobi (1 session)
    if a.full:
        from nextgem.gurobi_pool import hold_slot
        slot = hold_slot(timeout=3600)
        full = []
        try:
            for p in ("human1", "human2"):
                m = load_parent(p, solver="gurobi")
                m.objective = BIOMASS_REACTION
                R = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
                lab = labels[(labels.parent == p) & (labels.variant == "M2") & (labels.medium_key == "RPMI")
                             & (labels.f == 0.1)]
                sid = lab.sweep_id.iloc[0]
                bounds = by_id[sid].bounds
                t0 = time.perf_counter()
                with m:
                    m.medium = bounds
                    wt, _ = solve(m)
                    for g in R.columns:
                        with m:
                            m.genes.get_by_id(g).knock_out()
                            kg, kst = solve(m)
                        full.append(dict(parent=p, sweep_id=sid, gene=g, ratio_r_review=rr(kg / wt),
                                         ratio_r_w2=float(R.loc[sid, g]), status=kst))
                df = pd.DataFrame([r for r in full if r["parent"] == p])
                d = (df.ratio_r_review - df.ratio_r_w2).abs()
                summary.append(dict(check=f"full_sweep_gurobi_{p}", value=bool((d == 0).all()),
                                    detail=f"{sid} (RPMI M2 f=0.1): {len(df)} genes, n_diff {(d > 0).sum()}, "
                                           f"max |diff| {d.max():.2e}, {time.perf_counter() - t0:.0f} s"))
                del m
        finally:
            slot.release()
        fd = pd.DataFrame(full)
        fd[fd.ratio_r_review != fd.ratio_r_w2].to_csv(OUT / "ko_verify_full.tsv", sep="\t", index=False)
    s = pd.DataFrame(summary)
    s.to_csv(OUT / "ko_verify_summary.tsv", sep="\t", index=False)
    print(s.to_string())


if __name__ == "__main__":
    main()
