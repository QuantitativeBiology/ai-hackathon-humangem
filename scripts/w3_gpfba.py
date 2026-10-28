"""A8 exploratory score: gene-pFBA KO-cost tie-break on the generic parents (DESIGN §17 A8;
notes/methods/GPRTransform_assessment.md §7 B). Label-free; GLPK; M0.

Transform: results/gprt/proto/gpr_transform.py (reused as is, imported by path), exact "arm" isozyme encoding.
  w_g = 1 (generic parents);
  μ* = WT max biomass; C* = min Σ_g w_g·u_g subject to biomass >= 0.9·μ*;
  μ_g = KO max biomass (u_g = 0); C_g = min Σ w·u subject to u_g = 0 and biomass >= 0.9·μ_g;
  score_g = round(ratio_r − 1e-3·min(1, (C_g/μ_g)/(C*/μ*) − 1), 6); lethal KOs (μ_g ≈ 0) keep score = ratio_r.
Exact shortcut: a gene with u_g = 0 in the WT min-cost solution and ratio_r = 1 has C_g = C* (that solution stays
feasible and optimal), so its score is ratio_r. Every other gene with μ_g > 0 is solved.
ratio_r is taken from the transformed model's own KO (equal to results/baseline by A8 / the assessment; checked).

Output: results/w3/gpfba/gpfba_<parent>.parquet (gene, entrez, ratio_r, mu_ko, cost_ko, score, solved) and
gpfba_summary.json. The score is "exploratory": it never replaces ratio_r and enters no gate or selection.

Usage: .venv/bin/python scripts/w3_gpfba.py --workers 4
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
from w3_common import MAX_WORKERS, OUT_DIR, PROJECT, atomic_write_parquet, current, run_main  # noqa: E402

GDIR = OUT_DIR / "gpfba"
PROTO = PROJECT / "results" / "gprt" / "proto"
GROWTH_FRACTION = 0.9
COST_STEP = 1e-3
_W: dict = {}


def _build(parent):
    import cobra

    cobra.Configuration().solver = "glpk"
    sys.path.insert(0, str(PROTO))
    import gpr_transform as gt
    from nextgem.model import hams_medium, load_model

    m = load_model(parent)
    m.medium = hams_medium(m)
    m.objective = "MAR13082"
    tm, info = gt.transform(m, split="arm")
    tm.objective = "MAR13082"
    return tm, info, gt


def _cost_objective(tm, info, gt):
    return tm.problem.Objective(sum(tm.reactions.get_by_id(gt.usage(g)).forward_variable for g in info["genes"]),
                                direction="min")


def _init(parent):
    tm, info, gt = _build(parent)
    _W.update(tm=tm, info=info, gt=gt)


def _growth_and_cost(gene):
    """(mu, cost) for KO of ``gene`` (None = WT); cost at biomass >= 0.9 mu; cost NaN if mu ~ 0."""
    from nextgem import w3

    tm, info, gt = _W["tm"], _W["info"], _W["gt"]
    t0 = time.process_time()
    with tm:
        if gene is not None:
            gt.knock_out(tm, gene)
        tm.objective = "MAR13082"
        st = w3.solve_status(tm)
        mu = float(tm.solver.objective.value) if st == "optimal" else 0.0
        cost, used = np.nan, None
        if mu > 1e-6:
            tm.reactions.get_by_id("MAR13082").lower_bound = GROWTH_FRACTION * mu
            tm.objective = _cost_objective(tm, info, gt)
            st2 = w3.solve_status(tm)
            if st2 != "optimal":
                raise w3.SolverInfraError(f"min-cost LP {st2} for {gene}")
            cost = float(tm.solver.objective.value)
            if gene is None:
                pv = tm.solver.primal_values
                used = [g for g in info["genes"] if pv[tm.reactions.get_by_id(gt.usage(g)).forward_variable.name] > 1e-9]
    return gene, mu, cost, used, time.process_time() - t0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--model", nargs="*", default=["human1", "human2"])
    a = ap.parse_args(argv)
    if a.workers > MAX_WORKERS:
        raise SystemExit("at most 8 worker processes")
    from nextgem import paths

    t0 = time.time()
    summary = {}
    for p in a.model:
        tp = time.time()
        base = pd.read_parquet(paths.baseline_dir() / f"{p}_ko.parquet", columns=["gene", "entrez", "symbol", "ratio_r"])
        ctx = multiprocessing.get_context("spawn")
        with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init, initargs=(p,)) as pool:
            _, mu_wt, c_wt, used, cpu_wt = pool.submit(_growth_and_cost, None).result()
            used = set(used)
            todo = [g for g, r in zip(base.gene, base.ratio_r) if (g in used) or (0 < r < 1)]
            res = {}
            cpu = cpu_wt
            for n, (g, mu, c, _, cp) in enumerate(pool.map(_growth_and_cost, todo, chunksize=4), 1):
                res[g] = (mu, c)
                cpu += cp
                if n % 100 == 0:
                    print(f"[{p}] {n}/{len(todo)} genes, {time.time() - tp:.0f} s", flush=True)
        rows = []
        unit_wt = c_wt / mu_wt
        mism = 0
        for g, e, s, r in zip(base.gene, base.entrez, base.symbol, base.ratio_r):
            if g in res:
                mu, c = res[g]
                rr = float(np.clip(np.round(mu / mu_wt, 6), 0, 1)) + 0.0
                mism += int(abs(rr - r) > 1e-6)
                if mu > 1e-6 and np.isfinite(c):
                    x = (c / mu) / unit_wt - 1.0
                    score = round(r - COST_STEP * min(1.0, x), 6)
                else:
                    x, score = np.nan, r
                rows.append({"gene": g, "entrez": e, "symbol": s, "ratio_r": r, "mu_ko": mu, "cost_ko": c,
                             "rel_cost_minus1": x, "score": score, "solved": True})
            else:
                rows.append({"gene": g, "entrez": e, "symbol": s, "ratio_r": r, "mu_ko": np.nan, "cost_ko": np.nan,
                             "rel_cost_minus1": 0.0 if r == 1 else np.nan, "score": r, "solved": False})
        out = pd.DataFrame(rows)
        f = atomic_write_parquet(out, GDIR / f"gpfba_{p}.parquet")
        n_dist_r = out.ratio_r.round(6).nunique()
        n_dist_s = out.score.round(6).nunique()
        summary[p] = {"mu_wt": mu_wt, "cost_wt": c_wt, "n_used_genes_wt": len(used), "n_solved": len(todo),
                      "ratio_mismatch_vs_baseline": mism, "distinct_ratio_r": int(n_dist_r), "distinct_score": int(n_dist_s),
                      "n_genes_score_ne_ratio": int((out.score != out.ratio_r).sum()),
                      "n_negative_rel_cost": int((out.rel_cost_minus1 < -1e-9).sum()),
                      "tied_at_1_ratio": int((out.ratio_r == 1).sum()), "tied_at_1_score": int((out.score == 1).sum()),
                      "cpu_s_lp": round(cpu, 1), "wall_s": round(time.time() - tp, 1)}
        if mism:
            raise RuntimeError(f"{p}: {mism} transformed KO ratios differ from results/baseline")
        current().input(paths.baseline_dir() / f"{p}_ko.parquet", PROTO / "gpr_transform.py")
        current().output(f)
        print(json.dumps({p: summary[p]}, indent=1), flush=True)
    (GDIR / "gpfba_summary.json").write_text(json.dumps(summary, indent=1))
    current().output(GDIR / "gpfba_summary.json")
    return 0


if __name__ == "__main__":
    w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 4
    raise SystemExit(run_main(main, solver="glpk", workers=w))
