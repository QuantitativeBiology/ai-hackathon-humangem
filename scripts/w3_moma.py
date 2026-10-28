"""W3 O-MOMA (DESIGN §9 W3): QP MOMA built once with cobra ``add_moma`` per (parent, pFBA-WT reference) on M0, KO
bounds toggled. Gurobi only, through ``nextgem.gurobi_pool`` (one slot per worker process; W3 holds <= 3).

Score: ratio_r of the biomass flux of the MOMA solution (KO biomass flux / WT biomass flux, rounded to 6 decimals,
clipped to [0, 1]); DESIGN §6.1. Exact shortcut: a KO that closes only reactions with zero flux in the reference
leaves the (unique) MOMA optimum at the reference, so ratio 1. Genes whose KOs close the identical reaction set
share one QP.

Reference: pFBA on M0 computed once on GLPK in the parent process (so every worker uses the identical vector) and
written to results/w3/moma/reference_<parent>.parquet.

--pilot N: time N QPs per parent (KO sets drawn with SEEDS[0]) and extrapolate to all KO sets that touch the
reference support; also the §12 check ``add_moma`` (toggled bounds) == ``cobra.flux_analysis.moma`` on 5 genes.
The frozen rule: O-MOMA is dropped if the projected cost exceeds 1 CPU-h per model.
Without --pilot: the full sweep (only run if the pilot passes); QPs that fail count as NaN (>0.1% -> arm dropped).

Usage: .venv/bin/python scripts/w3_moma.py --pilot 30 --workers 3
       .venv/bin/python scripts/w3_moma.py --workers 3 [--model human2]
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
from w3_common import OUT_DIR, W3_GUROBI_SLOTS, atomic_write_tsv, current, run_main  # noqa: E402

MDIR = OUT_DIR / "moma"
QP_TIME_LIMIT = 300.0
_W: dict = {}


def _init(parent):
    from nextgem import w3
    from nextgem.config import SEEDS
    from nextgem.w2_runtime import load_parent  # holds one pool slot for this process (Gurobi)
    from cobra.flux_analysis import add_moma
    from nextgem.model import hams_medium

    m = load_parent(parent, solver="gurobi")
    p = m.solver.problem.Params
    p.Threads, p.Seed, p.TimeLimit = 1, SEEDS[0], QP_TIME_LIMIT
    m.medium = hams_medium(m)
    m.objective = w3.BIOMASS
    ref = pd.read_parquet(MDIR / f"reference_{parent}.parquet").set_index("reaction").flux
    import cobra

    sol = cobra.Solution(objective_value=float(ref[w3.BIOMASS]), status="optimal", fluxes=ref)
    _W["ref"] = sol
    _W["wt_bio"] = float(ref[w3.BIOMASS])
    add_moma(m, solution=sol, linear=False)
    _W["model"] = m


def _qp(ko_rxns):
    m = _W["model"]
    t0 = time.process_time()
    with m:
        for rid in ko_rxns:
            m.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        m.slim_optimize()
        st = m.solver.status
        bio = m.reactions.get_by_id("MAR13082").flux if st == "optimal" else np.nan
    return list(ko_rxns), st, bio, time.process_time() - t0


def _check5(genes):
    """§12: add_moma with toggled bounds equals cobra.flux_analysis.moma (fresh build per gene) on 5 genes."""
    from cobra.flux_analysis import moma

    out = []
    m = _W["model"]
    for g in genes:
        rx = [r.id for r in m.genes.get_by_id(g).reactions if not r.gpr.eval([g])]
        _, st, bio_toggled, _ = _qp(rx)
        from nextgem.w2_runtime import load_parent
        from nextgem.model import hams_medium

        fresh = _W.get("fresh")
        if fresh is None:
            fresh = load_parent(m.notes.get("nextgem_release"), solver="gurobi")
            fresh.solver.problem.Params.Threads = 1
            fresh.medium = hams_medium(fresh)
            fresh.objective = "MAR13082"
            _W["fresh"] = fresh
        with fresh:
            fresh.genes.get_by_id(g).knock_out()
            s = moma(fresh, solution=_W["ref"], linear=False)
            bio_ref = s.fluxes["MAR13082"]
        out.append((g, bio_toggled, float(bio_ref)))
    return out


def reference(parent):
    """pFBA reference on M0 (GLPK, parent process)."""
    import cobra

    cobra.Configuration().solver = "glpk"
    from cobra.flux_analysis import pfba
    from nextgem import w3
    from nextgem.model import hams_medium, load_model

    m = load_model(parent)
    m.medium = hams_medium(m)
    m.objective = w3.BIOMASS
    sol = pfba(m)
    ref = sol.fluxes.rename("flux").rename_axis("reaction").reset_index()
    MDIR.mkdir(parents=True, exist_ok=True)
    f = MDIR / f"reference_{parent}.parquet"
    ref.to_parquet(f, index=False)
    gd = w3.gene_disabled_reactions(m)
    groups = w3.ko_groups(gd)
    support = set(ref.reaction[ref.flux.abs() > w3.ZERO_FLUX_TOL])
    touch = sorted([k for k in groups if k & support], key=lambda s: tuple(sorted(s)))
    return f, groups, touch, sorted(g.id for g in m.genes), gd


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", type=int, default=0)
    ap.add_argument("--workers", type=int, default=W3_GUROBI_SLOTS)
    ap.add_argument("--model", nargs="*", default=["human1", "human2"])
    a = ap.parse_args(argv)
    if a.workers > W3_GUROBI_SLOTS:
        raise SystemExit(f"W3 holds at most {W3_GUROBI_SLOTS} Gurobi slots")
    from nextgem.config import SEEDS

    t0 = time.time()
    summary = {}
    for p in a.model:
        f, groups, touch, genes, gd = reference(p)
        jobs = touch
        if a.pilot:
            rng = np.random.default_rng(SEEDS[0])
            jobs = [touch[i] for i in sorted(rng.choice(len(touch), size=min(a.pilot, len(touch)), replace=False))]
        res = []
        ctx = multiprocessing.get_context("spawn")
        tp = time.time()
        with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init, initargs=(p,)) as pool:
            for n, r in enumerate(pool.map(_qp, jobs, chunksize=1), 1):
                res.append(r)
                if n % 25 == 0:
                    print(f"[{p}] {n}/{len(jobs)} QPs, {time.time() - tp:.0f} s", flush=True)
            chk = list(pool.submit(_check5, sorted({g for k in jobs[:5] for g, rx in gd.items() if rx == k})[:5]).result()) if a.pilot else []
        wt = float(pd.read_parquet(f).set_index("reaction").flux["MAR13082"])
        cpu = [r[3] for r in res]
        st = pd.Series([r[1] for r in res]).value_counts().to_dict()
        rec = {"n_ko_sets_total": len(groups), "n_ko_sets_touching_reference": len(touch), "n_qp": len(res),
               "qp_cpu_s_mean": float(np.mean(cpu)), "qp_cpu_s_max": float(np.max(cpu)), "status_counts": st,
               "projected_cpu_h_per_model": float(np.mean(cpu)) * len(touch) / 3600.0,
               "wall_s": round(time.time() - tp, 1), "wt_biomass_ref": wt}
        if chk:
            rec["check_add_moma_vs_moma"] = [{"gene": g, "toggled": b1, "fresh": b2, "abs_diff": abs(b1 - b2)} for g, b1, b2 in chk]
            rec["check_pass"] = all(abs(b1 - b2) <= 1e-6 * max(1.0, abs(b2)) for _, b1, b2 in chk)
        rec["rule_drop_if_gt_1_cpu_h"] = rec["projected_cpu_h_per_model"] > 1.0
        summary[p] = rec
        print(json.dumps({p: rec}, indent=1, default=float), flush=True)
        if not a.pilot:
            ko_ratio = {}
            for rx, s, bio, _ in res:
                # infeasible KO network -> ratio 0 (as the FBA engine); any other status is a failed solve (NaN)
                ko_ratio[frozenset(rx)] = (float(np.clip(np.round(bio / wt, 6), 0, 1)) + 0.0 if s == "optimal"
                                           else (0.0 if s == "infeasible" else np.nan))
            rows = []
            for g in genes:
                k = gd.get(g, frozenset())
                rows.append({"gene": g, "ratio_r": ko_ratio.get(k, 1.0) if k else 1.0,
                             "solved": bool(k) and k in ko_ratio})
            out = pd.DataFrame(rows)
            nan_frac = float(out.ratio_r.isna().mean())
            rec["failed_fraction"] = nan_frac
            rec["rule_drop_if_failed_gt_0.1pct"] = nan_frac > 0.001
            out.to_parquet(MDIR / f"moma_{p}.parquet", index=False)
            current().output(MDIR / f"moma_{p}.parquet")
        current().output(f)
    tag = "pilot" if a.pilot else "full"
    (MDIR / f"moma_{tag}.json").write_text(json.dumps(summary, indent=1, default=float))
    current().output(MDIR / f"moma_{tag}.json")
    current().note(total_wall_s=round(time.time() - t0, 1))
    return 0


if __name__ == "__main__":
    w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else W3_GUROBI_SLOTS
    raise SystemExit(run_main(main, seeds=[11], solver="gurobi", workers=w))
