"""F1 step 2 (label-free): KO sweeps of the curated P1 arms on every Human1 M2(f=1) medium signature (A15, post hoc).

Arms (scripts/f1_common.py): TYMS only, GPR rule only, F1 = GPR rule + TYMS; all on Human1 + W3
B-loo:cofactor_pool_biomass (the W3 spec, ``nextgem.w3.Variant``) under each M2(f=1) bounds signature used by any
line (W2 ids from results/w2/controls/line_sweeps.tsv; bounds rebuilt with the M6 code, ``m6_sweeps.bounds_by_id``).
The P1 reference rows are the pinned M6/W3 sweeps (``m6_scores.sweep_ratios``), not recomputed.

Engine: ``nextgem.w3.ko_sweep`` (the W2/W3/M6 KO engine, zero-flux and identical-closure shortcuts), GLPK, A10
limits (300 s per LP + one retry inside ``w3.solve_status``; 1,800 s CPU per sweep job via ITIMER_PROF; a timeout,
a non-optimal LP after the retry or a failed wild type is a failed build, never re-run with relaxed limits).
- TYMS: full sweep (bounds change the LP, so every gene can move).
- GPR, F1: GPR edits change only the KO closures of the genes listed in results/f1/closure_changes_human1.tsv and
  leave the wild-type LP unchanged, so only those genes are re-solved; every other gene keeps the ratio of the
  arm's base sweep (P1 for GPR, TYMS for F1). Exactness is checked by full spot-check sweeps (bitwise).
- Spot checks: P1 full sweep on 2 signatures == the M6/W3 rows; F1 full sweep on the same 2 == composed F1;
  GPR full sweep on 1 == composed GPR.

Outputs (results/f1/sweeps/):
  ratio_r_delta_human1.parquet  long (arm, w2_sweep_id, gene, ratio_r, ratio_r_P1) for every (arm, signature,
                                gene) whose ratio differs from P1 (compact; :func:`f1_common`-style reconstruction
                                in :func:`arm_ratios`)
  sweeps.tsv (one row per job), spot_checks.tsv, sweeps.json

Usage: .venv/bin/python scripts/f1_sweeps.py --workers 8 [--allow-dirty]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import multiprocessing
import signal
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f1_common as fc  # noqa: E402
import m6_common as mc  # noqa: E402

W2 = mc.PROJECT / "results" / "w2"
SPOT_MEDIA = ("RPMI", "DMEM")  # first medium key of the spot-check signatures (the two largest line groups)
_W: dict = {}


class JobCPUTimeout(Exception):
    pass


def _on_prof(signum, frame):
    raise JobCPUTimeout()


def _init(arm_key: str):
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import w3
    from nextgem.model import load_model

    m = load_model(fc.PARENT)
    arm = fc.ARM[arm_key]
    fc.apply_arm(m, arm, fc.load_curated_rules() if arm.gpr else {})
    _W.update(arm=arm_key, model=m, genes=sorted(g.id for g in m.genes), closures=w3.gene_closures(m))
    signal.signal(signal.SIGPROF, _on_prof)


def _job(job: dict) -> dict:
    from nextgem import w3
    from nextgem.knockouts import WT_GROWTH_TOL

    assert job["arm"] == _W["arm"]
    v = w3.Variant.from_json(job["variant_json"])
    genes = _W["genes"] if job["mode"] == "full" else list(job["genes"])
    t0, c0 = time.time(), time.process_time()
    status, r = "ok", None
    signal.setitimer(signal.ITIMER_PROF, float(job["cpu_limit"]))
    try:
        r = w3.ko_sweep(_W["model"], genes, _W["closures"], job["bounds"], v)
    except JobCPUTimeout:
        status = "failed-timeout"
    except w3.SolverInfraError as e:
        status = "failed-timeout" if "time" in str(e).lower() else "failed-nonoptimal"
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0.0)
    out = {k: job[k] for k in ("arm", "mode", "purpose", "w2_sweep_id")}
    if r is not None:
        if r["wt_status"] != "optimal":
            status = "failed-infeasible"
        elif not r["wt_growth"] > WT_GROWTH_TOL:
            status = "failed-wt"
        out.update({k: r[k] for k in ("wt_growth", "wt_status", "n_solved", "n_genes_lp", "n_infeasible",
                                      "shortcut_verify_max_dev")})
        ratio = r["ratio_r"] if status == "ok" else np.full(len(genes), np.nan)
    else:
        ratio = np.full(len(genes), np.nan)
    out.update(status=status, ratio_r=np.asarray(ratio, dtype=float), genes=genes, n_genes=len(genes),
               cpu_s=time.process_time() - c0, wall_s=time.time() - t0)
    return out


def p1_rows(sigs) -> pd.DataFrame:
    """P1 ratio_r (W2 sweep id x model genes) from the pinned M6/W3 sweeps."""
    from m6_scores import sweep_ratios

    R = sweep_ratios()[0][fc.PARENT]
    P = R.loc[[mc.sweep_key(mc.BLOO, s) for s in sigs]]
    P.index = list(sigs)
    return P


def arm_ratios(arm_key: str, P1: pd.DataFrame, delta: pd.DataFrame, meta: pd.DataFrame | None = None) -> pd.DataFrame:
    """Reconstruct an arm's full ratio_r frame (W2 sweep id x genes) from P1 and the stored differences; rows of
    failed arm jobs (``meta`` status != ok) become NaN."""
    R = P1.copy()
    if arm_key == "P1":
        return R
    d = delta[delta.arm == arm_key]
    for sid, g in d.groupby("w2_sweep_id"):
        R.loc[sid, g.gene.to_numpy()] = g.ratio_r.to_numpy()
    if meta is not None:
        bad = meta[(meta.arm == arm_key) & (meta.purpose == "arm") & (meta.status != "ok")].w2_sweep_id
        R.loc[list(bad)] = np.nan
    return R


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    if a.workers > fc.MAX_WORKERS:
        raise SystemExit("at most 8 worker processes")
    inputs = [W2 / "controls" / "line_sweeps.tsv", fc.curated_rules_path(), fc.F1 / "closure_changes_human1.tsv",
              fc.F1 / "edits_human1.tsv", mc.SWDIR / "ratio_r_human1.parquet", mc.SWDIR / "catalog.tsv",
              mc.PROJECT / "results" / "w3" / "sweeps" / "ratio_r_human1.parquet",
              mc.PROJECT / "results" / "w3" / "structure" / "variants_human1.json"]
    with fc.provenance(inputs=inputs, solver="glpk", workers=a.workers, allow_dirty=a.allow_dirty) as prov:
        outs, summ = run(a)
        prov.output(*outs)
        prov.note(**{k: v for k, v in summ.items() if k not in ("a10_limits", "status")})
    return 0 if summ["spot_all_bitwise_equal"] else 3


def run(a):
    from m6_sweeps import bounds_by_id, needed, variant_json

    t_start = time.time()
    ls = pd.read_csv(W2 / "controls" / "line_sweeps.tsv", sep="\t")
    sigs = needed(ls)[fc.PARENT][mc.M2]
    bounds, info = bounds_by_id({fc.PARENT: {mc.M2: sigs}})
    vj = variant_json(fc.PARENT, mc.BLOO)
    changed = list(pd.read_csv(fc.F1 / "closure_changes_human1.tsv", sep="\t").gene)
    first_key = {s: info[s]["medium_keys"].split(";")[0] for s in sigs}
    spot = []
    for k in SPOT_MEDIA:
        hit = [s for s in sigs if first_key[s] == k or k in info[s]["medium_keys"].split(";")]
        spot.append(hit[0])
    jobs = []
    for s in sigs:
        jobs.append(dict(arm="TYMS", mode="full", purpose="arm", w2_sweep_id=s))
        jobs.append(dict(arm="GPR", mode="partial", purpose="arm", w2_sweep_id=s, genes=changed))
        jobs.append(dict(arm="F1", mode="partial", purpose="arm", w2_sweep_id=s, genes=changed))
    for s in spot:
        jobs.append(dict(arm="P1", mode="full", purpose="spot_check", w2_sweep_id=s))
        jobs.append(dict(arm="F1", mode="full", purpose="spot_check", w2_sweep_id=s))
    jobs.append(dict(arm="GPR", mode="full", purpose="spot_check", w2_sweep_id=spot[0]))
    for j in jobs:
        j.update(bounds=bounds[j["w2_sweep_id"]], variant_json=vj, cpu_limit=mc.JOB_CPU_LIMIT_S)
    print(f"{len(sigs)} M2(f=1) signatures, {len(changed)} genes with changed closures, {len(jobs)} jobs", flush=True)
    res = []
    ctx = multiprocessing.get_context("spawn")
    for arm in ("TYMS", "GPR", "F1", "P1"):
        aj = [j for j in jobs if j["arm"] == arm]
        with cf.ProcessPoolExecutor(max_workers=min(a.workers, len(aj)), mp_context=ctx, initializer=_init,
                                    initargs=(arm,)) as pool:
            for n, r in enumerate(pool.map(_job, aj, chunksize=1), 1):
                res.append(r)
                print(f"  [{arm}] {n}/{len(aj)} {r['mode']} {r['purpose']} {r['w2_sweep_id']}: {r['status']}, "
                      f"WT {r.get('wt_growth', np.nan):.6g}, {r.get('n_solved', 0)} LPs, {r['cpu_s']:.0f} s CPU, "
                      f"{time.time() - t_start:.0f} s", flush=True)
    P1 = p1_rows(sigs)
    genes = list(P1.columns)
    # compose arm frames
    full = {}
    for r in res:
        if r["purpose"] == "arm" and r["mode"] == "full":
            full[(r["arm"], r["w2_sweep_id"])] = pd.Series(r["ratio_r"], index=r["genes"]).reindex(genes)
    comp = {"TYMS": P1.copy() * np.nan}
    for s in sigs:
        comp["TYMS"].loc[s] = full[("TYMS", s)]
    comp["GPR"], comp["F1"] = P1.copy(), comp["TYMS"].copy()
    for r in res:
        if r["purpose"] == "arm" and r["mode"] == "partial":
            if r["status"] != "ok":
                comp[r["arm"]].loc[r["w2_sweep_id"]] = np.nan
            else:
                comp[r["arm"]].loc[r["w2_sweep_id"], r["genes"]] = r["ratio_r"]
    rows = []
    for arm in ("TYMS", "GPR", "F1"):
        X = comp[arm]
        diff = ~np.isclose(X.to_numpy(), P1.to_numpy(), rtol=0, atol=0, equal_nan=True)
        ii, jj = np.nonzero(diff)
        for i, j in zip(ii, jj):
            rows.append({"arm": arm, "w2_sweep_id": X.index[i], "gene": genes[j], "ratio_r": X.iat[i, j],
                         "ratio_r_P1": P1.iat[i, j]})
    delta = pd.DataFrame(rows, columns=["arm", "w2_sweep_id", "gene", "ratio_r", "ratio_r_P1"])
    meta = pd.DataFrame([{k: v for k, v in r.items() if k not in ("ratio_r", "genes")} for r in res])
    meta["medium_keys"] = meta.w2_sweep_id.map(lambda s: info[s]["medium_keys"])
    meta["serum"] = meta.w2_sweep_id.map(lambda s: info[s]["serum"])
    # reconstruction round trip
    for arm in ("TYMS", "GPR", "F1"):
        Y = arm_ratios(arm, P1, delta, meta)
        if not np.array_equal(Y.to_numpy(), comp[arm].to_numpy(), equal_nan=True):
            raise RuntimeError(f"{arm}: reconstruction from the stored differences is not exact")
    # spot checks (bitwise)
    spots = []
    for r in res:
        if r["purpose"] != "spot_check":
            continue
        s = r["w2_sweep_id"]
        ref = P1.loc[s] if r["arm"] == "P1" else comp[r["arm"]].loc[s]
        new = pd.Series(r["ratio_r"], index=r["genes"]).reindex(genes)
        spots.append({"arm": r["arm"], "w2_sweep_id": s, "medium_keys": info[s]["medium_keys"],
                      "reference": "M6/W3 pinned P1 sweep" if r["arm"] == "P1" else f"composed {r['arm']} (partial re-solve)",
                      "status": r["status"], "bitwise_equal": bool(np.array_equal(ref.to_numpy(), new.to_numpy(), equal_nan=True)),
                      "max_abs_diff": float(np.nanmax(np.abs(ref.to_numpy() - new.to_numpy()))),
                      "n_lethal": int((new < 0.01).sum())})
    spots = pd.DataFrame(spots)
    D = fc.SWDIR
    outs = [fc.atomic_write_parquet(delta, D / "ratio_r_delta_human1.parquet"), fc.atomic_write_tsv(meta, D / "sweeps.tsv"),
            fc.atomic_write_tsv(spots, D / "spot_checks.tsv")]
    nd = delta.groupby("arm").gene.nunique().to_dict()
    summ = {"n_signatures": len(sigs), "n_jobs": len(jobs), "status": meta.status.value_counts().to_dict(),
            "n_failed": int((meta.status != "ok").sum()),
            "genes_with_any_change_vs_P1": {k: int(v) for k, v in nd.items()},
            "job_cpu_s": round(float(meta.cpu_s.sum()), 1), "wall_s": round(time.time() - t_start, 1),
            "workers": a.workers, "solver": "glpk", "spot_all_bitwise_equal": bool(spots.bitwise_equal.all()),
            "max_shortcut_dev": float(np.nanmax(meta.shortcut_verify_max_dev)),
            "a10_limits": {"lp_s": 300, "lp_retries": 1, "job_cpu_s": mc.JOB_CPU_LIMIT_S}}
    outs.append(fc.atomic_write_text(json.dumps(summ, indent=1) + "\n", D / "sweeps.json"))
    print(json.dumps(summ, indent=1), flush=True)
    print(spots.to_string(index=False), flush=True)
    return outs, summ


if __name__ == "__main__":
    raise SystemExit(main())
