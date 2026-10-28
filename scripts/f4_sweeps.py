"""F4a step 2 (label-free): KO sweeps of the P1 parent with each line's accepted genomic losses (A15 post hoc; D19).

Input: results/f4/losses/line_losses.tsv (scripts/f4_losses.py). One sweep job per distinct (M2(f=1) signature,
accepted-loss set); lines without an accepted loss equal P1 (pinned M6/W3 sweeps, not recomputed).

Engine: the P1 engine itself, ``nextgem.w3.ko_sweep`` (W2/W3/M6; zero-flux and identical-closure shortcuts, B-loo
applied with ``nextgem.w3.Variant``, the signature's bounds), GLPK, A10 limits (300 s per LP + one retry inside
``w3.solve_status``; 1,800 s CPU per job via ITIMER_PROF; timeout / non-optimal / failed WT = failed build, never
re-run with relaxed limits). Genomic constraint (A4 convention): inside the job the accepted genes are knocked out
(cobra ``Gene.knock_out``: a gene is non-functional for every later GPR evaluation, so a KO of an isozyme partner
closes the shared reaction and the KO of a lost gene closes nothing); the reactions this closes are asserted equal to
``context.loss_closure`` from f4_losses (modulo reactions already closed in the model); the KO closures are recomputed
in that lost state (``w3.gene_closures``).
- 'full' jobs (>= 1 reaction closed): every gene is swept.
- 'partial' jobs (no reaction closed; the WT LP is P1's): only genes whose KO closure differs from the generic one are
  re-solved; every other gene keeps P1's ratio. Exactness: spot check (a partial job also swept in full, bitwise).
Spot checks: an empty-loss job on the 2 largest signatures == the pinned P1 sweep (bitwise); cross-engine check on
2 full jobs against the W4 engine (``context.finish`` with the accepted losses + ``context.ko_sweep`` on the
``f4_common.p1_parent``; reported as max |Δratio_r|).

Outputs (results/f4/sweeps/): ratio_r_delta_human1.parquet (job_key, gene, ratio_r, ratio_r_P1) for every (job, gene)
whose ratio differs from P1 (reconstruction: :func:`job_ratios`), jobs.tsv, spot_checks.tsv, sweeps.json.

Usage: .venv/bin/python scripts/f4_sweeps.py --workers 12 [--pilot N] [--allow-dirty]
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
import f4_common as fc  # noqa: E402

fc.blas_single_thread()
_W: dict = {}


class JobCPUTimeout(Exception):
    pass


def _on_prof(signum, frame):
    raise JobCPUTimeout()


def _init():
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import w3
    from nextgem.model import load_model

    m = load_model(fc.PARENT)
    _W.update(model=m, genes=sorted(g.id for g in m.genes), closures=w3.gene_closures(m),
              zero0={r.id for r in m.reactions if r.bounds == (0, 0)})
    signal.signal(signal.SIGPROF, _on_prof)


def _job(job: dict) -> dict:
    from nextgem import w3
    from nextgem.knockouts import WT_GROWTH_TOL

    m = _W["model"]
    v = w3.Variant.from_json(job["variant_json"])
    t0, c0 = time.time(), time.process_time()
    status, r, genes, closure_ok = "ok", None, [], None
    signal.setitimer(signal.ITIMER_PROF, float(job["cpu_limit"]))
    try:
        with m:
            for g in job["accepted"]:
                m.genes.get_by_id(g).knock_out()
            closed_now = {x.id for x in m.reactions if x.bounds == (0, 0)} - _W["zero0"]
            exp = set(job["expected_closed"])
            closure_ok = bool(closed_now <= exp and (exp - closed_now) <= _W["zero0"])
            if not closure_ok:
                raise RuntimeError(f"{job['job_key']}: knock-out closure != context.loss_closure")
            clos = w3.gene_closures(m) if job["accepted"] else _W["closures"]
            if job["mode"] == "full":
                genes = _W["genes"]
            else:
                genes = [g for g in _W["genes"] if set(clos[g]) != set(_W["closures"][g])]
            if genes:
                r = w3.ko_sweep(m, genes, clos, job["bounds"], v)
    except JobCPUTimeout:
        status = "failed-timeout"
    except w3.SolverInfraError as e:
        status = "failed-timeout" if "time" in str(e).lower() else "failed-nonoptimal"
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0.0)
    out = {k: job[k] for k in ("job_key", "w2_sweep_id", "mode", "purpose")}
    out["n_accepted"] = len(job["accepted"])
    out["n_closed"] = len(job["expected_closed"])
    if r is not None:
        if r["wt_status"] != "optimal":
            status = "failed-infeasible"
        elif not r["wt_growth"] > WT_GROWTH_TOL:
            status = "failed-wt"
        out.update({k: r[k] for k in ("wt_growth", "wt_status", "n_solved", "n_genes_lp", "n_infeasible",
                                      "shortcut_verify_max_dev")})
        ratio = r["ratio_r"] if status == "ok" else np.full(len(genes), np.nan)
    else:
        ratio = np.full(len(genes), np.nan) if status != "ok" else np.zeros(0)
    out.update(status=status, ratio_r=np.asarray(ratio, dtype=float), genes=list(genes), n_genes=len(genes),
               closure_ok=closure_ok, cpu_s=time.process_time() - c0, wall_s=time.time() - t0)
    return out


def _xengine(job: dict) -> dict:
    """W4-engine cross-check (fresh process): context.finish with the accepted losses + context.ko_sweep."""
    from nextgem import context as cx

    c0 = time.process_time()
    parent = fc.p1_parent(job["bounds"], job["w2_sweep_id"])
    lm = cx.LineModel(line=job["job_key"], family="F0", closed=set())
    lm.info["_s_r"] = pd.Series(dtype="float64")
    cx.finish(parent, lm, set(job["accepted"]))
    ko = cx.ko_sweep(parent, lm)
    return {"job_key": job["job_key"], "ratio": ko.set_index("gene").ratio_r, "n_reopened": len(lm.reopened),
            "n_lost_closed": len(lm.lost_closed), "wt": lm.info["wt_growth"], "cpu_s": time.process_time() - c0}


def p1_rows(sigs) -> pd.DataFrame:
    """P1 ratio_r (W2 sweep id x Human1 genes) from the pinned M6/W3 sweeps (as F1)."""
    from m6_scores import sweep_ratios

    R = sweep_ratios()[0][fc.PARENT]
    P = R.loc[[f"{fc.BLOO}|{s}" for s in sigs]]
    P.index = list(sigs)
    return P


def job_ratios(P1: pd.DataFrame, delta: pd.DataFrame, jobs: pd.DataFrame) -> pd.DataFrame:
    """job_key x genes ratio_r: the job's signature row of P1 with the stored differences; failed jobs are NaN."""
    J = jobs[jobs.purpose == "arm"].set_index("job_key")
    R = P1.loc[J.w2_sweep_id.to_numpy()].copy()
    R.index = J.index
    for k, g in delta.groupby("job_key"):
        R.loc[k, g.gene.to_numpy()] = g.ratio_r.to_numpy()
    bad = J.index[J.status != "ok"]
    R.loc[bad] = np.nan
    return R


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=fc.MAX_WORKERS)
    ap.add_argument("--pilot", type=int, default=0, help="first N arm jobs only; outputs under results/f4/tmp/")
    ap.add_argument("--losses", help="line_losses.tsv (default results/f4/losses/line_losses.tsv)")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    if not 1 <= a.workers <= fc.MAX_WORKERS:
        raise SystemExit(f"--workers must be in 1..{fc.MAX_WORKERS}")
    losses = Path(a.losses) if a.losses else fc.LOSSDIR / "line_losses.tsv"
    if not a.pilot and fc.in_tmp(losses):
        raise SystemExit("production sweeps use the production loss table")
    out = fc.SWDIR if not a.pilot else fc.TMP / "sweeps_pilot"
    inputs = [losses, fc.PROJECT / "results" / "m6" / "sweeps" / "ratio_r_human1.parquet",
              fc.PROJECT / "results" / "m6" / "sweeps" / "catalog.tsv",
              fc.PROJECT / "results" / "w3" / "sweeps" / "ratio_r_human1.parquet",
              fc.PROJECT / "results" / "w3" / "structure" / "variants_human1.json",
              fc.PROJECT / "results" / "w2" / "controls" / "line_sweeps.tsv"]
    manifest = fc.MANIFEST if not a.pilot else out / "MANIFEST.tsv"
    with fc.provenance(inputs=inputs, solver="glpk", workers=a.workers, allow_dirty=a.allow_dirty or bool(a.pilot),
                       manifest=manifest) as prov:
        outs, summ = run(a, losses, out)
        prov.output(*outs)
        prov.note(step="f4_sweeps (F4a, label-free, D19)", **{k: v for k, v in summ.items() if k != "a10_limits"})
    return 0 if summ["spot_all_pass"] and summ["n_failed"] == 0 else 3


def run(a, losses: Path, out: Path):
    t_start = time.time()
    L = pd.read_csv(losses, sep="\t", dtype={"accepted": str, "closed_rxns": str, "job_key": str})
    L = L.fillna({"accepted": "", "closed_rxns": "", "job_key": ""})
    arm = L[L.n_accepted > 0].drop_duplicates("job_key").sort_values("job_key")
    sigs = sorted(set(L.w2_sweep_id))
    bounds, info = fc.signature_bounds(sigs)
    vj = fc.variant_json()
    jobs = []
    for r in arm.itertuples():
        acc = fc.split_genes(r.accepted)
        if fc.job_key(r.w2_sweep_id, acc) != r.job_key:
            raise RuntimeError(f"{r.ModelID}: job key mismatch")
        jobs.append(dict(job_key=r.job_key, w2_sweep_id=r.w2_sweep_id, accepted=acc,
                         expected_closed=fc.split_genes(r.closed_rxns), mode=r.sweep_mode, purpose="arm"))
    if a.pilot:
        jobs = jobs[: a.pilot]
    # spot checks
    big = L.w2_sweep_id.value_counts().index[:2].tolist()
    for s in big:
        jobs.append(dict(job_key=f"{s}|P1", w2_sweep_id=s, accepted=[], expected_closed=[], mode="full",
                         purpose="spot_p1"))
    part = [j for j in jobs if j["purpose"] == "arm" and j["mode"] == "partial"]
    if part:
        jobs.append({**part[0], "mode": "full", "purpose": "spot_partial"})
    full = [j for j in jobs if j["purpose"] == "arm" and j["mode"] == "full"]
    xjobs = [dict(j) for j in sorted(full, key=lambda j: -len(j["expected_closed"]))[:2]]
    for j in jobs + xjobs:
        j.update(bounds=bounds[j["w2_sweep_id"]], variant_json=vj, cpu_limit=fc.JOB_CPU_LIMIT_S)
    n_arm = sum(j["purpose"] == "arm" for j in jobs)
    print(f"{len(L)} lines, {int((L.n_accepted > 0).sum())} with an accepted loss; {n_arm} arm jobs "
          f"({sum(j['mode'] == 'full' for j in jobs if j['purpose'] == 'arm')} full); {len(jobs) - n_arm} spot jobs; "
          f"{len(xjobs)} cross-engine jobs", flush=True)
    ctx = multiprocessing.get_context("spawn")
    res = []
    with cf.ProcessPoolExecutor(max_workers=min(a.workers, len(jobs)), mp_context=ctx, initializer=_init) as pool:
        for n, r in enumerate(pool.map(_job, jobs, chunksize=1), 1):
            res.append(r)
            if n % 25 == 0 or n == len(jobs) or r["status"] != "ok":
                print(f"  {n}/{len(jobs)} {r['job_key']} {r['mode']} {r['purpose']}: {r['status']}, "
                      f"{r.get('n_solved', 0)} LPs, {r['cpu_s']:.0f} s CPU, {time.time() - t_start:.0f} s", flush=True)
    xres = []
    if xjobs:
        with cf.ProcessPoolExecutor(max_workers=min(a.workers, len(xjobs)), mp_context=ctx) as pool:
            xres = list(pool.map(_xengine, xjobs))
    P1 = p1_rows(sigs)
    genes = list(P1.columns)
    rows = []
    comp = {}
    for r in res:
        base = P1.loc[r["w2_sweep_id"]].copy()
        if r["status"] != "ok":
            new = base * np.nan
        else:
            new = base.copy()
            if r["genes"]:
                new.loc[r["genes"]] = r["ratio_r"]
        comp[(r["job_key"], r["purpose"])] = new
        if r["purpose"] == "arm" and r["status"] == "ok":
            d = new.to_numpy() != base.to_numpy()
            for g in np.asarray(genes)[d]:
                rows.append({"job_key": r["job_key"], "gene": g, "ratio_r": float(new[g]), "ratio_r_P1": float(base[g])})
    delta = pd.DataFrame(rows, columns=["job_key", "gene", "ratio_r", "ratio_r_P1"])
    meta = pd.DataFrame([{k: v for k, v in r.items() if k not in ("ratio_r", "genes")} for r in res])
    meta["medium_keys"] = meta.w2_sweep_id.map(lambda s: info[s]["medium_keys"])
    meta["serum"] = meta.w2_sweep_id.map(lambda s: info[s]["serum"])
    meta["n_lines"] = meta.job_key.map(L[L.job_key != ""].job_key.value_counts()).fillna(0).astype(int)
    R = job_ratios(P1, delta, meta)
    for r in res:  # reconstruction round trip
        if r["purpose"] == "arm":
            if not np.array_equal(R.loc[r["job_key"]].to_numpy(), comp[(r["job_key"], "arm")].to_numpy(), equal_nan=True):
                raise RuntimeError(f"{r['job_key']}: reconstruction from the stored differences is not exact")
    spots = []
    for r in res:
        if r["purpose"] == "spot_p1":
            new = comp[(r["job_key"], "spot_p1")]
            ref, what = P1.loc[r["w2_sweep_id"]], "pinned P1 sweep (M6/W3)"
        elif r["purpose"] == "spot_partial":
            new = comp[(r["job_key"], "spot_partial")]
            ref, what = comp[(r["job_key"], "arm")], "composed partial job"
        else:
            continue
        spots.append({"check": r["purpose"], "job_key": r["job_key"], "reference": what, "status": r["status"],
                      "pass": bool(r["status"] == "ok" and np.array_equal(ref.to_numpy(), new.to_numpy(), equal_nan=True)),
                      "max_abs_diff": float(np.nanmax(np.abs(ref.to_numpy() - new.to_numpy())))})
    for x in xres:
        mine = comp[(x["job_key"], "arm")]
        other = x["ratio"].reindex(genes).astype(float)
        dd = np.abs(mine.to_numpy() - other.to_numpy())
        spots.append({"check": "cross_engine_w4", "job_key": x["job_key"],
                      "reference": "W4 engine (context.finish + context.ko_sweep)", "status": "ok",
                      "pass": bool(np.nanmax(dd) <= 1e-6 and x["n_reopened"] == 0),
                      "max_abs_diff": float(np.nanmax(dd)), "n_genes_diff_gt_1e-6": int((dd > 1e-6).sum()),
                      "n_lethal_w3_engine": int((mine < 0.01).sum()), "n_lethal_w4_engine": int((other < 0.01).sum())})
    spots = pd.DataFrame(spots)
    am = meta[meta.purpose == "arm"]
    outs = [fc.atomic_write_parquet(delta, out / "ratio_r_delta_human1.parquet"),
            fc.atomic_write_tsv(meta, out / "jobs.tsv"), fc.atomic_write_tsv(spots, out / "spot_checks.tsv")]
    summ = {"post_hoc": "A15 post hoc / exploratory", "n_lines": int(len(L)),
            "n_lines_with_accepted_loss": int((L.n_accepted > 0).sum()), "n_arm_jobs": int(len(am)),
            "n_full": int((am["mode"] == "full").sum()), "n_partial": int((am["mode"] == "partial").sum()),
            "n_failed": int((meta.status != "ok").sum()), "status": meta.status.value_counts().to_dict(),
            "n_jobs_changing_any_ratio": int(delta.job_key.nunique()),
            "n_genes_changed_any_job": int(delta.gene.nunique()),
            "n_lethal_flips_to_lethal": int(((delta.ratio_r < 0.01) & (delta.ratio_r_P1 >= 0.01)).sum()),
            "n_lethal_flips_to_viable": int(((delta.ratio_r >= 0.01) & (delta.ratio_r_P1 < 0.01)).sum()),
            "job_cpu_s": round(float(meta.cpu_s.sum()) + sum(x["cpu_s"] for x in xres), 1),
            "wall_s": round(time.time() - t_start, 1), "workers": a.workers, "solver": "glpk",
            "spot_all_pass": bool(spots["pass"].all()) if len(spots) else True,
            "a10_limits": {"lp_s": 300, "lp_retries": 1, "job_cpu_s": fc.JOB_CPU_LIMIT_S}}
    outs.append(fc.atomic_write_json(summ, out / "sweeps.json"))
    print(json.dumps(summ, indent=1, default=float), flush=True)
    print(spots.to_string(index=False), flush=True)
    return outs, summ


if __name__ == "__main__":
    raise SystemExit(main())
