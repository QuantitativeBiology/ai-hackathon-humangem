"""M6 label-free KO sweeps: the W3 B-loo:cofactor_pool_biomass parent on every distinct M2(f=1) medium, both parents
(DESIGN §10.1, §17 A10-A12; RESEARCH_TREE D10). No label is read.

Reuse (not re-run): plain parents on M0 and on every M2(f=1) signature = the cached W2 sweeps
(results/w2/ko/ratio_r_<parent>.parquet); B-loo on M0 and on the W3 RPMI proxy (= the W2 classic-RPMI M2(f=1)
signature) = the W3 sweeps (results/w3/sweeps/ratio_r_<parent>.parquet). New: B-loo on each remaining M2(f=1)
signature. Spot checks: a few reused sweeps are recomputed here and must be bitwise equal on ratio_r.

Bounds: rebuilt with ``nextgem.w2_engine.build_sweeps`` from the frozen media library and matched to the W2 sweep
ids by bounds signature (every needed id must be found). Variant: the W3 spec (results/w3/structure/
variants_<parent>.json), applied with ``nextgem.w3.Variant`` inside ``nextgem.w3.ko_sweep`` (the exact W2/W3 KO
engine with the zero-flux and identical-closure shortcuts).

Solver and limits (A10): GLPK; 300 s per LP + one retry (``w3.solve_status``); 1,800 s CPU per sweep job
(``signal.ITIMER_PROF``; checked between LPs, so an LP in progress finishes first, <= 2 x 300 s). A timeout, a
non-optimal LP after the retry, an infeasible wild type or WT growth <= tolerance is a failed build: recorded with
its status, never re-run with relaxed limits.

Outputs (results/m6/sweeps/): ratio_r_<parent>.parquet (new B-loo sweeps; rows ``sweep`` = "<variant>|<W2 sweep
id>", columns = model genes), sweeps.tsv (one row per job: purpose, status, WT growth, LPs, CPU), spot_checks.tsv,
catalog.tsv (every sweep any M6 cell uses: source, serum flag, media, lines), sweeps.json.

Usage:
  .venv/bin/python scripts/m6_sweeps.py --workers 8 [--pilot 2]   # --pilot: first N new jobs, outputs to tmp/
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
import m6_common as mc  # noqa: E402

W2 = mc.PROJECT / "results" / "w2"
W3 = mc.PROJECT / "results" / "w3"
SPOT = (  # (parent, variant, step, medium_key or None) -> reference source
    ("human1", mc.T0, mc.M2, "DMEM"),
    ("human2", mc.T0, mc.M2, "DMEM"),
    ("human1", mc.BLOO, mc.M2, "RPMI"),
    ("human2", mc.BLOO, mc.M0, None),
)
_W: dict = {}


class JobCPUTimeout(Exception):
    pass


def _on_prof(signum, frame):
    raise JobCPUTimeout()


def _init(parent: str):
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import w3
    from nextgem.model import load_model

    m = load_model(parent)
    _W.update(parent=parent, model=m, genes=sorted(g.id for g in m.genes), closures=w3.gene_closures(m))
    signal.signal(signal.SIGPROF, _on_prof)


def _job(job: dict) -> dict:
    from nextgem import w3
    from nextgem.knockouts import WT_GROWTH_TOL

    assert job["parent"] == _W["parent"]
    v = w3.Variant.from_json(job["variant_json"])
    t0, c0 = time.time(), time.process_time()
    status, r = "ok", None
    signal.setitimer(signal.ITIMER_PROF, float(job["cpu_limit"]))
    try:
        r = w3.ko_sweep(_W["model"], _W["genes"], _W["closures"], job["bounds"], v)
    except JobCPUTimeout:
        status = "failed-timeout"
    except w3.SolverInfraError as e:  # non-optimal LP status after the one retry (incl. the 300 s LP limit)
        status = "failed-timeout" if "time" in str(e).lower() else "failed-nonoptimal"
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0.0)
    out = {k: job[k] for k in ("parent", "key", "variant", "w2_sweep_id", "purpose")}
    if r is not None:
        if r["wt_status"] != "optimal":
            status = "failed-infeasible"
        elif not r["wt_growth"] > WT_GROWTH_TOL:
            status = "failed-wt"
        out.update({k: r[k] for k in ("wt_growth", "wt_status", "n_solved", "n_genes_lp", "n_infeasible",
                                      "shortcut_verify_max_dev")})
        ratio = r["ratio_r"] if status == "ok" else np.full(len(_W["genes"]), np.nan)
    else:
        ratio = np.full(len(_W["genes"]), np.nan)
    out.update(status=status, ratio_r=ratio, genes=_W["genes"], cpu_s=time.process_time() - c0,
               wall_s=time.time() - t0)
    return out


def needed(line_sweeps: pd.DataFrame) -> dict:
    """parent -> {step: sorted W2 sweep ids used by any line}."""
    return {p: {s: sorted(set(mc.line_sweep_ids(line_sweeps, p, s))) for s in mc.STEPS} for p in mc.PARENTS}


def bounds_by_id(need: dict) -> tuple[dict, dict]:
    """W2 sweep id -> bounds, and -> {serum, medium_keys}, rebuilt from the frozen library (GLPK, no Gurobi)."""
    from nextgem import w2_engine

    lines, sheets, doc, phys, sens, pcs, exchanges, hams = w2_engine.sweep_context()
    sweeps = w2_engine.build_sweeps(exchanges, hams, pcs, doc, {}, {})
    by = {s.sweep_id: s for s in sweeps}
    out, info = {}, {}
    for p, steps in need.items():
        for step, ids in steps.items():
            for sid in ids:
                if sid not in by:
                    raise RuntimeError(f"{sid} ({step}) is not rebuilt by w2_engine.build_sweeps")
                s = by[sid]
                labs = [l for l in s.labels if (step == mc.M0 and l["variant"] == "M0")
                        or (l["variant"] == "M2" and l["f"] == 1.0 and l["gate_open"])]
                if not labs:
                    raise RuntimeError(f"{sid}: no {step} label")
                keys = sorted({l["medium_key"] for l in labs})
                serum = step == mc.M2 and any(doc[k].serum for k in keys if k in doc)
                out[sid], info[sid] = dict(s.bounds), {"serum": bool(serum), "medium_keys": ";".join(keys)}
    return out, info


def variant_json(parent: str, variant: str) -> dict:
    from nextgem import w3

    if variant == mc.T0:
        return w3.Variant(mc.T0).to_json()
    vj = json.loads((W3 / "structure" / f"variants_{parent}.json").read_text())
    hit = [v for v in vj["bloo"] if v["name"] == variant]
    if len(hit) != 1:
        raise RuntimeError(f"{parent}: variant {variant} not found once in the W3 spec")
    return hit[0]


def reused_rows() -> dict:
    """parent -> ratio_r frame (sweep key x genes) of every reusable cached sweep (W2 plain, W3 B-loo M0 / RPMI)."""
    ss = json.loads((W3 / "structure" / "structure_summary.json").read_text())["w2_sweep_ids"]
    out = {}
    for p in mc.PARENTS:
        w2 = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
        w2.index = [mc.sweep_key(mc.T0, s) for s in w2.index]
        w3s = pd.read_parquet(W3 / "sweeps" / f"ratio_r_{p}.parquet").set_index("sweep")
        b = w3s.loc[[f"M0|{mc.BLOO}", f"RPMI_M2f1|{mc.BLOO}"]]
        b.index = [mc.sweep_key(mc.BLOO, ss[p]["M0"]), mc.sweep_key(mc.BLOO, ss[p]["RPMI_M2f1"])]
        out[p] = pd.concat([w2, b.reindex(columns=w2.columns)])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--pilot", type=int, default=0, help="run only the first N new jobs; outputs under results/m6/tmp/")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    if a.workers > mc.MAX_WORKERS:
        raise SystemExit("at most 8 worker processes")
    inputs = [W2 / "controls" / "line_sweeps.tsv", W2 / "ko" / "sweep_labels.tsv",
              W3 / "structure" / "structure_summary.json"]
    inputs += [W2 / "ko" / f"ratio_r_{p}.parquet" for p in mc.PARENTS]
    inputs += [W3 / "sweeps" / f"ratio_r_{p}.parquet" for p in mc.PARENTS]
    inputs += [W3 / "structure" / f"variants_{p}.json" for p in mc.PARENTS]
    with mc.provenance(inputs=inputs, solver="glpk", workers=a.workers, allow_dirty=a.allow_dirty or bool(a.pilot)) as prov:
        rc = run(a, prov)
    return rc


def run(a, prov) -> int:
    t_start = time.time()
    ls = pd.read_csv(W2 / "controls" / "line_sweeps.tsv", sep="\t")
    need = needed(ls)
    bounds, info = bounds_by_id(need)
    reused = reused_rows()
    jobs = []
    for p in mc.PARENTS:
        have = set(reused[p].index)
        for step in mc.STEPS:
            for sid in need[p][step]:
                k = mc.sweep_key(mc.BLOO, sid)
                if k not in have:
                    jobs.append(dict(parent=p, key=k, variant=mc.BLOO, w2_sweep_id=sid, purpose="new"))
                if mc.sweep_key(mc.T0, sid) not in have:
                    raise RuntimeError(f"plain-parent sweep {sid} missing from the W2 cache")
    n_new = len(jobs)
    for p, variant, step, key in SPOT:
        ids = need[p][step]
        sid = ids[0] if step == mc.M0 else [s for s in ids if info[s]["medium_keys"].split(";")[0] == key or
                                             key in info[s]["medium_keys"].split(";")][0]
        jobs.append(dict(parent=p, key=mc.sweep_key(variant, sid), variant=variant, w2_sweep_id=sid, purpose="spot_check"))
    if a.pilot:
        jobs = jobs[: a.pilot]
    for j in jobs:
        j.update(bounds=bounds[j["w2_sweep_id"]], variant_json=variant_json(j["parent"], j["variant"]),
                 cpu_limit=mc.JOB_CPU_LIMIT_S)
    print(f"{n_new} new B-loo sweeps, {len(jobs) - min(n_new, len(jobs))} spot checks; running {len(jobs)} jobs",
          flush=True)
    res = []
    ctx = multiprocessing.get_context("spawn")
    for p in mc.PARENTS:
        pj = [j for j in jobs if j["parent"] == p]
        if not pj:
            continue
        with cf.ProcessPoolExecutor(max_workers=min(a.workers, len(pj)), mp_context=ctx, initializer=_init,
                                    initargs=(p,)) as pool:
            for n, r in enumerate(pool.map(_job, pj, chunksize=1), 1):
                res.append(r)
                print(f"  [{p}] {n}/{len(pj)} {r['purpose']} {r['key']}: {r['status']}, WT {r.get('wt_growth', np.nan):.4g}, "
                      f"{r.get('n_solved', 0)} LPs, {r['cpu_s']:.0f} s CPU, {time.time() - t_start:.0f} s", flush=True)
    out_dir = mc.TMP / "pilot" if a.pilot else mc.SWDIR
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = pd.DataFrame([{k: v for k, v in r.items() if k not in ("ratio_r", "genes")} for r in res])
    meta["n_lethal"] = [int(np.nansum(r["ratio_r"] < 0.01)) for r in res]
    meta["serum"] = meta.w2_sweep_id.map(lambda s: info[s]["serum"])
    meta["medium_keys"] = meta.w2_sweep_id.map(lambda s: info[s]["medium_keys"])
    outs = [mc.atomic_write_tsv(meta, out_dir / "sweeps.tsv")]
    # spot checks: bitwise equality with the reused sweep
    spot = []
    for r in res:
        if r["purpose"] != "spot_check":
            continue
        ref = reused[r["parent"]].loc[r["key"], r["genes"]].to_numpy(dtype=float)
        new = np.asarray(r["ratio_r"], dtype=float)
        src = "W2 ko (cached)" if r["variant"] == mc.T0 else "W3 sweeps (cached)"
        spot.append({"parent": r["parent"], "key": r["key"], "reference": src, "status": r["status"],
                     "n_genes": len(new), "bitwise_equal": bool(np.array_equal(ref, new, equal_nan=True)),
                     "max_abs_diff": float(np.nanmax(np.abs(ref - new))) if r["status"] == "ok" else np.nan,
                     "n_lethal": int((new < 0.01).sum())})
    if spot:
        outs.append(mc.atomic_write_tsv(pd.DataFrame(spot), out_dir / "spot_checks.tsv"))
    for p in mc.PARENTS:
        rr = [r for r in res if r["parent"] == p and r["purpose"] == "new"]
        if not rr:
            continue
        R = pd.DataFrame(np.vstack([r["ratio_r"] for r in rr]), columns=rr[0]["genes"],
                         index=pd.Index([r["key"] for r in rr], name="sweep"))
        outs.append(mc.atomic_write_parquet(R.reset_index(), out_dir / f"ratio_r_{p}.parquet"))
    # catalog of every sweep an M6 cell uses (label-free; line counts over all screened lines)
    if not a.pilot:
        w2m = pd.read_csv(W2 / "ko" / "sweeps.tsv", sep="\t").set_index("sweep_id")
        w3m = {p: pd.read_csv(W3 / "sweeps" / f"sweeps_{p}.tsv", sep="\t") for p in mc.PARENTS}
        ss = json.loads((W3 / "structure" / "structure_summary.json").read_text())["w2_sweep_ids"]
        newm = meta[meta.purpose == "new"].set_index("key")
        cat = []
        for p in mc.PARENTS:
            for step in mc.STEPS:
                n_lines = mc.line_sweep_ids(ls, p, step).value_counts()
                for sid in need[p][step]:
                    for v in mc.VARIANTS:
                        k = mc.sweep_key(v, sid)
                        if v == mc.T0:
                            src, wt, st = "W2", float(w2m.loc[sid, "wt_growth"]), "ok"
                        elif k in newm.index:
                            src, wt, st = "M6", float(newm.loc[k, "wt_growth"]), newm.loc[k, "status"]
                        else:
                            tag = "M0" if sid == ss[p]["M0"] else "RPMI_M2f1"
                            row = w3m[p][(w3m[p].medium == tag) & (w3m[p].variant == mc.BLOO)].iloc[0]
                            src, wt, st = "W3", float(row.wt_growth), "ok"
                        cat.append({"parent": p, "variant": v, "step": step, "w2_sweep_id": sid, "key": k,
                                    "source": src, "status": st, "wt_growth": wt, "serum": info[sid]["serum"],
                                    "medium_keys": info[sid]["medium_keys"], "n_lines": int(n_lines.get(sid, 0))})
        outs.append(mc.atomic_write_tsv(pd.DataFrame(cat), out_dir / "catalog.tsv"))
    summ = {"n_new": int((meta.purpose == "new").sum()), "n_spot": int((meta.purpose == "spot_check").sum()),
            "status": meta.status.value_counts().to_dict(), "job_cpu_s": round(float(meta.cpu_s.sum()), 1),
            "wall_s": round(time.time() - t_start, 1), "workers": a.workers, "solver": "glpk",
            "spot_all_bitwise_equal": bool(all(s["bitwise_equal"] for s in spot)) if spot else None,
            "max_shortcut_dev": float(np.nanmax(meta.shortcut_verify_max_dev)) if "shortcut_verify_max_dev" in meta else None,
            "a10_limits": {"lp_s": 300, "lp_retries": 1, "job_cpu_s": mc.JOB_CPU_LIMIT_S}}
    outs.append(mc.atomic_write_text(json.dumps(summ, indent=1) + "\n", out_dir / "sweeps.json"))
    prov.output(*outs)
    prov.note(**{k: v for k, v in summ.items() if k != "a10_limits"}, pilot=a.pilot)
    print(json.dumps(summ, indent=1), flush=True)
    if spot and not summ["spot_all_bitwise_equal"]:
        print("SPOT CHECK FAILED: a reused sweep is not reproduced bitwise", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
