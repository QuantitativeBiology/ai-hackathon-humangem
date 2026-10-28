"""W7-03 label-free KO sweeps: the 'artefact groups only' B-loo variant on Human1 (DESIGN §W7, §17 A11e/A12a;
notes/reviews/M4_W3_review.md §6). Exploratory, post hoc. No label is read.

Variant (fixed in autoresearch/ideas.tsv W7-03 before any sweep or score): remove from the cofactor pool MAR10065 only
the four artefact groups of results/w3/structure/variants_human1.json: B-loo:retinoid_catabolites,
B-loo:vitaminD_catabolites, B-loo:tocopherol_catabolites and B-loo:carnitine (combo). Ubiquinol, heme, lipoate,
cobamide and the unassigned vitamins stay in the pool.

Media: the same W2 sweep signatures as the M6 cells: M0 (Ham's) and every M2(f=1) signature used by any line for
human1 (results/w2/controls/line_sweeps.tsv). Engine, solver and A10 limits: exactly scripts/m6_sweeps.py (GLPK, 300 s
per LP + one retry, 1,800 s CPU per sweep job). One spot check: T0 on M0 must equal the W2 cache bitwise.

Outputs (results/w7/sweeps/): ratio_r_human1.parquet (rows "W7:artefact_bloo|<W2 sweep id>"), sweeps.tsv,
mev.tsv (W3-MEV diagnostic: mevalonate-core KO ratio per sweep), sweeps.json.

Usage: .venv/bin/python scripts/w7_artefact_sweeps.py --workers 8
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
import m6_common as mc  # noqa: E402
import m6_sweeps as ms  # noqa: E402

PARENT = "human1"
GROUPS = ("B-loo:retinoid_catabolites", "B-loo:vitaminD_catabolites", "B-loo:tocopherol_catabolites", "B-loo:carnitine")
VARIANT = "W7:artefact_bloo"
OUT = mc.PROJECT / "results" / "w7" / "sweeps"
MANIFEST = mc.PROJECT / "results" / "w7" / "MANIFEST.tsv"


def artefact_variant_json(parent: str = PARENT) -> dict:
    vj = json.loads((ms.W3 / "structure" / f"variants_{parent}.json").read_text())
    by = {v["name"]: v for v in vj["bloo"]}
    parts = [by[g] for g in GROUPS]
    for p in parts:
        assert p["kind"] == "pool_remove" and p["params"].get("pool", mc_pool()) == mc_pool(), p["name"]
    return {"name": VARIANT, "kind": "combo", "params": {"variants": parts}}


def mc_pool() -> str:
    from nextgem import w3

    return w3.COFACTOR_POOL


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args(argv)
    if a.workers > mc.MAX_WORKERS:
        raise SystemExit("at most 8 worker processes")
    from nextgem.run import provenance

    inputs = [ms.W2 / "controls" / "line_sweeps.tsv", ms.W2 / "ko" / f"ratio_r_{PARENT}.parquet",
              ms.W3 / "structure" / f"variants_{PARENT}.json"]
    with provenance(MANIFEST, inputs=inputs, solver="glpk", workers=a.workers, allow_dirty=True,
                    author="W7 auto-research engineer (Claude Opus)") as prov:
        return run(a, prov)


def run(a, prov) -> int:
    t0 = time.time()
    ls = pd.read_csv(ms.W2 / "controls" / "line_sweeps.tsv", sep="\t")
    need = {PARENT: {s: sorted(set(mc.line_sweep_ids(ls, PARENT, s))) for s in mc.STEPS}}
    bounds, info = ms.bounds_by_id(need)
    vj = artefact_variant_json()
    jobs = []
    for step in mc.STEPS:
        for sid in need[PARENT][step]:
            jobs.append(dict(parent=PARENT, key=mc.sweep_key(VARIANT, sid), variant=VARIANT, w2_sweep_id=sid,
                             purpose="new", step=step, variant_json=vj))
    m0 = need[PARENT][mc.M0][0]
    jobs.append(dict(parent=PARENT, key=mc.sweep_key(mc.T0, m0), variant=mc.T0, w2_sweep_id=m0, purpose="spot_check",
                     step=mc.M0, variant_json=ms.variant_json(PARENT, mc.T0)))
    for j in jobs:
        j.update(bounds=bounds[j["w2_sweep_id"]], cpu_limit=mc.JOB_CPU_LIMIT_S)
    print(f"{len(jobs)} jobs ({len(jobs) - 1} new)", flush=True)
    ctx = multiprocessing.get_context("spawn")
    res = []
    with cf.ProcessPoolExecutor(max_workers=min(a.workers, len(jobs)), mp_context=ctx, initializer=ms._init,
                                initargs=(PARENT,)) as pool:
        for n, r in enumerate(pool.map(ms._job, jobs, chunksize=1), 1):
            r["step"] = jobs[n - 1]["step"]
            res.append(r)
            print(f"  {n}/{len(jobs)} {r['key']}: {r['status']} WT {r.get('wt_growth', np.nan):.4g} "
                  f"{r['cpu_s']:.0f} s CPU", flush=True)
    meta = pd.DataFrame([{k: v for k, v in r.items() if k not in ("ratio_r", "genes")} for r in res])
    meta["serum"] = meta.w2_sweep_id.map(lambda s: info[s]["serum"])
    meta["medium_keys"] = meta.w2_sweep_id.map(lambda s: info[s]["medium_keys"])
    meta["n_lethal"] = [int(np.nansum(np.asarray(r["ratio_r"]) < 0.01)) for r in res]
    outs = [mc.atomic_write_tsv(meta, OUT / "sweeps.tsv")]
    genes = res[0]["genes"]
    new = [r for r in res if r["purpose"] == "new"]
    R = pd.DataFrame(np.vstack([r["ratio_r"] for r in new]), columns=genes,
                     index=pd.Index([r["key"] for r in new], name="sweep"))
    outs.append(mc.atomic_write_parquet(R.reset_index(), OUT / f"ratio_r_{PARENT}.parquet"))
    spot = [r for r in res if r["purpose"] == "spot_check"][0]
    w2 = pd.read_parquet(ms.W2 / "ko" / f"ratio_r_{PARENT}.parquet")
    ref = w2.loc[m0, genes].to_numpy(dtype=float)
    spot_eq = bool(np.array_equal(ref, np.asarray(spot["ratio_r"], dtype=float), equal_nan=True))
    # W3-MEV diagnostic (label-free): mevalonate core KO ratio per sweep, both this variant and P1's B-loo
    from nextgem import paths

    bk = pd.read_parquet(paths.baseline_dir() / f"{PARENT}_ko.parquet", columns=["gene", "symbol"])
    sym = {g: s for g, s in zip(bk.gene, bk.symbol) if s in mc.MEVALONATE_CORE}
    core = [g for g in genes if g in sym]
    if sorted(sym[g] for g in core) != sorted(mc.MEVALONATE_CORE):
        raise RuntimeError(f"mevalonate core genes not all found: {sorted(sym.values())}")
    rows = []
    for r in new:
        rr = pd.Series(r["ratio_r"], index=genes)
        for g in core:
            rows.append({"key": r["key"], "step": r["step"], "serum": info[r["w2_sweep_id"]]["serum"],
                         "gene": sym[g], "ratio_r": float(rr[g]), "lethal": bool(rr[g] < 0.01)})
    mev = pd.DataFrame(rows)
    outs.append(mc.atomic_write_tsv(mev, OUT / "mev.tsv"))
    ser = mev[mev.serum]
    summ = {"variant": vj["name"], "groups": list(GROUPS), "n_new": len(new), "status": meta.status.value_counts().to_dict(),
            "spot_T0_M0_bitwise_equal_W2": spot_eq, "job_cpu_s": round(float(meta.cpu_s.sum()), 1),
            "wall_s": round(time.time() - t0, 1), "workers": a.workers, "solver": "glpk",
            "mev_serum_sweeps": int(ser.key.nunique()),
            "mev_serum_sweeps_core_all_lethal": int(ser.groupby("key").lethal.all().sum()),
            "mev_core_max_ratio_serum": float(ser.ratio_r.max()) if len(ser) else None}
    outs.append(mc.atomic_write_text(json.dumps(summ, indent=1) + "\n", OUT / "sweeps.json"))
    prov.output(*outs)
    prov.note(**{k: v for k, v in summ.items() if k != "groups"})
    print(json.dumps(summ, indent=1), flush=True)
    return 0 if spot_eq else 3


if __name__ == "__main__":
    raise SystemExit(main())
