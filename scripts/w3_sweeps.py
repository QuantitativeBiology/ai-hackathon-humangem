"""W3 biomass-objective KO sweeps (DESIGN §9 W3): T0 (re-check), D-BS, D-ER, O-NGAM_k, B-loo_j, B-add_c (1e-4; 1e-3
sensitivity on M0), and the final-list T10 biomass, on M0 and on the documented-variant proxy (RPMI M2(f = 1)).

Variant specifications come from ``scripts/w3_structure.py`` (``results/w3/structure/variants_<parent>.json``), media
bounds from ``media_<parent>.json``. Exact KO engine with the W2 shortcuts (``nextgem.w3.ko_sweep``). T0 on each
medium must equal the cached sweep bitwise (results/baseline for M0; results/w2/ko for RPMI M2(f = 1)); the job aborts
otherwise.

Outputs (results/w3/sweeps/): ``ratio_r_<parent>.parquet`` (rows = sweep id ``<medium>|<variant>``, columns = genes),
``sweeps_<parent>.tsv`` (WT growth, LPs, shortcut verification, CPU). Label-free.

Usage: .venv/bin/python scripts/w3_sweeps.py --model human2 --workers 7
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
from w3_common import MAX_WORKERS, OUT_DIR, atomic_write_parquet, atomic_write_tsv, current, run_main  # noqa: E402

SWDIR = OUT_DIR / "sweeps"
SDIR = OUT_DIR / "structure"
_W: dict = {}


def _init(parent):
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import w3
    from nextgem.model import load_model

    m = load_model(parent)
    _W.update(model=m, genes=sorted(g.id for g in m.genes), closures=w3.gene_closures(m),
              media=json.loads((SDIR / f"media_{parent}.json").read_text()))


def _job(job):
    from nextgem import w3

    medium, vjson = job
    t0 = time.process_time()
    v = w3.Variant.from_json(vjson)
    r = w3.ko_sweep(_W["model"], _W["genes"], _W["closures"], _W["media"][medium], v)
    r.update(medium=medium, variant=v.name, cpu_s=time.process_time() - t0)
    return r


def jobs_for(variants: dict, media=("M0", "RPMI_M2f1")) -> list[tuple[str, dict]]:
    from nextgem import w3

    jobs = []
    for med in media:
        jobs.append((med, w3.Variant("T0").to_json()))
        for key in ("bloo", "t10", "badd"):
            jobs += [(med, v) for v in variants.get(key, [])]
        if med == "M0":
            jobs += [(med, v) for v in variants.get("badd_sens", [])]
        for key in (f"D-BS|{med}", f"D-ER|{med}", f"ngam|{med}"):
            jobs += [(med, v) for v in variants.get(key, [])]
    return jobs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=["human1", "human2"])
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--only", nargs="*", default=None, help="variant names to run (pilot)")
    ap.add_argument("--merge", action="store_true",
                    help="with --only: add/replace these sweeps in the main outputs instead of writing a pilot")
    a = ap.parse_args(argv)
    if a.workers > MAX_WORKERS:
        raise SystemExit("at most 8 worker processes")
    from nextgem import paths

    t0 = time.time()
    variants = json.loads((SDIR / f"variants_{a.model}.json").read_text())
    jobs = jobs_for(variants)
    if a.only:
        jobs = [j for j in jobs if j[1]["name"] in a.only]
    print(f"[{a.model}] {len(jobs)} sweeps, {a.workers} workers", flush=True)
    res = []
    ctx = multiprocessing.get_context("spawn")
    with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init, initargs=(a.model,)) as pool:
        for n, r in enumerate(pool.map(_job, jobs, chunksize=1), 1):
            res.append(r)
            print(f"  {n}/{len(jobs)} {r['medium']}|{r['variant']}: WT {r['wt_growth']:.4g}, {r['n_solved']} LPs, "
                  f"{r['cpu_s']:.0f} s CPU, {time.time() - t0:.0f} s", flush=True)
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem.model import load_model

    genes = sorted(g.id for g in load_model(a.model).genes)
    R = pd.DataFrame(np.vstack([r["ratio_r"] for r in res]), columns=genes,
                     index=pd.Index([f"{r['medium']}|{r['variant']}" for r in res], name="sweep"))
    # T0 equality with the cached sweeps (bitwise)
    base = pd.read_parquet(paths.baseline_dir() / f"{a.model}_ko.parquet", columns=["gene", "ratio_r"]).set_index("gene").ratio_r
    w2 = pd.read_parquet(paths.results_dir() / "w2" / "ko" / f"ratio_r_{a.model}.parquet")
    sid = json.loads((SDIR / "structure_summary.json").read_text())["w2_sweep_ids"][a.model]["RPMI_M2f1"]
    checks = {}
    if "M0|T0" in R.index:
        checks["M0_T0_equals_baseline"] = bool(np.array_equal(R.loc["M0|T0", genes].to_numpy(), base.reindex(genes).to_numpy(), equal_nan=True))
    if "RPMI_M2f1|T0" in R.index:
        checks["RPMI_T0_equals_w2"] = bool(np.array_equal(R.loc["RPMI_M2f1|T0", genes].to_numpy(), w2.loc[sid, genes].to_numpy(), equal_nan=True))
    if not all(checks.values()):
        raise RuntimeError(f"T0 equality failed: {checks}")
    meta = pd.DataFrame([{k: v for k, v in r.items() if k != "ratio_r"} for r in res])
    meta["n_lethal"] = [(r["ratio_r"] < 0.01).sum() for r in res]
    meta["n_partial"] = [((r["ratio_r"] >= 0.01) & (r["ratio_r"] < 1)).sum() for r in res]
    if a.merge:
        if not a.only:
            raise SystemExit("--merge needs --only")
        old_R = pd.read_parquet(SWDIR / f"ratio_r_{a.model}.parquet").set_index("sweep")
        old_m = pd.read_csv(SWDIR / f"sweeps_{a.model}.tsv", sep="\t")
        R = pd.concat([old_R.drop(index=[i for i in R.index if i in old_R.index]), R])
        key = meta.medium + "|" + meta.variant
        old_m = old_m[~(old_m.medium + "|" + old_m.variant).isin(set(key))]
        meta = pd.concat([old_m, meta.assign(added_by_merge=True)], ignore_index=True)
    out_dir = SWDIR / ("pilot" if (a.only and not a.merge) else "")
    pq = out_dir / f"ratio_r_{a.model}.parquet"
    out_dir.mkdir(parents=True, exist_ok=True)
    R.reset_index().to_parquet(pq.with_name(pq.name + ".tmp"), index=False, compression="zstd")
    pq.with_name(pq.name + ".tmp").replace(pq)
    mt = atomic_write_tsv(meta, out_dir / f"sweeps_{a.model}.tsv")
    rec = {"model": a.model, "n_sweeps": len(res), "checks": checks, "job_cpu_s": round(float(meta.cpu_s.sum()), 1),
           "wall_s": round(time.time() - t0, 1), "workers": a.workers, "solver": "glpk",
           "max_shortcut_dev": float(np.nanmax(meta.shortcut_verify_max_dev))}
    (out_dir / f"sweeps_{a.model}.json").write_text(json.dumps(rec, indent=1))
    current().input(SDIR / f"variants_{a.model}.json", SDIR / f"media_{a.model}.json",
                    paths.baseline_dir() / f"{a.model}_ko.parquet")
    current().output(pq, mt, out_dir / f"sweeps_{a.model}.json")
    current().note(**{k: v for k, v in rec.items() if k in ("n_sweeps", "checks", "job_cpu_s")})
    print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 7
    raise SystemExit(run_main(main, solver="glpk", workers=w))
