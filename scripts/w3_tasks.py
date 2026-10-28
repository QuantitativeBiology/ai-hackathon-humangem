"""W3 capacity matrix (DESIGN §9 W3, compute SF3): gene × task feasibility and max-output capacity after KO.

For one parent and one task list, every task is tested with its own inputs and every other boundary closed (the
Human-GEM / P57 protocol of DESIGN §6.8, which makes the task part medium-independent; the medium enters W3 only
through the biomass column, the demands and NGAM, see scripts/w3_sweeps.py). Deduplicated KO sets: genes whose KO
disables the identical reaction set share one test. Exact flux-support shortcut (see ``nextgem.w3``).

Lists:
  essential  : the release's 57 Essential tasks (incl. GR; P57 equality test input)
  cellfie    : the 193 parsable CellFie consensus tasks (Human2 file, both parents); also task-essential reactions
               (single-reaction removal over the WT pFBA support), used by the CellFie-inspired T5 activity rule
  full_subset: Full-list tasks 42–50, 52–68, 77–99, 104–118 (final a-priori list §5 sets T5/T6/T8; should-fail
               tasks excluded, logged)
  pseudo     : pseudo-task twins of the 56 non-GR Essential tasks for SEEDS[0:5] (negative control); feasibility only

Outputs (results/w3/tasks/): ``<list>_<parent>.parquet`` (gene × task long: feasible, capacity),
``<list>_<parent>_tasks.tsv`` (task meta: status, λ_WT, support sizes, solves, failing-gene counts),
``pseudo_<parent>_defs.tsv`` (pseudo-task swaps; candidate FROZEN rows). Label-free.

Usage: .venv/bin/python scripts/w3_tasks.py --model human2 --list essential --workers 7
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

TDIR = OUT_DIR / "tasks"
FULL_SUBSET = list(range(42, 51)) + list(range(52, 69)) + list(range(77, 100)) + list(range(104, 119))  # 1-based IDs
_W: dict = {}


def _setup(parent: str, which: str):
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import w3
    from nextgem.model import load_model
    from raven_toolbox.tasks.check import task_name_maps

    m = load_model(parent)
    base = w3.task_base(m)
    maps = task_name_maps(base)
    saved = {x.id: (base.constraints[x.id].lb, base.constraints[x.id].ub) for x in base.metabolites}
    original = {r.id for r in base.reactions}
    gd = w3.gene_disabled_reactions(base)
    groups = w3.ko_groups(gd)
    kos = sorted(groups, key=lambda s: tuple(sorted(s)))
    lst = "essential" if which == "pseudo" else ("full" if which == "full_subset" else which)
    tasks = w3.load_tasks(parent, lst)
    return dict(base=base, maps=maps, saved=saved, original=original, gd=gd, groups=groups, kos=kos, tasks=tasks,
                parent=parent, which=which)


def _init(parent, which):
    _W.update(_setup(parent, which))


def _comp_tokens():
    if "comp" not in _W:
        from nextgem import w3

        _W["comp"] = w3.compartment_tokens(_W["base"])
    return _W["comp"]


CKPT = OUT_DIR / "tmp" / "task_ckpt"  # per-job checkpoints (git-ignored; a restarted run reuses finished jobs)


def _job(job):
    import pickle

    ck = CKPT / f"{_W['parent']}_{_W['which']}_{job[2]}_{job[1]}.pkl"
    if ck.exists():
        return pickle.loads(ck.read_bytes())
    res = _job_compute(job)
    CKPT.mkdir(parents=True, exist_ok=True)
    tmp = ck.with_suffix(".tmp")
    tmp.write_bytes(pickle.dumps(res))
    tmp.replace(ck)
    return res


def _job_compute(job):
    from nextgem import w3

    kind, idx, seed = job
    t0 = time.process_time()
    task = _W["tasks"].task.iloc[idx]
    meta = {"idx": idx, "seed": seed, "swaps": None, "tries": None}
    if kind == "pseudo":
        rng = np.random.default_rng([seed, idx])
        res = w3.make_pseudo_task(task, _W["base"], _W["maps"], _W["saved"], rng, _comp_tokens())
        if res is None:
            meta.update(status="pseudo_not_generated", cpu_s=time.process_time() - t0)
            return meta, None, None
        task, swaps, tries = res
        meta.update(swaps=json.dumps(swaps), tries=tries)
    r = w3.task_capacity(_W["base"], task, _W["maps"], _W["saved"], _W["kos"], _W["original"],
                         essential_reactions=(_W["which"] == "cellfie"),
                         capacity=_W["which"] in ("essential", "full_subset"))
    meta.update(status=r["status"], error=r["error"], lam_wt=r["lam_wt"], n_required=r.get("n_required"),
                fset_size=r.get("fset_size"), lset_size=r.get("lset_size"), n_solves=r["n_solves"],
                ess_rxns=None if r["ess_rxns"] is None else json.dumps(r["ess_rxns"]), cpu_s=time.process_time() - t0)
    if r["status"] != "ok":
        return meta, None, None
    feas = np.array([r["feasible"][k] for k in _W["kos"]], dtype=bool)
    cap = np.array([w3.capacity_ratio(r["lam"][k], r["lam_wt"]) for k in _W["kos"]], dtype=float)
    return meta, feas, cap


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=["human1", "human2"])
    ap.add_argument("--list", required=True, choices=["essential", "cellfie", "full_subset", "pseudo"])
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--seeds", type=int, nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None, help="pilot: first N jobs only (outputs under tmp/)")
    a = ap.parse_args(argv)
    if a.workers + 1 > MAX_WORKERS + 1:
        raise SystemExit("at most 8 worker processes")
    from nextgem.config import SEEDS
    from nextgem.model import GEMS
    from nextgem import w3

    t0 = time.time()
    S = _setup(a.model, a.list)
    tasks = S["tasks"]
    if a.list == "full_subset":
        sel = [i for i in range(len(tasks)) if int(tasks.category.iloc[i]) in FULL_SUBSET]
    elif a.list == "pseudo":
        sel = [i for i in range(len(tasks)) if tasks.category.iloc[i] != "GR"]
    else:
        sel = list(range(len(tasks)))
    skipped_sf = [tasks.key.iloc[i] for i in sel if tasks.should_fail.iloc[i]]
    sel = [i for i in sel if not tasks.should_fail.iloc[i]]
    seeds = (a.seeds or SEEDS[0:5]) if a.list == "pseudo" else [None]
    jobs = [("pseudo" if a.list == "pseudo" else "task", i, s) for s in seeds for i in sel]
    if a.limit:
        jobs = jobs[: a.limit]
    print(f"[{a.model}/{a.list}] {len(S['kos'])} distinct KOs, {len(jobs)} jobs, {a.workers} workers", flush=True)
    metas, F, C = [], {}, {}
    ctx = multiprocessing.get_context("spawn")
    with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init,
                                initargs=(a.model, a.list)) as pool:
        for n, (meta, feas, cap) in enumerate(pool.map(_job, jobs, chunksize=1), 1):
            metas.append(meta)
            if feas is not None:
                F[(meta["seed"], meta["idx"])] = feas
                C[(meta["seed"], meta["idx"])] = cap
            if n % 20 == 0:
                print(f"  {n}/{len(jobs)} jobs, {time.time() - t0:.0f} s", flush=True)
    # gene-level long table
    gene_ko = {g: rx for g, rx in S["gd"].items() if rx}
    ko_index = {k: i for i, k in enumerate(S["kos"])}
    from nextgem.paths import baseline_dir

    gm = pd.read_parquet(baseline_dir() / f"{a.model}_ko.parquet", columns=["gene", "entrez", "symbol"])
    genes = list(gm.gene)
    gi = np.array([ko_index.get(gene_ko.get(g), -1) for g in genes])
    rows = []
    for (seed, idx), feas in F.items():
        cap = C[(seed, idx)]
        f = np.where(gi >= 0, feas[np.maximum(gi, 0)], True)
        c = np.where(gi >= 0, cap[np.maximum(gi, 0)], 1.0)
        rows.append(pd.DataFrame({"task": tasks.key.iloc[idx], "seed": -1 if seed is None else seed, "gene": genes,
                                  "feasible": f, "capacity": c}))
    long = pd.concat(rows, ignore_index=True).merge(gm, on="gene", how="left")
    meta = pd.DataFrame(metas)
    meta["task"] = meta.idx.map(lambda i: tasks.key.iloc[i])
    meta["description"] = meta.idx.map(lambda i: tasks.description.iloc[i])
    meta["category"] = meta.idx.map(lambda i: tasks.category.iloc[i])
    meta["group"] = meta.idx.map(lambda i: tasks.group.iloc[i])
    nfail = long.groupby(["task", "seed"]).feasible.apply(lambda s: int((~s).sum()))
    meta["n_genes_failing"] = [nfail.get((t, -1 if s is None else s), np.nan) for t, s in zip(meta.task, meta.seed)]
    meta["seed"] = meta.seed.fillna(-1).astype(int)
    out_dir = (TDIR / "pilot") if a.limit else TDIR
    stem = f"{a.list}_{a.model}"
    pq = atomic_write_parquet(long, out_dir / f"{stem}.parquet")
    mt = atomic_write_tsv(meta.drop(columns=["swaps"]) if a.list != "pseudo" else meta, out_dir / f"{stem}_tasks.tsv")
    rec = {"model": a.model, "list": a.list, "n_kos": len(S["kos"]), "n_jobs": len(jobs),
           "status_counts": meta.status.value_counts().to_dict(), "skipped_should_fail": skipped_sf,
           "n_solves": int(meta.n_solves.fillna(0).sum()), "job_cpu_s": round(float(meta.cpu_s.sum()), 1),
           "wall_s": round(time.time() - t0, 1), "workers": a.workers, "solver": "glpk"}
    (out_dir / f"{stem}.json").write_text(json.dumps(rec, indent=1))
    current().input(GEMS[a.model].sbml, w3.task_file(a.model, "essential" if a.list == "pseudo" else
                                                       ("full" if a.list == "full_subset" else a.list)))
    current().output(pq, mt, out_dir / f"{stem}.json")
    current().note(**{k: v for k, v in rec.items() if k in ("n_kos", "n_jobs", "n_solves", "job_cpu_s")})
    print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    w = 7
    if "--workers" in sys.argv:
        w = int(sys.argv[sys.argv.index("--workers") + 1])
    raise SystemExit(run_main(main, seeds=[11, 23, 37, 41, 53] if "pseudo" in sys.argv else [], solver="glpk", workers=w))
