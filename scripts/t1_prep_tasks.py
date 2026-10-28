"""T1 prep (5): the W3 Essential-task capacity matrix re-run on P1's model, i.e. the parent with the W3 biomass
variant B-loo:cofactor_pool_biomass (DESIGN §6.8 P57/PGR, §9 W3 per-task table). Label-free.

Protocol = scripts/w3_tasks.py --list essential (Human-GEM taskEssentialGenes / P57): each of the 57 tasks of the
release's metabolicTasks_Essential.txt with its own inputs and every other boundary closed; deduplicated KO sets
(genes whose KO disables the identical reaction set share one test); exact flux-support shortcut. Only feasibility is
computed (``capacity=False``): P57/PGR and the per-task table need nothing else. The only change is that the B-loo
variant (``nextgem.w3.Variant`` from results/w3/structure/variants_<parent>.json) is applied to the model first.
Under the paper protocol the GR task is growth on Ham's, so the medium (M2) does not enter P57/PGR (M6 README §3).

Solver and limits (A10): GLPK; 300 s per LP + one retry (``nextgem.w3.solve_status``); 1,800 s CPU per job
(ITIMER_PROF, as scripts/m6_sweeps.py). A failed job aborts the run (no partial matrix is written).

Also compared with the generic W3 matrix (results/w3/tasks/essential_<parent>.parquet): per task, the number of genes
whose feasibility differs. Expected: 0 for the 56 non-GR tasks (they never carry biomass flux), changes only in GR.

Outputs (results/m6/t1/): tasks_bloo_<parent>.parquet (task, gene, entrez, feasible), tasks_bloo_<parent>_tasks.tsv
(per-task meta + n_genes_failing + n_genes_changed_vs_generic), tasks_bloo_<parent>.json (run summary).

Usage: .venv/bin/python scripts/t1_prep_tasks.py --model human1 --workers 8
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
import t1_prep_common as tc  # noqa: E402

W3 = tc.PROJECT / "results" / "w3"
JOB_CPU_LIMIT_S = 1800
_W: dict = {}


class JobCPUTimeout(Exception):
    pass


def _on_prof(signum, frame):
    raise JobCPUTimeout()


def bloo_variant(parent: str):
    from nextgem import w3

    vj = json.loads((W3 / "structure" / f"variants_{parent}.json").read_text())
    hit = [v for v in vj["bloo"] if v["name"] == mc.BLOO]
    if len(hit) != 1:
        raise RuntimeError(f"{mc.BLOO} not found once in variants_{parent}.json")
    return w3.Variant.from_json(hit[0])


def setup(parent: str) -> dict:
    import cobra

    cobra.Configuration().solver = "glpk"
    from raven_toolbox.tasks.check import task_name_maps

    from nextgem import w3
    from nextgem.model import load_model

    m = load_model(parent)
    v = bloo_variant(parent)
    v.apply(m)  # permanent on this process's copy
    b = m.reactions.get_by_id(w3.BIOMASS)
    removed = [x for x in v.params["mets"] if x in {mm.id for mm in b.metabolites}]
    if removed:
        raise RuntimeError(f"B-loo not applied: {removed} still in {w3.BIOMASS}")
    base = w3.task_base(m)
    maps = task_name_maps(base)
    saved = {x.id: (base.constraints[x.id].lb, base.constraints[x.id].ub) for x in base.metabolites}
    original = {r.id for r in base.reactions}
    gd = w3.gene_disabled_reactions(base)
    groups = w3.ko_groups(gd)
    kos = sorted(groups, key=lambda s: tuple(sorted(s)))
    tasks = w3.load_tasks(parent, "essential")
    return dict(base=base, maps=maps, saved=saved, original=original, gd=gd, groups=groups, kos=kos, tasks=tasks,
                parent=parent, variant=v.to_json())


def _init(parent):
    _W.update(setup(parent))
    signal.signal(signal.SIGPROF, _on_prof)


def _job(idx: int):
    from nextgem import w3

    t0 = time.process_time()
    task = _W["tasks"].task.iloc[idx]
    meta = {"idx": idx}
    signal.setitimer(signal.ITIMER_PROF, float(JOB_CPU_LIMIT_S))
    try:
        r = w3.task_capacity(_W["base"], task, _W["maps"], _W["saved"], _W["kos"], _W["original"],
                             essential_reactions=False, capacity=False)
    except JobCPUTimeout:
        meta.update(status="failed-timeout", cpu_s=time.process_time() - t0)
        return meta, None
    except w3.SolverInfraError as e:
        meta.update(status="failed-nonoptimal", error=str(e)[:200], cpu_s=time.process_time() - t0)
        return meta, None
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0.0)
    meta.update(status=r["status"], error=r["error"], n_required=r.get("n_required"), fset_size=r.get("fset_size"),
                n_solves=r["n_solves"], cpu_s=time.process_time() - t0)
    if r["status"] != "ok":
        return meta, None
    return meta, np.array([r["feasible"][k] for k in _W["kos"]], dtype=bool)


def gene_long(S: dict, F: dict, gm: pd.DataFrame) -> pd.DataFrame:
    """(task, gene, entrez, feasible) for every model gene in the baseline map; genes whose KO closes nothing: True."""
    gene_ko = {g: rx for g, rx in S["gd"].items() if rx}
    ko_index = {k: i for i, k in enumerate(S["kos"])}
    genes = list(gm.gene)
    gi = np.array([ko_index.get(gene_ko.get(g), -1) for g in genes])
    rows = []
    for idx, feas in sorted(F.items()):
        f = np.where(gi >= 0, feas[np.maximum(gi, 0)], True)
        rows.append(pd.DataFrame({"task": S["tasks"].key.iloc[idx], "gene": genes, "feasible": f}))
    long = pd.concat(rows, ignore_index=True).merge(gm[["gene", "entrez"]], on="gene", how="left")
    return long[["task", "gene", "entrez", "feasible"]]


def compare_generic(long: pd.DataFrame, generic: pd.DataFrame) -> pd.Series:
    """task -> number of genes whose feasibility differs from the generic W3 matrix (seed -1)."""
    g = generic[generic.seed == -1][["task", "gene", "feasible"]]
    j = long.merge(g, on=["task", "gene"], how="outer", suffixes=("", "_generic"), indicator=True)
    if (j._merge != "both").any():
        raise RuntimeError("B-loo and generic task matrices cover different (task, gene) pairs")
    return (j.feasible.astype(bool) != j.feasible_generic.astype(bool)).groupby(j.task).sum().astype(int)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=list(mc.PARENTS))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None, help="pilot: first N tasks only; outputs under results/m6/tmp/t1/")
    a = ap.parse_args(argv)
    if not 1 <= a.workers <= tc.MAX_WORKERS:
        raise SystemExit(f"workers must be 1..{tc.MAX_WORKERS}")
    from nextgem import w3
    from nextgem.model import GEMS
    from nextgem.paths import baseline_dir

    t0 = time.time()
    S = setup(a.model)
    tasks = S["tasks"]
    sel = [i for i in range(len(tasks)) if not tasks.should_fail.iloc[i]]
    if len(sel) != len(tasks) or len(tasks) != 57:
        raise RuntimeError(f"expected 57 non-should-fail Essential tasks, got {len(sel)}/{len(tasks)}")
    if a.limit:
        sel = sel[: a.limit]
    print(f"[{a.model} + {mc.BLOO}] {len(S['kos'])} distinct KOs, {len(sel)} tasks, {a.workers} workers", flush=True)
    out_dir = tc.T1_TMP / "pilot" if a.limit else tc.T1
    stem = f"tasks_bloo_{a.model}"
    gm = pd.read_parquet(baseline_dir() / f"{a.model}_ko.parquet", columns=["gene", "entrez"])
    generic_f = W3 / "tasks" / f"essential_{a.model}.parquet"
    inputs = [GEMS[a.model].sbml, w3.task_file(a.model, "essential"), W3 / "structure" / f"variants_{a.model}.json",
              baseline_dir() / f"{a.model}_ko.parquet", generic_f]
    with tc.provenance(inputs=inputs, solver="glpk", workers=a.workers) as prov:
        metas, F = [], {}
        ctx = multiprocessing.get_context("spawn")
        with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init,
                                    initargs=(a.model,)) as pool:
            for n, (meta, feas) in enumerate(pool.map(_job, sel, chunksize=1), 1):
                metas.append(meta)
                if feas is not None:
                    F[meta["idx"]] = feas
                if n % 10 == 0:
                    print(f"  {n}/{len(sel)} tasks, {time.time() - t0:.0f} s", flush=True)
        meta = pd.DataFrame(metas).sort_values("idx")
        bad = meta[meta.status != "ok"]
        if len(bad):
            raise RuntimeError(f"failed task jobs (no matrix written): {bad[['idx', 'status']].to_dict('records')}")
        long = gene_long(S, F, gm)
        generic = pd.read_parquet(generic_f, columns=["task", "seed", "gene", "feasible"])
        if a.limit:
            generic = generic[generic.task.isin(set(long.task))]
        changed = compare_generic(long, generic)
        meta["task"] = meta.idx.map(lambda i: tasks.key.iloc[i])
        meta["category"] = meta.idx.map(lambda i: tasks.category.iloc[i])
        meta["description"] = meta.idx.map(lambda i: tasks.description.iloc[i])
        nfail = long.groupby("task").feasible.apply(lambda s: int((~s.astype(bool)).sum()))
        meta["n_genes_failing"] = meta.task.map(nfail)
        meta["n_genes_changed_vs_generic"] = meta.task.map(changed).fillna(0).astype(int)
        gr = meta.task[meta.category == "GR"].tolist()
        pq = tc.atomic_write_parquet(long, out_dir / f"{stem}.parquet")
        mt = tc.atomic_write_tsv(meta, out_dir / f"{stem}_tasks.tsv")
        rec = {"model": a.model, "variant": S["variant"], "list": "essential", "n_tasks": len(sel),
               "n_kos": len(S["kos"]), "gr_task": gr, "status_counts": meta.status.value_counts().to_dict(),
               "n_solves": int(meta.n_solves.fillna(0).sum()), "job_cpu_s": round(float(meta.cpu_s.sum()), 1),
               "wall_s": round(time.time() - t0, 1), "workers": a.workers, "solver": "glpk",
               "lp_time_limit_s": w3.LP_TIME_LIMIT, "job_cpu_limit_s": JOB_CPU_LIMIT_S,
               "n_genes_changed_vs_generic_nonGR": int(meta.n_genes_changed_vs_generic[meta.category != "GR"].sum()),
               "n_genes_changed_vs_generic_GR": int(meta.n_genes_changed_vs_generic[meta.category == "GR"].sum()),
               "n_genes_failing_GR": int(meta.n_genes_failing[meta.category == "GR"].sum())}
        js = tc.atomic_write_json(rec, out_dir / f"{stem}.json")
        prov.output(pq, mt, js)
        prov.note(**{k: rec[k] for k in ("n_kos", "n_tasks", "n_solves", "job_cpu_s")})
    print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
