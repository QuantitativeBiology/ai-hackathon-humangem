"""W1: task-based knockout calls on the generic parents (DESIGN.md §6.8, call rules P57 and PGR).

P57: a gene is essential if its knockout makes >= 1 of the 57 tasks of the release's own
``metabolicTasks_Essential.txt`` infeasible (each task with its own inputs; boundaries closed).
PGR: essential if the ``GR`` task (growth on Ham's) becomes infeasible.

The algorithm is that of Human-GEM ``code/test/taskEssentialGenes.py`` (one pFBA flux set per task; a
knockout is re-tested only on tasks whose flux-carrying reactions it disables; genes that disable an
identical reaction set share one test), with three changes required by DESIGN §11.4:
- every task touched by a knockout is tested (no early stop), so the full broken-task set is known;
- solver status is strict: ``optimal`` = feasible, ``infeasible`` = broken; anything else is retried
  once with NumericFocus = 2 and otherwise aborts the job (no silent "essential");
- the wild type must pass every non-should-fail task, otherwise the job aborts.
Equality with ``taskEssentialGenes.find_task_essential_categories`` is checked by
``scripts/w1_verify.py`` (§12, "the P57 call equals taskEssentialGenes.py on the generic parents").

Gurobi sessions: the parent plus ``--workers`` spawn workers (default 7, i.e. 8 sessions in total;
the coordinator's cap). No flock slot pool exists yet (DESIGN §11.4 gap); a licence / environment error
raises in the worker and aborts the whole job.

Usage: .venv/bin/python scripts/w1_p57.py --model human2 [--workers 7]
Output: results/w1/p57_<model>.parquet (+ .json run record). Label-free.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import multiprocessing
import os
import resource
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w1_common import OUT_DIR, PROJECT, atomic_write_parquet, current, run_main  # noqa: E402

TOL = 1e-8
_W: dict = {}


class SolverInfraError(RuntimeError):
    """Licence, environment or persistent numerical failure: abort with no output (DESIGN §11.4)."""


def _pin(model, seed: int) -> None:
    p = model.solver.problem.Params
    p.Threads = 1
    p.Seed = seed


def _restore(base, met_ids, saved):
    from raven_toolbox.tasks.check import _set_constraint_bounds

    for mid in met_ids:
        if mid in base.constraints:
            _set_constraint_bounds(base.constraints[mid], *saved[mid])


def _status_feasible(base) -> bool:
    """optimal -> True, infeasible -> False, anything else -> one NumericFocus=2 retry, else raise."""
    base.slim_optimize()
    st = base.solver.status
    if st == "optimal":
        return True
    if st == "infeasible":
        return False
    p = base.solver.problem.Params
    old = p.NumericFocus
    p.NumericFocus = 2
    try:
        base.slim_optimize()
        st2 = base.solver.status
    finally:
        p.NumericFocus = old
    if st2 == "optimal":
        return True
    if st2 == "infeasible":
        return False
    raise SolverInfraError(f"non-optimal LP status {st!r} then {st2!r} after NumericFocus=2")


def _task_flux_set(base, task, maps, original_ids, saved):
    from cobra.exceptions import OptimizationError
    from cobra.flux_analysis import pfba
    from raven_toolbox.tasks.check import apply_task_constraints

    with base:
        task_mets, error = apply_task_constraints(base, task, *maps)
        try:
            if error is not None:
                return None, f"task error: {error}"
            if not _status_feasible(base):
                return None, "infeasible on the wild type"
            try:
                fl = pfba(base).fluxes
            except OptimizationError as e:  # feasible LP but pFBA failed: infrastructure problem
                raise SolverInfraError(f"pFBA failed on a feasible task: {e}") from e
            return {r for r in original_ids if abs(fl.get(r, 0.0)) > TOL}, None
        finally:
            _restore(base, task_mets or saved, saved)


def _feasible_without(base, task, maps, rxns, saved) -> bool:
    from raven_toolbox.tasks.check import apply_task_constraints

    with base:
        task_mets, error = apply_task_constraints(base, task, *maps)
        try:
            if error is not None:
                raise SolverInfraError(f"task error on a knockout: {error}")
            for rid in rxns:
                base.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            return _status_feasible(base)
        finally:
            _restore(base, task_mets or saved, saved)


def _init(base, passing, maps, saved, seed):
    _pin(base, seed)
    _W.update(base=base, passing=passing, maps=maps, saved=saved)


def _work(key):
    broken, solves = [], 0
    for idx, task, fset in _W["passing"]:
        if not (key & fset):
            continue
        solves += 1
        if not _feasible_without(_W["base"], task, _W["maps"], key, _W["saved"]):
            broken.append(idx)
    return key, broken, solves


def gene_disabled_reactions(model) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for rxn in model.reactions:
        ids = [g.id for g in rxn.genes]
        for gid in ids:
            if not rxn.gpr.eval([gid]):
                out[gid].add(rxn.id)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=["human1", "human2"])
    ap.add_argument("--workers", type=int, default=7)
    a = ap.parse_args(argv)
    if a.workers + 1 > 8:
        raise SystemExit("Gurobi cap: parent + workers must be <= 8 sessions")

    from nextgem.config import SEEDS
    from nextgem.model import GEMS, load_model
    from raven_toolbox.tasks.check import task_name_maps
    from raven_toolbox.tasks.tasklist import parse_task_list

    t0, c0 = time.time(), time.process_time()
    gem = GEMS[a.model]
    model = load_model(a.model)
    model.solver = "gurobi"
    tasks = parse_task_list(gem.essential_tasks)
    base = model.copy()
    for r in base.boundary:
        r.bounds = (0.0, 0.0)
    _pin(base, SEEDS[0])
    maps = task_name_maps(base)
    original = {r.id for r in base.reactions}
    saved = {m.id: (base.constraints[m.id].lb, base.constraints[m.id].ub) for m in base.metabolites}

    testable = [(i, t) for i, t in enumerate(tasks) if not t.should_fail]
    passing, failing = [], []
    for i, t in testable:
        fset, why = _task_flux_set(base, t, maps, original, saved)
        if fset is None:
            failing.append((i, t.id, t.description, why))
        else:
            passing.append((i, t, fset))
    print(f"[{a.model}] {len(tasks)} tasks, {len(testable)} testable, {len(passing)} pass on the wild type", flush=True)
    if failing:
        raise SolverInfraError(f"wild type fails tasks (DESIGN requires all to pass): {failing}")

    gd = gene_disabled_reactions(base)
    groups: dict[frozenset, list[str]] = defaultdict(list)
    for g, rx in gd.items():
        if rx:
            groups[frozenset(rx)].append(g)
    print(f"[{a.model}] {len(gd)} genes with reactions, {len(groups)} distinct knockouts", flush=True)

    results, solves = {}, 0
    ctx = multiprocessing.get_context("spawn")
    with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init,
                                initargs=(base, passing, maps, saved, SEEDS[0])) as pool:
        for n, (key, broken, s) in enumerate(pool.map(_work, list(groups), chunksize=4), 1):
            results[key] = broken
            solves += s
            if n % 200 == 0:
                print(f"[{a.model}] {n}/{len(groups)} knockouts, {solves} solves, {time.time() - t0:.0f} s", flush=True)

    task_by_idx = {i: t for i, t in enumerate(tasks)}
    gr_idx = [i for i, t in testable if t.id == "GR"]
    if len(gr_idx) != 1:
        raise RuntimeError(f"expected exactly one GR task, found {gr_idx}")
    ko = pd.read_parquet(PROJECT / "results" / "baseline" / f"{a.model}_ko.parquet", columns=["gene", "entrez", "symbol"])
    rows = []
    for g in ko.gene:
        rx = gd.get(g, set())
        broken = sorted(results.get(frozenset(rx), [])) if rx else []
        rows.append({"gene": g, "n_disabled_rxns": len(rx),
                     "broken_tasks": json.dumps(broken),
                     "broken_categories": ",".join(sorted({task_by_idx[i].id for i in broken})),
                     "p57": bool(broken), "pgr": gr_idx[0] in broken})
    out = ko.merge(pd.DataFrame(rows), on="gene", how="left", validate="one_to_one")
    if len(out) != len(ko) or out.p57.isna().any():
        raise RuntimeError("gene table mismatch")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"p57_{a.model}.parquet"
    atomic_write_parquet(out, path)
    wall, cpu = time.time() - t0, time.process_time() - c0
    rec = {"model": a.model, "version": gem.version, "tasks_file": str(gem.essential_tasks), "n_tasks": len(tasks),
           "n_testable": len(testable), "n_passing_wt": len(passing), "n_distinct_knockouts": len(groups),
           "n_solves": solves, "p57_essential": int(out.p57.sum()), "pgr_essential": int(out.pgr.sum()),
           "workers": a.workers, "gurobi_sessions": a.workers + 1, "solver": "gurobi",
           "wall_s": round(wall, 1), "parent_cpu_s": round(cpu, 1),
           "peak_rss_parent_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)}
    (OUT_DIR / f"p57_{a.model}.json").write_text(json.dumps(rec, indent=1))
    current().input(gem.sbml, gem.essential_tasks, PROJECT / "results" / "baseline" / f"{a.model}_ko.parquet")
    current().output(path, OUT_DIR / f"p57_{a.model}.json")
    print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(run_main(main, seeds=[11], solver="gurobi", workers=7))
