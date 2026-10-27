"""W4 context-specific model extraction, genomic constraint, repair rule and KO sweep (DESIGN.md §9 W4).

Families (LP only; ftINIT/F1 excluded by D2):

* **F0** GPR pruning + repair (naive reference): close every reaction whose GPR is False when the
  line's *off* genes (signed score ``s_g < 0``) are non-functional.
* **F2** GIMME (Becker & Palsson 2008): biomass >= 0.9 x max on the line's medium; minimise
  ``sum_r c_r |v_r|`` with ``c_r = max(0, -s_r)``; keep reactions with ``s_r >= 0``, without
  evidence, or carrying flux; close the rest. troppo does not install (README), so this is ours.
* **F3** FASTCORE (Vlassis 2014): core = reactions with ``s_r >= 0``; consistency from cobra
  ``fastcc`` on the parent + medium (cached per parent/medium); FASTCORE's LP7/LP9 loop on top of
  cobra's LP7 helpers. Everything outside the returned support is closed.
* **F3b** CORDA (``corda`` 0.5.1, package defaults n=3, penalty 100, support 5); gene confidence from
  ``s_g`` (>=1: 3, [0,1): 2, [-1,0): 1, < -1: -1, no evidence: 0), reaction confidence with the
  package's own min/max GPR rule (= our presence convention).
* **F4** E-Flux (Colijn 2009): internal reactions with a GPR get ``|bound| <= max(10, 1000 cap_r /
  max cap)``, ``cap`` = GPR with AND = min, OR = sum of linear values ``2**s_g`` (0 for off genes);
  GPR-less reactions and exchanges keep their bounds. Not repaired (DESIGN: repair applies to F0,
  F2, F3, F3b, F5).
* **F5** LCP builder only (label-calibrated grid is fitted elsewhere): pruning from a given off-set.

Signed scores and their GPR aggregation (presence: AND = min, OR = max) treat "no evidence" (NaN)
so that ``s_r < 0`` iff the GPR is False with the off genes non-functional (see ``reaction_scores``).

Gene state in a line model (KO semantics; implementation choice mirroring ftINIT ``removeGenes =
true`` in Human-GEM's DepMap pipeline, README §Open items): genes lost under the genomic constraint
are non-functional everywhere. In presence families, an *open* reaction whose GPR is still True
with the off genes removed uses that pruned GPR (an off isozyme does not buffer); an open
reaction whose GPR is False under the off genes (kept by the method or re-opened by repair) keeps
its full GPR. E-Flux keeps full GPRs; a KO sets the gene's linear value to 0 and re-derives the bound.

Repair rule (fixed, DESIGN W4): one L1 LP ``min sum_r w_r |v_r|`` over the closed (repairable)
reactions with ``w_r = 1 + max(0, -s_r)``, s.t. biomass >= 1% of the generic WT on the medium; one
L1 LP per viability task (GR + the 9 ER tasks of the essential task list, closed boundaries as in
Human-GEM's task code), s.t. the task is feasible; re-open the union of reactions with ``|v| >
1e-6``. Reactions closed by an accepted loss are never re-opened.

Genomic constraint: ``accept_losses`` implements the order-independent drop on the parent + medium.
"""

from __future__ import annotations

import ast
import math
import time
import types
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from optlang.interface import INFEASIBLE, OPTIMAL
from optlang.symbolics import Zero

from nextgem.config import SEEDS
from nextgem.knockouts import rounded_ratio, solver_name
from nextgem.model import BIOMASS_REACTION

FAMILIES = ("F0", "F2", "F3", "F3b", "F4", "F5")
PRESENCE_FAMILIES = {"F0", "F2", "F3", "F3b", "F5"}
REPAIRED_FAMILIES = {"F0", "F2", "F3", "F3b", "F5"}
GIMME_FRACTION = 0.9
EFLUX_FLOOR = 10.0
EFLUX_SCALE = 1000.0
REPAIR_GROWTH_FRACTION = 0.01
LETHAL_FRACTION = 0.01  # a loss "kills WT growth" if growth < 1% of the generic WT (same floor as repair)
FLUX_TOL = 1e-6
FASTCORE_EPS = 1.0
VIABILITY_TASK_TYPES = ("GR", "ER")
LP_TIME_LIMIT_S = 300.0  # per LP solve, both solvers (review B2: nothing bounded a GLPK solve before)
CORDA_PARAMS = dict(n=3, penalty_factor=100, support=5)
# Amendment A4 (c): one declared solver per family, fixed for all reruns. F5 is not named in A4; it uses the
# F0 machinery (pruning + repair), whose GLPK = Gurobi check passed, so it is declared GLPK (flagged to the PI).
DECLARED_SOLVER = {"F0": "glpk", "F2": "glpk", "F4": "glpk", "F5": "glpk", "F3": "gurobi", "F3b": "gurobi"}
LOSS_RULE_SOLVER = "glpk"  # Stage-A phase 0 (family-independent accepted losses); fixed for all reruns


class BuildFailed(RuntimeError):
    """A failed build in the §11.4 sense (infeasible / no growth after the method's repair)."""


# ============================================================================= GPR machinery
def _compile(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.BoolOp):
        return ("and" if isinstance(node.op, ast.And) else "or", tuple(_compile(v) for v in node.values))
    raise TypeError(f"unsupported GPR node {ast.dump(node)}")


def compile_gprs(model) -> dict:
    """reaction id -> nested ('and'|'or', children) tuples / gene-id leaves (reactions with genes only)."""
    return {r.id: _compile(r.gpr.body) for r in model.reactions if r.genes and r.gpr.body is not None}


def gpr_true(node, false: set) -> bool:
    if isinstance(node, str):
        return node not in false
    op, ch = node
    if op == "and":
        return all(gpr_true(c, false) for c in ch)
    return any(gpr_true(c, false) for c in ch)


def gpr_genes(node) -> set:
    if isinstance(node, str):
        return {node}
    return set().union(*(gpr_genes(c) for c in node[1]))


def presence_score(node, s: dict) -> float:
    """AND = min, OR = max, with 'no evidence' (NaN) handled so that score < 0 <=> GPR False under
    the off genes {g: s_g < 0}: AND = min over known children (NaN if none); OR = max over known if it
    is >= 0 or no child is unknown, else NaN."""
    if isinstance(node, str):
        v = s.get(node, math.nan)
        return math.nan if v is None else v
    op, ch = node
    vals = [presence_score(c, s) for c in ch]
    known = [v for v in vals if not math.isnan(v)]
    if op == "and":
        return min(known) if known else math.nan
    if not known:
        return math.nan
    m = max(known)
    return m if (m >= 0 or len(known) == len(vals)) else math.nan


def capacity(node, lin: dict) -> float:
    """E-Flux capacity: AND = min, OR = sum (NaN = unknown: AND ignores it, OR becomes unknown)."""
    if isinstance(node, str):
        return lin.get(node, math.nan)
    op, ch = node
    vals = [capacity(c, lin) for c in ch]
    if op == "and":
        known = [v for v in vals if not math.isnan(v)]
        return min(known) if known else math.nan
    return math.nan if any(math.isnan(v) for v in vals) else float(sum(vals))


def reaction_scores(gprs: dict, s_g: pd.Series) -> pd.Series:
    s = {k: float(v) for k, v in s_g.items() if not pd.isna(v)}
    return pd.Series({rid: presence_score(n, s) for rid, n in gprs.items()}, dtype="float64")


# ============================================================================= parent context
@dataclass
class Parent:
    """Per-process parent state: configured model + caches shared by every line build."""

    model: object
    name: str
    medium_name: str
    gprs: dict
    wt_generic: float
    rxn_ids: list
    genes: list
    gene_rxns: dict
    entrez: dict
    tasks: list  # viability tasks (raven Task objects)
    all_tasks: list
    name_maps: tuple
    consistent: set | None = None  # FASTCORE consistent set (cobra fastcc), lazily loaded
    internal_with_gpr: set = field(default_factory=set)

    @property
    def solver(self) -> str:
        return solver_name(self.model)


def set_solver_params(model, seed: int = SEEDS[0]) -> None:
    """Gurobi: Threads = 1, Seed = SEEDS[0] (§7), TimeLimit; GLPK: optlang timeout. A solve that hits the
    limit is non-optimal and goes through ``_status_retry`` (one retry, then a failed solve, §11.4)."""
    if solver_name(model) == "gurobi":
        p = model.solver.problem.Params
        p.Threads = 1
        p.Seed = seed
        p.TimeLimit = LP_TIME_LIMIT_S
    elif solver_name(model) == "glpk":
        model.solver.configuration.timeout = int(LP_TIME_LIMIT_S)


def make_parent(model, name: str = "human2", medium=None, medium_name: str = "M0", tasks_file=None, genes_tsv=None) -> Parent:
    """Configure ``model`` in place (medium, biomass objective, solver params) and build the caches."""
    from raven_toolbox.tasks.tasklist import parse_task_list

    from nextgem.context_inputs import gene_table
    from nextgem.knockouts import configure
    from nextgem.model import tasks_file_for

    configure(model, medium)
    gt = gene_table(model, genes_tsv)
    tf = tasks_file_for(model, tasks_file)
    all_tasks = parse_task_list(str(tf))
    viability = [t for t in all_tasks if t.id in VIABILITY_TASK_TYPES]
    return parent_from_model(model, name, medium_name, gt.entrez.to_dict(), viability, all_tasks)


def parent_from_model(model, name: str = "toy", medium_name: str = "toy", entrez: dict | None = None,
                      tasks=(), all_tasks=()) -> Parent:
    """Parent for an already configured model (medium + objective set); used directly by toy tests."""
    from raven_toolbox.tasks.check import task_name_maps

    set_solver_params(model)
    if BIOMASS_REACTION in model.reactions and model.reactions.get_by_id(BIOMASS_REACTION).objective_coefficient == 1:
        _biomass_objective(model)
    wt = model.slim_optimize()
    if model.solver.status != OPTIMAL or wt <= 0:
        raise RuntimeError(f"generic parent does not grow on {medium_name}: {model.solver.status}")
    internal = {r.id for r in model.reactions if r.genes and not r.boundary}
    return Parent(
        model=model, name=name, medium_name=medium_name, gprs=compile_gprs(model), wt_generic=float(wt),
        rxn_ids=[r.id for r in model.reactions], genes=sorted(g.id for g in model.genes),
        gene_rxns={g.id: sorted(r.id for r in g.reactions) for g in model.genes},
        entrez=dict(entrez or {}), tasks=list(tasks), all_tasks=list(all_tasks),
        name_maps=task_name_maps(model), internal_with_gpr=internal,
    )


def consistent_reactions(parent: Parent, cache_path=None) -> set:
    """FASTCC-consistent reaction set of the parent on its medium (cobra ``fastcc``; cached)."""
    import json
    from pathlib import Path

    from cobra.flux_analysis import fastcc

    if parent.consistent is not None:
        return parent.consistent
    if cache_path is not None and Path(cache_path).exists():
        d = json.loads(Path(cache_path).read_text())
        assert d["parent"] == parent.name and d["medium"] == parent.medium_name and d["n_reactions"] == len(parent.rxn_ids)
        parent.consistent = set(d["consistent"])
        return parent.consistent
    t0 = time.process_time()
    cm = fastcc(parent.model, flux_threshold=FASTCORE_EPS)
    parent.consistent = {r.id for r in cm.reactions}
    del cm
    if cache_path is not None:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        Path(cache_path).write_text(json.dumps({
            "parent": parent.name, "medium": parent.medium_name, "n_reactions": len(parent.rxn_ids),
            "solver": parent.solver, "cpu_s": time.process_time() - t0, "consistent": sorted(parent.consistent)}))
    return parent.consistent


# ============================================================================= LP helpers
def _status_retry(model) -> str:
    """slim_optimize; on a non-optimal status (numerical or time limit) retry once (§11.4): Gurobi with
    NumericFocus = 2, GLPK from the standard (slack) basis. A second failure is returned as is."""
    model.slim_optimize()
    st = model.solver.status
    if st in (OPTIMAL, INFEASIBLE):
        return st
    if solver_name(model) == "gurobi":
        p = model.solver.problem.Params
        old = p.NumericFocus
        p.NumericFocus = 2
        model.slim_optimize()
        st = model.solver.status
        p.NumericFocus = old
    elif solver_name(model) == "glpk":
        import swiglpk

        swiglpk.glp_std_basis(model.solver.problem)
        model.slim_optimize()
        st = model.solver.status
    return st


def _fluxes(model, rxn_ids) -> dict:
    pv = model.solver.primal_values
    out = {}
    for rid in rxn_ids:
        r = model.reactions.get_by_id(rid)
        out[rid] = pv[r.forward_variable.name] - pv[r.reverse_variable.name]
    return out


def _set_objective(model, coefs: dict, direction: str) -> None:
    """Install a fresh sparse objective at the *solver* level (not tracked by cobra contexts).

    Performance rule for this module: never change the objective through cobra inside a ``with
    model:`` block while an attached objective is current; cobra's reset re-assigns the old object,
    which rebuilds its sympy expression from all 25k columns (~0.1 s per reset on Human-GEM). Every
    auxiliary LP installs a fresh objective here and restores the biomass objective with
    ``_biomass_objective`` in a ``finally``.
    """
    obj = model.problem.Objective(Zero, direction=direction, sloppy=True)
    model.solver.objective = obj
    if coefs:
        obj.set_linear_coefficients(coefs)


def _biomass_objective(model) -> None:
    r = model.reactions.get_by_id(BIOMASS_REACTION)
    _set_objective(model, {r.forward_variable: 1.0, r.reverse_variable: -1.0}, "max")


def _l1_objective(model, weights: dict) -> None:
    coefs = {}
    for rid, w in weights.items():
        r = model.reactions.get_by_id(rid)
        coefs[r.forward_variable] = w
        coefs[r.reverse_variable] = w
    _set_objective(model, coefs, "min")


def _apply_bounds(model, bounds: dict) -> None:
    get = model.reactions.get_by_id
    for rid, b in bounds.items():
        get(rid).bounds = b


# ============================================================================= families
def f0_closed(parent: Parent, off: set) -> set:
    return {rid for rid, n in parent.gprs.items() if not gpr_true(n, off)}


def gimme_closed(parent: Parent, s_r: pd.Series, fraction: float = GIMME_FRACTION,
                 noise_seed: int | None = None, noise_scale: float = 1e-3) -> tuple[set, dict]:
    """GIMME closed set. ``noise_seed`` adds ``noise_scale * U(0,1)`` to every reaction's cost
    (alternate-optima sampling for the Stage-A A7 stability study)."""
    m = parent.model
    cost = {rid: -v for rid, v in s_r.items() if not math.isnan(v) and v < 0}
    if noise_seed is not None:
        rng = np.random.default_rng(noise_seed)
        u = rng.random(len(parent.rxn_ids))
        cost = {rid: cost.get(rid, 0.0) + noise_scale * u[i] for i, rid in enumerate(parent.rxn_ids)}
    with m:
        m.reactions.get_by_id(BIOMASS_REACTION).lower_bound = fraction * parent.wt_generic
        try:
            _l1_objective(m, cost)
            st = _status_retry(m)
            if st != OPTIMAL:
                raise BuildFailed(f"GIMME LP {st}")
            candidates = [rid for rid, v in s_r.items() if not math.isnan(v) and v < 0]
            fl = _fluxes(m, candidates)
            obj = m.solver.objective.value
        finally:
            _biomass_objective(m)
    closed = {rid for rid in candidates if abs(fl[rid]) <= FLUX_TOL}
    return closed, {"gimme_objective": obj, "gimme_low_used": len(candidates) - len(closed)}


def fastcore_support(parent: Parent, core: set, eps: float = FASTCORE_EPS) -> tuple[set, dict]:
    """FASTCORE (Vlassis 2014, Algorithm 1) on the fastcc-consistent subnetwork. Returns support A.

    LP7 uses cobra's fastcc auxiliary variables (z in [0, eps], z <= sigma*v toggled per call);
    LP9 minimises sum |v_p| over the non-core reactions P not yet in A, with sigma*v_k >= 0.99 eps on K.
    Deviation from the MATLAB original (logged, not fatal): a core reaction that cannot be activated
    even as a flipped singleton is dropped from J ('unreachable_core') instead of aborting.
    """
    from cobra.flux_analysis.fastcc import LARGE_VALUE as LARGE
    from cobra.flux_analysis.fastcc import _add_lp7_vars, _flip_coefficients
    m = parent.model
    consistent = parent.consistent
    assert consistent is not None, "call consistent_reactions() first"
    N = [rid for rid in parent.rxn_ids if rid in consistent]
    irrev = {rid for rid in N if not m.reactions.get_by_id(rid).reversibility}
    C = set(core) & consistent
    stats = {"n_core": len(C), "n_lp": 0, "unreachable_core": 0, "singleton_iters": 0}
    inconsistent = [rid for rid in parent.rxn_ids if rid not in consistent]
    with m:
        _apply_bounds(m, {rid: (0.0, 0.0) for rid in inconsistent})
        _add_lp7_vars(m, [m.reactions.get_by_id(r) for r in sorted(C)], eps)
        cons = {r: m.constraints.get(f"aux_constraint_{r}") for r in C}
        auxv = {r: m.variables.get(f"auxiliary_{r}") for r in C}
        sign = {r: 1 for r in C}

        def lp7(J):
            for r in J:
                cons[r].lb = 0.0
            _set_objective(m, {auxv[r]: 1.0 for r in J}, "max")
            m.slim_optimize()
            stats["n_lp"] += 1
            st = m.solver.status
            fl = _fluxes(m, J) if st == OPTIMAL else None  # read before the model is modified again
            for r in J:
                cons[r].lb = -LARGE
            if st != OPTIMAL:
                raise BuildFailed(f"FASTCORE LP7 {st}")
            return {r for r in J if sign[r] * fl[r] >= 0.99 * eps}

        def lp9(K, P):
            with m:
                for r in K:
                    rx = m.reactions.get_by_id(r)
                    if sign[r] > 0:
                        rx.lower_bound = max(rx.lower_bound, 0.99 * eps)
                    else:
                        rx.upper_bound = min(rx.upper_bound, -0.99 * eps)
                _l1_objective(m, {p: 1.0 for p in P})
                st = _status_retry(m)
                stats["n_lp"] += 1
                if st != OPTIMAL:
                    raise BuildFailed(f"FASTCORE LP9 {st}")
                fl = _fluxes(m, N)
            return {r for r in N if abs(fl[r]) > FLUX_TOL}

        def find_sparse_mode(J, P, singleton):
            if not J:
                return set()
            Jl = sorted(J)
            K = lp7([Jl[0]] if singleton else Jl)
            if not K:
                return set()
            return lp9(K, P)

        P = set(N) - C
        A = find_sparse_mode(C & irrev, P, False)
        J = C - A
        flipped = singleton = False
        while J:
            P = P - A
            supp = find_sparse_mode(J, P, singleton)
            A |= supp
            if J & A:
                J = J - A
                flipped = False
            else:
                Jl = sorted(J)
                jirev = ({Jl[0]} if singleton else set(J)) - irrev
                if flipped or not jirev:
                    if singleton:
                        J = J - {Jl[0]}  # see docstring: logged instead of aborting
                        stats["unreachable_core"] += 1
                        flipped = False
                    else:
                        flipped, singleton = False, True
                else:
                    _flip_coefficients(m, jirev)
                    for r in jirev:
                        sign[r] = -sign[r]
                    flipped = True
                if singleton:
                    stats["singleton_iters"] += 1
        _biomass_objective(m)
    return A, stats


def corda_confidence_genes(s_g: pd.Series) -> dict:
    out = {}
    for g, v in s_g.items():
        if pd.isna(v):
            continue
        out[g] = 3 if v >= 1 else 2 if v >= 0 else 1 if v >= -1 else -1
    return out


def _corda_numpy_shim():
    """corda 0.5.1 calls np.in1d, removed in NumPy 2.4; give the corda module a numpy proxy with it."""
    import corda.corda as cc

    if getattr(cc.np, "_nextgem_shim", False):
        return
    proxy = types.SimpleNamespace(**{k: getattr(np, k) for k in dir(np) if not k.startswith("__")})
    proxy.in1d = lambda a, b, assume_unique=False, invert=False: np.isin(
        np.ravel(a), b, assume_unique=assume_unique, invert=invert)
    proxy._nextgem_shim = True
    cc.np = proxy


def corda_closed(parent: Parent, s_g: pd.Series) -> tuple[set, dict]:
    from corda import CORDA
    from corda.util import reaction_confidence

    _corda_numpy_shim()
    m = parent.model
    conf_g = corda_confidence_genes(s_g)
    conf = {r.id: int(reaction_confidence(r, conf_g)) if r.genes else 0 for r in m.reactions}
    t0 = time.process_time()
    c = CORDA(m, conf, **CORDA_PARAMS)
    set_solver_params(c.model)
    t_setup = time.process_time() - t0
    c.build()
    inc = c.included
    closed = {rid for rid in parent.rxn_ids if not inc.get(rid, False)}
    counts = pd.Series(list(conf.values())).value_counts().to_dict()
    del c
    return closed, {"corda_setup_cpu_s": t_setup, **{f"corda_conf_{k}": int(v) for k, v in counts.items()}}


def eflux_bounds(parent: Parent, s_g: pd.Series, lost: set = frozenset()) -> tuple[dict, dict, float]:
    """E-Flux bounds for internal GPR reactions. Returns (bounds, linear values, cap_max)."""
    lin = {g: (2.0 ** v if v >= 0 else 0.0) for g, v in s_g.items() if not pd.isna(v)}
    for g in lost:
        lin[g] = 0.0
    caps = {rid: capacity(parent.gprs[rid], lin) for rid in parent.internal_with_gpr}
    finite = [c for c in caps.values() if not math.isnan(c)]
    cap_max = max(finite) if finite else 1.0
    bounds = {}
    for rid, c in caps.items():
        if math.isnan(c):
            continue
        bounds[rid] = _eflux_bound(parent.model.reactions.get_by_id(rid), c, cap_max)
    return bounds, lin, cap_max


def _eflux_bound(rxn, cap, cap_max):
    b = max(EFLUX_FLOOR, EFLUX_SCALE * cap / cap_max)
    lb, ub = rxn.lower_bound, rxn.upper_bound  # parent bounds (called outside any line context)
    return (max(lb, -b) if lb < 0 else lb, min(ub, b) if ub > 0 else ub)


# ============================================================================= genomic constraint
def loss_closure(parent: Parent, genes: set) -> set:
    """Reactions whose (full) GPR is False when ``genes`` are non-functional."""
    if not genes:
        return set()
    cand = set().union(*(parent.gene_rxns.get(g, []) for g in genes))
    return {rid for rid in cand if rid in parent.gprs and not gpr_true(parent.gprs[rid], genes)}


def _growth_with_closed(parent: Parent, closed: set) -> float:
    m = parent.model
    with m:
        _apply_bounds(m, {rid: (0.0, 0.0) for rid in closed})
        st = _status_retry(m)
        if st == INFEASIBLE:
            return 0.0
        if st != OPTIMAL:
            raise BuildFailed(f"loss-rule growth LP {st}")
        return float(m.solver.objective.value)


def accept_losses(parent: Parent, lost, expr: pd.Series | None) -> dict:
    """Order-independent drop (DESIGN W4 genomic constraint), on the parent with the line's medium.

    1. each loss alone; 2. drop those that alone kill WT growth (model-contradicted); 3. apply the
    rest jointly; 4. while the joint set kills growth, drop the loss with the highest expression
    (ties: ascending Entrez ID; missing expression ranks lowest), re-test. Input order is irrelevant:
    the set is canonicalised first.
    """
    lost = sorted(set(lost) & set(parent.genes))
    floor = LETHAL_FRACTION * parent.wt_generic
    contradicted = []
    for g in lost:
        cl = loss_closure(parent, {g})
        if cl and _growth_with_closed(parent, cl) < floor:
            contradicted.append(g)
    remaining = [g for g in lost if g not in set(contradicted)]

    def key(g):
        e = expr.get(g, np.nan) if expr is not None else np.nan
        ez = parent.entrez.get(g)
        ez = int(ez) if ez is not None and str(ez).isdigit() else 10**12
        return (-(e if not pd.isna(e) else -np.inf), ez, g)

    order = sorted(remaining, key=key)
    dropped = []
    growth = _growth_with_closed(parent, loss_closure(parent, set(remaining)))
    while remaining and growth < floor:
        g = order.pop(0)
        dropped.append(g)
        remaining = [x for x in remaining if x != g]
        growth = _growth_with_closed(parent, loss_closure(parent, set(remaining)))
    return {"lost": lost, "contradicted": contradicted, "dropped_joint": dropped,
            "accepted": sorted(remaining), "joint_growth": growth}


# ============================================================================= repair rule
def _apply_task(m, task, name_maps):
    """Close boundaries, apply a raven task; returns the metabolite-constraint snapshot to restore."""
    from raven_toolbox.tasks.check import _metabolite_bounds, apply_task_constraints

    bounds, missing = _metabolite_bounds(task, *name_maps)
    if missing:
        return None, f"unknown metabolites {missing}"
    saved = {mid: (m.constraints[mid].lb, m.constraints[mid].ub) for mid in bounds}
    for r in m.boundary:
        r.bounds = (0.0, 0.0)
    _, err = apply_task_constraints(m, task, *name_maps)
    return saved, err


def _restore_mets(m, saved):
    from raven_toolbox.tasks.check import _set_constraint_bounds

    for mid, (lb, ub) in (saved or {}).items():
        _set_constraint_bounds(m.constraints[mid], lb, ub)


def repair(parent: Parent, fixed_bounds: dict, repairable: set, weights: dict) -> dict:
    """The fixed repair rule. ``fixed_bounds``: bounds of the line model that are NOT repairable
    (accepted-loss closures, E-Flux bounds); ``repairable``: closed reactions the LPs may use."""
    m = parent.model
    w = {rid: weights.get(rid, 1.0) for rid in repairable}
    out = {"growth_lp": None, "task_lps": {}, "reopen": set()}
    if not repairable:
        out["growth_lp"] = "nothing_to_repair"
        return out
    with m:
        _apply_bounds(m, fixed_bounds)
        with m:
            m.reactions.get_by_id(BIOMASS_REACTION).lower_bound = REPAIR_GROWTH_FRACTION * parent.wt_generic
            _l1_objective(m, w)
            st = _status_retry(m)
            out["growth_lp"] = st
            if st == OPTIMAL:
                fl = _fluxes(m, list(repairable))
                got = {r for r, v in fl.items() if abs(v) > FLUX_TOL}
                out["reopen"] |= got
                out["n_growth"] = len(got)
        _set_objective(m, {}, "max")  # cached zero objective: cheap cobra resets inside the task loop
        try:
            _repair_tasks(parent, m, w, out)
        finally:
            _biomass_objective(m)
    return out


def _repair_tasks(parent, m, w, out):
    for task in parent.tasks:
        saved = None
        try:
            with m:
                saved, err = _apply_task(m, task, parent.name_maps)
                if err:
                    out["task_lps"][task.description] = f"error: {err}"
                    continue
                _l1_objective(m, w)
                st = _status_retry(m)
                if st == OPTIMAL:
                    fl = _fluxes(m, list(w))
                    got = {r for r, v in fl.items() if abs(v) > FLUX_TOL}
                    out["reopen"] |= got
                    out["task_lps"][task.description] = len(got)
                else:
                    out["task_lps"][task.description] = st
        finally:
            _restore_mets(m, saved)


# ============================================================================= line model
@dataclass
class LineModel:
    line: str
    family: str
    closed: set  # closed by extraction (incl. accepted-loss closures if genomic)
    lost_closed: set = field(default_factory=set)
    reopened: set = field(default_factory=set)
    eflux: dict = field(default_factory=dict)  # rid -> (lb, ub) for F4
    off: set = field(default_factory=set)
    lost: set = field(default_factory=set)
    lin: dict = field(default_factory=dict)
    cap_max: float = 1.0
    info: dict = field(default_factory=dict)

    @property
    def final_closed(self) -> set:
        return (self.closed - self.reopened) | self.lost_closed

    def bounds(self) -> dict:
        b = dict(self.eflux)
        for rid in self.final_closed:
            b[rid] = (0.0, 0.0)
        return b


def extract(parent: Parent, family: str, s_g: pd.Series, *, off_override: set | None = None,
            gimme_noise_seed: int | None = None) -> LineModel:
    """Run one family's extraction on a line's signed gene scores (no genomic constraint, no repair)."""
    t0 = time.process_time()
    s_r = reaction_scores(parent.gprs, s_g)
    off = set(s_g.index[(s_g < 0).to_numpy()]) if off_override is None else set(off_override)
    lm = LineModel(line="", family=family, closed=set(), off=off)
    info = {}
    if family in ("F0", "F5"):
        lm.closed = f0_closed(parent, off)
    elif family == "F2":
        lm.closed, info = gimme_closed(parent, s_r, noise_seed=gimme_noise_seed)
    elif family == "F3":
        consistent_reactions(parent)
        core = {rid for rid, v in s_r.items() if not math.isnan(v) and v >= 0}
        A, info = fastcore_support(parent, core)
        lm.closed = set(parent.rxn_ids) - A
    elif family == "F3b":
        lm.closed, info = corda_closed(parent, s_g)
    elif family == "F4":
        lm.eflux, lm.lin, lm.cap_max = eflux_bounds(parent, s_g)
    else:
        raise ValueError(family)
    lm.info = {"extract_cpu_s": time.process_time() - t0, "n_closed_extract": len(lm.closed), **info}
    lm.info["_s_r"] = s_r
    return lm


def finish(parent: Parent, lm: LineModel, accepted_losses: set | None = None) -> LineModel:
    """Genomic constraint (optional) + repair (presence families) + final WT growth check."""
    t0 = time.process_time()
    lm.lost = set(accepted_losses or ())
    lm.lost_closed = loss_closure(parent, lm.lost)
    if lm.family == "F4" and lm.lost:
        # E-Flux: lost genes contribute 0 capacity; recompute bounds, closures stay closed
        lin = dict(lm.lin)
        lin.update({g: 0.0 for g in lm.lost})
        lm.lin = lin
        for rid in parent.internal_with_gpr:
            c = capacity(parent.gprs[rid], lin)
            if not math.isnan(c):
                lm.eflux[rid] = _eflux_bound(parent.model.reactions.get_by_id(rid), c, lm.cap_max)
    rep = {"growth_lp": "not_applicable", "task_lps": {}, "reopen": set()}
    if lm.family in REPAIRED_FAMILIES:
        s_r = lm.info["_s_r"]
        repairable = lm.closed - lm.lost_closed
        weights = {}
        for rid in repairable:
            v = s_r.get(rid, math.nan)
            weights[rid] = 1.0 if math.isnan(v) else 1.0 + max(0.0, -v)
        fixed = {rid: (0.0, 0.0) for rid in lm.lost_closed}
        rep = repair(parent, fixed, repairable, weights)
        lm.reopened = set(rep["reopen"])
    m = parent.model
    with m:
        _apply_bounds(m, lm.bounds())
        st = _status_retry(m)
        wt = float(m.solver.objective.value) if st == OPTIMAL else 0.0
    n_rxn = len(parent.rxn_ids)
    lm.info.update({
        "repair_cpu_s": time.process_time() - t0, "repair_growth_lp": rep["growth_lp"],
        "repair_task_lps": rep["task_lps"], "n_reopened": len(lm.reopened),
        "repair_frac_reactions": len(lm.reopened) / n_rxn, "n_lost_accepted": len(lm.lost),
        "n_lost_closed": len(lm.lost_closed), "n_final_closed": len(lm.final_closed),
        "n_open": n_rxn - len(lm.final_closed), "wt_status": st, "wt_growth": wt,
        "wt_ratio_generic": wt / parent.wt_generic,
        # failed build (§11.4): repaired families must reach the repair floor (1% of generic WT); F4 has no
        # repair rule and scaled bounds (growth in its own units), so for it "no growth" means <= 1e-6
        "failed": (st != OPTIMAL) or (wt < REPAIR_GROWTH_FRACTION * parent.wt_generic * (1 - 1e-6)
                                      if lm.family in REPAIRED_FAMILIES else wt <= 1e-6),
    })
    return lm


# ============================================================================= KO sweep
def closure_map(parent: Parent, lm: LineModel) -> tuple[dict, dict]:
    """gene -> reactions closed by its KO in the line model; gene -> {rid: new bounds} (E-Flux)."""
    final_closed = lm.final_closed
    close, rebound = {}, {}
    presence = lm.family in PRESENCE_FAMILIES
    base_false = set(lm.lost)
    pruned_false = base_false | lm.off
    for rid, node in parent.gprs.items():
        if rid in final_closed:
            continue
        F = pruned_false if (presence and gpr_true(node, pruned_false)) else base_false
        for g in gpr_genes(node):
            if g in F:
                continue
            if not gpr_true(node, F | {g}):
                close.setdefault(g, set()).add(rid)
            elif lm.family == "F4" and rid in lm.eflux and lm.lin.get(g, 0.0) > 0:
                lin = dict(lm.lin)
                lin[g] = 0.0
                c = capacity(node, {k: lin.get(k, math.nan) for k in gpr_genes(node)})
                if not math.isnan(c):
                    nb = _eflux_bound(parent.model.reactions.get_by_id(rid), c, lm.cap_max)
                    if nb != lm.eflux[rid]:
                        rebound.setdefault(g, {})[rid] = nb
    return close, rebound


def _pfba_fluxes(parent: Parent, m, wt: float) -> dict:
    """Minimal-L1 flux at WT growth (biomass >= wt (1 - 1e-9)); falls back to the FBA flux. Any optimal
    flux keeps the zero-flux shortcut exact; minimising the support only reduces the number of solves."""
    with m:
        m.reactions.get_by_id(BIOMASS_REACTION).lower_bound = wt * (1 - 1e-9)
        try:
            _l1_objective(m, {rid: 1.0 for rid in parent.rxn_ids})
            st = _status_retry(m)
            if st == OPTIMAL:
                return _fluxes(m, parent.rxn_ids)
        finally:
            _biomass_objective(m)
    m.slim_optimize()
    return _fluxes(m, parent.rxn_ids)


def ko_sweep(parent: Parent, lm: LineModel, genes=None) -> pd.DataFrame:
    """Single-gene KO ratios of a line model (zero-flux shortcut is exact: the WT pFBA solution stays
    feasible when every affected reaction carries zero flux / stays within its new bound)."""
    genes = parent.genes if genes is None else list(genes)
    m = parent.model
    t0 = time.process_time()
    close, rebound = closure_map(parent, lm)
    rows = []
    with m:
        _apply_bounds(m, lm.bounds())
        st = _status_retry(m)
        if st != OPTIMAL:
            raise BuildFailed(f"sweep WT {st}")
        wt = float(m.solver.objective.value)
        fl = _pfba_fluxes(parent, m, wt)
        n_solved = 0
        for g in genes:
            cl = close.get(g, set())
            rb = rebound.get(g, {})
            if not cl and not rb:
                rows.append((g, wt, "no_effect"))
                continue
            if all(abs(fl[r]) <= FLUX_TOL for r in cl) and all(
                    b[0] - FLUX_TOL <= fl[r] <= b[1] + FLUX_TOL for r, b in rb.items()):
                rows.append((g, wt, "zero_flux"))
                continue
            with m:
                _apply_bounds(m, {**rb, **{r: (0.0, 0.0) for r in cl}})
                s = _status_retry(m)
                n_solved += 1
                if s == OPTIMAL:
                    rows.append((g, float(m.solver.objective.value), s))
                elif s == INFEASIBLE:
                    rows.append((g, 0.0, s))
                else:
                    rows.append((g, math.nan, s))
    df = pd.DataFrame(rows, columns=["gene", "ko_growth", "status"])
    df["wt_growth"] = wt
    df["ratio"] = df.ko_growth / wt
    df["ratio_r"] = rounded_ratio(df.ratio)
    df["solver"] = solver_name(m)
    lm.info.update({"sweep_cpu_s": time.process_time() - t0, "sweep_n_solved": n_solved,
                    "sweep_n_nonoptimal": int(df.ratio.isna().sum()), "sweep_n_lethal": int((df.ratio_r < 0.01).sum())})
    return df


def task_profile(parent: Parent, lm: LineModel) -> dict:
    """Pass/fail of every task of the essential list (57) on the line model (closed boundaries)."""
    m = parent.model
    out = {}
    with m:
        _apply_bounds(m, lm.bounds())
        _set_objective(m, {}, "max")
        try:
            _profile_tasks(parent, m, out)
        finally:
            _biomass_objective(m)
    return out


def _profile_tasks(parent, m, out):
    for i, task in enumerate(parent.all_tasks):
        saved = None
        try:
            with m:
                saved, err = _apply_task(m, task, parent.name_maps)
                if err:
                    out[f"{i:02d}_{task.description}"] = False
                    continue
                m.slim_optimize()
                out[f"{i:02d}_{task.description}"] = (m.solver.status == OPTIMAL) != task.should_fail
        finally:
            _restore_mets(m, saved)


def summary_row(lm: LineModel) -> dict:
    keep = {k: v for k, v in lm.info.items() if not k.startswith("_") and not isinstance(v, (dict, set, list))}
    return {"line": lm.line, "family": lm.family, **keep}
