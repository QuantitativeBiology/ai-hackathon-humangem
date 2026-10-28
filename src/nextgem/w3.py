"""W3 — objectives, metabolic tasks and biomass (DESIGN.md §9 W3; amendments A2, A8, A9). Label-free library.

Nothing here opens a label file. Training-tier scoring lives in ``scripts/w3_score_train.py`` and reuses
``nextgem.w1`` / ``nextgem.metrics`` for every metric.

Contents
--------
1. Task lists: Human-GEM ``metabolicTasks_*.txt`` with the 11 functional groups of the Essential list (``#`` header
   rows), the CellFie consensus list and the Full list.
2. Task capacity engine (the W3 "capacity matrix", compute SF3). For every task t and every *distinct* gene
   knockout (genes whose KO disables the identical reaction set share one test, as in Human-GEM
   ``taskEssentialGenes.py``):
     - ``feasible``: the task as written is feasible after the KO (each task with its own inputs, every other
       boundary closed; identical to the P57 call rule of DESIGN §6.8);
     - ``capacity``: λ_KO / λ_WT, where λ = max uniform scaling of the task's *required* items (outputs with
       LB > 0 and task equations with LB > 0) subject to the task's other constraints, capped at LAMBDA_MAX
       (*ours*: DESIGN says "max output capacity" without a formula; λ ≥ 1 ⟺ the required outputs can be met).
   Exact shortcut (as taskEssentialGenes): a KO is re-tested only if it disables a reaction carrying flux in a
   feasible WT solution of that LP (pFBA solution for feasibility, the λ-optimum for capacity).
3. Biomass-objective KO engine for the W3 variants (D-BS/D-ER demands, O-NGAM_k, B-loo_j, B-add_c, T10 biomass),
   with the same exact shortcuts and status policy as ``nextgem.w2_engine.ko_task``.
4. Pseudo-task generator (negative control): required outputs swapped for random producible metabolites of the
   same compartment, SEEDS[0:5].
5. CellFie-inspired task activity (T5) from expression (local min–max–mean thresholds, 25th/75th percentiles).
6. Mode score assembly (pure functions over gene tables).
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from nextgem import paths

BIOMASS = "MAR13082"
COFACTOR_POOL = "MAR10065"
NGAM_RXN = "MAR03964"
BIOMASS_EXCLUDED = ("MAM01371c", "MAM02040c")  # ATP, H2O: never left out (DESIGN W3 B-loo)
LAMBDA_MAX = 1000.0
FEAS_TOL = 1e-8
ZERO_FLUX_TOL = 1e-9
RATIO_DECIMALS = 6
BADD_COEF = 1e-4  # mmol/gDW per unit biomass (DESIGN W3 B-add); sensitivity 1e-3
BADD_COEF_SENS = 1e-3
DEMAND_FRACTION = 0.01  # D-BS / D-ER: lb = 1% of the max under the same medium with biomass free
NGAM_LEVELS = (0.0, 0.1, 0.3)  # O-NGAM_k: fraction of the maximum MAR03964 flux on the medium

# Cofactor-pool groups (DESIGN W3 B-loo; final list §5). Matched on metabolite names of MAR10065 species.
# Retinoid and vitamin-D catabolites are kept as two groups (ours; see README "ambiguities").
POOL_GROUPS = {
    "retinoid_catabolites": re.compile(r"retin|vitamin a", re.I),
    "vitaminD_catabolites": re.compile(r"vitamin d", re.I),
    "tocopherol_catabolites": re.compile(r"chromanol|cehc|tocopher", re.I),
    "carnitine": re.compile(r"^l-carnitine$", re.I),
    "heme": re.compile(r"^heme$", re.I),
    "ubiquinol": re.compile(r"^ubiquinol$", re.I),
    "lipoate": re.compile(r"^lipoic acid$|^lipoate$", re.I),
    "cobamide": re.compile(r"cobamide|cobalamin", re.I),
}

# Single-metabolite B-loo components of other pools, added from label-free structural findings (not in the frozen
# DESIGN text): (S)-dihydroorotate in the small-metabolite pool MAR10064 makes CAD lethal even with uridine
# (notes/reviews/M3_W5_review.md §4, coordinator 2026-10-27; added before any W3 scoring).
EXTRA_BLOO = {"metpool:dihydroorotate": ("MAR10064", "(S)-dihydroorotate", "c")}

# B-add components (final list §5, verbatim names) -> how they enter MAR13082 (ours: metabolite choice per
# parent is resolved by name/compartment in scripts/w3_structure.py and frozen as candidate rows).
#   kind "sink": the metabolite is a biomass reactant; "turnover": reactants -> products as written.
BADD_SPECS = {
    "hemeA": {"kind": "sink", "names": [("heme a", "m"), ("heme a", "c")]},
    "dolichol_P": {"kind": "sink", "names": [("dolichyl-phosphate", "r")]},
    "GPI_anchor_precursor": {"kind": "sink",
                             "names": [("Mannosyl-3-(Phosphoethanolaminyl-Mannosyl)-Glucosaminyl-Acylphosphatidylinositol (M4A))", "r")]},
    "FeS_cluster": {"kind": "sink", "names": [("4Fe4S iron-sulfur cluster", "m")]},
    "SAM_sink": {"kind": "sink", "names": [("SAM", "c")]},
    "spermidine": {"kind": "sink", "names": [("spermidine", "c")]},
    "NAD_turnover": {"kind": "turnover", "reactants": [("NAD+", "c"), ("H2O", "c")],
                     "products": [("nicotinamide", "c"), ("ADP-ribose", "c"), ("H+", "c")]},
    # a turnover sink, not static presence: GSH (like NAD+) is in the metabolite pool only statically (final list §5)
    "GSH_turnover": {"kind": "turnover", "reactants": [("GSH", "c")], "products": []},
}

# D-ER demands (ours: "ER as ATP and NAD(P)H demand"; cytosolic)
DER_SPECS = {
    "ATP_hydrolysis": ([("ATP", "c"), ("H2O", "c")], [("ADP", "c"), ("Pi", "c"), ("H+", "c")]),
    "NADH_oxidation": ([("NADH", "c"), ("H+", "c")], [("NAD+", "c")]),
    "NADPH_oxidation": ([("NADPH", "c"), ("H+", "c")], [("NADP+", "c")]),
}


class SolverInfraError(RuntimeError):
    """Solver / licence / persistent numerical failure: abort the job with no output (DESIGN §11.4)."""


# ============================================================================= 1. task lists
TASK_FILES = {
    "essential": "metabolicTasks_Essential.txt",
    "cellfie": "metabolicTasks_CellfieConsensus.txt",
    "full": "metabolicTasks_Full.txt",
}


def task_file(parent: str, which: str) -> Path:
    """The release's own task file; CellFie exists only in Human2 (v2.1.1) and is used for both parents (ours)."""
    from nextgem.model import GEMS

    f = GEMS[parent].root / "data" / "metabolicTasks" / TASK_FILES[which]
    if not f.exists() and which == "cellfie":
        f = GEMS["human2"].root / "data" / "metabolicTasks" / TASK_FILES[which]
    return f


def essential_groups(path: Path) -> list[str]:
    """Functional group of each task of an Essential-style file (the ``#`` header row above it), in task order."""
    groups, cur = [], None
    for ln in Path(path).read_text().splitlines()[1:]:
        f = ln.split("\t")
        if f[0].strip() == "#":
            if len(f) > 1 and f[1].strip() and not f[1].startswith("ER: "):
                cur = f[1].strip()
            continue
        if len(f) > 1 and f[1].strip():
            groups.append(cur)
    return groups


def load_tasks(parent: str, which: str) -> pd.DataFrame:
    """One row per task: key (unique), idx (0-based file order), category/id, description, group, should_fail,
    task (raven Task object)."""
    from raven_toolbox.tasks.tasklist import parse_task_list

    f = task_file(parent, which)
    tasks = parse_task_list(str(f))
    groups = essential_groups(f) if which == "essential" else [None] * len(tasks)
    if len(groups) != len(tasks):
        raise RuntimeError(f"{f}: {len(groups)} group rows for {len(tasks)} tasks")
    rows = []
    for i, (t, g) in enumerate(zip(tasks, groups)):
        key = f"{which[:3].upper()}{i + 1:03d}"
        rows.append({"key": key, "idx": i, "list": which, "category": t.id, "description": t.description.strip(),
                     "group": g, "should_fail": bool(t.should_fail), "task": t})
    return pd.DataFrame(rows)


# ============================================================================= 2. task capacity engine
def gene_disabled_reactions(model) -> dict[str, frozenset]:
    """gene -> reactions whose GPR is False when only that gene is off (Human-GEM taskEssentialGenes logic)."""
    out: dict[str, set] = defaultdict(set)
    for rxn in model.reactions:
        for g in rxn.genes:
            if not rxn.gpr.eval([g.id]):
                out[g.id].add(rxn.id)
    return {g: frozenset(r) for g, r in out.items()}


def ko_groups(gd: dict[str, frozenset]) -> dict[frozenset, list[str]]:
    groups: dict[frozenset, list[str]] = defaultdict(list)
    for g, rx in gd.items():
        if rx:
            groups[rx].append(g)
    return dict(groups)


def task_base(model):
    """Copy with every boundary closed (as raven check_tasks / taskEssentialGenes)."""
    base = model.copy()
    for r in base.boundary:
        r.bounds = (0.0, 0.0)
    return base


LP_TIME_LIMIT = 300  # s per LP (as A10's LP limit); GLPK occasionally stalls on a degenerate LP


def solve_status(model) -> str:
    """slim_optimize; 'optimal' / 'infeasible'; anything else (incl. the 300 s time limit) retried once with presolve
    (GLPK: a fresh solve from scratch) / NumericFocus = 2 (Gurobi), else SolverInfraError."""
    cfg = model.solver.configuration
    if getattr(cfg, "timeout", None) != LP_TIME_LIMIT:
        cfg.timeout = LP_TIME_LIMIT
    model.slim_optimize()
    st = model.solver.status
    if st in ("optimal", "infeasible"):
        return st
    try:  # GLPK: retry once with presolve (no NumericFocus); Gurobi: NumericFocus = 2
        if "gurobi" in model.solver.interface.__name__:
            p = model.solver.problem.Params
            old = p.NumericFocus
            p.NumericFocus = 2
            try:
                model.slim_optimize()
            finally:
                p.NumericFocus = old
        else:
            cfg = model.solver.configuration
            old = cfg.presolve
            cfg.presolve = True
            try:
                model.slim_optimize()
            finally:
                cfg.presolve = old
    except Exception as e:  # noqa: BLE001
        raise SolverInfraError(f"retry failed: {e}") from e
    st2 = model.solver.status
    if st2 in ("optimal", "infeasible"):
        return st2
    raise SolverInfraError(f"non-optimal LP status {st!r} then {st2!r}")


def _support(model, original_ids: set[str]) -> set[str]:
    pv = model.solver.primal_values
    out = set()
    for r in model.reactions:
        if r.id in original_ids and abs(pv[r.forward_variable.name]) + abs(pv[r.reverse_variable.name]) > FEAS_TOL:
            out.add(r.id)
    return out


class TaskContext:
    """Applies one task to ``base`` (inside the caller's ``with base:``) and offers feasibility and λ solves.

    Metabolite-constraint bounds are not tracked by cobra contexts, so ``close()`` restores them explicitly
    (as Human-GEM taskEssentialGenes does).
    """

    def __init__(self, base, task, maps, saved):
        from raven_toolbox.tasks.check import apply_task_constraints

        self.base, self.saved = base, saved
        self.task_mets, self.error = apply_task_constraints(base, task, *maps)
        self.required: list[tuple[str, str, float]] = []  # (kind 'met'|'rxn', id, LB)
        if self.error is None:
            for mid in self.task_mets:
                c = base.constraints[mid]
                if c.lb is not None and c.lb > 0:
                    self.required.append(("met", mid, float(c.lb)))
            for r in base.reactions:
                if r.id.startswith("TASK_TMP_") and r.lower_bound > 0:
                    self.required.append(("rxn", r.id, float(r.lower_bound)))
        self._lam = None

    def close(self):
        from raven_toolbox.tasks.check import _set_constraint_bounds

        for mid in (self.task_mets or self.saved):
            if mid in self.base.constraints:
                _set_constraint_bounds(self.base.constraints[mid], *self.saved[mid])

    def feasible(self) -> bool:
        return solve_status(self.base) == "optimal"

    def lambda_setup(self):
        """Inside a nested ``with base:``: relax required LBs to 0, add λ and its rows, maximise λ."""
        from raven_toolbox.tasks.check import _set_constraint_bounds

        m = self.base
        lam = m.problem.Variable("w3_lambda", lb=0.0, ub=LAMBDA_MAX)
        cons = []
        self._relaxed = []
        for kind, ident, lb in self.required:
            if kind == "met":
                c = m.constraints[ident]
                self._relaxed.append((ident, c.lb, c.ub))
                _set_constraint_bounds(c, 0.0, c.ub)
                cons.append(m.problem.Constraint(c.expression - lb * lam, lb=0.0, name=f"w3_req_{ident}"))
            else:
                r = m.reactions.get_by_id(ident)
                r.lower_bound = 0.0
                cons.append(m.problem.Constraint(r.flux_expression - lb * lam, lb=0.0, name=f"w3_req_{ident}"))
        m.add_cons_vars([lam, *cons])
        m.objective = m.problem.Objective(lam, direction="max")
        self._lam = lam

    def lambda_teardown(self):
        from raven_toolbox.tasks.check import _set_constraint_bounds

        for ident, lb, ub in self._relaxed:
            _set_constraint_bounds(self.base.constraints[ident], lb, ub)

    def lam_value(self) -> float:
        st = solve_status(self.base)
        return float(self._lam.primal) if st == "optimal" else 0.0


def task_capacity(base, task, maps, saved, ko_list: list[frozenset], original_ids: set[str],
                  essential_reactions: bool = False, capacity: bool = True) -> dict:
    """Feasibility and λ capacity of one task for every KO set in ``ko_list`` (exact flux-support shortcut).

    Returns dict(status, error, lam_wt, n_required, fset_size, lset_size, feasible{ko: bool}, lam{ko: float},
    n_solves, ess_rxns (optional: reactions whose single removal makes the task infeasible)).
    """
    from cobra.flux_analysis import pfba
    from cobra.exceptions import OptimizationError

    out = {"status": "ok", "error": None, "lam_wt": np.nan, "feasible": {}, "lam": {}, "n_solves": 0,
           "ess_rxns": None}
    with base:
        tc = TaskContext(base, task, maps, saved)
        try:
            if tc.error is not None:
                out.update(status="task_error", error=tc.error)
                return out
            out["n_required"] = len(tc.required)
            if not tc.feasible():
                out.update(status="wt_infeasible")
                return out
            sol = pfba(base)
            fset = {rid for rid, v in sol.fluxes.items() if rid in original_ids and abs(v) > FEAS_TOL}
            out["fset_size"] = len(fset)
            # feasibility per KO
            for ko in ko_list:
                if ko & fset:
                    with base:
                        for rid in ko:
                            base.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                        out["feasible"][ko] = tc.feasible()
                    out["n_solves"] += 1
                else:
                    out["feasible"][ko] = True
            if essential_reactions:
                ess = []
                for rid in sorted(fset):
                    with base:
                        base.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                        if not tc.feasible():
                            ess.append(rid)
                    out["n_solves"] += 1
                out["ess_rxns"] = ess
            # λ capacity
            if not capacity:
                out["lam_wt"] = np.nan
                out["lam"] = {ko: np.nan for ko in ko_list}
            elif tc.required:
                with base:
                    tc.lambda_setup()
                    try:
                        lam_wt = tc.lam_value()
                        if lam_wt <= 0:
                            raise SolverInfraError(f"λ_WT = {lam_wt} for a feasible task")
                        lset = _support(base, original_ids)
                        out["lset_size"] = len(lset)
                        for ko in ko_list:
                            if ko & lset:
                                with base:
                                    for rid in ko:
                                        base.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                                    out["lam"][ko] = tc.lam_value()
                                out["n_solves"] += 1
                            else:
                                out["lam"][ko] = lam_wt
                    finally:
                        tc.lambda_teardown()
                out["lam_wt"] = lam_wt
            else:
                out["lam_wt"] = 1.0
                out["lam"] = {ko: (1.0 if out["feasible"][ko] else 0.0) for ko in ko_list}
        except OptimizationError as e:
            raise SolverInfraError(f"pFBA failed on a feasible task: {e}") from e
        finally:
            tc.close()
    return out


def capacity_ratio(lam_ko: float, lam_wt: float) -> float:
    if not lam_wt > 0:
        return np.nan
    return float(np.clip(np.round(min(lam_ko, lam_wt) / lam_wt, RATIO_DECIMALS), 0.0, 1.0)) + 0.0


# ============================================================================= 3. biomass-objective KO engine
@dataclass
class Variant:
    """A biomass-objective W3 variant applied inside ``with model:``.

    kind: 'none' | 'remove' (biomass reactant ids) | 'pool_remove' (MAR10065 species ids) | 'add' (stoich dict
    added to MAR13082) | 'fix' (reaction id -> fixed flux) | 'demands' (list of (rxn_id, {met: coef}, lb)) |
    'combo' (list of Variants).
    """

    name: str
    kind: str = "none"
    params: dict = field(default_factory=dict)

    def apply(self, model) -> None:
        import cobra

        k, p = self.kind, self.params
        if k == "none":
            return
        if k == "remove":
            b = model.reactions.get_by_id(BIOMASS)
            b.subtract_metabolites({model.metabolites.get_by_id(m): b.get_coefficient(m) for m in p["mets"]})
        elif k == "pool_remove":
            r = model.reactions.get_by_id(p.get("pool", COFACTOR_POOL))
            r.subtract_metabolites({model.metabolites.get_by_id(m): r.get_coefficient(m) for m in p["mets"]})
        elif k == "add":
            b = model.reactions.get_by_id(BIOMASS)
            b.add_metabolites({model.metabolites.get_by_id(m): c for m, c in p["stoich"].items()}, combine=True)
        elif k == "fix":
            for rid, v in p["fluxes"].items():
                model.reactions.get_by_id(rid).bounds = (v, v)
        elif k == "demands":
            rx = []
            for rid, st, lb in p["demands"]:
                r = cobra.Reaction(rid, lower_bound=lb, upper_bound=1000.0)
                r.add_metabolites({model.metabolites.get_by_id(m): c for m, c in st.items()})
                rx.append(r)
            model.add_reactions(rx)
        elif k == "combo":
            for v in p["variants"]:
                v.apply(model)
        else:
            raise ValueError(k)

    def to_json(self) -> dict:
        if self.kind == "combo":
            return {"name": self.name, "kind": "combo", "params": {"variants": [v.to_json() for v in self.params["variants"]]}}
        return {"name": self.name, "kind": self.kind, "params": self.params}

    @staticmethod
    def from_json(d: dict) -> "Variant":
        if d["kind"] == "combo":
            return Variant(d["name"], "combo", {"variants": [Variant.from_json(x) for x in d["params"]["variants"]]})
        return Variant(d["name"], d["kind"], d["params"])


def gene_closures(model) -> dict[str, list[str]]:
    """gene -> reactions newly closed by ``Gene.knock_out()`` (same as nextgem.w2_engine.gene_closures)."""
    from nextgem.w2_engine import gene_closures as _gc

    return _gc(model)


def solve_growth(model) -> tuple[float, str]:
    st = solve_status(model)
    return (model.solver.objective.value if st == "optimal" else 0.0), st


def ko_sweep(model, genes: list[str], closures: dict, bounds: dict, variant: Variant, *, verify_n: int = 10,
             seed: int = 11, extra_ids=()) -> dict:
    """Single-gene KO sweep (biomass objective) under ``bounds`` + ``variant``; exact shortcuts as w2_engine."""
    from nextgem.knockouts import WT_GROWTH_TOL

    ko = np.ones(len(genes))
    status = np.array(["shortcut"] * len(genes), dtype=object)
    names = [(r.id, r.forward_variable.name, r.reverse_variable.name) for r in model.reactions]
    with model:
        model.medium = bounds
        variant.apply(model)
        model.objective = BIOMASS
        wt, wst = solve_growth(model)
        if wst != "optimal" or not wt > WT_GROWTH_TOL:
            return {"wt_growth": wt, "wt_status": wst, "ratio_r": np.full(len(genes), np.nan), "n_solved": 0,
                    "n_genes_lp": 0, "shortcut_verify_max_dev": np.nan, "n_infeasible": 0}
        pv = model.solver.primal_values
        carrying = {rid for rid, f, r in names if rid in model.reactions and abs(pv.get(f, 0.0)) + abs(pv.get(r, 0.0)) > ZERO_FLUX_TOL}
        idx = [i for i, g in enumerate(genes) if closures[g] and any(r in carrying for r in closures[g])]
        groups: dict[frozenset, list[int]] = {}
        for i in idx:
            groups.setdefault(frozenset(closures[genes[i]]), []).append(i)
        for members in groups.values():
            with model:
                model.genes.get_by_id(genes[members[0]]).knock_out()
                g, st = solve_growth(model)
            for i in members:
                ko[i] = g / wt
                status[i] = st
        rng = np.random.default_rng(seed)
        short = [i for i in range(len(genes)) if status[i] == "shortcut" and closures[genes[i]]]
        dev = 0.0
        for i in (rng.choice(short, size=min(verify_n, len(short)), replace=False) if short else []):
            with model:
                model.genes.get_by_id(genes[i]).knock_out()
                g, _ = solve_growth(model)
            dev = max(dev, abs(g / wt - 1.0))
    r = np.clip(np.round(ko, RATIO_DECIMALS), 0.0, 1.0) + 0.0
    return {"wt_growth": wt, "wt_status": wst, "ratio_r": r, "n_solved": len(groups), "n_genes_lp": len(idx),
            "shortcut_verify_max_dev": dev, "n_infeasible": int((status == "infeasible").sum())}


def met_by_name(model, name: str, comp: str):
    hits = [m for m in model.metabolites if m.name == name and m.compartment == comp]
    return hits[0] if len(hits) == 1 else (None if not hits else hits)


# ============================================================================= 4. pseudo-tasks
def make_pseudo_task(task, base, maps, saved, rng: np.random.Generator, comp_mets: dict[str, list],
                     max_tries: int = 200):
    """Copy of ``task`` whose *required* items are replaced by random producible metabolites of the same
    compartment (same LB/UB); task equations with LB > 0 become a required output in the compartment of the
    equation's first reactant. Inputs, non-required outputs and CHANGED rows are kept. Producible = the pseudo
    task is feasible on the wild type. Returns (pseudo Task, list of (original, replacement), n_tries) or None."""
    import copy

    from raven_toolbox.tasks.check import apply_task_constraints

    used = {t.upper() for t, _, _ in task.inputs} | {t.upper() for t, _, _ in task.outputs}
    req_out = [(i, t, lb, ub) for i, (t, lb, ub) in enumerate(task.outputs) if lb > 0 and "[" in t]
    req_eq = [(i, e, lb, ub) for i, (e, lb, ub) in enumerate(task.equations) if lb > 0]
    slots = []  # (kind, index, compartment, lb, ub, original)
    for i, t, lb, ub in req_out:
        slots.append(("out", i, t.rsplit("[", 1)[1].rstrip("]"), lb, ub, t))
    for i, e, lb, ub in req_eq:
        lhs = e.split("=>")[0].split("<=")[0]
        first = lhs.split("+")[0].strip()
        first = re.sub(r"^[0-9.]+\s+", "", first)
        comp = first.rsplit("[", 1)[1].rstrip("]") if "[" in first else "c"
        slots.append(("equ", i, comp, lb, ub, e))
    if not slots:
        return None
    tries = 0
    new = copy.deepcopy(task)
    swaps = []
    outputs = list(new.outputs)
    eq_keep = [x for j, x in enumerate(new.equations) if j not in {s[1] for s in slots if s[0] == "equ"}]
    extra_out = []
    for kind, i, comp, lb, ub, orig in slots:
        cands = comp_mets.get(comp, [])
        order = rng.permutation(len(cands))
        ok = False
        for j in order[:max_tries]:
            tok = cands[j]
            if tok.upper() in used:
                continue
            tries += 1
            trial_out = list(outputs)
            trial_extra = list(extra_out)
            if kind == "out":
                trial_out[i] = (tok, lb, ub)
            else:
                trial_extra.append((tok, lb, ub))
            trial = copy.deepcopy(new)
            trial.outputs = trial_out + trial_extra
            trial.equations = eq_keep
            with base:
                tm, err = apply_task_constraints(base, trial, *maps)
                try:
                    feas = err is None and solve_status(base) == "optimal"
                finally:
                    from raven_toolbox.tasks.check import _set_constraint_bounds

                    for mid in (tm or saved):
                        if mid in base.constraints:
                            _set_constraint_bounds(base.constraints[mid], *saved[mid])
            if feas:
                outputs, extra_out = trial_out, trial_extra
                used.add(tok.upper())
                swaps.append((orig, tok))
                ok = True
                break
        if not ok:
            return None
    new.outputs = outputs + extra_out
    new.equations = eq_keep
    new.description = f"PSEUDO of: {task.description}"
    return new, swaps, tries


def compartment_tokens(model, exclude_names=()) -> dict[str, list[str]]:
    """``name[comp]`` tokens per compartment for pseudo-task sampling, sorted (deterministic). Excludes
    pool/pseudo metabolites (names ending in '_biomass' or containing 'pool') and generic tokens."""
    out: dict[str, set] = defaultdict(set)
    ex = {n.lower() for n in exclude_names}
    for m in model.metabolites:
        if not m.name or "biomass" in m.name.lower() or "pool" in m.name.lower() or m.name.lower() in ex:
            continue
        if not m.reactions:
            continue
        out[m.compartment].add(f"{m.name}[{m.compartment}]")
    return {c: sorted(v) for c, v in out.items()}


# ============================================================================= 5. CellFie-inspired activity
CELLFIE_LOW_Q, CELLFIE_HIGH_Q = 0.25, 0.75


def cellfie_thresholds(expr_log2: pd.DataFrame, fit_lines) -> pd.Series:
    """Local min–max–mean threshold per gene (CellFie 'local' rule): the gene's mean linear TPM over the fitting
    lines, clipped to the [25th, 75th] percentiles of all fitting-line gene values."""
    X = np.power(2.0, expr_log2.reindex(index=list(fit_lines)).astype(float)) - 1.0
    vals = X.to_numpy().ravel()
    vals = vals[np.isfinite(vals)]
    lo, hi = np.quantile(vals, CELLFIE_LOW_Q), np.quantile(vals, CELLFIE_HIGH_Q)
    return X.mean(axis=0).clip(lower=lo, upper=hi)


def cellfie_gene_scores(expr_log2: pd.DataFrame, lines, thr: pd.Series) -> pd.DataFrame:
    X = np.power(2.0, expr_log2.reindex(index=list(lines), columns=thr.index).astype(float)) - 1.0
    return 5.0 * np.log1p(X / thr.to_numpy()[None, :])


def eval_gpr_score(node, gs: dict) -> float:
    """AND = min, OR = max over gene scores; genes without a score are ignored (NaN if none scored)."""
    if isinstance(node, str):
        return gs.get(node, np.nan)
    op, ch = node
    vals = [eval_gpr_score(c, gs) for c in ch]
    vals = [v for v in vals if not np.isnan(v)]
    if not vals:
        return np.nan
    return min(vals) if op == "and" else max(vals)


CELLFIE_ACTIVE = 5.0 * math.log(2.0)  # gene at its threshold scores 5·ln 2; task active if mean ≥ this


# ============================================================================= 6. mode scores
def mode_e(ratio: pd.Series, fails: pd.Series | None = None) -> pd.Series:
    """e = 1 if the KO fails a required task, else 1 − ratio_r (DESIGN §6.1). NaN ratio -> 0 (absent gene)."""
    e = (1.0 - ratio).fillna(0.0)
    if fails is not None:
        e = e.where(~fails.reindex(e.index, fill_value=False).astype(bool), 1.0)
    return e.round(6)
