"""W2 sweep definitions and the per-medium generic KO engine (label-free).

Sweeps run per distinct (parent x bounds signature), not per line (DESIGN W2 "Computation"). A sweep definition
carries every ladder label that maps to its signature (e.g. M2(f) of a serum-free medium equals its M1).

KO engine (same semantics and status policy as ``nextgem.knockouts``; equality is tested on M0 against
results/baseline): ``Gene.knock_out()`` in a context, FBA, ratio = KO / WT growth, ``ratio_r`` rounded to 6
decimals and clipped to [0, 1]. Three exact shortcuts skip LPs whose answer is known:
  (a) the KO closes no reaction (GPR-buffered) -> the LP is unchanged -> ratio 1;
  (b) every reaction the KO closes carries zero flux (|v| <= ZERO_FLUX_TOL) in the WT optimum -> that optimum stays
      feasible, so the KO optimum equals WT -> ratio 1 (DESIGN §12 "zero-flux shortcut is exact");
  (c) genes whose KOs close the identical reaction set (e.g. subunits of one complex) share one LP solve.
Workers are ``spawn`` processes holding both parents on Gurobi, one thread each; the pool size is the Gurobi session
budget (the parent process builds bounds on GLPK-loaded models and holds no Gurobi session).
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from nextgem import media
from nextgem.knockouts import RATIO_DECIMALS, WT_GROWTH_TOL
from nextgem.model import BIOMASS_REACTION, hams_medium
from nextgem.w2_runtime import PARENTS, load_parent, w2_dir

ZERO_FLUX_TOL = 1e-9


@dataclass
class SweepDef:
    parent: str
    signature: str
    bounds: dict[str, float]
    labels: list[dict] = field(default_factory=list)  # each: variant, medium_key, f, gate_open, role

    @property
    def sweep_id(self) -> str:
        return f"{self.parent}:{self.signature}"


def gate_matters(bounds: dict[str, float]) -> bool:
    return any(media.EXCHANGES[c][0] in bounds for c in media.GATED)


def build_sweeps(exchanges: dict[str, list[str]], hams: dict[str, dict[str, float]],
                 pcs: dict[str, media.ParentConstraints], doc_media: dict[str, media.Medium],
                 phys_media: dict[str, media.Medium], sens_media: dict[str, media.Medium],
                 f_levels=media.SERUM_F_LEVELS) -> list[SweepDef]:
    """All ladder sweeps per parent, deduplicated by bounds signature.

    ``exchanges``/``hams``: per parent, the exchange ids (M-open) and Ham's medium (M0) — passed in so that no
    model object is needed here.
    """

    out: dict[tuple[str, str], SweepDef] = {}

    def add(parent, bounds, **label):
        sig = media.bounds_signature(bounds)
        sd = out.setdefault((parent, sig), SweepDef(parent, sig, bounds))
        sd.labels.append(label)

    for p in exchanges:
        pc = pcs[p]
        add(p, {e: media.UPTAKE_OPEN for e in exchanges[p]}, variant="M-open", medium_key="", f=None,
            gate_open=True, role="negative_control")
        add(p, dict(hams[p]), variant="M0", medium_key="Hams_GR", f=None, gate_open=True, role="baseline")
        groups = [("documented", doc_media), ("physiological", phys_media), ("sensitivity", sens_media)]
        for role, meds in groups:
            for key, med in meds.items():
                # physiological counterfactuals (and their sensitivities) are M2-only (DESIGN: + the M2 serum model)
                variants = ("M2",) if "physiological_counterfactual" in med.flags else ("M1", "M2")
                for v in variants:
                    for f in ([None] if v == "M1" else list(f_levels)):
                        b = _bounds(v, med, pc, f, True)
                        add(p, b, variant=v, medium_key=key, f=f, gate_open=True, role=role)
                        if gate_matters(b):
                            add(p, _bounds(v, med, pc, f, False), variant=v + "g", medium_key=key, f=f,
                                gate_open=False, role=role)
    return list(out.values())


def _bounds(variant, med, pc, f, gate_open):
    return media.uptake_bounds(None, variant, med, pc, f=f, gate_open=gate_open)


def labels_frame(sweeps: list[SweepDef]) -> pd.DataFrame:
    rows = []
    for s in sweeps:
        for lab in s.labels:
            rows.append(dict(sweep_id=s.sweep_id, parent=s.parent, signature=s.signature,
                             n_open=len(s.bounds), **lab))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------- worker side
_MODELS: dict = {}
_CLOSURES: dict = {}
_GENES: dict = {}
_VARNAMES: dict = {}


_SOLVER = {"name": "gurobi"}


def _init_worker(parents: tuple[str, ...], solver: str = "gurobi") -> None:
    from nextgem.w2_runtime import load_parent

    # Gurobi: load_parent holds one shared-pool slot per worker process until it exits. GLPK: no slot.
    _SOLVER["name"] = solver
    for p in parents:
        m = load_parent(p, solver=solver)
        m.objective = BIOMASS_REACTION
        _MODELS[p] = m
        _GENES[p] = sorted(g.id for g in m.genes)
        _CLOSURES[p] = gene_closures(m)
        _VARNAMES[p] = [(r.id, r.forward_variable.name, r.reverse_variable.name) for r in m.reactions]


def gene_closures(model) -> dict[str, list[str]]:
    """gene -> reactions newly closed by its KO (medium-independent)."""
    out = {}
    for g in model.genes:
        before = {r.id: r.bounds for r in g.reactions}
        with model:
            g.knock_out()
            out[g.id] = [r.id for r in g.reactions if r.bounds == (0, 0) and before[r.id] != (0, 0)]
    return out


def _solve(model) -> tuple[float, str]:
    from optlang.interface import INFEASIBLE, OPTIMAL

    model.slim_optimize(error_value=math.nan)
    st = model.solver.status
    if st == OPTIMAL:
        return model.solver.objective.value, st
    return (0.0 if st == INFEASIBLE else math.nan), st


def wt_task(args) -> dict:
    parent, sweep_id, bounds = args
    m = _MODELS[parent]
    t0 = time.perf_counter()
    with m:
        m.medium = bounds
        g, st = _solve(m)
    return dict(sweep_id=sweep_id, parent=parent, wt_growth=g, wt_status=st, runtime_s=time.perf_counter() - t0)


def wt_with_override_task(args) -> dict:
    """WT growth for an arbitrary bounds dict (calibration grid, floor non-binding check)."""
    parent, tag, bounds = args
    m = _MODELS[parent]
    with m:
        m.medium = bounds
        g, st = _solve(m)
    return dict(tag=tag, parent=parent, growth=g, status=st)


def ko_task(args) -> dict:
    """Full single-gene KO sweep for one (parent, bounds) with the two exact shortcuts."""
    from optlang.interface import INFEASIBLE, OPTIMAL

    parent, sweep_id, bounds, verify_n = args
    m = _MODELS[parent]
    genes = _GENES[parent]
    clos = _CLOSURES[parent]
    t0 = time.perf_counter()
    ko = np.ones(len(genes))
    status = np.array(["shortcut"] * len(genes), dtype=object)
    with m:
        m.medium = bounds
        wt, wst = _solve(m)
        if wst != OPTIMAL:
            raise RuntimeError(f"{sweep_id}: WT status {wst}")
        pv = m.solver.primal_values  # keyed by variable name: forward = reaction id, reverse = id_reverse_<hash>
        carrying = {rid for rid, fwd, rev in _VARNAMES[parent] if abs(pv[fwd]) + abs(pv[rev]) > ZERO_FLUX_TOL}
        idx_solve = [i for i, g in enumerate(genes) if clos[g] and any(r in carrying for r in clos[g])]
        # (c) genes whose KOs close the identical reaction set give the identical LP: solve once per set
        groups: dict[frozenset, list[int]] = {}
        for i in idx_solve:
            groups.setdefault(frozenset(clos[genes[i]]), []).append(i)
        for members in groups.values():
            with m:
                m.genes.get_by_id(genes[members[0]]).knock_out()
                g, st = _solve(m)
            if st not in (OPTIMAL, INFEASIBLE):
                raise RuntimeError(f"{sweep_id}: KO {genes[members[0]]} status {st}")
            for i in members:
                ko[i] = g / wt if wt > WT_GROWTH_TOL else math.nan
                status[i] = st
        # verification sample: solve some shortcut genes explicitly; must equal WT growth
        rng = np.random.default_rng(11)
        short = [i for i in range(len(genes)) if status[i] == "shortcut" and clos[genes[i]]]
        max_dev = 0.0
        for i in rng.choice(short, size=min(verify_n, len(short)), replace=False) if short else []:
            with m:
                m.genes.get_by_id(genes[i]).knock_out()
                g, _ = _solve(m)
            max_dev = max(max_dev, abs(g / wt - 1.0))
    ratio_r = np.clip(np.round(ko, RATIO_DECIMALS), 0.0, 1.0) + 0.0
    if not wt > WT_GROWTH_TOL:
        ratio_r[:] = np.nan
    return dict(sweep_id=sweep_id, parent=parent, wt_growth=wt, n_solved=len(groups), n_genes_lp=len(idx_solve),
                solver=_SOLVER["name"],
                n_infeasible=int((status == INFEASIBLE).sum()), ratio_r=ratio_r.astype(np.float64),
                runtime_s=time.perf_counter() - t0, shortcut_verify_max_dev=max_dev, n_verified=min(verify_n, len(short)))


def make_pool(processes: int, parents=("human1", "human2"), solver: str = "gurobi"):
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context

    return ProcessPoolExecutor(processes, mp_context=get_context("spawn"), initializer=_init_worker,
                               initargs=(tuple(parents), solver))


def sweep_context(use_mult: bool = True):
    """Everything needed to define sweeps (parents loaded on GLPK: no Gurobi session in this process)."""
    lib = w2_dir() / "media_library"
    lines = pd.read_csv(lib / "line_media.tsv", sep="\t")
    sheets = media.load_sheets()
    doc = media.documented_media(set(lines.medium_text_used.dropna()), sheets)
    assert set(doc) == set(lines.medium_key), set(lines.medium_key) ^ set(doc)
    phys = {f"M-{n}": media.physiological_medium(n, sheets) for n in ("HPLM", "Plasmax")}
    sens = {**media.sensitivity_media(sheets), **media.physiological_sensitivity_media(sheets)}
    pcs = media.floor_constraints(w2_dir() / "floor" / "floor.tsv", w2_dir() / "floor" / "requirements.tsv",
                                  w2_dir() / "floor" / "floor_mult.tsv" if use_mult else None)
    exchanges, hams = {}, {}
    for p in PARENTS:
        m = load_parent(p, solver="glpk")
        exchanges[p] = [r.id for r in m.exchanges]
        hams[p] = hams_medium(m)
        del m
    return lines, sheets, doc, phys, sens, pcs, exchanges, hams


def genes_task(parent: str) -> list[str]:
    """Sorted gene ids of a parent as held by the worker (column order of the ratio_r matrices)."""
    return _GENES[parent]
