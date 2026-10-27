"""Reviewer's own minimal implementations (independent of nextgem.context except where stated).

GPR logic uses cobra's own GPR AST (``reaction.gpr.body``) and ``GPR.eval(knockouts)``; LPs are built with
plain cobra/optlang calls. Semantics follow DESIGN W4 + README §5 (ratified by amendment A4):
  presence: AND = min, OR = max; "no evidence" genes are ignored by AND and make an OR unknown unless a known
  child is >= 0.   E-Flux capacity: AND = min, OR = sum, linear value 2**s for s >= 0, 0 for s < 0.
"""

from __future__ import annotations

import ast
import math

import numpy as np

BIO = "MAR13082"
TOL = 1e-6


# ------------------------------------------------------------------ GPR evaluation
def presence(node, s):
    """Return (value, known). Unknown = no evidence."""
    if isinstance(node, ast.Name):
        v = s.get(node.id)
        return (v, True) if v is not None and not math.isnan(v) else (math.nan, False)
    vals = [presence(c, s) for c in node.values]
    known = [v for v, k in vals if k]
    if isinstance(node.op, ast.And):
        return (min(known), True) if known else (math.nan, False)
    # OR
    if not known:
        return (math.nan, False)
    mx = max(known)
    if mx >= 0 or len(known) == len(vals):
        return (mx, True)
    return (math.nan, False)


def capacity(node, lin):
    if isinstance(node, ast.Name):
        v = lin.get(node.id)
        return v if v is not None else math.nan
    vals = [capacity(c, lin) for c in node.values]
    if isinstance(node.op, ast.And):
        k = [v for v in vals if not math.isnan(v)]
        return min(k) if k else math.nan
    return math.nan if any(math.isnan(v) for v in vals) else float(sum(vals))


def reaction_scores(model, s):
    out = {}
    for r in model.reactions:
        if r.genes and r.gpr.body is not None:
            out[r.id] = presence(r.gpr.body, s)[0]
    return out


# ------------------------------------------------------------------ GIMME (Becker & Palsson 2008)
def gimme(model, s_r, wt_max, fraction=0.9):
    """min sum_r c_r |v_r|, c_r = threshold - expression (= -s_r) for s_r < 0; biomass >= fraction * wt_max.
    Keep: s_r >= 0, unknown, or |v_r| > TOL. Returns (closed, objective)."""
    low = {rid: -v for rid, v in s_r.items() if not math.isnan(v) and v < 0}
    with model:
        model.reactions.get_by_id(BIO).lower_bound = fraction * wt_max
        terms = []
        for rid, c in low.items():
            r = model.reactions.get_by_id(rid)
            terms += [c * r.forward_variable, c * r.reverse_variable]
        model.objective = model.problem.Objective(sum(terms), direction="min")
        model.slim_optimize()
        assert model.solver.status == "optimal", model.solver.status
        obj = model.solver.objective.value
        flux = {rid: model.reactions.get_by_id(rid).flux for rid in low}
    closed = {rid for rid, v in flux.items() if abs(v) <= TOL}
    return closed, obj


# ------------------------------------------------------------------ E-Flux (Colijn 2009), DESIGN F4
def eflux_lin(s):
    return {g: (2.0 ** v if v >= 0 else 0.0) for g, v in s.items() if v is not None and not math.isnan(v)}


def eflux_bounds(model, lin, cap_max=None, floor=10.0, scale=1000.0):
    caps = {}
    for r in model.reactions:
        if r.genes and not r.boundary and r.gpr.body is not None:
            caps[r.id] = capacity(r.gpr.body, lin)
    if cap_max is None:
        cap_max = max(c for c in caps.values() if not math.isnan(c))
    b = {}
    for rid, c in caps.items():
        if math.isnan(c):
            continue
        r = model.reactions.get_by_id(rid)
        bb = max(floor, scale * c / cap_max)
        lb = -bb if r.lower_bound < 0 else r.lower_bound
        ub = bb if r.upper_bound > 0 else r.upper_bound
        lb = max(lb, r.lower_bound)
        ub = min(ub, r.upper_bound)
        b[rid] = (lb, ub)
    return b, cap_max


# ------------------------------------------------------------------ KO sweep (brute force; A4 (a) convention)
def ko_presence(model, final_closed, off, lost=frozenset(), genes=None):
    """Brute-force single-gene KO ratios of a presence-family line model. An open reaction whose GPR is True
    with off+lost genes removed uses the pruned GPR (KO of an off gene has no effect there); otherwise the
    full GPR (lost genes removed)."""
    eff = {}
    for r in model.reactions:
        if r.id in final_closed or not r.genes:
            continue
        ko_set = set(lost) | set(off)
        eff[r.id] = ko_set if r.gpr.eval(knockouts=ko_set) else set(lost)
    gene_rxns = {}
    for rid in eff:
        for g in model.reactions.get_by_id(rid).genes:
            gene_rxns.setdefault(g.id, []).append(rid)
    return _sweep(model, {rid: (0.0, 0.0) for rid in final_closed}, genes or [g.id for g in model.genes],
                  lambda g: {rid: (0.0, 0.0) for rid in gene_rxns.get(g, [])
                             if not model.reactions.get_by_id(rid).gpr.eval(knockouts=eff[rid] | {g})})


def ko_eflux(model, base_bounds, lin, cap_max, lost_closed=frozenset(), lost=frozenset(), literal=False, genes=None):
    """E-Flux KO: the gene's linear value is set to 0 and the bounds re-derived. Reactions whose GPR is False
    after the KO are closed (author's reading) unless ``literal`` (then they get the floor bound like any
    zero-capacity reaction)."""
    pb = {r.id: (r.lower_bound, r.upper_bound) for r in model.reactions}  # parent bounds, before any line bounds
    fixed = dict(base_bounds)
    fixed.update({rid: (0.0, 0.0) for rid in lost_closed})
    gene_rxns = {}
    for r in model.reactions:
        if r.genes and not r.boundary and r.id not in lost_closed:
            for g in r.genes:
                gene_rxns.setdefault(g.id, []).append(r.id)
    gene_rxns_all = {}
    for r in model.reactions:
        if r.genes and r.id not in lost_closed:
            for g in r.genes:
                gene_rxns_all.setdefault(g.id, []).append(r.id)

    def change(g):
        out = {}
        l2 = dict(lin)
        l2[g] = 0.0
        for rid in gene_rxns_all.get(g, []):
            r = model.reactions.get_by_id(rid)
            if not literal and not r.gpr.eval(knockouts=set(lost) | {g}):
                out[rid] = (0.0, 0.0)
                continue
            if rid not in base_bounds:
                continue  # unknown capacity or boundary: keeps its default bound
            c = capacity(r.gpr.body, l2)
            if math.isnan(c):
                continue
            bb = max(10.0, 1000.0 * c / cap_max)
            plb, pub = pb[rid]
            lb = max(-bb if plb < 0 else plb, plb)
            ub = min(bb if pub > 0 else pub, pub)
            if (lb, ub) != base_bounds[rid]:
                out[rid] = (lb, ub)
        return out

    return _sweep(model, fixed, genes or [g.id for g in model.genes], change)


def _sweep(model, line_bounds, genes, change):
    res = {}
    with model:
        for rid, b in line_bounds.items():
            model.reactions.get_by_id(rid).bounds = b
        model.objective = BIO
        wt = model.slim_optimize()
        assert model.solver.status == "optimal"
        for g in genes:
            ch = change(g)
            if not ch:
                res[g] = 1.0
                continue
            with model:
                for rid, b in ch.items():
                    model.reactions.get_by_id(rid).bounds = b
                v = model.slim_optimize()
                st = model.solver.status
            res[g] = (v / wt) if st == "optimal" else (0.0 if st == "infeasible" else np.nan)
    return wt, res
