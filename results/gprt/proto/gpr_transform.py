"""Stoichiometric GPR transformation (Machado, Herrgard & Rocha 2016, PLoS Comput Biol e1005140) for cobra models.

Own reimplementation (prototype, label-free); no code from GPRTransform or framed/reframed is imported or run.
Reference used: the paper's Methods ("Model transformation") and the structure of reframed's gpr_transform
(read only).

Transformation (paper, Fig 1 / Methods):
  1. every GPR is put in disjunctive normal form (DNF); each DNF term is one enzyme ("protein"):
     a set of genes that must all be present (complex subunits);
  2. every reaction with a GPR is made irreversible (forward / backward copies, each with bounds [0, |b|]);
  3. reactions with > 1 enzyme are split into one copy per enzyme (isozyme copies);
  4. each gene becomes a pseudo-metabolite ``gpr__<gene>``; every enzyme copy consumes 1 unit of each of
     its subunit genes per unit flux (S_gpr, coefficient -1);
  5. one enzyme-usage reaction ``u__<gene>``: -> gpr__<gene>, bounds [0, usage_ub] (paper: unbounded),
     whose flux is the total flux carried by the gene's enzyme (subunit) over all reactions it catalyses.

Gene knockout = ``u__<gene>`` bounded to 0.

Two isozyme encodings:
  * ``split="paper"``: as in the paper and in reframed: each isozyme copy inherits the original bounds,
    so a reaction with k enzymes can carry up to k x its original bound (capacity relaxation when the
    bound is binding).
  * ``split="arm"`` (default; exact): the irreversible copies keep the original bounds and consume one
    unit of a per-reaction catalysis token ``cat__<rxn>``; each enzyme is an "arm" reaction
    (subunit genes -> cat__<rxn>, unbounded). Same flux space as the original model projected on the
    original reactions, and the same gene-usage semantics as the paper.

Non-paper simplifications (documented in the assessment): reactions without a GPR (exchanges,
spontaneous, transport without genes) are copied unchanged (the paper makes every reaction irreversible;
only gene-associated reactions need it); DNF terms are made minimal by absorption (A or (A and B) -> A),
which leaves every Boolean GPR evaluation, hence every KO, unchanged.
"""

from __future__ import annotations

import math
import time
from collections import Counter
from itertools import product

import cobra

GENE_COMPARTMENT = "gpr"


def gene_met(g: str) -> str:
    return f"gpr__{g}"


def usage(g: str) -> str:
    return f"u__{g}"


# ----------------------------------------------------------------------------- GPR -> DNF
def compile_gpr(rxn) -> object | None:
    """Nested ('and'|'or', children) tuples with gene-id leaves (same format as nextgem.context)."""
    import ast

    def _c(node):
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.BoolOp):
            return ("and" if isinstance(node.op, ast.And) else "or", tuple(_c(v) for v in node.values))
        raise TypeError(ast.dump(node))

    if not rxn.genes or rxn.gpr.body is None:
        return None
    return _c(rxn.gpr.body)


def dnf(node, max_terms: int = 100_000) -> list[frozenset]:
    """Minimal DNF terms of a compiled GPR (absorption applied); raises if it exceeds ``max_terms``."""
    if isinstance(node, str):
        return [frozenset([node])]
    op, ch = node
    parts = [dnf(c, max_terms) for c in ch]
    if op == "or":
        terms = [t for p in parts for t in p]
    else:
        n = math.prod(len(p) for p in parts)
        if n > max_terms:
            raise OverflowError(f"DNF expansion {n} > {max_terms}")
        terms = [frozenset().union(*combo) for combo in product(*parts)]
    return minimal(terms)


def minimal(terms) -> list[frozenset]:
    uniq = sorted(set(terms), key=lambda t: (len(t), sorted(t)))
    out = []
    for t in uniq:
        if not any(s <= t for s in out):
            out.append(t)
    return out


def gpr_class(node) -> str:
    """'single' (one gene), 'or' (isozymes only), 'and' (one complex), 'nested' (mixed AND/OR)."""
    if node is None:
        return "none"
    if isinstance(node, str):
        return "single"
    ops = set()

    def walk(n):
        if isinstance(n, str):
            return
        ops.add(n[0])
        for c in n[1]:
            walk(c)

    walk(node)
    if ops == {"or"}:
        return "or"
    if ops == {"and"}:
        return "and"
    return "nested"


def gpr_stats(model) -> dict:
    """Counts of GPR types and DNF sizes for a cobra model."""
    cls, n_terms, max_sub, dnf_blowup, depth_not_dnf, n_complex = Counter(), [], 0, 0, 0, 0
    for r in model.reactions:
        node = compile_gpr(r)
        c = gpr_class(node)
        cls[c] += 1
        if node is None:
            continue
        terms = dnf(node)
        n_terms.append(len(terms))
        n_complex += any(len(t) > 1 for t in terms)
        max_sub = max(max_sub, max(len(t) for t in terms))
        if sum(len(t) for t in terms) > _n_leaves(node):
            dnf_blowup += 1  # DNF has more gene occurrences than the written rule (nested AND of ORs)
        if c == "nested" and not _is_dnf(node):
            depth_not_dnf += 1
    return {
        "reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes),
        "gpr_none": cls["none"], "gpr_single": cls["single"], "gpr_or_only": cls["or"], "gpr_and_only": cls["and"],
        "gpr_nested_mixed": cls["nested"], "gpr_nested_not_dnf": depth_not_dnf,
        "gpr_with_isozymes(dnf_terms>1)": sum(1 for k in n_terms if k > 1),
        "gpr_with_complex(any term>1 gene)": n_complex,
        "dnf_terms_max": max(n_terms), "dnf_terms_total": sum(n_terms), "dnf_subunits_max": max_sub,
        "gpr_dnf_expanded_more_than_written": dnf_blowup,
    }


def _n_leaves(node) -> int:
    return 1 if isinstance(node, str) else sum(_n_leaves(c) for c in node[1])


def _is_dnf(node) -> bool:
    """OR of ANDs of genes (or simpler)."""
    if isinstance(node, str):
        return True
    op, ch = node
    if op == "and":
        return all(isinstance(c, str) for c in ch)
    return all(isinstance(c, str) or (c[0] == "and" and all(isinstance(x, str) for x in c[1])) for c in ch)


# ----------------------------------------------------------------------------- transform
def transform(model, split: str = "arm", usage_ub: float = math.inf, max_terms: int = 100_000):
    """Return (transformed cobra.Model, info). ``model`` should already carry the medium and objective.

    info: {"map": orig rxn id -> [(new id, sign)], "terms": rxn id -> list of DNF terms, counts...}.
    Objective coefficients of gene-associated reactions are carried to forward (+c) / backward (-c) copies.
    """
    if split not in ("paper", "arm"):
        raise ValueError(split)
    t0 = time.process_time()
    tm = cobra.Model(f"{model.id}_gprt_{split}")
    mets = {m.id: cobra.Metabolite(m.id, formula=m.formula, name=m.name, compartment=m.compartment, charge=m.charge)
            for m in model.metabolites}
    genes = sorted(g.id for g in model.genes)
    gmets = {g: cobra.Metabolite(gene_met(g), name=f"enzyme {g}", compartment=GENE_COMPARTMENT) for g in genes}
    new_mets = list(mets.values()) + list(gmets.values())
    rxns, mapping, terms_of, obj = [], {}, {}, {}
    counts = Counter()

    def new_rxn(rid, name, lb, ub, stoich, subsystem=""):
        r = cobra.Reaction(rid, name=name, lower_bound=lb, upper_bound=ub, subsystem=subsystem)
        r.add_metabolites(stoich)
        rxns.append(r)
        return r

    for r in model.reactions:
        base = {mets[m.id]: c for m, c in r.metabolites.items()}
        node = compile_gpr(r)
        c_obj = r.objective_coefficient
        if node is None:
            new_rxn(r.id, r.name, r.lower_bound, r.upper_bound, base, r.subsystem)
            mapping[r.id] = [(r.id, 1)]
            if c_obj:
                obj[r.id] = c_obj
            counts["copied_no_gpr"] += 1
            continue
        terms = dnf(node, max_terms)
        terms_of[r.id] = terms
        lb, ub = r.lower_bound, r.upper_bound
        dirs = []
        if ub > 0:
            dirs.append(("f", 1, max(lb, 0.0), ub))
        if lb < 0:
            dirs.append(("b", -1, max(-ub, 0.0), -lb))
        if not dirs:  # lb = ub = 0: blocked; keep a blocked forward copy so ids stay mappable
            dirs.append(("f", 1, 0.0, 0.0))
        if any(lo > 0 for _, _, lo, _ in dirs):
            counts["forced_flux_gpr_reactions"] += 1
        mapping[r.id] = []
        multi = len(terms) > 1
        counts["isozyme_reactions"] += multi
        counts["complex_reactions"] += any(len(t) > 1 for t in terms)
        if multi and split == "arm":
            cat = cobra.Metabolite(f"cat__{r.id}", name=f"catalysis {r.id}", compartment=GENE_COMPARTMENT)
            new_mets.append(cat)
            for d, s, lo, hi in dirs:
                rid = f"{r.id}__{d}"
                st = {m: s * c for m, c in base.items()}
                st[cat] = st.get(cat, 0) - 1.0
                new_rxn(rid, r.name, lo, hi, st, r.subsystem)
                mapping[r.id].append((rid, s))
                if c_obj:
                    obj[rid] = s * c_obj
            for j, t in enumerate(terms, 1):
                st = {gmets[g]: -1.0 for g in sorted(t)}
                st[cat] = 1.0
                new_rxn(f"{r.id}__arm{j}", f"enzyme {j} of {r.id}", 0.0, math.inf, st, "GPR arm")
                counts["arm_reactions"] += 1
        else:
            for d, s, lo, hi in dirs:
                for j, t in enumerate(terms, 1):
                    rid = f"{r.id}__{d}" if not multi else f"{r.id}__{d}__iso{j}"
                    st = {m: s * c for m, c in base.items()}
                    for g in t:
                        st[gmets[g]] = -1.0
                    # paper: each isozyme copy inherits the bounds; a forced lower bound would be multiplied,
                    # so isozyme copies get lower bound 0 (forced flux reactions are counted above)
                    new_rxn(rid, r.name, lo if not multi else 0.0, hi, st, r.subsystem)
                    mapping[r.id].append((rid, s))
                    if c_obj:
                        obj[rid] = s * c_obj
    for g in genes:
        new_rxn(usage(g), f"usage of {g}", 0.0, usage_ub, {gmets[g]: 1.0}, "GPR usage")
    t1 = time.process_time()
    tm.add_metabolites(new_mets)
    tm.add_reactions(rxns)
    tm.objective = {tm.reactions.get_by_id(k): v for k, v in obj.items()}
    tm.objective_direction = model.objective_direction
    info = {
        "split": split, "map": mapping, "terms": terms_of, "genes": genes,
        "reactions_before": len(model.reactions), "metabolites_before": len(model.metabolites),
        "reactions_after": len(tm.reactions), "metabolites_after": len(tm.metabolites),
        "usage_reactions": len(genes), "build_cpu_s": t1 - t0, "add_cpu_s": time.process_time() - t1,
        **counts,
    }
    return tm, info


def knock_out(tm, gene: str) -> None:
    """Gene KO in the transformed model (use inside ``with tm:``)."""
    tm.reactions.get_by_id(usage(gene)).bounds = (0.0, 0.0)


def merged_fluxes(tm, info, primal=None) -> dict:
    """Net flux of each original reaction from a transformed solution."""
    if primal is None:
        primal = {r.id: r.flux for r in tm.reactions}
    return {rid: sum(s * primal[n] for n, s in parts) for rid, parts in info["map"].items()}


def usage_fluxes(tm, info) -> dict:
    pv = tm.solver.primal_values
    return {g: pv[tm.reactions.get_by_id(usage(g)).forward_variable.name] for g in info["genes"]}
