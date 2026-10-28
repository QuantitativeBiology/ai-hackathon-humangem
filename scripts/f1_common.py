"""F1 (post hoc / exploratory, DESIGN §17 A15; RESEARCH_TREE D17): GPR curation with the stricter paralog rule on top
of P1. Shared definitions and pure helpers. Label-free: nothing here opens a label file.

The D17 rule (fixed in the tree before any F1 scoring), as implemented here:

* Edit type: AND -> OR merges only (the W5 AND -> OR edit: AND-ed *leaf* sibling genes are merged into one OR; a
  connected component of qualifying pairs is merged only if it is a clique, otherwise it is logged as blocked).
  W5's other edit types (OR -> AND, pan-cancer-absent OR-branch removal) are not part of D17 and are not applied.
* A pair of AND-ed leaf siblings (a, b) qualifies iff all of
    1. Ensembl 116 paralog pair (any homology type) with min(%id a->b, %id b->a) >= 30 (``MIN_PID``);
    2. same compartment: both genes have a non-empty Human-GEM ``genes.tsv`` ``compartments`` annotation of the
       edited parent and the two compartment *sets* are equal (order-insensitive; missing annotation = not same);
    3. no shared CORUM 5.1 / Complex Portal complex (the W5 AND -> OR leg, kept: D17 makes the paralog rule stricter,
       it does not drop W5's complex check).
* Reaction-level exclusions (no merge is applied in an excluded reaction; the edit unit is the single reaction, not
  W5's GPR group, so a GPR shared by an eligible and an excluded reaction is edited only in the eligible one):
    E1 exchange: ``reaction.boundary`` or subsystem "Exchange/demand reactions";
    E2 transport: subsystem "Transport reactions" or metabolites in more than one compartment;
    E3 pseudo/lumped (label-free, from Human-GEM annotations; :func:`lumped_reason`):
         (a) subsystem "Artificial reactions" or "Pool reactions" (Human-GEM's own pseudo-reaction subsystems), or
         (b) composition-weighted lumping: one side of the reaction has >= 2 metabolites with non-integer
             coefficients that sum to 1 (a mixture drawn from a composition, e.g. MAR07161
             0.18 ATP + 0.3 CTP + 0.34 GTP + 0.18 UTP -> RNA + PPi; MAR07160 dNTPs -> DNA; MAR07162-07164);
    E4 large complex: some AND node of the reaction's canonical GPR spans >= 6 distinct genes (``MAX_AND_GENES``).
* Plus the H1 TYMS irreversibility fix (DESIGN §3.1): MAR04644 bounds (0, 1000), as Human2.

Arms (all on P1's parent and medium: Human1 + W3 B-loo:cofactor_pool_biomass on M2(f=1), line-specific medium as in
results/w2/controls/line_sweeps.tsv): P1 (reference; the M6 cell H1+Bloo/M2), TYMS only, GPR rule only, and the
curated P1 = GPR rule + TYMS (the F1 candidate).
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
F1 = PROJECT / "results" / "f1"
MANIFEST = F1 / "MANIFEST.tsv"
SWDIR = F1 / "sweeps"
TRDIR = F1 / "train"
TMP = F1 / "tmp"  # git-ignored
AUTHOR = "F1 engineer agent (Claude Opus; A15 post hoc)"
MAX_WORKERS = 8
PARENT = "human1"

# ---- the D17 rule
MIN_PID = 30.0
MAX_AND_GENES = 6  # a reaction with an AND node spanning >= 6 genes is excluded
LUMPED_SUBSYSTEMS = ("Artificial reactions", "Pool reactions")
TRANSPORT_SUBSYSTEM = "Transport reactions"
EXCHANGE_SUBSYSTEM = "Exchange/demand reactions"
COMP_TOL = 1e-6
TYMS_REACTION = "MAR04644"
TYMS_BOUNDS = (0.0, 1000.0)

# ---- arms (method names are the score-file / look-spec column names)
P1_METHOD = "M6:H1+Bloo/M2(f=1)"  # = M6 cell H1+Bloo/M2 (P1, D11)


@dataclass(frozen=True)
class Arm:
    key: str
    tyms: bool
    gpr: bool
    method: str
    role: str


ARMS = (
    Arm("P1", False, False, P1_METHOD, "reference"),
    Arm("TYMS", True, False, "F1:H1+Bloo+TYMS/M2(f=1)", "component"),
    Arm("GPR", False, True, "F1:H1+Bloo+GPR/M2(f=1)", "component"),
    Arm("F1", True, True, "F1:H1+Bloo+GPR+TYMS/M2(f=1)", "candidate"),
)
ARM = {a.key: a for a in ARMS}
CANDIDATE = ARM["F1"]

TG_MIN_DELTA = 0.005
TG_MIN_FOLDS = 9


# ============================================================================= rule helpers (pure)
def is_exchange(boundary: bool, subsystem: str) -> bool:
    return bool(boundary) or (subsystem or "").strip() == EXCHANGE_SUBSYSTEM


def is_transport(subsystem: str, compartments: set) -> bool:
    return (subsystem or "").strip() == TRANSPORT_SUBSYSTEM or len(compartments) > 1


def composition_lumped(stoich: dict) -> bool:
    """True if one side has >= 2 non-integer coefficients summing to 1 (composition-weighted mixture)."""
    for sign in (-1, 1):
        fr = [abs(c) for c in stoich.values() if c * sign > 0 and abs(abs(c) - round(abs(c))) > 1e-9]
        if len(fr) >= 2 and abs(sum(fr) - 1.0) <= COMP_TOL:
            return True
    return False


def lumped_reason(subsystem: str, stoich: dict) -> str:
    """'' if not pseudo/lumped, else the reason (E3a subsystem / E3b composition)."""
    s = (subsystem or "").strip()
    if s in LUMPED_SUBSYSTEMS:
        return f"subsystem:{s}"
    if composition_lumped(stoich):
        return "composition-weighted lumped stoichiometry"
    return ""


def _genes(node) -> set:
    if isinstance(node, str):
        return {node}
    return set().union(*(_genes(c) for c in node[1]))


def max_and_genes(tree) -> int:
    """Largest number of distinct genes under any AND node of a canonical GPR tree (0 if no AND node)."""
    if tree is None or isinstance(tree, str):
        return 0
    op, ch = tree
    here = len(_genes(tree)) if op == "and" else 0
    return max([here] + [max_and_genes(c) for c in ch])


def comp_set(s) -> frozenset:
    if s is None or (isinstance(s, float) and math.isnan(s)):
        return frozenset()
    return frozenset(t.strip() for t in str(s).split(";") if t.strip())


def same_compartment(ca, cb) -> bool:
    a, b = comp_set(ca), comp_set(cb)
    return bool(a) and a == b


def reaction_exclusions(*, boundary: bool, subsystem: str, stoich: dict, met_comps: set, tree) -> list[str]:
    """Every reaction-level exclusion that applies (E1-E4), in that order; [] = eligible."""
    out = []
    if is_exchange(boundary, subsystem):
        out.append("E1 exchange")
    if is_transport(subsystem, met_comps):
        out.append("E2 transport")
    lr = lumped_reason(subsystem, stoich)
    if lr:
        out.append(f"E3 pseudo/lumped ({lr})")
    if max_and_genes(tree) >= MAX_AND_GENES:
        out.append(f"E4 AND complex >= {MAX_AND_GENES} genes")
    return out


@dataclass(frozen=True)
class PairFacts:
    paralog: bool  # Ensembl 116 pair of any type
    ptype: str | None
    pid: float  # min %id of both directions (NaN if not a pair)
    shared_complexes: tuple
    comp_a: frozenset
    comp_b: frozenset

    @property
    def pid_ok(self) -> bool:
        return self.paralog and self.pid >= MIN_PID

    @property
    def comp_ok(self) -> bool:
        return bool(self.comp_a) and self.comp_a == self.comp_b

    @property
    def complex_ok(self) -> bool:
        return not self.shared_complexes

    @property
    def qualifies(self) -> bool:
        return self.pid_ok and self.comp_ok and self.complex_ok

    def fail_reasons(self) -> list[str]:
        r = []
        if not self.paralog:
            r.append("not an Ensembl paralog pair")
        elif not self.pid_ok:
            r.append(f"min %id {self.pid:.1f} < {MIN_PID:g}")
        if not self.comp_ok:
            r.append("compartment differs or missing")
        if not self.complex_ok:
            r.append("shared complex")
        return r


def curate_tree(tree, facts, rec: list, ctx: dict, apply: bool = True):
    """D17 AND -> OR merges on a canonical tree, bottom-up. ``facts(a, b)`` -> PairFacts. Appends one record per
    examined paralog pair (kind 'pair'), per merged clique (kind 'merge') and per blocked non-clique component
    (kind 'blocked') to ``rec``. With ``apply=False`` (reaction excluded) nothing is merged but records are written
    with ``applied=False``. Returns the new tree (canonicalised with ``nextgem.w5.canon``)."""
    from nextgem import w5

    if tree is None or isinstance(tree, str):
        return tree
    op, ch = tree
    ch = tuple(curate_tree(c, facts, rec, ctx, apply) for c in ch)
    node = (op, ch)
    if op != "and":
        return w5.canon(node)
    leaves = sorted(c for c in ch if isinstance(c, str))
    ok = {}
    for a, b in combinations(leaves, 2):
        f = facts(a, b)
        if f.paralog:
            rec.append({**ctx, "kind": "pair", "a": a, "b": b, "paralog_type": f.ptype, "min_pid": f.pid,
                        "comp_a": ";".join(sorted(f.comp_a)), "comp_b": ";".join(sorted(f.comp_b)),
                        "shared_complexes": ";".join(f.shared_complexes), "pid_ok": f.pid_ok, "comp_ok": f.comp_ok,
                        "complex_ok": f.complex_ok, "pair_qualifies": f.qualifies,
                        "fail": "; ".join(f.fail_reasons()), "applied": bool(apply and f.qualifies)})
        if f.qualifies:
            ok[(a, b)] = f
    if not ok:
        return w5.canon(node)
    adj: dict[str, set] = {}
    for a, b in ok:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    seen, comps = set(), []
    for v in sorted(adj):
        if v in seen:
            continue
        stack, comp = [v], set()
        while stack:
            x = stack.pop()
            if x not in comp:
                comp.add(x)
                stack.extend(adj[x] - comp)
        seen |= comp
        comps.append(sorted(comp))
    rest, merged = list(ch), []
    for comp in comps:
        pairs = list(combinations(comp, 2))
        if not all(p in ok for p in pairs):
            rec.append({**ctx, "kind": "blocked", "genes": ";".join(comp), "reason": "component not a clique",
                        "applied": False})
            continue
        rec.append({**ctx, "kind": "merge", "genes": ";".join(comp), "old_node": w5.to_rule(node),
                    "pairs": ";".join(f"{a}|{b}:{ok[(a, b)].ptype},pid={ok[(a, b)].pid:.1f}" for a, b in pairs),
                    "applied": bool(apply)})
        if apply:
            for g in comp:
                rest.remove(g)
            merged.append(("or", tuple(comp)))
    return w5.canon(("and", tuple(rest + merged)))


# ============================================================================= io helpers
def atomic_write_tsv(df: pd.DataFrame, path: Path, float_format: str | None = "%.6g") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, float_format=float_format)
    os.replace(tmp, path)
    return path


def atomic_write_parquet(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_parquet(tmp, index=False, compression="zstd")
    os.replace(tmp, path)
    return path


def atomic_write_text(text: str, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)
    return path


def provenance(*, inputs=(), seeds=(), solver="", workers=1, allow_dirty=False):
    from nextgem.run import provenance as _prov

    F1.mkdir(parents=True, exist_ok=True)
    return _prov(MANIFEST, inputs=inputs, seeds=seeds, solver=solver, workers=workers, allow_dirty=allow_dirty,
                 author=AUTHOR)


def curated_rules_path() -> Path:
    return F1 / f"curated_rules_{PARENT}.tsv"


def load_curated_rules() -> dict:
    """reaction id -> curated gene_reaction_rule (only edited reactions)."""
    df = pd.read_csv(curated_rules_path(), sep="\t", dtype=str)
    return dict(zip(df.rxn, df.new_rule))


def apply_arm(model, arm: Arm, rules: dict) -> None:
    """Make ``model`` (a loaded Human1, not inside a context) the arm's parent: TYMS bounds and/or curated GPRs."""
    if arm.tyms:
        model.reactions.get_by_id(TYMS_REACTION).bounds = TYMS_BOUNDS
    if arm.gpr:
        for rid, rule in rules.items():
            r = model.reactions.get_by_id(rid)
            before = {g.id for g in r.genes}
            r.gene_reaction_rule = rule
            if {g.id for g in r.genes} != before:
                raise RuntimeError(f"{rid}: the curated rule changes the gene set")
