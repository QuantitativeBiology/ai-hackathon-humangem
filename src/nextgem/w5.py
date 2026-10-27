"""W5 GPR curation: GPR trees, single-edit enumeration, protected isozyme relations, local KO engine (label-free).

DESIGN.md §9 W5. Everything here is label-free; training labels enter only in ``scripts/w5_score_train.py``
(through ``nextgem.cv.lolo`` / ``labels.fold``).

GPR trees
---------
A GPR is a nested tuple ``(op, children)`` with ``op`` in {"and", "or"} or a gene-id string (the format of
``nextgem.context.compile_gprs``). :func:`canon` flattens nested same-op nodes, removes duplicate children and
sorts children, so equal Boolean structures written differently get one key. Reactions are grouped by the
canonical GPR (a "GPR group"): an edit acts on the gene logic of an enzyme, so it is applied to every reaction
that carries that GPR (*ours*; see README "Deviations": the literal per-reaction reading would score 0 flips
for enzymes that catalyse several individually essential reactions, e.g. the 14 RNR reactions).

Single edits (enumerated per GPR group; node addressed by its path of child indices):
- ``OR2AND``   an OR node becomes an AND of the same children;
- ``AND2OR``   two children i < j of an AND node are merged into ``(ci or cj)`` (the minimal AND -> OR edit;
               for a two-child AND this equals the node-level edit);
- ``RMBRANCH`` child i of an OR node with >= 2 children is removed.

KO semantics (identical to ``nextgem.knockouts``/``cobra.Gene.knock_out``): KO of gene g closes every reaction
whose GPR is False with g off. Only genes of an edited GPR can change their closed-reaction set, so an edit is
re-evaluated on those genes only (DESIGN compute SF4), and identical closed sets share one LP.
"""

from __future__ import annotations

import math
import os
import re
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from nextgem import paths

PARENTS = ("human1", "human2")
SHORT = {"human1": "H1", "human2": "H2"}
LETHAL = 0.01  # binary call: ratio_r < 0.01 (metrics.FBA_CUTOFF, DESIGN §6.1)
ZERO_FLUX_TOL = 1e-9  # zero-flux shortcut (exact; same tolerance as nextgem.w2_engine)
LP_TIMEOUT_S = 60  # GLPK per-LP limit on the warm-started basis (a stalled warm start was observed once)
LP_RETRY_TIMEOUT_S = 300  # retry from scratch (presolve), = nextgem.context.LP_TIME_LIMIT_S
N_TARGETS = 200  # DESIGN W5: curate the top 200
TYMS_REACTION = "MAR04644"  # DESIGN §3.1: made irreversible in optimised H1 pipelines only
AUTHOR = "W5 engineer (Claude Opus)"


def w5_dir() -> Path:
    return paths.results_dir() / "w5"


def tmp_dir() -> Path:
    return w5_dir() / "tmp"


# ============================================================================= GPR trees
def canon(node):
    """Canonical form: flatten nested same-op nodes, drop duplicate children, sort children; collapse 1-child nodes."""
    if isinstance(node, str):
        return node
    op, ch = node
    flat = []
    for c in ch:
        c = canon(c)
        if not isinstance(c, str) and c[0] == op:
            flat.extend(c[1])
        else:
            flat.append(c)
    uniq = sorted(set(flat), key=to_rule)
    if len(uniq) == 1:
        return uniq[0]
    return (op, tuple(uniq))


def to_rule(node, top: bool = True) -> str:
    if isinstance(node, str):
        return node
    op, ch = node
    s = f" {op} ".join(to_rule(c, False) for c in ch)
    return s if top else f"({s})"


def genes_of(node) -> set:
    if isinstance(node, str):
        return {node}
    return set().union(*(genes_of(c) for c in node[1]))


def gpr_true(node, off) -> bool:
    if isinstance(node, str):
        return node not in off
    op, ch = node
    if op == "and":
        return all(gpr_true(c, off) for c in ch)
    return any(gpr_true(c, off) for c in ch)


def parse_rule(rule: str):
    """Parse a cobra gene_reaction_rule string into a canonical tree (None for an empty rule)."""
    rule = rule.strip()
    if not rule:
        return None
    toks = re.findall(r"\(|\)|[^\s()]+", rule)
    pos = 0

    def expr():
        nonlocal pos
        left = term()
        items = [left]
        while pos < len(toks) and toks[pos].lower() == "or":
            pos += 1
            items.append(term())
        return items[0] if len(items) == 1 else ("or", tuple(items))

    def term():
        nonlocal pos
        items = [atom()]
        while pos < len(toks) and toks[pos].lower() == "and":
            pos += 1
            items.append(atom())
        return items[0] if len(items) == 1 else ("and", tuple(items))

    def atom():
        nonlocal pos
        if pos >= len(toks):
            raise ValueError(f"truncated GPR {rule!r}")
        t = toks[pos]
        pos += 1
        if t == "(":
            e = expr()
            if pos >= len(toks) or toks[pos] != ")":
                raise ValueError(f"unbalanced GPR {rule!r}")
            pos += 1
            return e
        if t.lower() in ("and", "or", ")"):
            raise ValueError(f"bad GPR token {t!r} in {rule!r}")
        return t

    tree = expr()
    if pos != len(toks):
        raise ValueError(f"trailing tokens in GPR {rule!r}")
    return canon(tree)


def get_node(tree, path):
    for i in path:
        tree = tree[1][i]
    return tree


def replace_node(tree, path, new):
    if not path:
        return new
    op, ch = tree
    ch = list(ch)
    ch[path[0]] = replace_node(ch[path[0]], path[1:], new)
    return (op, tuple(ch))


def iter_nodes(tree, path=()):
    """Yield (path, node) for every internal node, pre-order."""
    if isinstance(tree, str):
        return
    yield path, tree
    for i, c in enumerate(tree[1]):
        yield from iter_nodes(c, path + (i,))


@dataclass(frozen=True)
class Edit:
    kind: str  # OR2AND | AND2OR | RMBRANCH
    path: tuple
    idx: tuple = ()  # (i, j) for AND2OR, (i,) for RMBRANCH

    def key(self) -> str:
        return f"{self.kind}@{'.'.join(map(str, self.path)) or 'root'}:{','.join(map(str, self.idx))}"


def apply_edit(tree, e: Edit):
    """Apply one edit to a canonical tree; returns the canonical result."""
    node = get_node(tree, e.path)
    op, ch = node
    if e.kind == "OR2AND":
        assert op == "or"
        new = ("and", ch)
    elif e.kind == "AND2OR":
        assert op == "and"
        i, j = e.idx
        merged = ("or", (ch[i], ch[j]))
        rest = [c for k, c in enumerate(ch) if k not in (i, j)]
        new = ("and", tuple(rest + [merged])) if rest else merged
    elif e.kind == "RMBRANCH":
        assert op == "or" and len(ch) >= 2
        (i,) = e.idx
        new = ("or", tuple(c for k, c in enumerate(ch) if k != i))
    else:
        raise ValueError(e.kind)
    return canon(replace_node(tree, e.path, new))


def enumerate_edits(tree) -> list[Edit]:
    """Every single edit of a canonical tree (see module docstring)."""
    out = []
    for path, (op, ch) in iter_nodes(tree):
        if op == "or":
            out.append(Edit("OR2AND", path))
            for i in range(len(ch)):
                out.append(Edit("RMBRANCH", path, (i,)))
        else:
            for i in range(len(ch)):
                for j in range(i + 1, len(ch)):
                    out.append(Edit("AND2OR", path, (i, j)))
    return out


def or_branch_unique_genes(node) -> list[set]:
    """For an OR node: per branch, the genes that are not in every branch (the genes that make branches differ)."""
    sets = [genes_of(c) for c in node[1]]
    common = set.intersection(*sets)
    return [s - common for s in sets]


# ============================================================================= model / gene maps
def load_parent(parent: str, *, tyms_irreversible: bool = False):
    """Generic parent on GLPK (the caller must have set ``cobra.Configuration().solver`` before), M0 medium,
    biomass objective. ``tyms_irreversible`` (H1 only, DESIGN §3.1): MAR04644 bounds set to Human2's (0, 1000)."""
    from nextgem.knockouts import configure
    from nextgem.model import load_model

    m = load_model(parent)
    configure(m)
    if tyms_irreversible:
        if parent != "human1":
            raise ValueError("the TYMS edit is H1-only")
        r = m.reactions.get_by_id(TYMS_REACTION)
        r.bounds = (0.0, 1000.0)
    return m


def gene_table(parent: str) -> pd.DataFrame:
    """Model genes (Ensembl) with symbol, primary Entrez (int or NA), UniProt list, compartments (genes.tsv)."""
    from nextgem.model import GEMS

    g = pd.read_csv(GEMS[parent].genes_tsv, sep="\t", dtype=str).rename(columns={"genes": "gene"})
    g["entrez"] = pd.to_numeric(g.geneEntrezID.str.split(";").str[0], errors="coerce").astype("Int64")
    g["uniprot"] = g.geneUniProtID.fillna("").str.split(";").apply(lambda v: [x for x in v if x])
    g["symbol"] = g.geneSymbols
    return g[["gene", "symbol", "entrez", "uniprot", "compartments", "geneAliases"]].set_index("gene")


# ============================================================================= protected isozyme relations
# notes/apriori_lists_final.md §3, used verbatim (A2). Each entry quotes its text fragment from §3; the test
# tests/test_w5.py::test_protected_tokens_verbatim checks that every quoted fragment occurs in §3.
PROTECTED_S3 = [
    {"text": "`MAR04363` ENO1/ENO3/ENO2", "reactions": ["MAR04363"], "genes": ["ENO1", "ENO3", "ENO2"]},
    {"text": "`MAR04375`, `MAR04356`, `MAR08762` ALDOC/ALDOB/ALDOA", "reactions": ["MAR04375", "MAR04356", "MAR08762"],
     "genes": ["ALDOC", "ALDOB", "ALDOA"]},
    {"text": "`MAR04281`, `MAR04280`, `MAR01453` LDHB/LDHA, LDHB/LDHAL6B", "reactions": ["MAR04281", "MAR04280", "MAR01453"],
     "genes": ["LDHB", "LDHA", "LDHAL6B"]},
    {"text": "`MAR04319` HKDC1/HK1/HK2/HK3", "reactions": ["MAR04319"], "genes": ["HKDC1", "HK1", "HK2", "HK3"]},
    {"text": "`MAR04379`, `MAR04301`, `MAR04567` PFKP/PFKL/PFKM", "reactions": ["MAR04379", "MAR04301", "MAR04567"],
     "genes": ["PFKP", "PFKL", "PFKM"]},
    {"text": "`MAR04171`, `MAR04193`, `MAR04210` PKM/PKLR", "reactions": ["MAR04171", "MAR04193", "MAR04210"],
     "genes": ["PKM", "PKLR"]},
    {"text": "`MAR04371` PGAM2/PGAM1/BPGM/PGAM4", "reactions": ["MAR04371"], "genes": ["PGAM2", "PGAM1", "BPGM", "PGAM4"]},
    {"text": "`MAR03829` GOT1/GOT1L1", "reactions": ["MAR03829"], "genes": ["GOT1", "GOT1L1"]},
    {"text": "`MAR03912`, `MAR06725` GOT1/GOT2", "reactions": ["MAR03912", "MAR06725"], "genes": ["GOT1", "GOT2"]},
    {"text": "`MAR04340` MTHFD2/MTHFD1L/MTHFD2L", "reactions": ["MAR04340"], "genes": ["MTHFD2", "MTHFD1L", "MTHFD2L"]},
    {"text": "IMPDH1/IMPDH2", "reactions": [], "genes": ["IMPDH1", "IMPDH2"]},
    {"text": "CTPS1/CTPS2", "reactions": [], "genes": ["CTPS1", "CTPS2"]},
    {"text": "NME1/NME2", "reactions": [], "genes": ["NME1", "NME2"]},
    {"text": "ADSS1/ADSS2", "reactions": [], "genes": ["ADSS1", "ADSS2"]},
    {"text": "PRPS1/PRPS2", "reactions": [], "genes": ["PRPS1", "PRPS2"]},
    {"text": "CKB/CKM (and CKMT1A/CKMT1B)", "reactions": [], "genes": ["CKB", "CKM", "CKMT1A", "CKMT1B"]},
    {"text": "`(COX4I1 or COX4I2)`", "reactions": [], "genes": ["COX4I1", "COX4I2"]},
    {"text": "`(COX6A1 or COX6A2)`", "reactions": [], "genes": ["COX6A1", "COX6A2"]},
    {"text": "`(COX6B1 or COX6B2)`", "reactions": [], "genes": ["COX6B1", "COX6B2"]},
    {"text": "`(COX7A1 or COX7A2)`", "reactions": [], "genes": ["COX7A1", "COX7A2"]},
    {"text": "ATP6V1B1/ATP6V1B2", "reactions": [], "genes": ["ATP6V1B1", "ATP6V1B2"]},
    # "ATP6V0A1-A4": A3 is HGNC TCIRG1 (alias ATP6V0A3); resolved through genes.tsv aliases
    {"text": "ATP6V0A1-A4", "reactions": [], "genes": ["ATP6V0A1", "ATP6V0A2", "ATP6V0A3", "ATP6V0A4"]},
    {"text": "GLUD1/GLUD2", "reactions": [], "genes": ["GLUD1", "GLUD2"]},
    {"text": "ME1/ME2/ME3", "reactions": [], "genes": ["ME1", "ME2", "ME3"]},
    {"text": "ACSS3/ACSS2/ACSS1", "reactions": [], "genes": ["ACSS3", "ACSS2", "ACSS1"]},
    {"text": "SULT2B1/SULT1A3 (`MAR01919`)", "reactions": ["MAR01919"], "genes": ["SULT2B1", "SULT1A3"]},
]
APRIORI_FINAL = paths.PROJECT_ROOT / "notes" / "apriori_lists_final.md"


def section3_text() -> str:
    t = APRIORI_FINAL.read_text()
    a = t.index("# 3. Final isozyme exceptions")
    b = t.index("# 4. Final metabolic-drug list")
    return t[a:b]


@dataclass
class Protected:
    reactions: set
    groups: list  # list of sets of Ensembl gene ids
    unresolved: list

    def or_node_protected(self, reactions, node) -> str | None:
        """Reason string if this OR node is protected (in a listed reaction, or >= 2 branches carry genes of one
        protected group), else None."""
        hit = sorted(set(reactions) & self.reactions)
        if hit:
            return f"listed reaction {hit[0]}"
        br = [genes_of(c) for c in node[1]]
        for k, grp in enumerate(self.groups):
            if sum(1 for b in br if b & grp) >= 2:
                return f"protected group {k}"
        return None


def resolve_protected(parent: str) -> Protected:
    gt = gene_table(parent)
    by_sym: dict[str, set] = {}
    for gid, row in gt.iterrows():
        for s in [row.symbol] + (str(row.geneAliases).split(";") if isinstance(row.geneAliases, str) else []):
            if isinstance(s, str) and s:
                by_sym.setdefault(s.strip(), set()).add(gid)
    groups, unresolved, rx = [], [], set()
    for ent in PROTECTED_S3:
        rx |= set(ent["reactions"])
        grp = set()
        for s in ent["genes"]:
            ids = by_sym.get(s, set())
            # a symbol match wins over alias matches
            exact = {g for g in ids if gt.loc[g, "symbol"] == s}
            ids = exact or ids
            if not ids:
                unresolved.append(s)
            grp |= ids
        groups.append(grp)
    return Protected(rx, groups, unresolved)


# ============================================================================= local KO engine
class NonOptimal(RuntimeError):
    pass


def rounded(r: float) -> float:
    return float(min(1.0, max(0.0, round(r, 6))) + 0.0)


@dataclass
class KOEngine:
    """Ratios of closed-reaction sets on one configured model (GLPK), with caching and the exact zero-flux shortcut."""

    model: object
    wt: float = 0.0
    flux: dict = field(default_factory=dict)
    cache: dict = field(default_factory=dict)
    n_lp: int = 0
    n_retry: int = 0

    def __post_init__(self):
        sol = self.model.optimize()
        if sol.status != "optimal" or sol.objective_value <= 1e-6:
            raise NonOptimal(f"wild type: {sol.status} {sol.objective_value}")
        self.wt = float(sol.objective_value)
        self.flux = sol.fluxes.to_dict()

    def ratio(self, closed) -> float:
        closed = frozenset(closed)
        if not closed:
            return 1.0
        if all(abs(self.flux[r]) <= ZERO_FLUX_TOL for r in closed):
            return 1.0
        if closed in self.cache:
            return self.cache[closed]
        v = self._solve(closed)
        self.cache[closed] = v
        return v

    def _solve(self, closed) -> float:
        """One KO LP. GLPK time limit LP_TIMEOUT_S on the warm-started basis; a non-optimal, non-infeasible status is
        retried once from scratch (presolve on, limit LP_RETRY_TIMEOUT_S), then raises (never silently ratio 1)."""
        from optlang.interface import INFEASIBLE, OPTIMAL

        m = self.model
        cfg = m.solver.configuration
        st = None
        for attempt, (tmo, pre) in enumerate(((LP_TIMEOUT_S, "auto"), (LP_RETRY_TIMEOUT_S, True))):
            cfg.timeout = tmo
            cfg.presolve = pre
            with m:
                for r in closed:
                    m.reactions.get_by_id(r).bounds = (0.0, 0.0)
                m.slim_optimize(error_value=math.nan)
                st = m.solver.status
                val = m.solver.objective.value if st == OPTIMAL else None
            if st in (OPTIMAL, INFEASIBLE):
                break
            self.n_retry += 1
        cfg.presolve = "auto"
        cfg.timeout = LP_TIMEOUT_S
        if st == OPTIMAL:
            g = val
        elif st == INFEASIBLE:
            g = 0.0
        else:
            raise NonOptimal(f"KO status {st} after retry for {sorted(closed)[:5]}")
        self.n_lp += 1
        return rounded(g / self.wt)


def closures(gprs: dict) -> dict:
    """gene -> frozenset of reactions its KO closes, for reaction -> tree ``gprs``."""
    out: dict[str, set] = {}
    for rid, t in gprs.items():
        for g in genes_of(t):
            if not gpr_true(t, {g}):
                out.setdefault(g, set()).add(rid)
    return {g: frozenset(v) for g, v in out.items()}


def model_gprs(model) -> dict:
    """reaction -> canonical tree, for reactions with a GPR."""
    out = {}
    for r in model.reactions:
        t = parse_rule(r.gene_reaction_rule)
        if t is not None:
            out[r.id] = t
    return out


def groups_of(gprs: dict) -> dict:
    """canonical rule -> sorted reaction ids."""
    g: dict[str, list] = {}
    for rid, t in gprs.items():
        g.setdefault(to_rule(t), []).append(rid)
    return {k: sorted(v) for k, v in sorted(g.items())}


def edited_closure(base: frozenset, rxns: frozenset, new_tree, gene) -> frozenset:
    """Closed set of ``gene`` after the reactions ``rxns`` get ``new_tree`` (None: gene not in the new tree)."""
    closes = new_tree is not None and gene in genes_of(new_tree) and not gpr_true(new_tree, {gene})
    return (base - rxns) | (rxns if closes else frozenset())


def apply_tree_edits(gprs: dict, new_trees: dict) -> dict:
    """Copy of reaction->tree with ``new_trees`` (reaction -> canonical tree) substituted."""
    out = dict(gprs)
    out.update(new_trees)
    return out


def full_sweep(engine: KOEngine, gprs: dict, genes) -> pd.Series:
    """ratio_r for every gene in ``genes`` under ``gprs`` (genes in no GPR close nothing -> 1)."""
    cl = closures(gprs)
    return pd.Series({g: engine.ratio(cl.get(g, frozenset())) for g in genes}, name="ratio_r")


_POOL_ENGINE: list = [None]


def _pool_solve(cs):
    return cs, _POOL_ENGINE[0].ratio(cs)


def local_sweep(engine: KOEngine, gprs: dict, base: pd.Series, new_trees: dict, workers: int = 1) -> pd.Series:
    """ratio_r of every gene in ``base`` after substituting ``new_trees`` (reaction -> tree) into ``gprs``; only the
    genes of edited GPRs are re-evaluated (exact: no other gene's closed set changes), the rest keep ``base``."""
    import multiprocessing as mp

    out = base.copy()
    if not new_trees:
        return out
    gp2 = apply_tree_edits(gprs, new_trees)
    touched = set().union(*(genes_of(gprs[r]) for r in new_trees))
    cl = closures(gp2)
    need = {g: cl.get(g, frozenset()) for g in touched}
    todo = sorted({cs for cs in need.values() if cs and cs not in engine.cache
                   and not all(abs(engine.flux[r]) <= ZERO_FLUX_TOL for r in cs)}, key=lambda s: (len(s), sorted(s)))
    if workers > 1 and len(todo) > 20:
        _POOL_ENGINE[0] = engine
        with mp.get_context("fork").Pool(workers) as pool:
            for cs, v in pool.imap_unordered(_pool_solve, todo, chunksize=4):
                engine.cache[cs] = v
    for g, cs in need.items():
        if g in out.index:
            out[g] = engine.ratio(cs)
    return out


def e_by_entrez(gene_ids, entrez, ratio, U) -> pd.Series:
    """Per-U-gene score e = 1 - ratio_r, Entrez dedup exactly as ``nextgem.w1.generic_table`` (first row per Entrez in
    results/baseline order); genes absent from the parent score 0 (§3.3)."""
    df = pd.DataFrame({"entrez": entrez, "ratio_r": ratio}).dropna(subset=["entrez"])
    df = df.assign(entrez=df.entrez.astype(int)).drop_duplicates("entrez").set_index("entrez").ratio_r
    return (1.0 - df.reindex(sorted(U))).fillna(0.0).round(6)


# ============================================================================= 2 x 2 factorial (DESIGN W5 "Ablation")
FACTORIAL_ARMS = [("generic_gpr", "generic_constraints"), ("curated_gpr", "generic_constraints"),
                  ("generic_gpr", "best_constraints"), ("curated_gpr", "best_constraints")]


def curated_for_extracted(extracted_model, curated_trees_: dict):
    """Apply curated GPRs to an *already-extracted* model (DESIGN W5: no re-extraction): every reaction present in the
    extracted model whose generic GPR was curated gets the curated tree; reactions removed by extraction stay removed.
    Returns reaction -> tree for the extracted model (feed to ``local_sweep`` / ``full_sweep``)."""
    gp = model_gprs(extracted_model)
    return apply_tree_edits(gp, {r: t for r, t in curated_trees_.items() if r in gp})


def best_constraint_pipeline(*_a, **_k):
    """Stub: the 'best constraint pipeline' arm needs the W2/W3/W4 shortlists (not yet available)."""
    raise NotImplementedError("constraint arm of the W5 2x2 factorial: awaiting the W2/W3/W4 shortlists (DESIGN §10.1)")


# ============================================================================= provenance helper
def run_main(main, *, inputs=(), seeds=(), solver="glpk", workers=1, reads_labels=False, manifest=None) -> int:
    """Run ``main(prov)`` under ``nextgem.run.provenance`` (DESIGN §11.4). The W5 code is run uncommitted
    (the brief forbids commits), so ``allow_dirty=True`` and the SHA-256 of every W5 code file is recorded."""
    import hashlib

    from nextgem.run import provenance

    code = sorted(list((paths.SRC_DIR).glob("w5*.py")) + list(paths.SCRIPTS_DIR.glob("w5_*.py")))
    code_sha = {str(p.relative_to(paths.PROJECT_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in code}
    with provenance(manifest or (w5_dir() / "MANIFEST.tsv"), inputs=inputs, seeds=seeds, solver=solver,
                    workers=workers, allow_dirty=True, author=AUTHOR) as prov:
        prov.note(w5_code_sha256=code_sha)
        if reads_labels:
            prov.note(label_access="labels.train()/labels.fold(k) via nextgem.cv.lolo only", **_frozen_label_hashes())
        rc = main(prov)
        if rc:
            raise SystemExit(rc)
    return 0


def _frozen_label_hashes() -> dict:
    out = {}
    f = paths.frozen_tsv()
    for ln in f.read_text().splitlines()[1:]:
        p = ln.split("\t")
        if len(p) >= 4 and p[1] == "label_file":
            out[f"frozen_{p[0]}"] = p[3]
    return out


def atomic_tsv(df: pd.DataFrame, path: Path, float_format="%.6g") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, float_format=float_format)
    with open(tmp, "rb+") as f:
        os.fsync(f.fileno())
    os.replace(tmp, path)
    return path


def atomic_parquet(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_parquet(tmp, index=False, compression="zstd")
    os.replace(tmp, path)
    return path


class Timer:
    def __init__(self):
        self.t0, self.c0 = time.time(), time.process_time()

    def lap(self) -> dict:
        return {"wall_s": round(time.time() - self.t0, 1), "cpu_s_self": round(time.process_time() - self.c0, 1)}
