"""W5 evidence (label-free): Ensembl 116 paralogs, CORUM 5.1 / Complex Portal complexes, ProCan / Nusinow protein
co-abundance and detection, DepMap 26Q3 expression; the DESIGN W5 edit rules and the evidence-ID validator.

Rules (DESIGN §9 W5 table, applied to the top-200 target GPR groups only; protected OR relations of
notes/apriori_lists_final.md §3 are never edited):

- OR -> AND: every cross-branch gene pair (genes that distinguish the branches) is **not** an Ensembl paralog pair,
  **and** co-members of >= 1 CORUM or Complex Portal complex, **and** has protein co-abundance Pearson r >= 0.5 over
  >= 30 lines (ProCan if it has >= 30 shared lines, else Nusinow).
- AND -> OR: AND-ed sibling genes that are Ensembl paralogs **and** share no complex are merged into one OR (a
  clique of such pairs is merged; non-clique components are logged, not merged).
- Remove OR branch (pan-cancer absent; the UniProt "unrelated function" path needs a UniProt file that is not on
  disk, see notes/USER_DOWNLOADS.md): a leaf branch gene that is a soluble enzyme (proxy, see
  :func:`EvidenceDB.soluble`), never detected in either proteomics set, and has TPM < 1 in >= 95% of lines.

Evidence IDs (each checked against the raw local files by :func:`validate`, an independent code path):
``CORUM:<complex_id>``, ``CPX:<CPX-n>``, ``ENSEMBL116_PARALOG:<ENSG>|<ENSG>``, ``NOT_ENSEMBL116_PARALOG:<ENSG>|<ENSG>``,
``NO_SHARED_COMPLEX:<ENSG>|<ENSG>``, ``PROCAN_COAB:<acc>|<acc>`` / ``NUSINOW_COAB:<acc>|<acc>``,
``NOT_DETECTED_PROCAN:<acc>``, ``NOT_DETECTED_NUSINOW:<acc>``, ``DEPMAP26Q3_EXPR_TPMlt1:<entrez>``.
"""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from nextgem import paths
from nextgem import w5

COAB_R = 0.5
COAB_MIN_N = 30
COESS_R = 0.4  # W5b
COESS_MIN_N = 30  # *ours*: same minimum as co-abundance
ABSENT_TPM_FRAC = 0.95
EXPRESSED_LOG = 1.0  # log2(TPM+1) < 1  <=>  TPM < 1

PARALOG_PID_SENSITIVITY = 30.0  # W5a_pid30 sensitivity: W1's paralog convention (%id >= 30), declared before scoring
PARALOGS = "paralogs/ensembl_human_paralogs.tsv"
CORUM = "complexes/corum_humanComplexes.txt"
CPX = "complexes/complexportal_9606.tsv"
PROCAN = "proteomics_goncalves/Protein_matrix_averaged_20250211.tsv"
NUSINOW = "proteomics_nusinow/protein_quant_current_normalized.csv.gz"
EXPR_RAW = "depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
DENYLIST = ("1-s2.0-S240547122100329X-mmc2.xlsx", "41588_2021_967_MOESM2_ESM.xlsx")  # De Kegel 2021, Ito 2021 (§11.2)


def data(rel: str) -> Path:
    return paths.data_root() / rel


def evidence_inputs() -> list[Path]:
    return [data(x) for x in (PARALOGS, CORUM, CPX, PROCAN, NUSINOW)] + [paths.results_dir() / "w1" / "tmp" / "expr.parquet"]


_ACC = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})")


def acc(x: str) -> str | None:
    m = _ACC.match(str(x).strip())
    return m.group(1) if m else None


def pair_key(a: str, b: str) -> str:
    return "|".join(sorted((a, b)))


# ============================================================================= evidence database
@dataclass
class EvidenceDB:
    genes: pd.DataFrame  # union of H1/H2 gene tables (index Ensembl)
    paralog: dict = field(default_factory=dict)  # frozenset(g1, g2) -> (type, min %id)
    complexes: dict = field(default_factory=dict)  # gene -> set of complex evidence ids
    procan: pd.DataFrame | None = None  # lines x acc
    nusinow: pd.DataFrame | None = None  # lines x acc
    expr_frac_low: pd.Series | None = None  # entrez -> fraction of lines with TPM < 1
    expr_n: pd.Series | None = None
    transport_genes: set = field(default_factory=set)
    unmapped_complex_members: int = 0
    min_pid: float = 0.0  # paralog = Ensembl pair with min(%id both directions) >= min_pid (0 = literal DESIGN rule)

    # ---------------------------------------------------------------- pair facts
    def is_paralog(self, a, b) -> bool:
        v = self.paralog.get(frozenset((a, b)))
        return v is not None and v[1] >= self.min_pid

    def shared_complexes(self, a, b) -> list[str]:
        return sorted(self.complexes.get(a, set()) & self.complexes.get(b, set()))

    def protein(self, gene, src) -> str | None:
        mat = self.procan if src == "PROCAN" else self.nusinow
        for u in self.genes.uniprot.get(gene, []) or []:
            if u in mat.columns and mat[u].notna().any():
                return u
        return None

    def coab(self, a, b) -> dict:
        """ProCan if >= 30 shared lines, else Nusinow; dict(source, r, n, id) or source None."""
        for src in ("PROCAN", "NUSINOW"):
            mat = self.procan if src == "PROCAN" else self.nusinow
            pa, pb = self.protein(a, src), self.protein(b, src)
            if pa is None or pb is None:
                continue
            x = mat[[pa, pb]].dropna()
            if len(x) >= COAB_MIN_N:
                r = float(np.corrcoef(x[pa], x[pb])[0, 1])
                return {"source": src, "r": r, "n": len(x), "id": f"{src}_COAB:{pair_key(pa, pb)}"}
        return {"source": None, "r": np.nan, "n": 0, "id": None}

    # ---------------------------------------------------------------- gene facts
    def detected(self, gene) -> dict:
        out = {}
        for src in ("PROCAN", "NUSINOW"):
            out[src] = self.protein(gene, src) is not None
        return out

    def soluble(self, gene) -> bool:
        """Label-free proxy (no UniProt topology on disk): Human-GEM genes.tsv location has no 'Cell membrane' and the
        gene catalyses no transport reaction (metabolites in > 1 compartment) in either parent."""
        comp = str(self.genes.compartments.get(gene, "") or "")
        return "Cell membrane" not in comp and gene not in self.transport_genes

    def pan_cancer_absent(self, gene) -> dict:
        ent = self.genes.entrez.get(gene)
        det = self.detected(gene)
        frac = float(self.expr_frac_low.get(int(ent), np.nan)) if pd.notna(ent) else np.nan
        sol = self.soluble(gene)
        ok = sol and not det["PROCAN"] and not det["NUSINOW"] and frac >= ABSENT_TPM_FRAC
        ids = []
        if ok:
            accs = self.genes.uniprot.get(gene, []) or []
            ids = [f"NOT_DETECTED_{src}:{a}" for a in accs for src in ("PROCAN", "NUSINOW")]
            ids.append(f"DEPMAP26Q3_EXPR_TPMlt1:{int(ent)}")
        return {"absent": bool(ok), "soluble_proxy": sol, "procan": det["PROCAN"], "nusinow": det["NUSINOW"],
                "frac_tpm_lt1": frac, "n_expr_lines": int(self.expr_n.get(int(ent), 0)) if pd.notna(ent) else 0,
                "ids": ids}


def _read_paralogs(model_genes: set) -> dict:
    cols = ["Gene stable ID", "Human paralogue gene stable ID", "Paralogue %id. target Human gene identical to query gene",
            "Paralogue %id. query gene identical to target Human gene", "Human paralogue homology type"]
    out = {}
    for ch in pd.read_csv(data(PARALOGS), sep="\t", usecols=cols, chunksize=500_000, dtype=str):
        ch = ch[ch[cols[0]].isin(model_genes) & ch[cols[1]].isin(model_genes)]
        for a, b, p1, p2, t in ch.itertuples(index=False):
            if a == b:
                continue
            k = frozenset((a, b))
            pid = min(float(p1), float(p2))
            if k not in out or pid > out[k][1]:
                out[k] = (t, pid)
    return out


def _read_complexes(genes: pd.DataFrame) -> tuple[dict, int]:
    by_acc: dict[str, set] = {}
    by_sym: dict[str, set] = {}
    for g, row in genes.iterrows():
        for u in row.uniprot or []:
            by_acc.setdefault(u, set()).add(g)
        if isinstance(row.symbol, str):
            by_sym.setdefault(row.symbol, set()).add(g)
    out: dict[str, set] = {}
    unmapped = 0
    c = pd.read_csv(data(CORUM), sep="\t", dtype=str, usecols=["complex_id", "organism", "subunits_uniprot_id",
                                                                "subunits_gene_name"])
    c = c[c.organism.str.contains("Human", na=False)]
    for cid, _, us, gs in c.itertuples(index=False):
        accs = [acc(u) for u in str(us).split(";")]
        syms = str(gs).split(";")
        for i, a in enumerate(accs):
            gg = by_acc.get(a, set()) if a else set()
            if not gg and i < len(syms):
                gg = by_sym.get(syms[i].strip(), set())
            if not gg:
                unmapped += 1
            for g in gg:
                out.setdefault(g, set()).add(f"CORUM:{cid}")
    p = pd.read_csv(data(CPX), sep="\t", dtype=str, usecols=["#Complex ac", "Expanded participant list"])
    for cid, lst in p.itertuples(index=False):
        for tok in str(lst).split("|"):
            a = acc(tok)
            for g in by_acc.get(a, set()) if a else set():
                out.setdefault(g, set()).add(f"CPX:{cid}")
    return out, unmapped


def _read_procan(accs: set) -> pd.DataFrame:
    hdr = pd.read_csv(data(PROCAN), sep="\t", nrows=0).columns
    keep = [c for c in hdr if c in accs]
    df = pd.read_csv(data(PROCAN), sep="\t", skiprows=[1, 2], usecols=[hdr[0], hdr[1]] + keep)
    df = df.rename(columns={hdr[0]: "model_name", hdr[1]: "model_id"}).set_index("model_id").drop(columns="model_name")
    return df.astype(float)


def _read_nusinow(accs: set) -> pd.DataFrame:
    df = pd.read_csv(data(NUSINOW))
    df = df[df.Uniprot_Acc.isin(accs)]
    ch = [c for c in df.columns if "_TenPx" in c and not c.endswith("_Peptides")]
    m = df.set_index("Uniprot_Acc")[ch]
    m = m[~m.index.duplicated()]
    line = ["_".join(c.split("_")[:-1]) for c in ch]  # CCLE name (replicate channels averaged)
    m.columns = line
    m = m.T.groupby(level=0).mean()
    return m.astype(float)


def _expression():
    e = pd.read_parquet(paths.results_dir() / "w1" / "tmp" / "expr.parquet").set_index("ModelID")
    e.columns = e.columns.astype(int)
    n = e.notna().sum()
    frac = (e < EXPRESSED_LOG).sum() / n.replace(0, np.nan)
    return frac, n


def transport_genes(model) -> set:
    out = set()
    for r in model.reactions:
        if r.genes and len({m.compartment for m in r.metabolites}) > 1:
            out |= {g.id for g in r.genes}
    return out


def build_db(models: dict) -> EvidenceDB:
    """``models``: parent -> cobra model (for gene tables and transport-reaction genes)."""
    gts = [w5.gene_table(p) for p in models]
    genes = pd.concat(gts)
    genes = genes[~genes.index.duplicated(keep="last")]  # H2 rows win (loaded last)
    mg = set(genes.index)
    db = EvidenceDB(genes=genes)
    db.paralog = _read_paralogs(mg)
    db.complexes, db.unmapped_complex_members = _read_complexes(genes)
    accs = {u for v in genes.uniprot for u in (v or [])}
    db.procan = _read_procan(accs)
    db.nusinow = _read_nusinow(accs)
    db.expr_frac_low, db.expr_n = _expression()
    for m in models.values():
        db.transport_genes |= transport_genes(m)
    return db


# ============================================================================= rules
def or_pairs(node) -> list[tuple]:
    """Cross-branch gene pairs of an OR node (genes that distinguish the branches); [] if some branch has none."""
    uniq = w5.or_branch_unique_genes(node)
    if any(not u for u in uniq):
        return []
    pairs = set()
    for i, j in combinations(range(len(uniq)), 2):
        for a in uniq[i]:
            for b in uniq[j]:
                if a != b:
                    pairs.add(tuple(sorted((a, b))))
    return sorted(pairs)


def pair_facts(db: EvidenceDB, a, b) -> dict:
    co = db.coab(a, b)
    sc = db.shared_complexes(a, b)
    par = db.is_paralog(a, b)
    ptype, pid = db.paralog.get(frozenset((a, b)), (None, np.nan))
    return {"a": a, "b": b, "paralog": par, "paralog_type": ptype, "paralog_min_pid": pid, "shared_complexes": sc,
            "coab_source": co["source"], "coab_r": co["r"], "coab_n": co["n"], "coab_id": co["id"]}


def or2and_ok(facts: list[dict]) -> bool:
    return bool(facts) and all((not f["paralog"]) and f["shared_complexes"] and f["coab_source"] is not None
                               and f["coab_r"] >= COAB_R for f in facts)


def or2and_ids(facts) -> list[str]:
    ids = []
    for f in facts:
        pk = pair_key(f["a"], f["b"])
        ids += [f"NOT_ENSEMBL116_PARALOG:{pk}", f"{f['shared_complexes'][0]}@{pk}", f["coab_id"]]
    return ids


def and2or_ok(f: dict) -> bool:
    return f["paralog"] and not f["shared_complexes"]


@dataclass
class Curation:
    edits: list = field(default_factory=list)  # applied edit records
    or_nodes: list = field(default_factory=list)  # every OR node examined (for W5b and the report)
    and_pairs: list = field(default_factory=list)  # every AND leaf-sibling pair examined
    removals_blocked: list = field(default_factory=list)


def curate_tree(tree, rxns, db: EvidenceDB, prot: w5.Protected, cur: Curation, ctx: dict):
    """W5a on one canonical tree (bottom-up). Returns the new canonical tree; appends records to ``cur``."""
    if isinstance(tree, str):
        return tree
    op, ch = tree
    ch = tuple(curate_tree(c, rxns, db, prot, cur, ctx) for c in ch)
    node = (op, ch)
    if op == "or":
        reason = prot.or_node_protected(rxns, node)
        # remove OR branch: pan-cancer absent leaf genes
        keep, removed = list(ch), []
        for c in ch:
            if isinstance(c, str):
                pa = db.pan_cancer_absent(c)
                if pa["absent"]:
                    if reason:
                        cur.removals_blocked.append({**ctx, "gene": c, "protected": reason, **{k: v for k, v in pa.items() if k != "ids"}})
                    else:
                        removed.append((c, pa))
        if removed and len(removed) < len(ch):
            for c, pa in removed:
                keep.remove(c)
                cur.edits.append({**ctx, "kind": "RMBRANCH", "genes": c, "evidence_ids": ";".join(pa["ids"]),
                                  "detail": f"soluble_proxy={pa['soluble_proxy']};procan={pa['procan']};nusinow={pa['nusinow']};"
                                            f"frac_tpm_lt1={pa['frac_tpm_lt1']:.4f};n_lines={pa['n_expr_lines']}",
                                  "old_node": w5.to_rule(node), "rule": "pan-cancer absent"})
        elif removed:
            cur.removals_blocked.append({**ctx, "gene": ";".join(c for c, _ in removed), "protected": "all branches absent"})
        node = w5.canon(("or", tuple(keep)))
        if isinstance(node, str) or node[0] != "or":
            return node
        pairs = or_pairs(node)
        facts = [pair_facts(db, a, b) for a, b in pairs]
        ok = or2and_ok(facts)
        rec = {**ctx, "node": w5.to_rule(node), "n_pairs": len(pairs), "protected": reason or "",
               "all_nonparalog": bool(facts) and all(not f["paralog"] for f in facts),
               "all_cocomplex": bool(facts) and all(bool(f["shared_complexes"]) for f in facts),
               "all_coab": bool(facts) and all(f["coab_source"] is not None and f["coab_r"] >= COAB_R for f in facts),
               "min_coab_r": min((f["coab_r"] for f in facts), default=np.nan),
               "pairs": ";".join(f"{f['a']}|{f['b']}" for f in facts), "W5a_OR2AND": bool(ok and not reason)}
        cur.or_nodes.append(rec)
        if ok and not reason:
            cur.edits.append({**ctx, "kind": "OR2AND", "genes": ";".join(sorted(w5.genes_of(node))),
                              "evidence_ids": ";".join(or2and_ids(facts)),
                              "detail": ";".join(f"{f['a']}|{f['b']}:r={f['coab_r']:.3f},n={f['coab_n']},{f['coab_source']}" for f in facts),
                              "old_node": w5.to_rule(node), "rule": "complex+non-paralog+co-abundance"})
            return w5.canon(("and", node[1]))
        return node
    # AND node: merge cliques of AND-ed leaf siblings that are paralogs and share no complex
    leaves = [c for c in ch if isinstance(c, str)]
    ok_edges = {}
    for a, b in combinations(sorted(leaves), 2):
        f = pair_facts(db, a, b)
        if f["paralog"]:
            cur.and_pairs.append({**ctx, "a": a, "b": b, "paralog_type": f["paralog_type"], "paralog_min_pid": f["paralog_min_pid"],
                                  "shared_complexes": ";".join(f["shared_complexes"]), "W5a_AND2OR": and2or_ok(f)})
        if and2or_ok(f):
            ok_edges[(a, b)] = f
    if not ok_edges:
        return node
    adj: dict[str, set] = {}
    for a, b in ok_edges:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    seen, comps = set(), []
    for v in sorted(adj):
        if v in seen:
            continue
        stack, comp = [v], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(adj[x] - comp)
        seen |= comp
        comps.append(sorted(comp))
    rest = [c for c in ch]
    merged = []
    for comp in comps:
        clique = all((a, b) in ok_edges for a, b in combinations(comp, 2))
        if not clique:
            cur.removals_blocked.append({**ctx, "gene": ";".join(comp), "protected": "AND2OR component not a clique"})
            continue
        for g in comp:
            rest.remove(g)
        merged.append(("or", tuple(comp)))
        ids = []
        for a, b in combinations(comp, 2):
            ids += [f"ENSEMBL116_PARALOG:{pair_key(a, b)}", f"NO_SHARED_COMPLEX:{pair_key(a, b)}"]
        cur.edits.append({**ctx, "kind": "AND2OR", "genes": ";".join(comp), "evidence_ids": ";".join(ids),
                          "detail": ";".join(f"{a}|{b}:{ok_edges[(a, b)]['paralog_type']},pid={ok_edges[(a, b)]['paralog_min_pid']:.1f}"
                                             for a, b in combinations(comp, 2)),
                          "old_node": w5.to_rule(node), "rule": "paralogs+no shared complex"})
    return w5.canon(("and", tuple(rest + merged)))


def w5b_candidates(tree, rxns, db: EvidenceDB, prot: w5.Protected, path=()) -> list[dict]:
    """OR nodes of a (W5a-curated) tree eligible for the W5b co-essentiality rule: not protected, with cross-branch
    pairs. ``anchored`` = all pairs non-paralog and co-complex (primary W5b); ``loose`` = all pairs non-paralog."""
    out = []
    for p, node in w5.iter_nodes(tree, path):
        if node[0] != "or" or prot.or_node_protected(rxns, node):
            continue
        pairs = or_pairs(node)
        if not pairs:
            continue
        facts = [pair_facts(db, a, b) for a, b in pairs]
        nonpar = all(not f["paralog"] for f in facts)
        cocx = all(bool(f["shared_complexes"]) for f in facts)
        if not nonpar:
            continue
        out.append({"path": p, "pairs": pairs, "anchored": cocx, "loose": True,
                    "coab_ok": or2and_ok(facts)})
    return out


def convert_or_nodes(tree, paths_):
    """OR -> AND at every path in ``paths_`` (applied on the raw tree, then canonicalised once)."""
    paths_ = set(map(tuple, paths_))

    def rec(node, p):
        if isinstance(node, str):
            return node
        op, ch = node
        ch = tuple(rec(c, p + (i,)) for i, c in enumerate(ch))
        if p in paths_:
            assert op == "or"
            op = "and"
        return (op, ch)

    return w5.canon(rec(tree, ()))


# ============================================================================= validator (independent code path)
def validate(edits: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """Check every evidence ID of ``edits`` against the raw local files by plain text scanning.

    Returns one row per (edit, evidence id) with ``ok`` and a reason. Never uses EvidenceDB."""
    ids = sorted({i for v in edits.evidence_ids.fillna("") for i in v.split(";") if i})
    for f in DENYLIST:
        if any(f in i for i in ids):
            raise AssertionError(f"denylisted source {f} cited")
    need_par = {i.split(":", 1)[1] for i in ids if "PARALOG:" in i}
    par_rows = set()
    with open(data(PARALOGS)) as fh:
        next(fh)
        for ln in fh:
            a, b = ln.split("\t", 2)[:2]
            k = "|".join(sorted((a, b)))
            if k in need_par:
                par_rows.add(k)
    corum = {}
    cr = pd.read_csv(data(CORUM), sep="\t", dtype=str, usecols=["complex_id", "subunits_uniprot_id", "subunits_gene_name"])
    for cid, us, gs in cr.fillna("").itertuples(index=False):
        corum[cid] = (set(x.strip() for x in us.split(";")), set(x.strip() for x in gs.split(";")))
    cpx = {}
    cp = pd.read_csv(data(CPX), sep="\t", dtype=str, usecols=["#Complex ac", "Expanded participant list"])
    for cid, lst in cp.fillna("").itertuples(index=False):
        cpx[cid] = {acc(t) for t in lst.split("|")} - {None}
    procan_cols = set(open(data(PROCAN)).readline().rstrip("\n").split("\t"))
    nus_acc = set()
    with gzip.open(data(NUSINOW), "rt", newline="") as fh:
        rd = csv.reader(fh)
        ia = next(rd).index("Uniprot_Acc")
        for row in rd:
            nus_acc.add(row[ia].strip())
    expr_hdr = open(data(EXPR_RAW)).readline()
    expr_entrez = set(re.findall(r"\((\d+)\)", expr_hdr))

    def accs_of(g):
        return set(genes.uniprot.get(g, []) or []), genes.symbol.get(g)

    def complex_contains(cid, g):
        a, s = accs_of(g)
        if cid.startswith("CORUM:"):
            us, gs = corum.get(cid[6:], (None, None))
            return us is not None and bool(a & us or (s in gs))
        if cid.startswith("CPX:"):
            mem = cpx.get(cid[4:])
            return mem is not None and bool(a & mem)
        return False

    def all_complex_ids_with(g):
        a, s = accs_of(g)
        out = {f"CORUM:{k}" for k, (us, gs) in corum.items() if a & us or s in gs}
        return out | {f"CPX:{k}" for k, mem in cpx.items() if a & mem}

    rows = []
    for _, e in edits.iterrows():
        for i in (x for x in str(e.evidence_ids).split(";") if x):
            kind, val = i.split(":", 1)
            ok, why = False, ""
            if kind == "ENSEMBL116_PARALOG":
                ok, why = val in par_rows, "pair row in Ensembl 116 paralog file"
            elif kind == "NOT_ENSEMBL116_PARALOG":
                ok, why = val not in par_rows, "pair absent from Ensembl 116 paralog file"
            elif kind in ("CORUM", "CPX"):
                cid, pair = val.split("@")
                ok = cid in (corum if kind == "CORUM" else cpx) and all(complex_contains(f"{kind}:{cid}", g) for g in pair.split("|"))
                why = "complex exists and lists both genes of the pair"
            elif kind == "NO_SHARED_COMPLEX":
                a, b = val.split("|")
                ok, why = not (all_complex_ids_with(a) & all_complex_ids_with(b)), "no CORUM/CPX complex lists both"
            elif kind == "PROCAN_COAB":
                ok, why = all(x in procan_cols for x in val.split("|")), "both accessions are ProCan columns"
            elif kind == "NUSINOW_COAB":
                ok, why = all(x in nus_acc for x in val.split("|")), "both accessions are Nusinow rows"
            elif kind == "NOT_DETECTED_PROCAN":
                ok, why = val not in procan_cols, "accession absent from ProCan matrix"
            elif kind == "NOT_DETECTED_NUSINOW":
                ok, why = val not in nus_acc, "accession absent from Nusinow"
            elif kind == "DEPMAP26Q3_EXPR_TPMlt1":
                ok, why = val in expr_entrez, "Entrez column in DepMap 26Q3 expression"
            else:
                why = "unknown evidence kind"
            rows.append({"edit_id": e.edit_id, "evidence_id": i, "ok": bool(ok), "check": why})
    return pd.DataFrame(rows)
