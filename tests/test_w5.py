"""W5 GPR curation (DESIGN §12 W5 row): toy edit flips the KO; evidence rules; evidence validator; protected list.

All label-free. Test names are stable (DESIGN §11.3 test-run check).
"""

from __future__ import annotations

import gzip

import cobra
import pandas as pd
import pytest

from nextgem import w5
from nextgem import w5_evidence as ev


# ----------------------------------------------------------------------------- GPR trees
def test_w5_parse_canon_roundtrip():
    t = w5.parse_rule("(B and A) or (A and B) or C")
    assert t == ("or", (("and", ("A", "B")), "C")) or t == ("or", ("C", ("and", ("A", "B"))))
    assert w5.parse_rule(w5.to_rule(t)) == t
    assert w5.parse_rule("A or (B or C)") == ("or", ("A", "B", "C"))
    assert w5.parse_rule("") is None
    with pytest.raises(ValueError):
        w5.parse_rule("A and (B or C")


def test_w5_enumerate_edit_kinds_and_counts():
    t = w5.parse_rule("(A and B and C) or D")
    es = w5.enumerate_edits(t)
    kinds = pd.Series([e.kind for e in es]).value_counts().to_dict()
    # root OR: 1 OR2AND + 2 RMBRANCH; inner AND of 3: 3 pairwise AND2OR
    assert kinds == {"OR2AND": 1, "RMBRANCH": 2, "AND2OR": 3}
    res = {w5.to_rule(w5.apply_edit(t, e)) for e in es}
    assert "A and B and C and D" in res  # OR2AND flattens into one AND
    assert "D" in res and "A and B and C" in res  # branch removals
    assert any("(A or B)" in r for r in res)


def test_w5_gpr_eval_matches_cobra():
    rule = "(G1 and G2) or (G3 and (G4 or G5))"
    t = w5.parse_rule(rule)
    m = cobra.Model("t")
    r = cobra.Reaction("R")
    m.add_reactions([r])
    r.gene_reaction_rule = rule
    for g in ["G1", "G2", "G3", "G4", "G5"]:
        assert w5.gpr_true(t, {g}) == r.gpr.eval({g})


def _toy_model():
    """A -> B -> biomass; two parallel B-producing reactions R1 (G1) and R2 (G2), plus R3 (G3 and G4)."""
    m = cobra.Model("toy")
    a, b, c = (cobra.Metabolite(x, compartment="c") for x in ("a", "b", "c"))
    ex = cobra.Reaction("EX_a", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({a: -1})
    r1 = cobra.Reaction("R1", upper_bound=1000)
    r1.add_metabolites({a: -1, b: 1})
    r2 = cobra.Reaction("R2", upper_bound=1000)
    r2.add_metabolites({b: -1, c: 1})
    bio = cobra.Reaction("BIO", upper_bound=1000)
    bio.add_metabolites({c: -1})
    m.add_reactions([ex, r1, r2, bio])
    r1.gene_reaction_rule = "G1 or G2"
    r2.gene_reaction_rule = "G3 and G4"
    m.objective = "BIO"
    return m


def test_w5_toy_edit_flips_ko():
    """DESIGN §12 W5 known answer: OR->AND makes the isozyme lethal; AND->OR rescues the subunit; the local
    re-evaluation equals a full recomputation with cobra's own knock_out on the edited model."""
    cobra.Configuration().solver = "glpk"
    m = _toy_model()
    eng = w5.KOEngine(m)
    gp = w5.model_gprs(m)
    base = w5.full_sweep(eng, gp, ["G1", "G2", "G3", "G4"])
    assert base.to_dict() == {"G1": 1.0, "G2": 1.0, "G3": 0.0, "G4": 0.0}
    cl = w5.closures(gp)
    t1 = gp["R1"]
    new1 = w5.apply_edit(t1, w5.Edit("OR2AND", ()))
    new2 = w5.apply_edit(gp["R2"], w5.Edit("AND2OR", (), (0, 1)))
    got = {}
    for g in ("G1", "G2"):
        got[g] = eng.ratio(w5.edited_closure(cl.get(g, frozenset()), frozenset({"R1"}), new1, g))
    for g in ("G3", "G4"):
        got[g] = eng.ratio(w5.edited_closure(cl.get(g, frozenset()), frozenset({"R2"}), new2, g))
    assert got == {"G1": 0.0, "G2": 0.0, "G3": 1.0, "G4": 1.0}
    # full recomputation with cobra on an edited copy
    m2 = m.copy()
    m2.reactions.R1.gene_reaction_rule = w5.to_rule(new1)
    m2.reactions.R2.gene_reaction_rule = w5.to_rule(new2)
    for g, v in got.items():
        with m2:
            m2.genes.get_by_id(g).knock_out()
            assert round((m2.slim_optimize(error_value=0.0) or 0.0) / eng.wt, 6) == v


def test_w5_zero_flux_shortcut_is_exact():
    cobra.Configuration().solver = "glpk"
    m = _toy_model()
    extra = cobra.Reaction("R9", upper_bound=1000)
    extra.add_metabolites({m.metabolites.a: -1, m.metabolites.c: 1})
    m.add_reactions([extra])
    extra.bounds = (0, 0)
    eng = w5.KOEngine(m)
    assert eng.ratio({"R9"}) == 1.0 and eng.n_lp == 0


# ----------------------------------------------------------------------------- protected list (§3, verbatim)
def test_w5_protected_tokens_verbatim():
    s3 = w5.section3_text()
    for ent in w5.PROTECTED_S3:
        assert ent["text"] in s3, ent["text"]
        for r in ent["reactions"]:
            assert r in s3


def test_w5_protected_or_node():
    p = w5.Protected(reactions={"MARX"}, groups=[{"A", "B"}], unresolved=[])
    assert p.or_node_protected(["MARX"], ("or", ("Q", "Z")))
    assert p.or_node_protected(["MARY"], ("or", ("A", "B", "C")))
    assert p.or_node_protected(["MARY"], ("or", (("and", ("A", "X")), "B")))
    assert p.or_node_protected(["MARY"], ("or", ("A", "C"))) is None


def test_w5_resolve_protected_human2():
    p = w5.resolve_protected("human2")
    assert "MAR04363" in p.reactions
    assert not p.unresolved, p.unresolved


# ----------------------------------------------------------------------------- evidence rules
class _DB(ev.EvidenceDB):
    """EvidenceDB with injected facts (no files)."""

    def __init__(self, paralog=(), complexes=None, coab=None, absent=()):
        super().__init__(genes=pd.DataFrame(columns=["uniprot", "symbol", "entrez", "compartments"]))
        self.paralog = {frozenset(p): ("within_species_paralog", 80.0) for p in paralog}
        self.complexes = complexes or {}
        self._coab = coab or {}
        self._absent = set(absent)

    def coab(self, a, b):
        r = self._coab.get(frozenset((a, b)))
        if r is None:
            return {"source": None, "r": float("nan"), "n": 0, "id": None}
        return {"source": "PROCAN", "r": r, "n": 40, "id": f"PROCAN_COAB:{a}|{b}"}

    def pan_cancer_absent(self, gene):
        ok = gene in self._absent
        return {"absent": ok, "soluble_proxy": True, "procan": not ok, "nusinow": not ok, "frac_tpm_lt1": 1.0 if ok else 0.0,
                "n_expr_lines": 100, "ids": [f"NOT_DETECTED_PROCAN:{gene}"] if ok else []}


def _prot(groups=(), reactions=()):
    return w5.Protected(reactions=set(reactions), groups=[set(g) for g in groups], unresolved=[])


def test_w5_rule_or2and_requires_all_three():
    t = w5.parse_rule("A or B")
    cx = {"A": {"CORUM:1"}, "B": {"CORUM:1"}}
    for kw, expect in [({"complexes": cx, "coab": {frozenset("AB"): 0.6}}, "A and B"),
                       ({"complexes": cx, "coab": {frozenset("AB"): 0.4}}, "A or B"),
                       ({"complexes": {}, "coab": {frozenset("AB"): 0.9}}, "A or B"),
                       ({"complexes": cx, "coab": {frozenset("AB"): 0.9}, "paralog": [("A", "B")]}, "A or B")]:
        cur = ev.Curation()
        new = ev.curate_tree(t, ["R"], _DB(**kw), _prot(), cur, {})
        assert w5.to_rule(new) == expect, kw


def test_w5_rule_and2or_paralogs_without_shared_complex():
    t = w5.parse_rule("R1 and R2 and R2B and T")
    cur = ev.Curation()
    new = ev.curate_tree(t, ["R"], _DB(paralog=[("R2", "R2B")]), _prot(), cur, {})
    assert new == w5.parse_rule("R1 and T and (R2 or R2B)")
    cur = ev.Curation()
    new = ev.curate_tree(t, ["R"], _DB(paralog=[("R2", "R2B")], complexes={"R2": {"CPX:9"}, "R2B": {"CPX:9"}}), _prot(), cur, {})
    assert new == t


def test_w5_rule_remove_absent_branch_and_protection():
    t = w5.parse_rule("A or B or C")
    cur = ev.Curation()
    assert w5.to_rule(ev.curate_tree(t, ["R"], _DB(absent={"C"}), _prot(), cur, {})) == "A or B"
    cur = ev.Curation()
    assert ev.curate_tree(t, ["R"], _DB(absent={"C"}), _prot(groups=[{"A", "C"}]), cur, {}) == t
    assert cur.removals_blocked and not cur.edits
    # never remove the last branch
    cur = ev.Curation()
    assert ev.curate_tree(w5.parse_rule("A or B"), ["R"], _DB(absent={"A", "B"}), _prot(), cur, {}) == w5.parse_rule("A or B")


def test_w5_protected_or_never_converted():
    t = w5.parse_rule("A or B")
    cx = {"A": {"CORUM:1"}, "B": {"CORUM:1"}}
    cur = ev.Curation()
    new = ev.curate_tree(t, ["MARP"], _DB(complexes=cx, coab={frozenset("AB"): 0.9}), _prot(reactions={"MARP"}), cur, {})
    assert new == t and not cur.edits


def test_w5_or_pairs_use_distinguishing_genes():
    node = w5.parse_rule("(A and B) or (A and C)")
    assert ev.or_pairs(node) == [("B", "C")]
    assert ev.or_pairs(w5.parse_rule("A or (A and B)")) == []


def test_w5_convert_or_nodes_multiple_paths():
    t = ("and", (("or", ("A", "B")), ("or", ("C", "D")), "E"))
    assert w5.to_rule(ev.convert_or_nodes(t, [(0,), (1,)])) == "A and B and C and D and E"


# ----------------------------------------------------------------------------- validator
def _toy_evidence_files(root):
    (root / "paralogs").mkdir()
    (root / "complexes").mkdir()
    (root / "proteomics_goncalves").mkdir()
    (root / "proteomics_nusinow").mkdir()
    (root / "depmap").mkdir()
    (root / "paralogs" / "ensembl_human_paralogs.tsv").write_text(
        "Gene stable ID\tHuman paralogue gene stable ID\tx\ty\tz\tw\tt\nENSG1\tENSG2\tP\t50\t50\tA\twithin_species_paralog\n")
    (root / "complexes" / "corum_humanComplexes.txt").write_text(
        "complex_id\tcomplex_name\torganism\tsubunits_uniprot_id\tsubunits_gene_name\n7\tX\tHuman\tP11111;P22222\tGA;GB\n")
    (root / "complexes" / "complexportal_9606.tsv").write_text(
        "#Complex ac\tRecommended name\tExpanded participant list\nCPX-5\tY\tP11111(1)|P33333(1)\n")
    (root / "proteomics_goncalves" / "Protein_matrix_averaged_20250211.tsv").write_text("uniprot_id\t\tP11111\tP22222\n")
    with gzip.open(root / "proteomics_nusinow" / "protein_quant_current_normalized.csv.gz", "wt") as f:
        f.write('"Protein_Id","Uniprot_Acc"\n"x","P33333"\n')
    (root / "depmap" / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv").write_text(",GA (101),GB (102)\n")


def test_w5_validator_checks_every_id(tmp_path, monkeypatch):
    _toy_evidence_files(tmp_path)
    monkeypatch.setenv("NEXTGEM_DATA_ROOT", str(tmp_path))
    genes = pd.DataFrame({"uniprot": [["P11111"], ["P22222"], ["P33333"], ["P44444"]], "symbol": ["GA", "GB", "GC", "GD"]},
                         index=["ENSG1", "ENSG2", "ENSG3", "ENSG4"])
    good = pd.DataFrame([
        {"edit_id": "e1", "genes": "ENSG1;ENSG2", "evidence_ids": "CORUM:7@ENSG1|ENSG2;PROCAN_COAB:P11111|P22222"},
        {"edit_id": "e2", "genes": "ENSG1;ENSG2", "evidence_ids": "ENSEMBL116_PARALOG:ENSG1|ENSG2"},
        {"edit_id": "e3", "genes": "ENSG4", "evidence_ids": "NOT_DETECTED_PROCAN:P44444;NOT_DETECTED_NUSINOW:P44444;DEPMAP26Q3_EXPR_TPMlt1:101"},
        {"edit_id": "e4", "genes": "ENSG2;ENSG3", "evidence_ids": "NO_SHARED_COMPLEX:ENSG2|ENSG3;NOT_ENSEMBL116_PARALOG:ENSG2|ENSG3"},
        {"edit_id": "e5", "genes": "ENSG1;ENSG3", "evidence_ids": "CPX:CPX-5@ENSG1|ENSG3;NUSINOW_COAB:P33333|P33333"},
    ])
    v = ev.validate(good, genes)
    assert v.ok.all(), v[~v.ok]
    bad = pd.DataFrame([
        {"edit_id": "b1", "genes": "ENSG1;ENSG3", "evidence_ids": "CORUM:7@ENSG1|ENSG3"},  # ENSG3 not in CORUM 7
        {"edit_id": "b2", "genes": "ENSG1;ENSG2", "evidence_ids": "CORUM:99@ENSG1|ENSG2"},  # no such complex
        {"edit_id": "b3", "genes": "ENSG1;ENSG3", "evidence_ids": "ENSEMBL116_PARALOG:ENSG1|ENSG3"},
        {"edit_id": "b4", "genes": "ENSG1", "evidence_ids": "NOT_DETECTED_PROCAN:P11111"},
        {"edit_id": "b5", "genes": "ENSG1;ENSG2", "evidence_ids": "NO_SHARED_COMPLEX:ENSG1|ENSG2"},
        {"edit_id": "b6", "genes": "ENSG1;ENSG2", "evidence_ids": "NOT_ENSEMBL116_PARALOG:ENSG1|ENSG2"},
        {"edit_id": "b7", "genes": "ENSG1", "evidence_ids": "DEPMAP26Q3_EXPR_TPMlt1:999"},
    ])
    v = ev.validate(bad, genes)
    assert (~v.ok).all(), v[v.ok]


def test_w5_validator_refuses_denylisted_source():
    e = pd.DataFrame([{"edit_id": "x", "genes": "A", "evidence_ids": "XLSX:41588_2021_967_MOESM2_ESM.xlsx"}])
    with pytest.raises(AssertionError):
        ev.validate(e, pd.DataFrame(columns=["uniprot", "symbol"]))
