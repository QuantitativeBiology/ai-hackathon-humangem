"""Label-free tests for the F1 follow-up (A15, post hoc; scripts/f1_*.py): the D17 rule, arms, reconstruction and the
validation spec shape."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import f1_common as fc  # noqa: E402

MAR07161 = {"MAM01371n": -0.18, "MAM01623n": -0.3, "MAM02034n": -0.34, "MAM03130n": -0.18, "MAM02759n": 1.0,
            "MAM02847n": 1.0}
CYOOM3I = {"MAM01826m": -4.0, "MAM02039m": -7.92, "MAM02630m": -1.0, "MAM01824m": 4.0, "MAM02039i": 4.0,
           "MAM02040m": 1.96, "MAM02631m": 0.02}
HALF_O2 = {"a": -1.0, "o2": -0.5, "b": 1.0, "h2o": 1.0}


def facts_table(table):
    def f(a, b):
        v = table.get(frozenset((a, b)))
        if v is None:
            return fc.PairFacts(False, None, float("nan"), (), frozenset({"Cytosol"}), frozenset({"Cytosol"}))
        return v
    return f


def pf(pid=50.0, ca="Cytosol", cb="Cytosol", sc=()):
    return fc.PairFacts(True, "within_species_paralog", pid, tuple(sc), fc.comp_set(ca), fc.comp_set(cb))


def test_rule_constants_fixed_by_d17():
    assert fc.MIN_PID == 30.0 and fc.MAX_AND_GENES == 6
    assert set(fc.LUMPED_SUBSYSTEMS) == {"Artificial reactions", "Pool reactions"}
    assert fc.TYMS_REACTION == "MAR04644" and fc.TYMS_BOUNDS == (0.0, 1000.0)


def test_composition_lumped_detects_mar07161_not_cyoom3i():
    assert fc.composition_lumped(MAR07161)
    assert fc.lumped_reason("Nucleotide metabolism", MAR07161).startswith("composition")
    assert not fc.composition_lumped(CYOOM3I)
    assert not fc.composition_lumped(HALF_O2)
    assert fc.lumped_reason("Pool reactions", HALF_O2) == "subsystem:Pool reactions"
    assert fc.lumped_reason("Purine metabolism", HALF_O2) == ""


def test_max_and_genes():
    assert fc.max_and_genes(None) == 0 and fc.max_and_genes("A") == 0
    t = ("or", (("and", ("A", "B")), ("and", ("C", "D", "E", "F", ("or", ("G", "H"))))))
    assert fc.max_and_genes(t) == 6
    assert fc.max_and_genes(("or", ("A", "B", "C", "D", "E", "F", "G"))) == 0


def test_reaction_exclusions():
    big = ("and", tuple("ABCDEF"))
    small = ("and", ("A", "B"))
    assert fc.reaction_exclusions(boundary=True, subsystem="x", stoich=HALF_O2, met_comps={"c"}, tree=small) == ["E1 exchange"]
    assert fc.reaction_exclusions(boundary=False, subsystem="Transport reactions", stoich=HALF_O2, met_comps={"c"},
                                  tree=small) == ["E2 transport"]
    assert fc.reaction_exclusions(boundary=False, subsystem="x", stoich=HALF_O2, met_comps={"c", "m"},
                                  tree=small) == ["E2 transport"]
    ex = fc.reaction_exclusions(boundary=False, subsystem="Nucleotide metabolism", stoich=MAR07161, met_comps={"n"},
                                tree=big)
    assert [e[:2] for e in ex] == ["E3", "E4"]
    assert fc.reaction_exclusions(boundary=False, subsystem="x", stoich=HALF_O2, met_comps={"c"}, tree=small) == []


def test_same_compartment_is_set_equality_and_fails_closed():
    assert fc.same_compartment("Cytosol;Nucleus", "Nucleus;Cytosol")
    assert not fc.same_compartment("Cytosol", "Cytosol;Nucleus")
    assert not fc.same_compartment("", "")
    assert not fc.same_compartment(float("nan"), "Cytosol")


def test_pair_facts_criteria():
    assert pf().qualifies
    assert not pf(pid=29.99).qualifies and "min %id" in pf(pid=29.99).fail_reasons()[0]
    assert pf(pid=30.0).qualifies
    assert not pf(cb="Mitochondria").qualifies
    assert not pf(ca="", cb="").qualifies
    assert not pf(sc=("CORUM:1",)).qualifies
    nonpar = fc.PairFacts(False, None, math.nan, (), frozenset({"x"}), frozenset({"x"}))
    assert not nonpar.qualifies and nonpar.fail_reasons()[0] == "not an Ensembl paralog pair"


def test_curate_tree_merges_qualifying_pair():
    from nextgem import w5

    tree = w5.canon(("and", ("A", "B", "C")))
    rec = []
    new = fc.curate_tree(tree, facts_table({frozenset("AB"): pf()}), rec, {"rxn": "R1"})
    assert new == w5.canon(("and", ("C", ("or", ("A", "B")))))
    kinds = [r["kind"] for r in rec]
    assert kinds.count("pair") == 1 and kinds.count("merge") == 1 and all(r["applied"] for r in rec)


def test_curate_tree_low_identity_or_excluded_reaction_unchanged():
    from nextgem import w5

    tree = w5.canon(("and", ("A", "B", "C")))
    rec = []
    assert fc.curate_tree(tree, facts_table({frozenset("AB"): pf(pid=25.0)}), rec, {}) == tree
    assert rec[0]["kind"] == "pair" and not rec[0]["pair_qualifies"]
    rec = []
    assert fc.curate_tree(tree, facts_table({frozenset("AB"): pf()}), rec, {}, apply=False) == tree
    assert [r["applied"] for r in rec] == [False, False]


def test_curate_tree_non_clique_blocked():
    from nextgem import w5

    tree = w5.canon(("and", ("A", "B", "C")))
    rec = []
    t = {frozenset("AB"): pf(), frozenset("BC"): pf()}  # A-C not a pair: component {A,B,C} is not a clique
    assert fc.curate_tree(tree, facts_table(t), rec, {}) == tree
    assert any(r["kind"] == "blocked" for r in rec)


def test_curate_tree_nested_and_nodes():
    from nextgem import w5

    tree = w5.parse_rule("(A and B and X) or (A and B and Y)")
    new = fc.curate_tree(tree, facts_table({frozenset("AB"): pf()}), [], {})
    assert new == w5.parse_rule("((A or B) and X) or ((A or B) and Y)")


def test_arms_and_methods_acceptable_to_harness():
    import evaluate_val as ev

    assert [a.key for a in fc.ARMS] == ["P1", "TYMS", "GPR", "F1"]
    assert (fc.CANDIDATE.tyms, fc.CANDIDATE.gpr) == (True, True) and fc.CANDIDATE.role == "candidate"
    assert fc.ARM["P1"].method == "M6:H1+Bloo/M2(f=1)"
    for a in fc.ARMS:
        assert not ev.LABEL_LIKE_RE.search(a.method) and not ev.ID_RE.match(a.method) and a.method.strip() == a.method
    assert len({a.method for a in fc.ARMS}) == 4


def test_contrasts_predeclared():
    import f1_score_train as st

    assert st.CONTRASTS[0] == ("F1", "P1", "primary")
    assert [c for c in st.CONTRASTS if c[2] == "primary"] == [("F1", "P1", "primary")]
    assert {("TYMS", "P1"), ("GPR", "P1"), ("F1", "TYMS")} <= {(a, b) for a, b, _ in st.CONTRASTS}
    assert st.YCELLS == ("L_prob", "L_prob_noCE", "L_prob_expr")


def test_tg_rule():
    import m6_common as mc

    assert mc.tg_pass(0.005, 0.001, 9) and not mc.tg_pass(0.0049, 0.001, 12)
    assert not mc.tg_pass(0.01, 0.0, 12) and not mc.tg_pass(0.01, 0.001, 8)


def test_arm_ratios_reconstruction():
    from f1_sweeps import arm_ratios

    P1 = pd.DataFrame([[1.0, 0.0, 0.5], [1.0, 0.0, 0.7]], index=["s1", "s2"], columns=["g1", "g2", "g3"])
    delta = pd.DataFrame({"arm": ["GPR", "GPR", "TYMS"], "w2_sweep_id": ["s1", "s2", "s2"], "gene": ["g2", "g2", "g1"],
                          "ratio_r": [1.0, 1.0, 0.0], "ratio_r_P1": [0.0, 0.0, 1.0]})
    meta = pd.DataFrame({"arm": ["TYMS", "TYMS"], "purpose": ["arm", "arm"], "w2_sweep_id": ["s1", "s2"],
                         "status": ["failed-timeout", "ok"]})
    g = arm_ratios("GPR", P1, delta)
    assert g.loc["s1", "g2"] == 1.0 and g.loc["s2", "g2"] == 1.0 and g.loc["s1", "g3"] == 0.5
    t = arm_ratios("TYMS", P1, delta, meta)
    assert t.loc["s1"].isna().all() and t.loc["s2", "g1"] == 0.0
    assert arm_ratios("P1", P1, delta).equals(P1)


def test_val_spec_passes_harness_static_checks():
    import evaluate_val as ev
    import f1_val_spec as vs

    train = {k: {"delta_lin": 0.01, "ci": [0.005, 0.015], "folds_positive": 12, "tg_pass": True}
             for k in ("F1_vs_P1", "TYMS_vs_P1", "GPR_vs_P1")}
    spec = vs.build_spec("0" * 64, {"scripts/f1_common.py": "1" * 64}, train)
    ev.check_spec(spec, vs.LOOK_ID, "m10")
    assert spec["workstream"] == spec["milestone"] == "M10" and spec["reference"] == "H1"
    assert spec["candidates"] == [fc.CANDIDATE.method] and fc.P1_METHOD in spec["methods"]


def test_edit_outputs_respect_rule():
    p = fc.F1 / "edits_human1.tsv"
    if not p.is_file():
        pytest.skip("results/f1/edits_human1.tsv not built")
    e = pd.read_csv(p, sep="\t")
    ex = pd.read_csv(fc.F1 / "excluded_reactions_human1.tsv", sep="\t")
    m = e[e.kind == "AND2OR"]
    assert not set(m.rxn) & set(ex.rxn)
    assert not m.symbols.str.contains("POLR").any()
    assert (e.kind == "TYMS_IRREVERSIBLE").sum() == 1
    assert "MAR07161" in set(pd.read_csv(fc.F1 / "lumped_reactions_human1.tsv", sep="\t").rxn)
    chk = pd.read_csv(fc.F1 / "evidence_check_human1.tsv", sep="\t")
    assert chk.ok.all()
