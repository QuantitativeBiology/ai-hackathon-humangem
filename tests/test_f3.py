"""Label-free tests for the F3 follow-up (A15, post hoc / descriptive; scripts/f3_trace.py): GPR structure,
flux diffs, the FN mechanism classes and the FP biomass-component check, on toy models with known answers."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import f3_trace as ft  # noqa: E402


def test_gpr_structure_kinds_and_alternatives():
    s = ft.gpr_structure("(a and b) or (a and c) or d", "a")
    assert s["kind"] == "OR-of-AND" and s["n_branches"] == 3 and s["n_genes"] == 4
    assert s["subunit"] and not s["closed_by_ko"] and s["alternatives"] == ["d"]
    s = ft.gpr_structure("a or b or c", "b")
    assert s["kind"] == "OR" and s["alternatives"] == ["a", "c"] and not s["closed_by_ko"]
    s = ft.gpr_structure("a and b", "b")
    assert s["kind"] == "AND" and s["closed_by_ko"] and s["alternatives"] == []
    s = ft.gpr_structure("a", "a")
    assert s["kind"] == "single" and s["closed_by_ko"]
    assert ft.gpr_structure("", "a")["kind"] == "none"
    # nested OR inside a branch is flattened to top-level branches
    assert len(ft.gpr_branches("a or (b or (c and d))")) == 3


def test_flux_diff_and_uptake_increase():
    wt = pd.Series({"R1": 10.0, "R2": 0.0, "EX_a": -10.0, "EX_b": 0.0, "Rx": 5.0})
    ko = pd.Series({"R1": 0.0, "R2": 10.0, "EX_a": -4.0, "EX_b": -6.0, "Rx": 5.0})
    d = ft.flux_diff(wt, ko, exclude=["R1"])
    assert d.index[0] == "R2" and set(d.index) == {"R2", "EX_a", "EX_b"}  # Rx unchanged, R1 excluded
    assert d.loc["R2", "delta"] == 10.0
    assert ft.uptake_increase(wt, ko, ["EX_a", "EX_b"]) == ["EX_b"]


@pytest.mark.parametrize("args,expected", [
    ((False, np.nan, 0.0, []), "outside_model"),
    ((True, 0.0, 5.0, []), "isozyme_in_gpr"),
    ((True, 1.0, 0.0, []), "function_not_represented"),
    ((True, 1.0, 5.0, ["EX_b"]), "medium_uptake_bypass"),
    ((True, 0.7, 5.0, []), "alternative_pathway"),
])
def test_classify_fn(args, expected):
    assert ft.classify_fn(*args) == expected


def test_modal_class_line_weighted():
    c = pd.Series(["alternative_pathway", "isozyme_in_gpr", "alternative_pathway"], index=["s1", "s2", "s3"])
    w = pd.Series([1, 10, 2], index=["s1", "s2", "s3"])
    assert ft.modal_class(c, w) == "isozyme_in_gpr"


def test_crosstab_counts():
    cells = pd.DataFrame({"mechanism_class": ["a", "a", "b"], "fa_x": [True, False, True]})
    t = ft.crosstab(cells, ["fa_x"]).set_index("mechanism_class")
    assert t.loc["a", "n_cells"] == 2 and t.loc["a", "fa_x"] == 1 and t.loc["all", "fa_x"] == 2


def _toy():
    """A -> B by R1 (g1 or g2) or by A -> C (R2, g3) -> B (R3, g4); B -> biomass; R4 (g5): D -> B uses EX_d; R5 (g6)
    is an orphan dead end (E -> F) never needed."""
    cobra = pytest.importorskip("cobra")
    cobra.Configuration().solver = "glpk"
    m = cobra.Model("toy")
    mets = {k: cobra.Metabolite(k, compartment="c") for k in "ABCDEF"}
    rx = []

    def R(rid, st, rule="", lb=0.0, ub=1000.0):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites({mets[k]: v for k, v in st.items()})
        r.gene_reaction_rule = rule
        rx.append(r)

    R("EX_a", {"A": -1}, lb=-10.0)
    R("EX_d", {"D": -1}, lb=-10.0)
    R("R1", {"A": -1, "B": 1}, "g1 or g2")
    R("R2", {"A": -1, "C": 1}, "g3")
    R("R3", {"C": -1, "B": 1}, "g4")
    R("R4", {"D": -2, "B": 1}, "g5")
    R("R5", {"E": -1, "F": 1}, "g6")
    R("BIO", {"B": -1})
    m.add_reactions(rx)
    m.objective = "BIO"
    return m


def test_toy_reaction_ko_and_bypass():
    m = _toy()
    wt = m.slim_optimize()
    assert wt == pytest.approx(15.0)  # 10 from A plus 5 from D
    with m:
        m.genes.get_by_id("g1").knock_out()
        assert m.reactions.R1.bounds != (0.0, 0.0)  # g2 remains: not addressable
        assert m.slim_optimize() == pytest.approx(wt)
    with m:  # all-reactions KO of R1: the A -> C -> B route carries the flux
        m.reactions.R1.bounds = (0.0, 0.0)
        assert m.slim_optimize() == pytest.approx(wt)
        with m:
            m.reactions.R2.bounds = (0.0, 0.0)
            assert m.slim_optimize() == pytest.approx(5.0)  # only D left
    from cobra.flux_analysis import pfba

    v = pfba(m).fluxes
    assert abs(v["R5"]) < 1e-9  # orphan reaction: function_not_represented
    assert ft.classify_fn(True, 1.0, abs(v["R5"]), []) == "function_not_represented"


def test_blocked_biomass_components_toy():
    cobra = pytest.importorskip("cobra")
    cobra.Configuration().solver = "glpk"
    m = cobra.Model("bio")
    a, b, c = (cobra.Metabolite(k, compartment="c") for k in ("a", "b", "c"))
    ex = cobra.Reaction("EX_a", lower_bound=-10.0, upper_bound=1000.0)
    ex.add_metabolites({a: -1})
    r1 = cobra.Reaction("R1", upper_bound=1000.0)
    r1.add_metabolites({a: -1, b: 1})
    r1.gene_reaction_rule = "g1"
    r2 = cobra.Reaction("R2", upper_bound=1000.0)
    r2.add_metabolites({a: -1, c: 1})
    r2.gene_reaction_rule = "g2"
    bio = cobra.Reaction("BIO", upper_bound=1000.0)
    bio.add_metabolites({b: -1, c: -1})
    m.add_reactions([ex, r1, r2, bio])
    m.objective = "BIO"
    wt = ft.blocked_biomass_components(m, "BIO", skip=())
    assert wt["b"] == pytest.approx(10.0) and wt["c"] == pytest.approx(10.0)
    with m:
        m.genes.get_by_id("g2").knock_out()
        ko = ft.blocked_biomass_components(m, "BIO", skip=())
    assert ko["c"] == 0.0 and ko["b"] == pytest.approx(10.0)
