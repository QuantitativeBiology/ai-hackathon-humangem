"""Label-free tests for the F4 follow-up (A15 post hoc; RESEARCH_TREE D19; scripts/f4_*.py): the rule constants, the
genomic constraint on the P1 engine (toy known answers), job keys, reconstruction and the pre-check verdicts."""

from __future__ import annotations

import sys
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import f4_common as fc  # noqa: E402

from nextgem import context as cx  # noqa: E402
from nextgem import context_inputs as ci  # noqa: E402
from nextgem.model import BIOMASS_REACTION as BIO  # noqa: E402

# glc -> a (gA or gB) -> b (gC) -> biomass; glc -> b (gD) is a second route
TOY = [
    ("EX_glc", "glc <=> ", "", -1000, 1000),
    ("R1", "glc --> a", "gA or gB", 0, 1000),
    ("R2", "a --> b", "gC", 0, 1000),
    ("R3", "glc --> b", "gD", 0, 5),
    (BIO, "b --> ", "", 0, 1000),
]


def toy_model():
    m = cobra.Model("toy")
    m.solver = "glpk"
    for rid, eq, gpr, lb, ub in TOY:
        r = cobra.Reaction(rid)
        m.add_reactions([r])
        r.build_reaction_from_string(eq, verbose=False)
        r.gene_reaction_rule = gpr
        r.bounds = (lb, ub)
    for met in m.metabolites:
        met.compartment = "c"
    m.objective = BIO
    return m


def test_rule_constants_fixed_by_d19():
    assert fc.LOSS_RULE == "r2-11"
    assert ci.CN_DEL == 0.3 and ci.CN_LOH_PROXY == 0.5 and ci.EXPR_SILENT == 1.0
    assert fc.F4B_CONFIG == ("F4", "R1", "rna", "real", "none")
    assert fc.TG_MIN_DELTA == 0.005 and fc.TG_MIN_FOLDS == 9 and fc.PANEL_MIN_DELTA == 0.005
    assert fc.MAX_WORKERS == 12 and fc.JOB_CPU_LIMIT_S == 1800
    assert fc.BLOO == "B-loo:cofactor_pool_biomass" and fc.P1_METHOD == "M6:H1+Bloo/M2(f=1)"
    assert cx.LP_TIME_LIMIT_S == 300.0 and cx.DECLARED_SOLVER["F4"] == "glpk" and cx.LOSS_RULE_SOLVER == "glpk"


def test_tg_and_panel_verdicts():
    assert fc.tg_pass(0.005, 0.001, 9)
    assert not fc.tg_pass(0.0049, 0.001, 12)
    assert not fc.tg_pass(0.01, 0.0, 12)
    assert not fc.tg_pass(0.01, 0.001, 8)
    assert not fc.tg_pass(float("nan"), 0.1, 12)
    assert fc.panel_proceeds(0.005) and fc.panel_proceeds(0.0051)
    assert not fc.panel_proceeds(0.0044) and not fc.panel_proceeds(float("nan"))


def test_job_key_is_order_independent_and_signature_specific():
    k = fc.job_key("human1:abc", ["gB", "gA"])
    assert k == fc.job_key("human1:abc", ["gA", "gB"]) == fc.job_key("human1:abc", fc.split_genes("gA;gB"))
    assert k != fc.job_key("human1:xyz", ["gA", "gB"]) and k != fc.job_key("human1:abc", ["gA"])
    assert fc.split_genes(fc.join_genes({"b", "a"})) == ["a", "b"]
    assert fc.split_genes("") == [] and fc.split_genes(float("nan")) == []
    assert fc.sweep_mode(0) == "partial" and fc.sweep_mode(3) == "full"


def test_loss_rule_known_answer_r2_11():
    lines = ["L1"]
    dmg = pd.DataFrame({"g2": [2.0], "g1lo": [1.0], "g1si": [1.0], "g1ok": [1.0], "gdel": [0.0], "gok": [0.0]},
                       index=lines)
    cn = pd.DataFrame({"g2": [1.0], "g1lo": [0.4], "g1si": [1.0], "g1ok": [1.0], "gdel": [0.29], "gok": [0.31]},
                      index=lines)
    rna = pd.DataFrame({"g2": [5.0], "g1lo": [5.0], "g1si": [0.5], "g1ok": [5.0], "gdel": [5.0], "gok": [5.0]},
                       index=lines)
    lost, flags = ci.lost_genes("L1", dmg, cn, rna, rule=fc.LOSS_RULE)
    assert lost == ["g1lo", "g1si", "g2", "gdel"] and not flags["cn_missing"]
    lost, flags = ci.lost_genes("L1", dmg, None, rna, rule=fc.LOSS_RULE)  # no CN: mutations only (flagged)
    assert lost == ["g1si", "g2"] and flags["cn_missing"]


def test_knockout_state_equals_loss_closure_and_a4_convention():
    """cobra knock-outs of the accepted genes close exactly context.loss_closure; in that state the KO of the lost
    gene's isozyme partner closes the shared reaction and the lost gene's own KO closes nothing (A4)."""
    from nextgem import w3

    m = toy_model()
    p = cx.parent_from_model(m.copy())
    assert cx.loss_closure(p, {"gA"}) == set()
    assert cx.loss_closure(p, {"gA", "gB"}) == {"R1"}
    generic = w3.gene_closures(m)
    assert generic["gB"] == [] and generic["gC"] == ["R2"]
    with m:
        m.genes.get_by_id("gA").knock_out()
        closed = {r.id for r in m.reactions if r.bounds == (0, 0)}
        assert closed == cx.loss_closure(p, {"gA"}) == set()
        lost_state = w3.gene_closures(m)
    assert lost_state["gB"] == ["R1"] and lost_state["gA"] == []
    assert w3.gene_closures(m)["gB"] == []  # context restored


def test_lost_state_sweep_known_answer_and_partial_shortcut():
    """P1 engine (w3.ko_sweep) in the lost state: gB becomes partially essential (only R3, capped at 5, remains)."""
    from nextgem import w3

    m = toy_model()
    v = w3.Variant("T0")
    genes = sorted(g.id for g in m.genes)
    bounds = {"EX_glc": 10.0}
    base = w3.ko_sweep(m, genes, w3.gene_closures(m), bounds, v)
    assert base["wt_growth"] == pytest.approx(10.0)
    assert dict(zip(genes, base["ratio_r"])) == {"gA": 1.0, "gB": 1.0, "gC": 0.5, "gD": 1.0}
    with m:
        m.genes.get_by_id("gA").knock_out()
        clos = w3.gene_closures(m)
        full = w3.ko_sweep(m, genes, clos, bounds, v)
        changed = [g for g in genes if set(clos[g]) != set(w3.gene_closures(toy_model())[g])]
        part = w3.ko_sweep(m, changed, clos, bounds, v)
    got = dict(zip(genes, full["ratio_r"]))
    assert got == {"gA": 1.0, "gB": 0.5, "gC": 0.5, "gD": 1.0}
    assert changed == ["gB"]  # gA: closure [] in both states (isozyme; then already lost)
    comp = dict(zip(genes, base["ratio_r"]))
    comp.update(dict(zip(changed, part["ratio_r"])))
    assert comp == got  # partial mode (no reaction closed) is exact


def test_job_ratios_reconstruction():
    from f4_sweeps import job_ratios

    P1 = pd.DataFrame([[1.0, 0.5, 0.0], [1.0, 1.0, 0.2]], index=["s1", "s2"], columns=["gA", "gB", "gC"])
    jobs = pd.DataFrame({"job_key": ["s1|k1", "s2|k2", "s1|k3", "s1|P1"], "w2_sweep_id": ["s1", "s2", "s1", "s1"],
                         "purpose": ["arm", "arm", "arm", "spot_p1"], "status": ["ok", "ok", "failed-timeout", "ok"]})
    delta = pd.DataFrame({"job_key": ["s1|k1", "s2|k2"], "gene": ["gB", "gC"], "ratio_r": [0.0, 0.0],
                          "ratio_r_P1": [0.5, 0.2]})
    R = job_ratios(P1, delta, jobs)
    assert list(R.index) == ["s1|k1", "s2|k2", "s1|k3"]
    assert R.loc["s1|k1"].tolist() == [1.0, 0.0, 0.0] and R.loc["s2|k2"].tolist() == [1.0, 1.0, 0.0]
    assert R.loc["s1|k3"].isna().all()


def test_eflux_build_and_score_ids_agree():
    import f4_eflux_build as eb
    import f4_eflux_score as es

    assert es.f4b_build_ids(["ACH-1"])["ACH-1"] == eb.build_id("ACH-1") + "|none"
    assert eb.build_id("ACH-1") == "F4b|ACH-1|F4|R1|rna|real"
    panel = pd.DataFrame({"ModelID": ["L1", "L2", "L3"]})
    specs = eb.enumerate_jobs(panel, pd.Series({"L1": "s2", "L2": "s1", "L3": "s2"}))
    assert [s["build_id"] for s in specs] == ["generic|s1", "generic|s2", "F4b|L1|F4|R1|rna|real",
                                              "F4b|L2|F4|R1|rna|real", "F4b|L3|F4|R1|rna|real"]
    assert [s["w2_sweep_id"] for s in specs[2:]] == ["s2", "s1", "s2"]


def test_line_signatures_match_w2_and_cover_training_lines():
    from nextgem import labels, w1

    sig = fc.line_signatures()
    split = labels.load_split().set_index("ModelID")
    E = fc.scored_training_lines(split, w1.line_features())
    assert len(E) == 686 and sig.reindex(E).notna().all()
    assert (split.loc[E, "tier"] == "train").all()
    assert sig.str.startswith("human1:").all()
