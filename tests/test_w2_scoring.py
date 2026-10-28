"""Toy tests for nextgem.w2_score (W2 M4 scoring helpers). No label file is read."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from nextgem import w2_score as ws


def test_parse_step_families_and_gate():
    assert ws.parse_step("M0").family == "reference"
    assert ws.parse_step("M-open").family == "negative_control"
    s = ws.parse_step("M2g(f=0.1)")
    assert s.gated and s.f == 0.1 and not s.eligible and s.gate_free_twin == "M2(f=0.1)"
    assert ws.parse_step("M1g").gate_free_twin == "M1"
    h = ws.parse_step("M-HPLM:M2g(f=1)")
    assert h.family == "physiological" and h.gated and h.gate_free_twin == "M-HPLM:M2(f=1)" and not h.eligible
    assert ws.parse_step("M2(f=1)").eligible and ws.parse_step("M1").eligible


def test_eligible_set_is_documented_gate_free_only():
    assert set(ws.ELIGIBLE) == {"M1", "M2(f=0.001)", "M2(f=0.01)", "M2(f=0.1)", "M2(f=1)"}
    assert all(not ws.parse_step(s).gated for s in ws.ELIGIBLE)
    assert ws.DOC_LADDER[0] == "M0"


def test_step_order():
    names = ["M2(f=1)", "M-Plasmax:M2(f=1)", "M1g", "M0", "M-open", "M1", "M2(f=0.001)", "M-HPLM:M2(f=1)"]
    o = ws.step_order(names)
    assert o[:4] == ["M-open", "M0", "M1", "M2(f=0.001)"]
    assert o.index("M1g") < o.index("M-HPLM:M2(f=1)") < o.index("M-Plasmax:M2(f=1)")


def test_sweep_e_matrix_and_line_matrix():
    ratio = pd.DataFrame({"g1": [1.0, 0.0], "g2": [0.5, 1.0], "g3": [0.0, 0.0]}, index=["s1", "s2"])
    emap = pd.Series([10, 20], index=["g1", "g2"])  # g3 has no Entrez -> dropped
    E = ws.sweep_e_matrix(ratio, emap, [10, 20, 30])
    assert E.loc["s1"].tolist() == [0.0, 0.5, 0.0] and E.loc["s2"].tolist() == [1.0, 0.0, 0.0]
    M = ws.line_matrix(pd.Series({"L1": "s2", "L2": "s1"}), E, ["L1", "L2"])
    assert M.loc["L1", 10] == 1.0 and M.loc["L2", 20] == 0.5
    with pytest.raises(KeyError):
        ws.line_matrix(pd.Series({"L1": "s2"}), E, ["L1", "L2"])


def test_gene_to_entrez_first_per_entrez():
    ko = pd.DataFrame({"gene": ["a", "b", "c", "d"], "entrez": [1.0, np.nan, 1.0, 2.0]})
    m = ws.gene_to_entrez(ko)
    assert m.to_dict() == {"a": 1, "d": 2}


def test_folds_positive_small_fold_not_positive():
    lin = pd.Series({f"L{i}": ("A" if i < 5 else "B" if i < 9 else "C") for i in range(12)})
    d = pd.Series(0.01, index=lin.index)
    r = ws.folds_positive(d, lin, ["A", "B", "C", "D"])
    # A: 5 lines positive; B: 4 lines (< 5) not positive; C: 3 lines not positive; D: absent
    assert r["folds_positive"] == 1 and r["n_folds"] == 4
    assert json.loads(r["fold_sizes"]) == {"A": 5, "B": 4, "C": 3, "D": 0}


def test_tg_pass_rule():
    assert ws.tg_pass(0.005, 0.001, 9)
    assert not ws.tg_pass(0.0049, 0.001, 12)
    assert not ws.tg_pass(0.01, 0.0, 12)
    assert not ws.tg_pass(0.01, 0.001, 8)
    assert not ws.tg_pass(np.nan, 0.1, 12)


def test_pick_shortlist_one_se_prefers_simplest():
    rows = pd.DataFrame({"step": ["M1", "M2(f=0.1)", "M2(f=1)"], "delta_lin": [0.010, 0.012, 0.011],
                         "tg_pass": [True, True, True]})
    se = {"M1": (0.002, 0.003), "M2(f=0.1)": (0.0, 0.0), "M2(f=1)": (0.001, 0.003)}
    d = ws.pick_shortlist(rows, se)
    assert d["shortlist"] == ["M1"] and d["best_by_delta"] == "M2(f=0.1)"
    se["M1"] = (0.002, 0.001)  # M1 now outside one SE -> best M2 within the band (larger delta first)
    assert ws.pick_shortlist(rows, se)["shortlist"] == ["M2(f=0.1)"]
    none = ws.pick_shortlist(rows.assign(tg_pass=False), se)
    assert none["shortlist"] == [] and "M0 is kept" in none["decision"]


def test_flip_counts():
    Y = pd.DataFrame([[1, 0, np.nan], [0, 1, 1]], index=["a", "b"], columns=[1, 2, 3], dtype=float)
    ref = pd.DataFrame([[False, True, True], [False, False, True]], index=Y.index, columns=Y.columns)
    new = pd.DataFrame([[True, True, False], [True, False, True]], index=Y.index, columns=Y.columns)
    fc = ws.flip_counts(ref, new, Y)
    # ref errors: a1 (FN), a2 (FP), b2 (FN) = 3; new fixes a1; breaks b1 (FP); a3 is NaN (ignored)
    assert fc["n_pairs"] == 5 and fc["errors_ref"] == 3 and fc["fp_ref"] == 1 and fc["fn_ref"] == 2
    assert fc["fixed"] == 1 and fc["broken"] == 1 and fc["fixed_fn"] == 1 and fc["fixed_fp"] == 0
    assert fc["net_flip_fraction"] == 0.0
    assert fc["per_line_fixed"].tolist() == [1, 0] and fc["per_line_broken"].tolist() == [0, 1]


def test_oracle_ratio_point_estimate():
    lin = pd.Series({"a": "X", "b": "X", "c": "Y"})
    d_doc = pd.Series({"a": 0.01, "b": 0.03, "c": 0.0})
    d_or = pd.Series({"a": 0.02, "b": 0.04, "c": 0.02})
    r = ws.oracle_ratio(d_doc, d_or, lin, B=200, seed=11)
    # lineage-equal: doc = (0.02 + 0) / 2 = 0.01; oracle = (0.03 + 0.02) / 2 = 0.025
    assert r["doc_minus_m0"] == pytest.approx(0.01) and r["oracle_minus_m0"] == pytest.approx(0.025)
    assert r["ratio"] == pytest.approx(0.4) and r["ratio_lo"] <= r["ratio"] <= r["ratio_hi"]


def test_shuffle_eligible_lineages():
    keys = pd.Series(["RPMI"] * 8 + ["DMEM+pyruvate"] * 8 + ["RPMI"] * 9 + ["DMEM"] * 7,
                     index=[f"L{i}" for i in range(32)])
    lin = pd.Series(["A"] * 16 + ["B"] * 16, index=keys.index)
    assert ws.shuffle_eligible_lineages(keys, lin) == ["A"]


def test_shuffled_sweeps_lookup():
    lab = pd.DataFrame({"parent": ["h", "h", "h", "h"], "variant": ["M1", "M2", "M2", "M2g"],
                        "f": [np.nan, 0.1, 1.0, 0.1], "gate_open": [True, True, True, False],
                        "medium_key": ["RPMI", "RPMI", "RPMI", "RPMI"], "sweep_id": ["s1", "s2", "s3", "s4"]})
    sk = pd.Series({"L1": "RPMI"})
    assert ws.shuffled_sweeps(sk, "M2(f=0.1)", "h", lab).tolist() == ["s2"]
    assert ws.shuffled_sweeps(sk, "M1", "h", lab).tolist() == ["s1"]
    with pytest.raises(KeyError):
        ws.shuffled_sweeps(pd.Series({"L1": "DMEM"}), "M1", "h", lab)


def test_addition_prune_safe_bound():
    # non-lethal genes have ratio >= 0.5 on the reference: safe while WT grows at most 50-fold
    assert ws.addition_prune_safe(10.0, 1.0, 0.5)
    assert ws.addition_prune_safe(50.0, 1.0, 0.5)
    assert not ws.addition_prune_safe(50.01, 1.0, 0.5)
    assert not ws.addition_prune_safe(1.0, 0.0, 0.5)
    # worst case at the bound: ratio_new >= 0.5 * 1 / 50 = 0.01 -> still non-lethal under the strict < 0.01 call
    assert 0.5 * 1.0 / 50.0 >= 0.01
