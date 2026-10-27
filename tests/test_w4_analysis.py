"""Known-answer tests for the Stage-A analysis library (src/nextgem/w4_analysis.py). No label file is read."""

import numpy as np
import pandas as pd
import pytest

from nextgem import w4_analysis as wa


def test_w4a_gates_known_answer():
    b = pd.DataFrame({"failed": [False, True, False, False], "wt_growth": [1.0, 0.0, 1.0, 1.0],
                      "task_profile": ["1111", "1111", "1011", "1111"],
                      "repair_frac_reactions": [0.01, 0.0, 0.0, 0.06], "cfg": ["a", "a", "b", "b"]})
    g = wa.build_gates(b, viability_idx=[0, 1])
    assert g.gate_all.tolist() == [True, False, False, False]
    assert g.gate_viability.tolist() == [True, True, False, True]
    c = wa.config_gates(g, ["cfg"]).set_index("cfg")
    assert c.loc["a", "frac_all"] == 0.5 and not c.loc["a", "config_gate_pass"]


def test_w4a_jaccard_and_r2_known_answer():
    X = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1]], bool)
    D2 = wa.jaccard_distance_sq(X)
    assert D2[0, 1] == 0 and D2[0, 2] == 1
    assert wa.r2(D2, ["a", "a", "b", "b"]) == pytest.approx(1.0)  # groups explain everything
    assert wa.r2(D2, ["a", "b", "a", "b"]) == pytest.approx(0.0)
    # partial overlap: rows {1,2} vs {2,3}: J = 1/3
    assert wa.jaccard_distance_sq(np.array([[1, 1, 0], [0, 1, 1]], bool))[0, 1] == pytest.approx(2 / 3)


def test_w4a_permanova_detects_structure_and_not_noise():
    rng = np.random.default_rng(3)
    base = rng.random((2, 200)) < 0.5
    X = np.vstack([base[i % 2] ^ (rng.random(200) < 0.05) for i in range(20)])
    meta = pd.DataFrame({"grp": [i % 2 for i in range(20)], "noise": rng.permutation([0, 1] * 10)})
    D2 = wa.jaccard_distance_sq(X)
    assert wa.permanova(D2, meta, "grp", perms=199)["p_perm"] < 0.01
    assert wa.permanova(D2, meta, "noise", perms=199)["p_perm"] > 0.05


def test_w4a_tg_panel_rule():
    lin = pd.Series([f"L{i // 3}" for i in range(30)], index=[f"l{i}" for i in range(30)])
    d = pd.Series(0.02, index=lin.index)
    t = wa.tg_panel(d, lin, B=200)
    assert t["tg_pass"] and t["folds_positive"] == 10
    d2 = d.copy()
    d2[lin.isin(["L0", "L1", "L2"])] = -0.001  # 7 of 10 folds positive -> fails the >= 8 bar
    assert not wa.tg_panel(d2, lin, B=200)["tg_pass"]
    d3 = d.copy()
    d3[["l0", "l1"]] = np.nan  # fold L0 has 1 estimable line -> "not positive"; 9 of 10 still pass
    t3 = wa.tg_panel(d3, lin, B=200)
    assert t3["folds_positive"] == 9 and t3["tg_pass"]
    assert not wa.tg_panel(pd.Series(0.004, index=lin.index), lin, B=200)["tg_pass"]  # below 0.005


def test_w4a_nested_selection_never_uses_heldout():
    lin = pd.Series(["A"] * 3 + ["B"] * 3 + ["C"] * 3, index=list("abcdefghi"))
    # config x is best everywhere except lineage A where y is best; A's choice must come from B, C only
    D = pd.DataFrame({"x": [0.0] * 3 + [0.05] * 6, "y": [0.9] * 3 + [0.01] * 6}, index=lin.index)
    outer, chosen = wa.nested_lolo_select(D, lin, order=["x", "y"])
    assert chosen["A"] == "x" and outer["a"] == 0.0
    assert chosen["B"] == "y"  # inner (A, C): y mean (0.9 + 0.01)/2 > x (0 + 0.05)/2
    tie = pd.DataFrame({"x": [0.1] * 9, "y": [0.1] * 9}, index=lin.index)
    assert set(wa.nested_lolo_select(tie, lin, order=["y", "x"])[1].values()) == {"y"}  # ties: simplicity order


def test_w4a_score_frame_convention():
    emap = pd.Series({1: "gA", 2: "gB"})
    ratio = pd.DataFrame({"gA": [0.25], "gB": [1.0]}, index=["L"])
    S = wa.score_frame(ratio, [1, 2, 3], emap)
    assert S.loc["L"].tolist() == [0.75, 0.0, 0.0]  # Entrez 3 absent from the parent -> 0
