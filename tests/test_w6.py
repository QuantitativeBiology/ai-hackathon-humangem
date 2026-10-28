"""Toy / known-answer tests for the W6 library (src/nextgem/w6.py). Label-free, GLPK, fast.

Known-answer controls (DESIGN §12 W6 row): synthetic information injection is recovered (positive), a block that is a
linear function of expression adds nothing (negative), permuting the flux within lineage removes the injected signal.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from nextgem import w6

ROOT = Path(__file__).resolve().parents[1]


# ----------------------------------------------------------------------------------------------- compound list
def test_w6_compound_table_matches_apriori_list():
    text = (ROOT / w6.APRIORI_FILE).read_text()
    t = w6.compound_table(text)
    assert len(t) == len(w6.COMPOUNDS) == 105
    assert t.off_model_control.sum() == 8
    assert set(t[t["class"] == "Statins and FDPS"].targets) == {"HMGCR", "FDPS"}
    assert "trifluridine+tipiracil" in set(t.compound)
    assert not any("KPT-9274" in c or "Negative" in c or "PMID" in c for c in t.compound)


def test_w6_parser_handles_groups_citations_and_prose():
    text = ("x\n# 4. Final metabolic-drug list\n\n| Class | Compounds |\n|---|---|\n"
            "| A | foo (G1; → G2), bar and baz (G3) [LIT X 2020, PMID 1]; qux (G4). Note excluded |\n"
            "| B | one/two, three. **Prose**: words, more |\n\n**Positive-control pairs**: G1–foo.\n\n# 5. next\n")
    p = w6.parse_apriori_compounds(text)
    assert list(p.token) == ["foo", "bar", "baz", "qux", "one/two", "three"]
    assert list(p["class"]) == ["A"] * 4 + ["B"] * 2
    with pytest.raises(ValueError):
        w6.compound_table(text)  # tokens differ from the transcription -> refuse


def test_w6_positive_control_pairs_reference_listed_compounds():
    text = (ROOT / w6.APRIORI_FILE).read_text()
    para = w6.positive_control_text(text)
    for genes, comps, _ in w6.POSITIVE_CONTROL_PAIRS:
        for c in comps:
            assert c in w6.COMPOUNDS, c
        for g in genes:
            stem = re.sub(r"\d+$", "", g)
            assert g in para or stem in para, g
    assert len(w6.positive_control_pairs_long()) == sum(len(g) * len(c) for g, c, _ in w6.POSITIVE_CONTROL_PAIRS)


def test_w6_name_matching_prefers_alias_and_most_lines():
    idx = w6.name_index([("A", "Fluorouracil, 5-FU"), ("B", "5-fluorouracil"), ("C", "methotrexate"), ("D", "Methotrexate")])
    m = w6.match_compounds(["5-FU", "methotrexate", "SHIN1"], idx, {"A": 10, "B": 50, "C": 5, "D": 7}).set_index("compound")
    assert m.loc["5-FU", "screen_id"] == "A"  # first alias 'fluorouracil' matches A only
    assert m.loc["methotrexate", "screen_id"] == "D" and m.loc["methotrexate", "n_ids"] == 2
    assert m.loc["SHIN1", "screen_id"] is None
    assert w6.norm_name("α-Lipoic acid") == "alphalipoicacid"


def test_w6_assay_medium_maps():
    lines = ["L1", "L2", "L3"]
    culture = pd.Series({"L1": "DMEM", "L2": "RPMI"})
    gd = pd.Series({"L1": "R", "L2": "D/F12", "L3": None})
    assert set(w6.assay_medium_map("prism", lines, culture)) == {"RPMI"}
    assert w6.assay_medium_map("gdsc", lines, culture, gd).to_dict() == {"L1": "RPMI", "L2": "DMEM:F12"}
    assert w6.assay_medium_map("culture", lines, culture).to_dict() == {"L1": "DMEM", "L2": "RPMI"}


# ----------------------------------------------------------------------------------------------- PCA engine
def test_w6_state_space_pca_equals_sklearn_on_repeated_rows():
    from sklearn.decomposition import PCA

    rng = np.random.default_rng(11)
    S = pd.DataFrame(rng.normal(size=(7, 30)), index=[f"s{i}" for i in range(7)])
    b = w6.Block.from_matrix("toy", S)
    tr = rng.integers(0, 7, size=40)
    allx = np.concatenate([tr, rng.integers(0, 7, size=10)])
    Z = w6.pca_fit_transform(b.G, tr, allx, 5)
    ref = PCA(5, svd_solver="full").fit(S.to_numpy()[tr]).transform(S.to_numpy()[allx])
    assert Z.shape == ref.shape
    for j in range(5):
        assert np.allclose(np.abs(Z[:, j]), np.abs(ref[:, j]), atol=1e-8)


def test_w6_constant_block_gives_no_components():
    S = pd.DataFrame(np.ones((4, 3)), index=list("abcd"))
    b = w6.Block.from_matrix("c", S)
    assert b.constant
    assert w6.pca_fit_transform(b.G, np.array([0, 1, 2]), np.array([0, 3]), 20).shape == (2, 0)
    one_state = w6.Block.from_matrix("x", pd.DataFrame(np.random.default_rng(11).normal(size=(3, 4)), index=list("abc")))
    assert w6.pca_fit_transform(one_state.G, np.array([1, 1, 1]), np.array([0, 1]), 20).shape == (2, 0)


def test_w6_permute_within_lineage_preserves_multisets():
    lines = [f"l{i}" for i in range(30)]
    lin = pd.Series(["A"] * 10 + ["B"] * 12 + ["C"] * 8, index=lines)
    st = pd.Series([f"s{i % 4}" for i in range(30)], index=lines)
    p = w6.permute_within_lineage(st, lin, 11)
    for L in "ABC":
        assert sorted(p[lin == L]) == sorted(st[lin == L])
    assert not p.equals(st)
    assert p.equals(w6.permute_within_lineage(st, lin, 11))


def test_w6_reaction_capacity_matches_context_capacity():
    from nextgem.context import capacity

    rng = np.random.default_rng(11)
    gprs = {"R1": ("and", ("g1", "g2")), "R2": ("or", ("g1", ("and", ("g2", "g3")))), "R3": ("or", ("g1", "g9")),
            "R4": ("and", ("g9", "g3")), "R5": "g2"}
    lin = pd.DataFrame(rng.uniform(0, 50, size=(6, 3)), columns=["g1", "g2", "g3"])
    cap = w6.reaction_capacity(gprs, lin)
    assert "R3" not in cap  # OR with an unknown child is unknown
    for i in range(6):
        d = lin.iloc[i].to_dict()
        for rid in cap:
            assert cap.loc[i, rid] == pytest.approx(capacity(gprs[rid], d))


# ----------------------------------------------------------------------------------------------- statistics
def test_w6_holm_and_wilcoxon_known_values():
    h = w6.holm({"a": 0.01, "b": 0.04, "c": 0.03})
    assert h == pytest.approx({"a": 0.03, "c": 0.06, "b": 0.06})
    assert w6.wilcoxon_greater([0, 0, 0]) == 1.0
    assert w6.wilcoxon_greater(np.linspace(0.1, 1, 12)) < 0.001
    assert w6.wilcoxon_greater(-np.linspace(0.1, 1, 12)) > 0.99


def test_w6_decision_rows_verdicts():
    rng = np.random.default_rng(11)
    pos = pd.Series(0.1 + 0.01 * rng.normal(size=40))
    null = pd.Series(0.01 * rng.normal(size=40))
    perm = {s: pd.Series(0.01 * rng.normal(size=40)) for s in w6.PERM_SEEDS}
    t = w6.decision_rows({"x": pos, "y": null, "z": pd.Series(np.zeros(30))},
                         {"x": perm, "y": perm, "z": perm}).set_index("analysis")
    assert t.loc["x", "verdict"] == "supported" and t.loc["x", "above_permuted"]
    assert t.loc["z", "verdict"] == "not supported" and t.loc["z", "p_wilcoxon_one_sided"] == 1.0
    assert t.loc["y", "p_holm"] >= t.loc["y", "p_wilcoxon_one_sided"]


def test_w6_fit_free_rho_constant_within_lineage_is_nan():
    lines = [f"l{i}" for i in range(40)]
    lin = pd.Series(["A"] * 20 + ["B"] * 20, index=lines)
    y = pd.Series(np.arange(40.0), index=lines)
    const = pd.Series([0.0] * 20 + [1.0] * 20, index=lines)  # varies across, not within, lineages
    assert np.isnan(w6.fit_free_rho(const, y, lin)[0])
    r, n = w6.fit_free_rho(y * 2, y, lin)
    assert r == pytest.approx(1.0, abs=1e-5) and n == 2


# ----------------------------------------------------------------------------------------------- known-answer controls
def _toy_world(seed=11, n_lin=6, per=40, n_states=8, p_rna=10):
    rng = np.random.default_rng(seed)
    lines = [f"ACHX-{i:04d}" for i in range(n_lin * per)]
    lineage = pd.Series(np.repeat([f"L{j}" for j in range(n_lin)], per), index=lines)
    rna = pd.DataFrame(rng.normal(size=(len(lines), p_rna)), index=lines)
    S = pd.DataFrame(rng.normal(size=(n_states, 25)), index=[f"st{i}" for i in range(n_states)])
    state_of = pd.Series(rng.choice(S.index, size=len(lines)), index=lines)  # medium-like, independent of RNA
    return rng, lines, lineage, rna, S, state_of


def _views(rna, S, state_of, extra=None):
    br = w6.Block.from_matrix("RNA", rna)
    bf = w6.Block.from_matrix("FLUX", S)
    v_r = w6.View(br, pd.Series(rna.index, index=rna.index), "RNA")
    v_f = w6.View(bf, state_of, "FLUX")
    uni = {"RNA": list(rna.index), "FLUX": list(state_of.index)}
    out = {"RNA": v_r, "FLUX": v_f}
    if extra is not None:
        bn = w6.Block.from_matrix("NEG", extra)
        out["NEG"] = w6.View(bn, pd.Series(extra.index, index=extra.index), "NEG")
        uni["NEG"] = list(extra.index)
    return out, uni


def test_w6_synthetic_positive_control_is_recovered_and_permutation_removes_it():
    rng, lines, lineage, rna, S, state_of = _toy_world()
    v, uni = _views(rna, S, state_of)
    cache = w6.PCACache()
    real, permd = [], []
    vp = w6.View(v["FLUX"].block, w6.permute_within_lineage(state_of, lineage, 23), "FLUXperm")
    uni["FLUXperm"] = list(lines)
    for i in range(8):
        z = w6.random_state_direction(S, 100 + i)
        y = w6.synthetic_positive_response(z, state_of, 200 + i) + 0.3 * rna.iloc[:, 0]
        base = w6.heldout_scores(y, w6.lolo_ridge(y, [v["RNA"]], uni, lineage, cache), lineage)["rho"]
        withf = w6.heldout_scores(y, w6.lolo_ridge(y, [v["RNA"], v["FLUX"]], uni, lineage, cache), lineage)["rho"]
        withp = w6.heldout_scores(y, w6.lolo_ridge(y, [v["RNA"], vp], uni, lineage, cache), lineage)["rho"]
        real.append(withf - base)
        permd.append(withp - base)
    real, permd = np.array(real), np.array(permd)
    assert np.median(real) > 0.2
    assert w6.wilcoxon_greater(real) < 0.01
    assert np.median(permd) < np.median(real) - 0.15


def test_w6_negative_control_linear_of_expression_adds_nothing():
    rng, lines, lineage, rna, S, state_of = _toy_world(p_rna=10)
    neg = w6.linear_of_expression(rna, 11, dim=40)  # rank <= 10 < PCA_K: carries only RNA information
    v, uni = _views(rna, S, state_of, extra=neg)
    cache = w6.PCACache()
    d = []
    for i in range(10):
        y = rna.to_numpy() @ rng.normal(size=10) + rng.normal(size=len(lines))
        y = pd.Series(y, index=lines)
        base = w6.heldout_scores(y, w6.lolo_ridge(y, [v["RNA"]], uni, lineage, cache), lineage)["rho"]
        withn = w6.heldout_scores(y, w6.lolo_ridge(y, [v["RNA"], v["NEG"]], uni, lineage, cache), lineage)["rho"]
        d.append(withn - base)
    assert np.max(np.abs(d)) < 0.02
    assert w6.wilcoxon_greater(d) > 0.01


def test_w6_line_constant_flux_adds_exactly_nothing():
    """A generic model under one assay medium (PRISM) gives one state for every line: Δρ is exactly 0."""
    rng, lines, lineage, rna, S, _ = _toy_world()
    v, uni = _views(rna, S, pd.Series("st0", index=lines))
    cache = w6.PCACache()
    y = pd.Series(rna.iloc[:, 0].to_numpy() + rng.normal(size=len(lines)), index=lines)
    a = w6.lolo_ridge(y, [v["RNA"]], uni, lineage, cache)
    b = w6.lolo_ridge(y, [v["RNA"], v["FLUX"]], uni, lineage, cache)
    assert np.allclose(a.to_numpy(), b.to_numpy())


def test_w6_subtype_scores_detect_flux_separable_classes():
    rng, lines, lineage, rna, S, state_of = _toy_world(n_lin=1, per=60, n_states=2)
    v, _ = _views(rna, S, state_of)
    lab = state_of.map({"st0": "a", "st1": "b"})
    f_flux = w6.subtype_scores(lab, [v["FLUX"]])
    f_rna = w6.subtype_scores(lab, [v["RNA"]])
    assert f_flux["knn_macro_f1"] == pytest.approx(1.0)
    assert f_rna["knn_macro_f1"] < 0.8
    assert f_flux["silhouette"] > f_rna["silhouette"]


# ----------------------------------------------------------------------------------------------- simulation
def _toy_model():
    import cobra

    cobra.Configuration().solver = "glpk"
    m = cobra.Model("toy")
    A_e, A, B, X = (cobra.Metabolite(i, compartment=c) for i, c in (("A_e", "e"), ("A_c", "c"), ("B_c", "c"), ("X_c", "c")))
    ex = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({A_e: -1})
    tr = cobra.Reaction("TR_A", lower_bound=0, upper_bound=1000)
    tr.add_metabolites({A_e: -1, A: 1})
    r1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    r1.gene_reaction_rule = "g1"
    r2 = cobra.Reaction("R2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({A: -2, B: 1})  # costlier bypass of R1
    r2.gene_reaction_rule = "g2"
    rx = cobra.Reaction("RX", lower_bound=0, upper_bound=1000)
    rx.add_metabolites({B: -1, X: 1})
    rx.gene_reaction_rule = "g3"
    bio = cobra.Reaction("MAR13082", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({B: -1, X: -1})
    m.add_reactions([ex, tr, r1, r2, rx, bio])
    m.objective = "MAR13082"
    return m


def test_w6_pfba_and_koflux_shortcut_are_exact():
    from nextgem import w3

    m = _toy_model()
    names = w6.reaction_var_names(m)
    genes = sorted(g.id for g in m.genes)
    clos = w3.gene_closures(m)
    sym = {"S1": "g1", "S2": "g2", "S3": "g3"}
    r = w6.simulate_state(m, names, genes, clos, {"EX_A": 10.0}, w3.Variant("T0"), ["S1", "S2", "S3"], sym, "MAR13082")
    rid = [n[0] for n in names]
    v = dict(zip(rid, r["flux"]))
    assert r["wt_growth"] == pytest.approx(5.0)
    assert v["R2"] == pytest.approx(0.0) and v["R1"] == pytest.approx(10.0)
    assert r["n_shortcut_koflux"] == 1  # g2 (R2 carries no flux): WT pFBA kept, exact
    assert np.allclose(r["koflux"]["S2"], r["flux"])
    with m:
        m.medium = {"EX_A": 10.0}
        m.genes.get_by_id("g1").knock_out()
        mu, vk, _ = w6.pfba(m, names, "MAR13082")
    assert np.allclose(r["koflux"]["S1"], vk) and r["ko_growth"]["S1"] == pytest.approx(mu) == pytest.approx(2.5)
    assert r["ko_growth"]["S3"] == 0.0 and not np.any(r["koflux"]["S3"])
    ratio = dict(zip(genes, r["ratio_r"]))
    assert ratio == pytest.approx({"g1": round(mu / 5.0, 6), "g2": 1.0, "g3": 0.0})


def test_w6_feature_transforms():
    v = np.array([0.0, 2.0, -4.0])
    assert np.allclose(w6.flux_features(v, 2.0), np.arcsinh([0, 1, -2]))
    assert np.allclose(w6.koflux_features(np.zeros(3), 2.0), 0.0)
    X = pd.DataFrame({"a": [1.0, np.nan, 3.0, 4.0], "b": [np.nan, 2.0, 2.0, 2.0]})
    Y = w6.impute_line_floor(X, q=0.0)
    assert Y.loc[0, "b"] == 1.0 and Y.loc[1, "a"] == 2.0


def test_w6_code_has_no_unseeded_randomness():
    for p in [ROOT / "src/nextgem/w6.py", *sorted((ROOT / "scripts").glob("w6_*.py"))]:
        s = p.read_text()
        assert "default_rng()" not in s and "random_state=None" not in s, p
        label_file_stem = "CRISPR" + "Gene"  # split so that shell tools quoting this file never name a label file
        assert label_file_stem not in s and "locked_dir" not in s and "nextgem.labels" not in s, p
        assert "import labels" not in s and "labels.train" not in s and "labels.fold" not in s, p


def test_w6_ridge_loo_matches_sklearn_ridgecv():
    from sklearn.linear_model import RidgeCV

    rng = np.random.default_rng(11)
    for n, p in ((80, 5), (200, 40), (30, 60)):
        X = rng.normal(size=(n, p))
        y = X[:, 0] - 0.5 * X[:, 1] + rng.normal(size=n)
        Xte = rng.normal(size=(17, p))
        pred, a = w6.ridge_loo_fit_predict(X, y, Xte)
        m = RidgeCV(alphas=w6.RIDGE_ALPHAS).fit(X, y)
        assert a == pytest.approx(m.alpha_)
        assert np.allclose(pred, m.predict(Xte), atol=1e-8)
