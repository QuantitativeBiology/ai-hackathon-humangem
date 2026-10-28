"""W7 library tests (toy data only; no label file is opened)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from nextgem import w7


def _toy(seed=0, n_genes=30):
    rng = np.random.default_rng(seed)
    genes = list(range(1, n_genes + 1))
    tiers = {"train": [f"ACH-{i:06d}" for i in range(1, 13)], "validation": [f"ACH-{i:06d}" for i in range(13, 17)],
             "test": [f"ACH-{i:06d}" for i in range(17, 21)]}
    lines = sum(tiers.values(), [])
    tier_of = pd.Series({l: t for t, ls in tiers.items() for l in ls})
    lineage = pd.Series({l: f"L{i % 3}" for i, l in enumerate(lines)})

    def mat():
        return pd.DataFrame(rng.normal(size=(len(lines), n_genes)).astype("float32"), index=lines, columns=genes)

    omics = {k: mat() for k in ("expr", "cn", "mut_dam", "mut_hot")}
    gfeat = pd.DataFrame({"n_isozyme_genes_h2": rng.integers(0, 3, n_genes)}, index=pd.Index(genes, name="entrez"))
    lf = pd.DataFrame({"base_medium_class": rng.choice(["RPMI", "DMEM"], len(lines)), "serum": 1.0}, index=lines)
    gem_line = {"e_P1": mat().clip(0, 1), "e_Gstar": mat().clip(0, 1)}
    gem_gene = pd.DataFrame({"p57_H1": rng.integers(0, 2, n_genes).astype("float32")}, index=pd.Index(genes, name="entrez"))
    inp = w7.Inputs(genes=genes, omics=omics, aggs=(mat(), mat()), gfeat=gfeat, lf=lf,
                    class_codes={"DMEM": 0, "RPMI": 1}, gem_line=gem_line, gem_gene=gem_gene, tier_of=tier_of,
                    lineage_of=lineage)
    base = rng.random(n_genes)
    rows = []
    for l in tiers["train"]:
        y = (rng.random(n_genes) < base).astype(float)
        for g, v in zip(genes, y):
            rows.append({"ModelID": l, "entrez": g, "lineage": lineage[l], "L_prob": v, "gene_effect": -v})
    return inp, pd.DataFrame(rows), tiers


IDEA = w7.Idea("T-1", extra=("e_P1", "p57_H1"), hgb=tuple(sorted({**w7.HGB_M2, "max_iter": 20}.items())))


def test_fit_rejects_non_training_lines():
    inp, fit_long, tiers = _toy()
    bad = fit_long[fit_long.ModelID == tiers["train"][0]].assign(ModelID=tiers["validation"][0])
    with pytest.raises(w7.LeakageError):
        w7.fit(inp, pd.concat([fit_long, bad]), IDEA, seed=11)
    with pytest.raises(w7.LeakageError):
        w7.assert_training_only([tiers["test"][0]], inp.tier_of)


def test_fitting_never_touches_non_training_lines():
    """Poison every feature row of validation and test lines (omics, aggregates, GEM): the fitted model and its
    predictions on held-out training lines must be bitwise unchanged."""
    inp, fit_long, tiers = _toy()
    held = tiers["train"][-4:]
    fit_rows = fit_long[~fit_long.ModelID.isin(held)]
    p0 = w7.predict(inp, w7.fit(inp, fit_rows, IDEA, seed=11), held)
    other = tiers["validation"] + tiers["test"]
    rng = np.random.default_rng(99)
    for d in (inp.omics, inp.gem_line):
        for k in d:
            d[k].loc[other] = (rng.normal(size=(len(other), len(inp.genes))) * 1e3).astype("float32")
    inp.aggs = tuple(a.copy() for a in inp.aggs)
    for a in inp.aggs:
        a.loc[other] = np.nan
    p1 = w7.predict(inp, w7.fit(inp, fit_rows, IDEA, seed=11), held)
    assert np.array_equal(p0.to_numpy(), p1.to_numpy())


def test_held_out_training_rows_do_not_enter_fit():
    """Labels of held-out lines are not in fit_long; changing their features does not change the fit either."""
    inp, fit_long, tiers = _toy()
    held = tiers["train"][-4:]
    fit_rows = fit_long[~fit_long.ModelID.isin(held)]
    f0 = w7.fit(inp, fit_rows, IDEA, seed=11)
    probe = tiers["train"][:2]
    p0 = w7.predict(inp, f0, probe)
    inp.omics["expr"].loc[held] = 1e6
    p1 = w7.predict(inp, w7.fit(inp, fit_rows, IDEA, seed=11), probe)
    assert np.array_equal(p0.to_numpy(), p1.to_numpy())
    assert set(f0.fit_lines) == set(fit_rows.ModelID)


def test_resid_and_mono_kinds_run():
    inp, fit_long, tiers = _toy()
    for idea in (w7.Idea("T-2", kind="resid", hgb=IDEA.hgb), w7.Idea("T-3", mono=True, hgb=IDEA.hgb)):
        p = w7.predict(inp, w7.fit(inp, fit_long, idea, seed=11), tiers["validation"])
        assert p.shape == (4, len(inp.genes)) and np.isfinite(p.to_numpy()).all()


def test_stack_scores():
    p = pd.DataFrame([[0.2, 0.8, 0.5]], columns=[1, 2, 3])
    e = pd.DataFrame([[1.0, 0.0, 0.0]], columns=[1, 2, 3])
    s0 = w7.stack_scores(p, e, 0.0)
    assert list(s0.iloc[0].rank()) == list(p.iloc[0].rank())
    s = w7.stack_scores(p, e, 10.0)
    assert s.iloc[0].idxmax() == 1


def test_registered_ideas_are_logged():
    from pathlib import Path

    f = Path(w7.__file__).resolve().parents[2] / "autoresearch" / "ideas.tsv"
    logged = set(pd.read_csv(f, sep="\t").id)
    assert set(w7.IDEAS) <= logged


def test_tg_rule():
    lin = pd.Series({f"l{i}": f"K{i % 12}" for i in range(120)})
    d = pd.Series(0.01, index=lin.index)
    r = w7.tg(d, lin, sorted(lin.unique()), B=200, seed=11)
    assert r["tg_pass"] and r["folds_positive"] == 12
    r = w7.tg(d * 0.4, lin, sorted(lin.unique()), B=200, seed=11)
    assert not r["tg_pass"]  # below the 0.005 floor


def test_deviation_features_use_training_norms_only():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(w7.__file__).resolve().parents[2] / "scripts"))
    import w7_common as wc

    inp, _, tiers = _toy()
    split = pd.DataFrame({"ModelID": inp.tier_of.index, "tier": inp.tier_of.values})
    d0 = wc.deviation_features(inp.omics, inp.aggs, inp.genes, split)
    other = tiers["validation"] + tiers["test"]
    inp.omics["expr"].loc[other] = 1e4
    d1 = wc.deviation_features(inp.omics, inp.aggs, inp.genes, split)
    tr = tiers["train"]
    assert np.allclose(d0["expr_dev"].loc[tr].to_numpy(), d1["expr_dev"].loc[tr].to_numpy())


def test_gene_groups_partition():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(w7.__file__).resolve().parents[2] / "scripts"))
    import w7_genecv

    genes = list(range(100, 160))
    g = w7_genecv.gene_groups(genes, 11)
    assert sorted(sum(g, [])) == genes and len(g) == 5 and g == w7_genecv.gene_groups(genes, 11)


@pytest.mark.parametrize("idea", [w7.Idea("T-4", mono=True, extra=("knn_prior",), mono_extra=("knn_prior",),
                                          hgb=IDEA.hgb), w7.Idea("T-5", kind="knn_only")])
def test_knn_prior_ideas_never_touch_non_training_lines(idea):
    inp, fit_long, tiers = _toy()
    held = tiers["train"][-4:]
    fit_rows = fit_long[~fit_long.ModelID.isin(held)]
    p0 = w7.predict(inp, w7.fit(inp, fit_rows, idea, seed=11), held)
    other = tiers["validation"] + tiers["test"]
    inp.omics["expr"].loc[other] = np.float32(1e3)
    inp._knn_Z = None
    p1 = w7.predict(inp, w7.fit(inp, fit_rows, idea, seed=11), held)
    assert np.array_equal(p0.to_numpy(), p1.to_numpy())


def test_knn_prior_excludes_self_and_own_lineage():
    inp, fit_long, tiers = _toy()
    from nextgem import w1

    lines = tiers["train"]
    Yw = w1.wide(fit_long, "L_prob", lines, inp.genes)
    lin = fit_long.drop_duplicates("ModelID").set_index("ModelID").lineage
    # with k = 1 and a lineage exclusion, each line's prior equals the label row of a line from another lineage
    kp = w7.knn_prior(inp, Yw, lines, exclude_lineage=lin, fit_lineage=lin, k=1)
    for l in lines:
        match = [m for m in lines if np.array_equal(Yw.loc[m].to_numpy(), kp.loc[l].to_numpy())]
        assert match and all(lin[m] != lin[l] for m in match)
