"""§12 "Labels": Entrez join on a toy matrix; NaN counting; fold scoping (a memorising method gets null AP
on the held-out fold); kappa(L-prob, L-ge05) on training (positive control); label variants incl. L-ge06;
per-tier files and their FROZEN.tsv hashes; rebuild reproducibility."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import cohen_kappa_score

import toy


@pytest.fixture()
def toy_env(built_toy, monkeypatch):
    monkeypatch.setenv("NEXTGEM_DATA_ROOT", str(built_toy.data))
    monkeypatch.setenv("NEXTGEM_STATE_ROOT", str(built_toy.state))
    return built_toy


def test_build_outputs_and_frozen_rows(toy_env):
    t = toy_env
    for f in (t.data / "labels" / "labels_train.parquet", t.data / "labels" / "label_derived_lists.parquet",
              t.data / "labels" / "label_ids.json", t.state / "results" / "labels" / "dropcounts.tsv",
              t.state / "results" / "labels" / "build_report.json"):
        assert f.exists(), f
    from nextgem import frozen, labels

    for tier in ("validation", "test"):  # sealed files exist (existence check only: os.stat, not open)
        assert labels.tier_file(tier).exists() and labels.tier_file(tier).parent.name == "locked"
    assert labels.fingerprints_file().exists() and (t.data / "locked" / "HOOK_PROBE").exists()
    rows = {r["item"]: r for r in frozen.read_rows()}
    rep = json.loads((t.state / "results" / "labels" / "build_report.json").read_text())
    for tier in ("train", "validation", "test"):
        assert rows[f"labels.{tier}"]["sha256_or_value"] == rep["tier_sha256"][tier]
    assert rows["label.primary"]["sha256_or_value"] == "L_prob"
    # check_frozen.py (permitted to hash label files) passes on the toy project
    r = t.run("check_frozen.py")
    assert r.returncode == 0, r.stdout + r.stderr


def test_entrez_join(toy_env):
    from nextgem import labels

    assert labels.ids()["U_entrez"] == toy.U_ENTREZ  # 1000..1069; 1070..1079 not in U; "G75"/5555 not joined
    df = labels.train()
    assert sorted(df.entrez.unique()) == toy.U_ENTREZ
    assert 1005 in set(df.entrez)  # model symbol "G5_ALIAS" differs; joined on Entrez
    assert len(df) == len(toy_env.ids_in("train")) * len(toy.U_ENTREZ)


def test_label_variants_known_cells(toy_env):
    from nextgem import labels

    df = labels.train().set_index(["ModelID", "entrez"])
    c = df.loc[("ACH-000001", 1040)]  # ge -0.55, dp 0.5
    assert (c.L_ge05, c.L_ge06, c.L_ge10, c.L_prob) == (1, 0, 0, 0)
    c = df.loc[("ACH-000001", 1041)]  # ge -1.2, dp 0.9
    assert (c.L_ge05, c.L_ge06, c.L_ge10, c.L_prob) == (1, 1, 1, 1)
    c = df.loc[("ACH-000001", 1042)]  # ge -0.6 (strict <), dp 0.51
    assert (c.L_ge05, c.L_ge06, c.L_ge10, c.L_prob) == (1, 0, 0, 1)
    # every label equals its definition on the toy matrix
    ge = toy_env.ge.stack(future_stack=True).rename("ge")
    dp = toy_env.dp.stack(future_stack=True).rename("dp")
    ref = pd.concat([ge, dp], axis=1).loc[df.index]
    for name, (src, op, thr) in {"L_prob": ("dp", ">", .5), "L_ge05": ("ge", "<", -.5), "L_ge06": ("ge", "<", -.6),
                                  "L_ge10": ("ge", "<", -1.0)}.items():
        v = ref[src]
        expect = ((v > thr) if op == ">" else (v < thr)).astype(float).where(v.notna())
        got = df[name].astype("Float64").to_numpy(dtype=float, na_value=np.nan)
        assert np.array_equal(np.isnan(got), np.isnan(expect.to_numpy())) and np.allclose(got[~np.isnan(got)], expect.dropna())
    # the §6.8 grid helper agrees with L_ge06 at -0.6
    assert labels.ge_threshold_label(df, -0.6).astype("Float64").equals(df.L_ge06.reset_index(drop=True).astype("Float64"))


def test_nan_accounting(toy_env):
    dc = pd.read_csv(toy_env.state / "results" / "labels" / "dropcounts.tsv", sep="\t", dtype={"key": str})
    gene = dc[dc.scope == "gene"].set_index(["source", "key"]).n_nan
    assert gene[("gene_effect", "1010")] == 2 and gene[("gene_effect", "1020")] == 1
    assert gene[("dep_prob", "1030")] == 2 and gene[("dep_prob", "1031")] == 1
    assert gene.xs("gene_effect").sum() == 3  # (25, 75) is outside U
    assert gene.xs("dep_prob").sum() == 3
    tot = dc[dc.scope == "tier_total"].set_index(["source", "tier"]).n_nan
    assert tot[("gene_effect", "train")] == 3 and tot[("dep_prob", "train")] == 1
    assert tot[("dep_prob", "test")] == 2  # rows 25 and 27 are LinT
    lib = dc[dc.scope == "library_class"].groupby(["source", "key"]).n_nan.sum()
    assert lib[("gene_effect", "KY-only")] == 1  # row 2 is among the first 4 (KY-only) lines
    from nextgem import labels

    df = labels.train().set_index(["ModelID", "entrez"])
    assert df.loc[("ACH-000003", 1010), ["L_ge05", "L_ge06", "L_ge10"]].isna().all()
    assert not pd.isna(df.loc[("ACH-000003", 1010), "L_prob"])


def test_kappa_rule_positive_control_and_ledger(toy_env):
    from nextgem import labels, ledger

    tr = labels.train()
    ok = tr.L_prob.notna() & tr.L_ge05.notna()
    expected = cohen_kappa_score(tr.L_prob[ok].astype(int), tr.L_ge05[ok].astype(int))
    rows = ledger.read("validation")
    assert [r["event"] for r in rows] == ["non_look"]
    d = json.loads(rows[0]["detail"])
    assert abs(d["kappa"] - expected) < 1e-6 and d["kappa"] >= 0.5 and d["primary_label"] == "L_prob"
    assert ledger.verify("validation") == []
    log = toy.git(toy_env.state, "log", "--format=%s").stdout
    assert "ledger(validation): non_look kappa_rule" in log


def test_kappa_rule_switches_when_low(tmp_path_factory):
    t = toy.make_toy(tmp_path_factory.mktemp("toylow"), kappa="low")
    r = toy.build(t)
    assert r.returncode == 0, r.stderr[-1500:]
    from nextgem import frozen, paths

    rows = {x["item"]: x for x in frozen.read_rows(t.state / "FROZEN.tsv")}
    assert rows["label.primary"]["sha256_or_value"] == "L_ge05"


def test_rebuild_reproduces_hashes_and_overwrites_nothing(toy_env):
    from nextgem import ledger

    before = toy_env.state.joinpath("FROZEN.tsv").read_text()
    r = toy.build(toy_env)
    assert r.returncode == 0, r.stderr[-1500:]
    assert "reproduced FROZEN.tsv hashes exactly" in r.stdout
    assert toy_env.state.joinpath("FROZEN.tsv").read_text() == before
    assert len(ledger.read("validation")) == 1


def test_label_derived_lists_need_purpose(toy_env):
    from nextgem import labels

    with pytest.raises(PermissionError):
        labels.label_derived_list("depmap_inferred_common_essentials", purpose="fitting")
    assert labels.label_derived_list("depmap_inferred_common_essentials", purpose="reporting") == {1000, 1001, 1075}
    rep = json.loads((toy_env.state / "results" / "labels" / "build_report.json").read_text())
    assert rep["label_derived_lists_in_U"]["depmap_inferred_common_essentials"] == 2
    assert rep["id_sets"]["rows_equal"] and rep["id_sets"]["gene_columns_equal"]


def test_sealed_tiers_unreachable_from_tests(toy_env):
    """Negative control: validation/test labels, raw files and fingerprints from a test module."""
    from nextgem import labels

    for fn in (labels.validation, labels.test, labels.load_fingerprints,
               lambda: labels._read_raw("gene_effect", ["G1 (1001)"]),
               lambda: open(toy_env.data / "depmap" / toy.RAW_GE),
               lambda: pd.read_parquet(labels.tier_file("test"))):
        with pytest.raises(PermissionError):
            fn()


def test_fold_scoping_memoriser_gets_null_ap(toy_env):
    """A memorising method (predicts each (line, gene) label it has seen) gets exactly the random-ranking
    AP on every held-out line; the same method given all training labels (leaky control) gets AP = 1."""
    from nextgem import cv, labels

    def memoriser(fit, heldout_ids):
        mem = {(m, e): v for m, e, v in zip(fit.ModelID, fit.entrez, fit.L_prob.fillna(0).astype(int))}
        genes = sorted(fit.entrez.unique())
        rows = [(m, g, mem.get((m, g), 0)) for m in heldout_ids for g in genes]
        return pd.DataFrame(rows, columns=["ModelID", "entrez", "score"])

    def per_line(preds, held):
        m = preds.merge(held[["ModelID", "entrez", "L_prob"]].dropna(), on=["ModelID", "entrez"])
        return {mid: (toy.tap(g.score, g.L_prob.astype(int)), toy.null_ap(g.L_prob.astype(int)))
                for mid, g in m.groupby("ModelID") if g.L_prob.sum() > 0}

    res = cv.lolo(memoriser, per_line, columns=["L_prob"])
    assert set(res) == {"LinA", "LinB", "LinC"}  # rare lineage is never a fold
    vals = [v for d in res.values() for v in d.values()]
    assert vals and all(abs(a - b) < 1e-12 for a, b in vals)

    def leaky(fit, heldout_ids):
        return memoriser(labels.train(columns=["L_prob"]), heldout_ids)

    res = cv.lolo(leaky, per_line, columns=["L_prob"])
    vals = [v for d in res.values() for v in d.values()]
    assert all(abs(a - 1.0) < 1e-12 and b < 1 for a, b in vals)


def test_nested_lolo_never_shows_outer_fold(toy_env):
    from nextgem import cv

    seen = []

    def fp(fit, ids):
        seen.append((set(fit.lineage), set(ids)))
        return pd.DataFrame(columns=["ModelID", "entrez", "score"])

    cv.lolo(fp, lambda p, h: set(h.lineage), outer="LinA")
    assert len(seen) == 2
    for lins, ids in seen:
        assert "LinA" not in lins and not (ids & set(toy_env.ids_in("train")[:6]))


def test_tap_helper_known_answers():
    # brute force over all orderings of a tied block
    import itertools

    s = [3, 2, 2, 2, 1]
    y = [1, 0, 1, 1, 0]
    aps = []
    for perm in itertools.permutations([1, 2, 3]):
        order = [0, *perm, 4]
        yy = [y[i] for i in order]
        hits = np.cumsum(yy)
        aps.append(sum(hits[k] / (k + 1) for k in range(5) if yy[k]) / sum(yy))
    assert abs(toy.tap(s, y) - np.mean(aps)) < 1e-12
    assert abs(toy.tap([1, 0, 0], [1, 0, 0]) - 1.0) < 1e-12
