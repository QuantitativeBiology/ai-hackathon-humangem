"""Tests of scripts/t2_prep_spec.py (P2 = W7-05 spec for T2). Label-free: toy frames only, plus static checks of
the written spec when it exists."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT = Path(__file__).resolve().parents[1]


def _load(name):
    s = importlib.util.spec_from_file_location(name, PROJECT / "scripts" / f"{name}.py")
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


tp = _load("t2_prep_spec")
et = tp._evaluate_test()  # one module instance: same PreflightError class


def _gen():
    idx = pd.Index([7, 3, 11, 5], name="entrez")
    return pd.DataFrame({"e_H1": [0.0, 1.0, 0.25, 0.5], "e_H2": [0.1, 0.0, 1.0, 0.0]}, index=idx)


def test_generic_scores_pass_harness_loader_and_equal_generic(tmp_path):
    from nextgem import w1

    gen = _gen()
    df = tp.generic_scores(gen)
    assert list(df.columns) == ["entrez", "H1", "H2"] and list(df.entrez) == [3, 5, 7, 11]
    p = tmp_path / "s.parquet"
    df.to_parquet(p, index=False)
    got = et.ev.load_scores(p, "gene", ["H1", "H2"])  # the P2 score_methods when P2 is refitted
    lines, genes = ["ACH-000001", "ACH-000002"], list(gen.index)
    for m in ("H1", "H2"):
        M = et._tile(got.set_index("entrez")[m].astype(float), lines, genes)
        assert np.allclose(M.to_numpy(), w1._const(gen[f"e_{m}"], lines, genes).to_numpy(), rtol=0, atol=1e-12)


def _perline(method, seeds, taps, lines=("A", "B", "C")):
    rows = [{"method": method, "seed": s, "ModelID": l, "tap": t + 0.01 * s} for s in seeds for l, t in zip(lines, taps)]
    return pd.DataFrame(rows)


def test_line_deltas_uses_seed_for_fitted_and_single_seed_comparator():
    tabs = {"p2": _perline("W7-05", [11, 23], [0.9, 0.8, 0.7]), "gstar": _perline("Gstar", [0], [0.3, 0.2, 0.1])}
    d = tp.line_deltas(tabs, "p2", "gstar", 11)
    assert np.allclose(d.loc[["A", "B", "C"]].to_numpy(), [0.6 + 0.11, 0.6 + 0.11, 0.6 + 0.11])
    tabs["gstar"] = _perline("Gstar", [0], [0.3, 0.2], lines=("A", "B"))
    with pytest.raises(ValueError, match="different lines"):
        tp.line_deltas(tabs, "p2", "gstar", 11)


def test_skew_decision_rule_and_disagreement_flag():
    rng = np.random.default_rng(0)
    sym = pd.Series(rng.normal(size=500))
    skewed = pd.Series(rng.exponential(size=500) ** 2)
    r = tp.skew_decision({"C4": sym, "C5": sym})
    assert r["bootstrap_method"] == "percentile" and r["disagrees"] == []
    r = tp.skew_decision({"C4": skewed, "C5": sym})
    assert r["bootstrap_method"] == "bca" and r["disagrees"] == ["C5"] and r["decided_on"] == "C4"
    r = tp.skew_decision({"C4": sym, "C5": skewed})
    assert r["bootstrap_method"] == "percentile" and r["disagrees"] == ["C5"]


def test_build_spec_passes_harness_checks():
    prov = {p: "0" * 64 for p in tp.PROVENANCE_FILES}
    entry = {"path": "results/m6/t2/p2_scores_test.parquet", "sha256": "a" * 64, "layout": "gene"}
    spec = tp.build_spec(entry, "percentile", prov)
    tp.check(spec)
    assert spec["fitted_on_labels"] is True and spec["refit"] == {"recipe": "w7_idea", "idea": "W7-05", "method": "P2"}
    assert set(et.REFIT_CODE) <= set(spec["provenance"]) and spec["failed_builds"] == {"test": None}
    tmpl = json.loads((PROJECT / "results/m6/p2_spec.template.json").read_text())
    assert set(spec) == set(tmpl) - {"_doc"} | {"notes"}
    with pytest.raises(et.PreflightError):
        tp.check({**spec, "refit": {**spec["refit"], "idea": "W7-02"}})


@pytest.mark.skipif(not tp.SPEC_PATH.is_file(), reason="results/m6/p2_spec.json not written yet")
def test_written_spec_static_checks():
    spec = json.loads(tp.SPEC_PATH.read_text())
    tp.check(spec)
    assert spec["spec_id"] == tp.SPEC_ID and spec["primary"] == tp.PRIMARY
    assert spec["scores"]["test"]["path"] == "results/m6/t2/p2_scores_test.parquet"
    assert spec["provenance"].get("results/m6/t2/skew.json") is not None
