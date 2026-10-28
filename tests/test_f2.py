"""Label-free tests for the F2 follow-up (A15 post hoc; scripts/f2_*.py): own tAP and bootstrap, score composition,
and the validation spec shape (evaluate_val.check_spec)."""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import f2_common as fc  # noqa: E402


def _ap(y_ordered) -> float:
    y = np.asarray(y_ordered, dtype=bool)
    hits = np.cumsum(y)
    return float(np.sum(hits[y] / (np.flatnonzero(y) + 1)) / y.sum())


def _tap_bruteforce(s, y) -> float:
    """Mean AP over all orderings that permute items within tied score blocks."""
    s, y = np.asarray(s, float), np.asarray(y, bool)
    vals = sorted(set(s), reverse=True)
    blocks = [list(np.flatnonzero(s == v)) for v in vals]
    aps = [_ap(y[list(itertools.chain(*p))]) for p in itertools.product(*(itertools.permutations(b) for b in blocks))]
    return float(np.mean(aps))


def test_f2_tap_equals_bruteforce_on_ties():
    rng = np.random.default_rng(0)
    for _ in range(30):
        n = int(rng.integers(3, 8))
        s = rng.integers(0, 3, n).astype(float)
        y = rng.random(n) < 0.4
        if not y.any():
            y[0] = True
        assert fc.tap(s, y) == pytest.approx(_tap_bruteforce(s, y), abs=1e-12)


def test_f2_tap_untied_is_plain_ap():
    s = np.array([0.9, 0.8, 0.7, 0.6, 0.5])
    y = np.array([1, 0, 1, 0, 1], bool)
    assert fc.tap(s, y) == pytest.approx((1 + 2 / 3 + 3 / 5) / 3)


def test_f2_tap_agrees_with_project_metrics():
    from nextgem import metrics

    rng = np.random.default_rng(1)
    for _ in range(20):
        n = 400
        s = np.round(rng.random(n), 1) * (rng.random(n) < 0.5)
        y = rng.random(n) < 0.15
        assert fc.tap(s, y) == pytest.approx(metrics.tap(s, y), abs=1e-12)


def test_f2_boot_lineage_deterministic_and_centred():
    lin = pd.Series(["a"] * 10 + ["b"] * 30, index=[f"l{i}" for i in range(40)])
    v = pd.Series(np.r_[np.ones(10), np.zeros(30)], index=lin.index)
    r1 = fc.boot_lineage(v, lin, 500, 11)
    assert np.array_equal(r1, fc.boot_lineage(v, lin, 500, 11))
    assert np.allclose(r1, 0.5)  # constant within lineage -> every replicate is the lineage-equal mean
    assert fc.lineage_equal_mean(v, lin) == 0.5


def test_f2_contrast_fold_rule_needs_five_lines():
    idx = [f"l{i}" for i in range(12)]
    lin = pd.Series(["A"] * 6 + ["B"] * 4 + ["C"] * 2, index=idx)
    base = pd.DataFrame({"estimable": True, "tap": 0.3}, index=idx)
    up = base.assign(tap=0.31)
    r = fc.contrast(up, base, lin, ["A", "B", "C"], B=200, seed=11)
    assert r["folds_positive"] == 1 and r["delta_lin"] == pytest.approx(0.01)


def test_f2_e_rows_rounding_and_absent_genes():
    gmap = pd.Series([10, 20], index=["g1", "g2"])
    R = pd.DataFrame({"g1": [0.1234567], "g2": [np.nan], "gx": [0.5]}, index=["s"])
    e = fc.e_rows(R, gmap, [10, 20, 30])
    assert e.loc["s", 10] == round(1 - 0.1234567, 6) and e.loc["s", 20] == 0.0 and e.loc["s", 30] == 0.0


def test_f2_method_names_are_fixed():
    assert fc.P1_METHOD == "M6:H1+Bloo/M2(f=1)" and fc.F2_METHOD == "F2:H1+artBloo/M2(f=1)"
    assert fc.F2_VARIANT == "W7:artefact_bloo"


@pytest.mark.skipif(not (fc.M10VAL / "look_spec_M10_F2.json").exists(), reason="spec not written yet")
def test_f2_spec_passes_check_spec():
    import evaluate_val as ev

    spec = json.loads((fc.M10VAL / "look_spec_M10_F2.json").read_text())
    ev.check_spec(spec, "M10_F2", "m10")
    assert spec["candidates"] == [fc.F2_METHOD] and spec["reference"] == "H1"
    assert spec["methods"] == ["H1", "H2", fc.P1_METHOD, fc.F2_METHOD]
    assert spec["workstream"] == spec["milestone"] == "M10"
    assert "A15" in spec["notes"] and "D18" in spec["notes"] and "F2 − P1" in spec["notes"]


@pytest.mark.skipif(not (fc.M10VAL / "f2_val_scores.parquet").exists(), reason="scores not built yet")
def test_f2_val_scores_shape():
    df = pd.read_parquet(fc.M10VAL / "f2_val_scores.parquet")
    assert list(df.columns) == ["ModelID", "entrez", "H1", "H2", fc.P1_METHOD, fc.F2_METHOD]
    assert df.ModelID.nunique() == 198 and not df.duplicated(["ModelID", "entrez"]).any()
    assert (df[fc.F2_METHOD] != df[fc.P1_METHOD]).any()
