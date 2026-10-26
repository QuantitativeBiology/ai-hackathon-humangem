"""Known-answer tests for the §6.8 binary endpoints and §6.1 binary calls (toy data only)."""

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, precision_score, recall_score

from nextgem import metrics as M

SEED = 11


@pytest.mark.parametrize(
    "tp,fp,tn,fn,expect",
    [
        # hand-computed
        (3, 1, 4, 2, dict(precision=0.75, recall=0.6, specificity=0.8, accuracy=0.7, f1=2 / 3, mcc=10 / np.sqrt(600))),
        (5, 0, 5, 0, dict(precision=1.0, recall=1.0, specificity=1.0, accuracy=1.0, f1=1.0, mcc=1.0)),
        (0, 5, 0, 5, dict(precision=0.0, recall=0.0, specificity=0.0, accuracy=0.0, f1=0.0, mcc=-1.0)),
        (10, 20, 60, 10, dict(precision=1 / 3, recall=0.5, specificity=0.75, accuracy=0.7, f1=0.4,
                              mcc=(600 - 200) / np.sqrt(30 * 20 * 80 * 70))),
    ],
)
def test_hand_computed(tp, fp, tn, fn, expect):
    got = M.binary_metrics(tp, fp, tn, fn)
    for k, v in expect.items():
        assert got[k] == pytest.approx(v, abs=1e-12), k


def test_zero_denominators_nan_and_sensitivity_row():
    # constant 'non-essential' call: TP + FP = 0 -> precision and MCC undefined
    got = M.binary_metrics(0, 0, 90, 10)
    assert np.isnan(got["mcc"]) and np.isnan(got["precision"])
    assert got["recall"] == 0 and got["f1"] == 0 and got["specificity"] == 1 and got["accuracy"] == 0.9
    zero = M.binary_metrics(0, 0, 90, 10, zero_division=0.0)
    assert zero["mcc"] == 0.0 and zero["precision"] == 0.0 and zero["recall"] == 0.0
    # no true essentials: recall undefined; F1 undefined only if also no predicted essentials
    g = M.binary_metrics(0, 3, 7, 0)
    assert np.isnan(g["recall"]) and np.isnan(g["mcc"]) and g["f1"] == 0.0
    g0 = M.binary_metrics(0, 0, 10, 0)
    assert np.isnan(g0["f1"]) and np.isnan(g0["precision"]) and np.isnan(g0["recall"])
    # vectorised
    v = M.binary_metrics(np.array([3, 0]), np.array([1, 0]), np.array([4, 9]), np.array([2, 1]))
    assert v["mcc"][0] == pytest.approx(10 / np.sqrt(600)) and np.isnan(v["mcc"][1])


def test_against_sklearn_where_defined():
    rng = np.random.default_rng(SEED)
    n_checked = 0
    for _ in range(300):
        n = int(rng.integers(5, 300))
        truth = rng.random(n) < rng.uniform(0.02, 0.6)
        calls = rng.random(n) < rng.uniform(0.02, 0.6)
        tp, fp, tn, fn = M.confusion_counts(calls, truth.astype(float))
        if min(tp + fp, tp + fn, tn + fp, tn + fn) == 0:
            continue
        m = M.binary_metrics(tp, fp, tn, fn)
        assert m["mcc"] == pytest.approx(matthews_corrcoef(truth, calls), abs=1e-12)
        assert m["f1"] == pytest.approx(f1_score(truth, calls), abs=1e-12)
        assert m["precision"] == pytest.approx(precision_score(truth, calls), abs=1e-12)
        assert m["recall"] == pytest.approx(recall_score(truth, calls), abs=1e-12)
        assert m["specificity"] == pytest.approx(recall_score(truth, calls, pos_label=False), abs=1e-12)
        assert m["accuracy"] == pytest.approx(accuracy_score(truth, calls), abs=1e-12)
        n_checked += 1
    assert n_checked > 250


def test_constant_call_negative_control():
    """§12: a constant-call method gives MCC 0 - realised as NaN (zero denominator, primary) -> 0 (sensitivity row)."""
    rng = np.random.default_rng(SEED)
    truth = (rng.random(500) < 0.1).astype(float)
    for const in (False, True):
        tp, fp, tn, fn = M.confusion_counts(np.full(500, const), truth)
        assert np.isnan(M.binary_metrics(tp, fp, tn, fn)["mcc"])
        assert M.binary_metrics(tp, fp, tn, fn, zero_division=0.0)["mcc"] == 0.0
        assert matthews_corrcoef(truth, np.full(500, const)) == 0.0  # sklearn convention = our sensitivity row
    # a random caller independent of truth has MCC ~ 0
    calls = rng.random((200, 500)) < 0.2
    mccs = [M.binary_metrics(*M.confusion_counts(c, truth))["mcc"] for c in calls]
    assert abs(np.mean(mccs)) < 0.01


def test_confusion_excludes_nan_truth_and_validates():
    truth = np.array([1, 0, np.nan, 1, 0])
    calls = np.array([True, True, True, False, False])
    assert M.confusion_counts(calls, truth) == (1, 1, 1, 1)
    with pytest.raises(ValueError):
        M.confusion_counts(np.array([1.0, np.nan]), np.array([1.0, 0.0]))
    with pytest.raises(M.PairingError):
        M.confusion_counts(calls[:4], truth)


def toy_frames(seed=SEED, n_lines=12, n_genes=80):
    rng = np.random.default_rng(seed)
    idx = [f"ACH-{i:06d}" for i in range(n_lines)]
    cols = [f"G{j} ({j})" for j in range(n_genes)]
    effect = pd.DataFrame(rng.normal(-0.3, 0.4, (n_lines, n_genes)), idx, cols)
    effect.iloc[0, :5] = np.nan
    calls = pd.DataFrame(rng.random((n_lines, n_genes)) < 0.15, idx, cols)
    lin = pd.Series(np.repeat(["A", "B", "C"], n_lines // 3), index=idx)
    return calls, effect, lin


def test_per_line_matches_scalar_path_and_sklearn():
    calls, effect, _ = toy_frames()
    truth = M.labels_from_effect(effect)  # L-ge06: effect < -0.6
    t = M.binary_per_line(calls, truth)
    for line in t.index:
        m = truth.loc[line].notna()
        y, c = truth.loc[line][m].astype(bool), calls.loc[line][m]
        assert tuple(t.loc[line, ["tp", "fp", "tn", "fn"]]) == M.confusion_counts(c.to_numpy(), y.to_numpy(float))
        if not np.isnan(t.loc[line, "mcc"]):
            assert t.loc[line, "mcc"] == pytest.approx(matthews_corrcoef(y, c), abs=1e-12)
    assert t.loc[calls.index[0], ["tp", "fp", "tn", "fn"]].sum() == 75  # 5 NaN-truth pairs excluded
    with pytest.raises(M.PairingError):
        M.binary_per_line(calls.drop(columns=calls.columns[0]), truth)


def test_labels_from_effect_strict_and_nan():
    eff = pd.DataFrame({"a": [-0.6, -0.61, np.nan, 0.0]})
    lab = M.labels_from_effect(eff, -0.6)
    assert lab["a"].tolist()[:2] == [0.0, 1.0] and np.isnan(lab["a"].iloc[2]) and lab["a"].iloc[3] == 0.0


def test_threshold_grid():
    assert len(M.THRESHOLD_GRID) == 21
    assert M.THRESHOLD_GRID[0] == -1.0 and M.THRESHOLD_GRID[-1] == 0.0
    np.testing.assert_allclose(np.diff(M.THRESHOLD_GRID), 0.05)
    assert -0.6 in M.THRESHOLD_GRID
    calls, effect, _ = toy_frames()
    g = M.binary_threshold_grid(calls, effect)
    assert len(g) == 21 * len(calls)
    at06 = g[g["threshold"] == -0.6].drop(columns="threshold")
    pd.testing.assert_frame_equal(at06, M.binary_per_line(calls, M.labels_from_effect(effect, -0.6)))
    # positives are monotone in the threshold
    pos = g.assign(P=g.tp + g.fn).groupby("threshold")["P"].sum()
    assert pos.is_monotonic_increasing


def test_summaries():
    calls, effect, lin = toy_frames()
    t = M.binary_per_line(calls, M.labels_from_effect(effect))
    t.loc[t.index[0], "mcc"] = np.nan
    s = M.summarise_binary(t, lin).set_index(["metric", "nan_policy"])
    ex, zr = s.loc[("mcc", "exclude")], s.loc[("mcc", "zero")]
    assert ex["n_nan"] == 1 and ex["n_lines_used"] == len(t) - 1 and zr["n_lines_used"] == len(t)
    assert ex["line_mean"] == pytest.approx(t["mcc"].mean())
    assert zr["line_mean"] == pytest.approx(t["mcc"].fillna(0).mean())
    assert ex["lineage_equal_mean"] == pytest.approx(t["mcc"].groupby(lin).mean().mean())
    assert ex["median"] == pytest.approx(t["mcc"].median())


def test_fba_and_rate_calls():
    r = np.array([0.0, 0.009999, 0.01, 0.5, 1.0])
    assert M.fba_calls(r).tolist() == [True, True, False, False, False]
    fit = np.arange(100.0)
    thr = M.rate_threshold(fit, 0.1)
    assert thr == 90.0 and M.calls_at_threshold(fit, thr).sum() == 10
    tied = np.r_[np.zeros(90), np.ones(10) * 5]
    assert M.calls_at_threshold(tied, M.rate_threshold(tied, 0.05)).sum() == 10  # ties: rate >= target
    with pytest.raises(ValueError):
        M.rate_threshold(fit, 0.0)
