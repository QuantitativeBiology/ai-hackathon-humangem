"""Known-answer tests for tAP, its null, units and inference (DESIGN.md v3 §5.2, §6, §11.4, §12 'tAP and metrics').

Toy / synthetic data only; no label file is read.
"""

from fractions import Fraction
from itertools import permutations, product
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats
from sklearn.metrics import average_precision_score

from nextgem import metrics as M

TOL = 1e-12
SEED = 11  # SEEDS[0] (DESIGN §7)


# ---------------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------------
def brute_force_tap(scores, labels) -> Fraction:
    """Exact E[AP] by enumerating every ordering within every tied block (Fractions)."""
    scores, labels = list(map(float, scores)), list(map(int, labels))
    levels = sorted(set(scores), reverse=True)
    blocks = [[labels[i] for i in range(len(scores)) if scores[i] == v] for v in levels]
    n_pos = sum(labels)
    total, count = Fraction(0), 0
    for combo in product(*[list(permutations(b)) for b in blocks]):
        order = [y for blk in combo for y in blk]
        hits, ap = 0, Fraction(0)
        for k, y in enumerate(order, start=1):
            if y:
                hits += 1
                ap += Fraction(hits, k)
        total += ap / n_pos
        count += 1
    return total / count


def step_ap_untied(scores, labels) -> float:
    """Plain AP of an untied ranking (independent of nextgem)."""
    order = np.argsort(-np.asarray(scores), kind="stable")
    y = np.asarray(labels, float)[order]
    prec = np.cumsum(y) / np.arange(1, y.size + 1)
    return float((prec * y).sum() / y.sum())


def toy_case(rng, n, levels=3, p=0.4):
    s = rng.integers(0, levels, n).astype(float)
    y = (rng.random(n) < p).astype(float)
    if y.sum() == 0:
        y[rng.integers(n)] = 1.0
    return s, y


def near_binary_case(rng, n=2866, top=150, partial=40, n_pos=320, enrich=0.5):
    """Generic-FBA-shaped scorer: a top block of e = 1, a few partial values, everything else tied at 0."""
    s = np.zeros(n)
    idx = rng.permutation(n)
    s[idx[:top]] = 1.0
    s[idx[top : top + partial]] = rng.choice([0.25, 0.5, 0.75], partial)
    y = np.zeros(n)
    n_top = int(enrich * top)
    y[idx[:n_top]] = 1.0
    rest = idx[top:]
    y[rng.choice(rest, n_pos - n_top, replace=False)] = 1.0
    return s, y


# ---------------------------------------------------------------------------------------------------
# §12: brute-force enumeration on toy rankings with ties
# ---------------------------------------------------------------------------------------------------
def test_hand_example():
    # block {pos, neg} then {pos}: orders (+ - | +) AP = 5/6 and (- + | +) AP = 7/12 -> 17/24
    assert M.tap([1, 1, 0], [1, 0, 1]) == pytest.approx(17 / 24, abs=TOL)
    assert brute_force_tap([1, 1, 0], [1, 0, 1]) == Fraction(17, 24)


def test_brute_force_enumeration_random_toys():
    rng = np.random.default_rng(SEED)
    n_cases = 0
    for _ in range(300):
        n = int(rng.integers(1, 8))
        s, y = toy_case(rng, n, levels=int(rng.integers(1, 4)))
        exact = float(brute_force_tap(s, y))
        assert M.tap(s, y) == pytest.approx(exact, abs=TOL)
        assert M.tap_hypergeometric(s, y) == pytest.approx(exact, abs=TOL)
        assert M.tap_weighted(s, y, np.ones(n)) == pytest.approx(exact, abs=TOL)
        n_cases += 1
    assert n_cases == 300


def test_brute_force_single_big_block():
    s, y = np.zeros(8), np.array([1, 0, 1, 0, 0, 1, 0, 0], float)
    exact = float(brute_force_tap(s, y))
    assert M.tap(s, y) == pytest.approx(exact, abs=TOL)
    assert M.expected_ap_random(8, 3) == pytest.approx(exact, abs=TOL)


# ---------------------------------------------------------------------------------------------------
# §12: O(n_g) form equals the hypergeometric form (incl. blocks whose binomials overflow)
# ---------------------------------------------------------------------------------------------------
def test_on_form_equals_hypergeometric_random():
    rng = np.random.default_rng(SEED + 1)
    for _ in range(200):
        n = int(rng.integers(2, 200))
        s, y = toy_case(rng, n, levels=int(rng.integers(1, 12)), p=rng.uniform(0.02, 0.9))
        assert M.tap(s, y) == pytest.approx(M.tap_hypergeometric(s, y), abs=TOL)


def test_on_form_equals_hypergeometric_large_blocks():
    rng = np.random.default_rng(SEED + 2)
    s, y = near_binary_case(rng)  # bottom block ~2,676 items with ~245 positives: C(n, p) overflows float
    assert M.tap(s, y) == pytest.approx(M.tap_hypergeometric(s, y), abs=TOL)
    s2 = np.zeros(2650)
    y2 = np.zeros(2650)
    y2[:250] = 1
    assert M.tap(s2, y2) == pytest.approx(M.tap_hypergeometric(s2, y2), abs=TOL)


def test_vectorised_and_weighted_forms_equal_tap():
    rng = np.random.default_rng(SEED + 3)
    S = rng.integers(0, 5, (60, 9)).astype(float)
    Y = (rng.random((60, 9)) < 0.3).astype(float)
    Y[rng.random((60, 9)) < 0.1] = np.nan
    t, n, p = M.tap_columns(S, Y)
    for j in range(9):
        m = ~np.isnan(Y[:, j])
        assert t[j] == pytest.approx(M.tap(S[m, j], Y[m, j]), abs=TOL)
        assert n[j] == m.sum() and p[j] == Y[m, j].sum()
    # integer weights == duplicated items
    for _ in range(50):
        s, y = toy_case(rng, 25, levels=4)
        w = rng.integers(0, 4, 25)
        if (w * y).sum() == 0:
            w[np.flatnonzero(y)[0]] = 1
        dup_s, dup_y = np.repeat(s, w), np.repeat(y, w)
        assert M.tap_weighted(s, y, w) == pytest.approx(M.tap(dup_s, dup_y), abs=TOL)


def test_stats_review_reference_value():
    # notes/design_review_r2/statistics.md: E[AP_random](192, 5) = 0.0507
    assert M.expected_ap_random(192, 5) == pytest.approx(0.0507, abs=5e-5)


# ---------------------------------------------------------------------------------------------------
# §12: constant scorer -> E[AP_random] exactly
# ---------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("n,p", [(1, 1), (2, 1), (2, 2), (10, 3), (192, 5), (2866, 320), (2450, 37), (5000, 4999)])
def test_constant_scorer_equals_null(n, p):
    y = np.zeros(n)
    y[:p] = 1
    for c in (0.0, 0.37, 1.0):
        assert M.tap(np.full(n, c), y) == pytest.approx(M.expected_ap_random(n, p), abs=TOL)
    assert M.tap_hypergeometric(np.zeros(n), y) == pytest.approx(M.expected_ap_random(n, p), abs=1e-10)


def test_null_vectorised_and_no_positive():
    n = np.array([10, 20, 30])
    p = np.array([1, 5, 30])
    e = M.expected_ap_random(n, p)
    assert e[2] == pytest.approx(1.0)
    assert e[0] == pytest.approx(M.harmonic(10) / 10)  # one positive: mean of 1/r
    assert np.isnan(M.expected_ap_random(10, 0))
    assert np.isnan(M.tap(np.arange(5.0), np.zeros(5)))


def test_line_constant_scorer_gene_tap_equals_null():
    """§6.2: for a line-constant scorer tAP_gene = E[AP_random] exactly."""
    rng = np.random.default_rng(SEED + 4)
    lines = [f"L{i}" for i in range(40)]
    genes = [f"g{j}" for j in range(30)]
    labels = pd.DataFrame((rng.random((40, 30)) < 0.4).astype(float), index=lines, columns=genes)
    scores = pd.DataFrame(np.tile(rng.random(30), (40, 1)), index=lines, columns=genes)
    g = M.per_gene(scores, labels)
    est = g[g["estimable"]]
    assert len(est) > 20
    np.testing.assert_allclose(est["tap"], est["e_random"], atol=TOL)
    assert M.m_delta_ap_gene(g) == pytest.approx(0.0, abs=TOL)


# ---------------------------------------------------------------------------------------------------
# §12: Monte-Carlo average of score + 1e-9 noise over 1,000 seeds equals tAP within MC error
# ---------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("case", ["toy", "near_binary"])
def test_noise_average_matches_tap(case):
    rng = np.random.default_rng(SEED + 5)
    if case == "toy":
        s, y = toy_case(rng, 40, levels=3)
    else:
        s, y = near_binary_case(rng, n=600, top=60, partial=15, n_pos=60)
    target = M.tap(s, y)
    aps = []
    for child in np.random.SeedSequence(SEED).spawn(1000):
        noisy = s + 1e-9 * np.random.default_rng(child).random(s.size)
        aps.append(step_ap_untied(noisy, y))
    aps = np.asarray(aps)
    se = aps.std(ddof=1) / np.sqrt(aps.size)
    assert abs(aps.mean() - target) < 4 * se + 1e-12
    # the untied AP equals sklearn's on the jittered ranking (independent implementation)
    noisy = s + 1e-9 * np.random.default_rng(0).random(s.size)
    assert step_ap_untied(noisy, y) == pytest.approx(average_precision_score(y, noisy), abs=1e-12)


# ---------------------------------------------------------------------------------------------------
# §12: perfect and reversed rankings
# ---------------------------------------------------------------------------------------------------
def test_perfect_and_reversed():
    n, p = 50, 7
    y = np.zeros(n)
    y[:p] = 1
    perfect = np.linspace(1, 0, n)
    assert M.tap(perfect, y) == pytest.approx(1.0, abs=TOL)
    assert M.tap(np.r_[np.ones(p), np.zeros(n - p)], y) == pytest.approx(1.0, abs=TOL)  # tied positives on top
    expected_rev = np.mean([j / (n - p + j) for j in range(1, p + 1)])
    assert M.tap(-perfect, y) == pytest.approx(expected_rev, abs=TOL)
    assert M.tap(-perfect, y) < M.expected_ap_random(n, p) < M.tap(perfect, y)
    # untied rankings: tAP == sklearn step AP
    rng = np.random.default_rng(SEED + 6)
    s = rng.random(n)
    assert M.tap(s, y) == pytest.approx(average_precision_score(y, s), abs=TOL)
    assert M.ap_step_sklearn(s, y) == pytest.approx(M.tap(s, y), abs=TOL)


def test_sklearn_sensitivity_row_penalises_ties():
    rng = np.random.default_rng(SEED + 7)
    s, y = near_binary_case(rng)
    assert M.ap_step_sklearn(s, y) < M.tap(s, y)
    assert np.isnan(M.ap_step_sklearn(s, np.zeros_like(y)))


def test_auroc_and_top_k():
    y = np.array([1, 0, 1, 0, 0], float)
    assert M.auroc([5, 4, 3, 2, 1], y) == pytest.approx(5 / 6)
    assert np.isnan(M.auroc([1, 2, 3, 4, 5], np.zeros(5)))
    assert M.expected_tp_at_k([5, 4, 3, 2, 1], y) == 1.0  # top 2 = {pos, neg}
    # top block of 4 tied items with 2 positives, k = 2 -> expected 1 positive
    assert M.expected_tp_at_k([1, 1, 1, 1, 0], [1, 1, 0, 0, 0]) == pytest.approx(1.0)
    assert M.expected_tp_at_k([1, 1, 1, 1, 0], [1, 1, 0, 0, 0], k=4) == pytest.approx(2.0)


def test_input_validation():
    with pytest.raises(ValueError):
        M.tap([1, 2], [1, np.nan])
    with pytest.raises(ValueError):
        M.tap([1, np.nan], [1, 0])
    with pytest.raises(ValueError):
        M.tap([1, 2], [1, 2])


# ---------------------------------------------------------------------------------------------------
# §5.2: the >= 5-positive NA rule (lines), >= 5/5 (genes), grid cell rule
# ---------------------------------------------------------------------------------------------------
def make_frames(Y, S=None, seed=SEED):
    rng = np.random.default_rng(seed)
    Y = np.asarray(Y, float)
    S = rng.random(Y.shape) if S is None else S
    idx = [f"ACH-{i:06d}" for i in range(Y.shape[0])]
    cols = [f"G{j} ({j})" for j in range(Y.shape[1])]
    return pd.DataFrame(S, idx, cols), pd.DataFrame(Y, idx, cols)


def test_line_na_rule():
    Y = np.zeros((3, 40))
    Y[0, :4] = 1  # 4 positives -> NA
    Y[1, :5] = 1  # 5 positives -> estimable
    Y[2, :6] = 1
    Y[2, 5] = np.nan  # NaN pair excluded: 5 positives remain -> estimable
    S, L = make_frames(Y)
    t = M.per_line(S, L)
    assert list(t["estimable"]) == [False, True, True]
    assert np.isnan(t["tap"].iloc[0]) and np.isnan(t["e_random"].iloc[0])
    assert t["n_pos"].tolist() == [4, 5, 5] and t["n"].tolist() == [40, 40, 39]
    assert t["lift"].iloc[1] == pytest.approx(t["tap"].iloc[1] / M.expected_ap_random(40, 5))
    assert t["delta_null"].iloc[1] == pytest.approx(t["tap"].iloc[1] - M.expected_ap_random(40, 5))


def test_gene_na_rule():
    Y = np.zeros((12, 3))
    Y[:5, 0] = 1  # 5 pos, 7 neg -> estimable
    Y[:8, 1] = 1  # 8 pos, 4 neg -> NA
    Y[:4, 2] = 1  # 4 pos -> NA
    S, L = make_frames(Y)
    g = M.per_gene(S, L)
    assert list(g["estimable"]) == [True, False, False]
    assert g["n_neg"].tolist() == [7, 4, 8]


def test_cell_rule():
    lt = pd.DataFrame({"estimable": [True, False, False, True]})
    assert not M.cell_is_na(lt)
    assert M.cell_is_na(pd.DataFrame({"estimable": [True, False, False, False]}))
    gt = pd.DataFrame({"estimable": [True] * 29 + [False] * 10})
    assert M.cell_is_na(lt, gt)
    assert not M.cell_is_na(lt, pd.DataFrame({"estimable": [True] * 30}))


# ---------------------------------------------------------------------------------------------------
# §6.5: set-mismatch / pairing assertions
# ---------------------------------------------------------------------------------------------------
def test_set_mismatch_assertions():
    rng = np.random.default_rng(SEED)
    S, L = make_frames((rng.random((10, 20)) < 0.4).astype(float))
    with pytest.raises(M.PairingError):
        M.per_line(S.drop(columns=S.columns[0]), L)
    with pytest.raises(M.PairingError):
        M.per_line(S.drop(index=S.index[0]), L)
    S_nan = S.copy()
    S_nan.iloc[0, 0] = np.nan
    with pytest.raises(M.PairingError):
        M.per_line(S_nan, L)
    # a NaN score on a NaN-label pair is fine (pair is excluded for every method)
    L2 = L.copy()
    L2.iloc[0, 0] = np.nan
    M.per_line(S_nan, L2)
    # reordered rows/columns are aligned, not mismatched
    pd.testing.assert_frame_equal(M.per_line(S.iloc[::-1, ::-1], L), M.per_line(S, L))
    # paired deltas need identical pairs
    a = M.per_line(S, L)
    L3 = L.copy()
    L3.iloc[1, 1] = np.nan
    b = M.per_line(S, L3)
    with pytest.raises(M.PairingError):
        M.paired_line_delta(a, b)
    with pytest.raises(M.PairingError):
        M.assert_same_pairs(S, L.rename(columns={L.columns[0]: "X (0)"}))
    with pytest.raises(M.PairingError):
        M.lineage_equal_mean(pd.Series([1.0, 2.0], index=["a", "b"]), pd.Series(["k"], index=["a"]))


# ---------------------------------------------------------------------------------------------------
# §11.4 failed-line no-information rule (per line: all tied; per gene: per-gene median, as a tie)
# ---------------------------------------------------------------------------------------------------
def test_failed_line_rules():
    rng = np.random.default_rng(SEED + 8)
    S, L = make_frames((rng.random((12, 50)) < 0.3).astype(float), S=rng.integers(0, 4, (12, 50)) / 4)
    failed = [S.index[2], S.index[7]]
    sl = M.no_information_scores(S, failed, "line")
    t = M.per_line(sl, L)
    for f in failed:
        n, p = t.loc[f, "n"], t.loc[f, "n_pos"]
        if p >= 5:
            assert t.loc[f, "tap"] == pytest.approx(M.expected_ap_random(n, p), abs=TOL)
    sg = M.no_information_scores(S, failed, "gene")
    med = S.drop(index=failed).median()
    for f in failed:
        np.testing.assert_allclose(sg.loc[f].to_numpy(), med.to_numpy())
    pd.testing.assert_frame_equal(sg.drop(index=failed), S.drop(index=failed))
    # the two failed lines are tied with each other for every gene
    assert (sg.loc[failed[0]] == sg.loc[failed[1]]).all()
    with pytest.raises(ValueError):
        M.no_information_scores(S, failed, "pooled")


# ---------------------------------------------------------------------------------------------------
# summaries: lineage-equal vs line-weighted, per-lineage, pooled
# ---------------------------------------------------------------------------------------------------
def test_lineage_equal_mean():
    v = pd.Series([1.0, 1.0, 1.0, 4.0, np.nan], index=list("abcde"))
    lin = pd.Series(["x", "x", "x", "y", "y"], index=list("abcde"))
    assert M.lineage_equal_mean(v, lin) == pytest.approx(2.5)
    assert M.line_weighted_mean(v) == pytest.approx(7 / 4)
    w = M.lineage_weights(lin.loc[list("abcd")])
    assert w.sum() == pytest.approx(1.0)


def test_per_lineage_and_pooled():
    rng = np.random.default_rng(SEED + 9)
    S, L = make_frames((rng.random((20, 60)) < 0.3).astype(float))
    lin = pd.Series(np.repeat(["A", "B", "C", "D"], 5), index=S.index)
    t = M.per_line(S, L)
    ref = M.per_line(pd.DataFrame(0.0, S.index, S.columns), L)
    pl = M.per_lineage_summary(t, lin, reference=ref)
    assert pl.loc["A", "mean_tap"] == pytest.approx(t.loc[lin == "A", "tap"].mean())
    assert pl.loc["A", "lift"] == pytest.approx(pl.loc["A", "mean_tap"] / pl.loc["A", "mean_e_random"])
    d = M.paired_line_delta(t, ref)
    assert pl.loc["B", "delta_vs_ref"] == pytest.approx(d[lin.loc[d.index] == "B"].mean())
    po = M.pooled(S, L)
    m = L.notna().to_numpy()
    assert po["tap"] == pytest.approx(M.tap(S.to_numpy()[m], L.to_numpy()[m]))
    assert po["lift"] == pytest.approx(po["tap"] / M.expected_ap_random(po["n"], po["n_pos"]))
    assert M.pooled(S, L, lines=S.index[:5])["n"] == 5 * 60


# ---------------------------------------------------------------------------------------------------
# inference: bootstrap, sign-flip, lineage t-interval, transfer, claims, gatekeeping, SESOIs
# ---------------------------------------------------------------------------------------------------
def lineage_series(sizes, seed=SEED):
    idx = [f"L{i}" for i in range(sum(sizes))]
    return pd.Series(np.repeat([f"K{k}" for k in range(len(sizes))], sizes), index=idx)


def test_bootstrap_counts_and_determinism():
    lin = lineage_series([5, 10, 3])
    C = M.stratified_bootstrap_counts(lin, 1000, SEED)
    for k, n in (("K0", 5), ("K1", 10), ("K2", 3)):
        assert (C[:, (lin == k).to_numpy()].sum(1) == n).all()
    np.testing.assert_array_equal(C, M.stratified_bootstrap_counts(lin, 1000, SEED))
    with pytest.raises(ValueError):
        M.stratified_bootstrap_counts(lin, 10, None)


def test_bootstrap_mean():
    lin = lineage_series([30, 10, 20])
    rng = np.random.default_rng(SEED)
    v = pd.Series(rng.normal(0.01, 0.02, len(lin)), index=lin.index)
    r = M.bootstrap_mean(v, lin, B=4000, seed=SEED)
    assert r.estimate == pytest.approx(M.lineage_equal_mean(v, lin))
    assert r.ci_low < r.estimate < r.ci_high
    assert r.replicates.mean() == pytest.approx(r.estimate, abs=3 * r.replicates.std() / 60)
    # stratified bootstrap SE ~ analytic SE of the lineage-equal mean (with (n-1)/n shrinkage)
    se = np.sqrt(sum(v[lin == k].var(ddof=0) / (lin == k).sum() for k in lin.unique())) / 3
    assert r.replicates.std() == pytest.approx(se, rel=0.08)
    rl = M.bootstrap_mean(v, lin, B=2000, seed=SEED, weighting="line")
    assert rl.estimate == pytest.approx(v.mean())
    rb = M.bootstrap_mean(v, lin, B=4000, seed=SEED, method="bca")
    assert rb.ci_low == pytest.approx(r.ci_low, abs=0.003) and rb.ci_high == pytest.approx(r.ci_high, abs=0.003)
    assert not M.skew_triggers_bca(v)
    assert M.skew_triggers_bca(np.r_[np.zeros(50), 1.0, 2.0, 5.0])


def test_sign_flip():
    lin = lineage_series([10, 10, 10])
    big = pd.Series(np.full(30, 0.05), index=lin.index)
    obs, p = M.sign_flip_test(big, lin, n_flip=999, seed=SEED)
    assert obs == pytest.approx(0.05) and p == pytest.approx(1 / 1000)
    rng = np.random.default_rng(SEED)
    ps = []
    for i in range(200):
        d = pd.Series(rng.normal(0, 1, 30), index=lin.index)
        ps.append(M.sign_flip_test(d, lin, n_flip=200, seed=i)[1])
    ps = np.asarray(ps)
    assert 0.0 < np.mean(ps <= 0.05) < 0.11  # approximately uniform under the null


def test_lineage_t_interval_and_exact_sign_flip():
    m = np.array([0.01, 0.02, -0.005, 0.015, 0.03, 0.0, 0.012])
    mean, lo, hi = M.lineage_t_interval(m)
    ref = stats.t.interval(0.95, 6, loc=m.mean(), scale=stats.sem(m))
    assert (lo, hi) == pytest.approx(ref)
    assert M.exact_sign_flip_lineage(np.full(7, 0.01)) == pytest.approx(1 / 128)
    assert M.exact_sign_flip_lineage(np.array([1, 1, 1, 1, 1, 1, -1.0])) == pytest.approx(8 / 128)  # <= 1 flip
    assert M.exact_sign_flip_lineage(np.array([0.05, 0.04, 0.03, 0.02, 0.01, 0.005, -0.001])) == pytest.approx(2 / 128)


def test_transfer_rule():
    lin = lineage_series([3] * 7)
    vals = np.repeat([0.02, 0.03, 0.01, 0.02, 0.025, -0.001, -0.002], 3)
    tr = M.transfer_summary(pd.Series(vals, index=lin.index), lin)
    assert tr["K"] == 7 and tr["k_required"] == 5 and tr["k_positive"] == 5
    assert tr["generalises"] == (tr["t_low"] > 0)
    vals2 = np.repeat([0.02, 0.03, 0.01, 0.02, -0.025, -0.001, -0.002], 3)
    assert not M.transfer_summary(pd.Series(vals2, index=lin.index), lin)["generalises"]
    # ceil(0.7 K) without float artefacts
    for K, req in ((5, 4), (7, 5), (10, 7), (12, 9)):
        lk = lineage_series([2] * K)
        assert M.transfer_summary(pd.Series(1.0, index=lk.index), lk)["k_required"] == req


def test_claim_tiers_and_sesoi():
    ok = dict(gate_status="confirmed")
    assert M.claim_tier(False, 0.05, 0.03, 0.01, gate_status="not confirmed") == "no detectable improvement"
    assert M.claim_tier(True, 0.005, 0.001, 0.01, **ok) == "detectable improvement"
    assert M.claim_tier(True, 0.012, 0.004, 0.01, **ok) == "improves"
    assert M.claim_tier(True, 0.03, 0.011, 0.01, **ok) == "meaningfully improves"
    assert M.claim_tier(True, 0.03, 0.011, None, **ok) == "detectable improvement"  # C3: no SESOI
    assert M.sesoi_c1b([0.03, 0.05], [0.10, 0.10]) == pytest.approx(0.004)


def test_claim_tier_respects_gatekeeping_state():
    """Review should-fix 2: after a failed gate a passing test gets no claim tier (old API ignored the gate)."""
    with pytest.raises(TypeError):
        M.claim_tier(True, 0.03, 0.011, 0.01)  # gate state is now required
    assert M.claim_tier(True, 0.03, 0.011, 0.01, gate_status="estimate only") == "estimate only (earlier gate failed)"
    assert M.claim_tier(True, 0.03, 0.011, 0.01, gate_status="pending") == "pending"
    with pytest.raises(ValueError):
        M.claim_tier(True, 0.03, 0.011, 0.01, gate_status="not confirmed")  # inconsistent
    with pytest.raises(ValueError):
        M.claim_tier(False, 0.03, 0.011, 0.01, gate_status="confirmed")
    with pytest.raises(ValueError):
        M.claim_tier(True, 0.03, 0.011, 0.01, gate_status="maybe")


def test_fixed_sequence_gatekeeping():
    t = M.fixed_sequence([("C1", True), ("C1b", True), ("C3", False), ("C4", True), ("C5", True)])
    assert t["status"].tolist() == ["confirmed", "confirmed", "not confirmed", "estimate only", "estimate only"]
    assert t["confirmatory"].tolist() == [True, True, True, False, False]
    t1 = M.fixed_sequence([("C1", True), ("C1b", True), ("C3", True), ("C4", None), ("C5", None)])
    assert t1["status"].tolist()[3:] == ["pending", "estimate only"]
    assert tuple(t["contrast"]) == M.GATE_ORDER


def test_fixed_sequence_enforces_gate_order():
    """Review F1: the frozen order C1 -> C1b -> C3 -> C4 -> C5 is enforced (old code accepted any order)."""
    with pytest.raises(ValueError):
        M.fixed_sequence([("C5", True), ("C1", True)])
    with pytest.raises(ValueError):
        M.fixed_sequence([("C1", True), ("C3", True)])  # C1b skipped
    with pytest.raises(ValueError):
        M.fixed_sequence([("C1", True), ("C6", True)])  # unknown contrast
    with pytest.raises(ValueError):
        M.fixed_sequence({"C1": True, "C2": True})
    with pytest.raises(ValueError):
        M.fixed_sequence({"C1": True, "C3": True})  # mapping must also be a prefix
    # a mapping in any key order is put into the frozen order; testing stops at the first failure
    t = M.fixed_sequence({"C5": True, "C3": False, "C1": True, "C4": True, "C1b": True})
    assert t["contrast"].tolist() == list(M.GATE_ORDER)
    assert t["status"].tolist() == ["confirmed", "confirmed", "not confirmed", "estimate only", "estimate only"]
    assert M.fixed_sequence([("C1", False)])["status"].tolist() == ["not confirmed"]


def test_confirmatory_contrast_end_to_end():
    lin = lineage_series([40, 30, 20, 15, 15, 10, 10])
    rng = np.random.default_rng(SEED)
    d = pd.Series(rng.normal(0.03, 0.02, len(lin)), index=lin.index)
    r = M.confirmatory_contrast(d, lin, contrast="C1", sesoi=M.SESOI, seed=SEED, B=2000, n_flip=2000)
    assert r["passed"] and r["p_sign_flip"] == pytest.approx(1 / 2001)
    assert "claim" not in r  # wording only through gatekept_claims
    assert r["transfer"]["generalises"]
    r2 = M.confirmatory_contrast(d, lin, contrast="C1", sesoi=M.SESOI, seed=SEED, B=2000, n_flip=2000)
    assert r2["ci_low"] == r["ci_low"]  # deterministic given the seed
    claims = M.gatekept_claims({"C1": r})
    assert claims["claim"].tolist() == ["meaningfully improves"]


def test_confirmatory_contrast_requires_contrast_specific_sesoi():
    """Review F2: no SESOI default; C1/C4/C5 use the frozen 0.01, C1b must pass its own lift-scaled SESOI."""
    lin = lineage_series([10, 10, 10])
    d = pd.Series(np.random.default_rng(SEED).normal(0.02, 0.01, 30), index=lin.index)
    with pytest.raises(TypeError):
        M.confirmatory_contrast(d, lin, seed=SEED, B=100, n_flip=100)  # old call: silent default 0.01
    with pytest.raises(TypeError):
        M.confirmatory_contrast(d, lin, contrast="C1b", seed=SEED, B=100, n_flip=100)
    with pytest.raises(ValueError):
        M.confirmatory_contrast(d, lin, contrast="C1b", sesoi=M.SESOI, seed=SEED, B=100, n_flip=100)
    with pytest.raises(ValueError):
        M.confirmatory_contrast(d, lin, contrast="C1", sesoi=0.004, seed=SEED, B=100, n_flip=100)
    with pytest.raises(ValueError):
        M.confirmatory_contrast(d, lin, contrast="C3", sesoi=M.SESOI, seed=SEED, B=100, n_flip=100)
    s_c1b = M.sesoi_c1b([0.03, 0.05], [0.10, 0.10])
    r = M.confirmatory_contrast(d, lin, contrast="C1b", sesoi=s_c1b, seed=SEED, B=500, n_flip=500)
    assert r["sesoi"] == pytest.approx(0.004) and r["contrast"] == "C1b"
    for c in ("C1", "C4", "C5"):
        assert M.confirmatory_contrast(d, lin, contrast=c, sesoi=0.01, seed=SEED, B=100, n_flip=100)["sesoi"] == 0.01


def test_gatekept_claims_downgrade_after_failure():
    """Review should-fix 2: a C4 that passes its own test after a C3 failure is 'estimate only'."""
    lin = lineage_series([20] * 7)
    rng = np.random.default_rng(SEED)
    good = pd.Series(rng.normal(0.03, 0.02, len(lin)), index=lin.index)
    kw = dict(seed=SEED, B=500, n_flip=500)
    res = {
        "C1": M.confirmatory_contrast(good, lin, contrast="C1", sesoi=0.01, **kw),
        "C1b": M.confirmatory_contrast(good, lin, contrast="C1b", sesoi=0.004, **kw),
        "C3": {"contrast": "C3", "passed": False, "estimate": 0.001, "ci_low": -0.002},
        "C4": M.confirmatory_contrast(good, lin, contrast="C4", sesoi=0.01, **kw),
        "C5": None,
    }
    assert res["C4"]["passed"]
    t = M.gatekept_claims(res).set_index("contrast")
    assert t.loc["C1", "claim"] == "meaningfully improves"
    assert t.loc["C3", "claim"] == "no detectable improvement"
    assert t.loc["C4", "claim"] == "estimate only (earlier gate failed)"
    assert t.loc["C5", "status"] == "estimate only"
    # C3 confirmed: only 'detectable' wording (no SESOI); a not-yet-run C4 is pending
    res["C3"] = {"contrast": "C3", "passed": True, "estimate": 0.05, "ci_low": 0.03}
    res["C4"] = None
    t2 = M.gatekept_claims(res).set_index("contrast")
    assert t2.loc["C3", "claim"] == "detectable improvement"
    assert t2.loc["C4", "claim"] == "pending" and t2.loc["C5", "claim"] == "estimate only (earlier gate failed)"
    with pytest.raises(ValueError):
        M.gatekept_claims({"C1": res["C1b"]})  # mislabelled result


def test_headroom_closed():
    lin = lineage_series([20, 20, 20])
    rng = np.random.default_rng(SEED)
    b2 = pd.Series(rng.normal(0.10, 0.02, 60), index=lin.index)
    p1 = 0.3 * b2 + rng.normal(0, 0.002, 60)
    h = M.headroom_closed(p1, b2, lin, B=2000, seed=SEED)
    assert h["reportable"] and h["estimate"] == pytest.approx(0.3, abs=0.01)
    assert h["ci_low"] < h["estimate"] < h["ci_high"]
    noisy = pd.Series(rng.normal(0, 0.05, 60), index=lin.index)
    assert not M.headroom_closed(p1, noisy, lin, B=2000, seed=SEED)["reportable"]
    with pytest.raises(M.PairingError):
        M.headroom_closed(p1.iloc[:-1], b2, lin, B=10, seed=SEED)


# ---------------------------------------------------------------------------------------------------
# C3: per-gene line bootstrap (vs brute force duplication) and within-lineage permutation null
# ---------------------------------------------------------------------------------------------------
def c3_frames(seed=SEED, n_lines=36, n_genes=25):
    rng = np.random.default_rng(seed)
    lin = lineage_series([12, 12, 12])
    Y = (rng.random((n_lines, n_genes)) < 0.35).astype(float)
    Y[rng.random(Y.shape) < 0.05] = np.nan
    Sa = np.round(np.clip(np.nan_to_num(Y) * 0.5 + rng.random(Y.shape) * 0.7, 0, 1), 1)  # informative, tied
    Sb = rng.integers(0, 3, Y.shape) / 2.0
    cols = [f"g{j}" for j in range(n_genes)]
    f = lambda a: pd.DataFrame(a, index=lin.index, columns=cols)  # noqa: E731
    return f(Sa), f(Sb), f(Y), lin


def test_gene_bootstrap_matches_duplication():
    Sa, Sb, L, lin = c3_frames()
    res = M.gene_bootstrap(Sa, Sb, L, lin, B=40, seed=SEED)
    est = M.m_delta_ap_gene(M.per_gene(Sa, L)) - M.m_delta_ap_gene(M.per_gene(Sb, L))
    assert res["estimate"] == pytest.approx(est)
    C = M.stratified_bootstrap_counts(lin, 40, SEED)
    elig = M.per_gene(Sa, L).query("estimable").index
    for b in range(40):
        rows = np.repeat(np.arange(len(lin)), C[b])
        vals = []
        for S in (Sa, Sb):
            s_dup = S[elig].iloc[rows].reset_index(drop=True)
            l_dup = L[elig].iloc[rows].reset_index(drop=True)
            g = M.per_gene(s_dup, l_dup, min_pos=1, min_neg=1)  # replicate rule: drop 0 pos / 0 neg only
            vals.append(g.loc[g["estimable"], "delta_null"].mean())
        assert res["replicates"][b] == pytest.approx(vals[0] - vals[1], abs=1e-10)


def test_c3_decision_combines_bootstrap_and_null():
    """Review should-fix 1: C3 = bootstrap LB > 0 (P1 - B3) AND P1 above the 97.5th permutation percentile."""
    Sa, Sb, L, lin = c3_frames()
    perfect = L.fillna(0.0)
    r = M.c3_decision(perfect, Sb, L, lin, seed=SEED, B=200, n_perm=100)
    assert r["bootstrap_passed"] and r["permutation_passed"] and r["passed"]
    assert r["estimate"] == pytest.approx(
        M.m_delta_ap_gene(M.per_gene(perfect, L)) - M.m_delta_ap_gene(M.per_gene(Sb, L)))
    r_again = M.c3_decision(perfect, Sb, L, lin, seed=SEED, B=200, n_perm=100)
    assert r_again["ci_low"] == r["ci_low"] and r_again["null_threshold"] == r["null_threshold"]
    # identical P1 and B3: the paired contrast is exactly 0, so the bootstrap arm fails even though P1 beats its null
    same = M.c3_decision(perfect, perfect, L, lin, seed=SEED, B=200, n_perm=100)
    assert same["permutation_passed"] and not same["bootstrap_passed"] and not same["passed"]
    # B3 better than P1: the permutation arm alone would not be enough
    worse = M.c3_decision(Sb, perfect, L, lin, seed=SEED, B=200, n_perm=100)
    assert not worse["passed"]
    # pairing vs B3 is asserted
    with pytest.raises(M.PairingError):
        M.c3_decision(perfect, Sb.drop(columns=Sb.columns[0]), L, lin, seed=SEED, B=10, n_perm=10)
    Sb_nan = Sb.copy()
    Sb_nan.iloc[0, 0] = np.nan
    if not np.isnan(L.iloc[0, 0]):
        with pytest.raises(M.PairingError):
            M.c3_decision(perfect, Sb_nan, L, lin, seed=SEED, B=10, n_perm=10)
    # §6.3 claim condition: the expressed-pairs version must agree in sign with CI LB > 0
    assert M.c3_claim_condition(r, {"estimate": 0.02, "ci_low": 0.001})
    assert not M.c3_claim_condition(r, {"estimate": 0.02, "ci_low": -0.001})
    assert not M.c3_claim_condition(same, {"estimate": 0.02, "ci_low": 0.001})
    # feeds the gate
    t = M.gatekept_claims({"C1": {"passed": True, "estimate": 0.02, "ci_low": 0.012, "sesoi": 0.01},
                           "C1b": {"passed": True, "estimate": 0.01, "ci_low": 0.002, "sesoi": 0.004},
                           "C3": r})
    assert t["claim"].tolist() == ["meaningfully improves", "improves", "detectable improvement"]


def test_permute_within_lineage_preserves_structure():
    rng = np.random.default_rng(SEED)
    _, _, L, lin = c3_frames()
    Y = L.to_numpy()
    P = M.permute_within_lineage(Y, lin.to_numpy(), rng)
    np.testing.assert_array_equal(np.isnan(P), np.isnan(Y))
    for k in lin.unique():
        m = (lin == k).to_numpy()
        np.testing.assert_array_equal(np.nansum(P[m], 0), np.nansum(Y[m], 0))
    assert not np.array_equal(np.nan_to_num(P), np.nan_to_num(Y))


def test_c3_permutation_null():
    Sa, _, L, lin = c3_frames()
    perfect = L.fillna(0.0)
    r = M.c3_permutation_null(perfect, L, lin, n_perm=200, seed=SEED)
    assert r["passed"] and r["p"] == pytest.approx(1 / 201)
    # a line-constant scorer carries no line-specific information: mΔ = 0 = every null draw -> no pass
    const = pd.DataFrame(np.tile(np.arange(L.shape[1], dtype=float), (len(L), 1)), index=L.index, columns=L.columns)
    rc = M.c3_permutation_null(const, L, lin, n_perm=50, seed=SEED)
    assert rc["observed"] == pytest.approx(0.0, abs=TOL) and not rc["passed"]
    np.testing.assert_allclose(rc["null"], 0.0, atol=1e-12)
    # a lineage-only scorer: labels depend on lineage only through the permutation-invariant counts
    lin_score = pd.DataFrame(
        np.repeat(np.arange(3.0)[:, None], 12, 0).repeat(1, 0) @ np.ones((1, L.shape[1])), index=L.index, columns=L.columns
    )
    rl = M.c3_permutation_null(lin_score, L, lin, n_perm=50, seed=SEED)
    np.testing.assert_allclose(rl["null"], rl["observed"], atol=1e-12)  # within-lineage permutation is invariant
    assert not rl["passed"]


# ---------------------------------------------------------------------------------------------------
# seeds (DESIGN §7): no unseeded randomness in metrics.py
# ---------------------------------------------------------------------------------------------------
def test_no_unseeded_randomness():
    src = Path(M.__file__).read_text()
    assert "default_rng()" not in src
    assert "random_state=None" not in src
    assert "np.random.seed" not in src
