"""Metrics and statistics of DESIGN.md v3 (frozen): §5.2, §6.1-6.6, §6.8, §11.4 (failed builds), §12.

Pure functions only: inputs are arrays / Series / DataFrames, nothing is read from disk, and every
random draw takes an explicit ``seed`` (int, ``np.random.SeedSequence`` or ``np.random.Generator``;
``None`` is refused, DESIGN §7).

Conventions
-----------
* Score matrices and label matrices are DataFrames **lines x genes** (rows = ModelID, columns = genes).
  Labels are {0, 1} with NaN for pairs excluded under §5.1; scores are finite on every non-NaN label pair.
* Higher score = more essential (e = 1 - ratio_r). Ties are exact float equality; scores are expected to
  be rounded upstream (ratio_r = clip(round(KO/WT, 6), 0, 1), §6.1).
* ``tap`` = tie-aware expected AP (McSherry & Najork 2008): the exact expectation of (uninterpolated) AP
  over uniformly random orderings within each tied score block. NaN when there is no positive.

Three algebraically identical tAP implementations are kept and cross-checked in the tests:
  ``tap_hypergeometric`` - the defining hypergeometric order-statistic sum (log-binomials, no overflow);
  ``tap``                - the O(n_g) per-item form of §6.1 (primary implementation);
  ``_SortedColumns``     - the O(n_g) form summed in closed form per block with harmonic numbers, with
                           integer item weights (= duplicated items); used by the bootstrap / permutation
                           code where the same ranking is re-scored thousands of times.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import product

import numpy as np
import pandas as pd
from scipy import special, stats

# ----------------------------------------------------------------------------------------------------
# Frozen constants (DESIGN.md v3)
# ----------------------------------------------------------------------------------------------------
MIN_POS_LINE = 5  # §5.2: AP_line needs >= 5 positives
MIN_POS_GENE = 5  # §5.2: AP_gene needs >= 5 positive ...
MIN_NEG_GENE = 5  # ... and >= 5 negative lines
CELL_MIN_LINE_FRAC = 0.5  # §5.2: grid cell NA if < 50% of its lines estimable ...
CELL_MIN_GENES = 30  # ... or < 30 of its genes estimable
ALPHA = 0.025  # §6.3: one-sided alpha of every confirmatory contrast
SESOI = 0.01  # §6.4: C1, C4, C5 (absolute, lineage-equal Delta)
GATE_ORDER = ("C1", "C1b", "C3", "C4", "C5")  # §6.3 fixed sequence
C3_NULL_PERCENTILE = 97.5  # §6.3
TRANSFER_FRACTION = 0.7  # §6.4: Delta_k > 0 in >= ceil(0.7 K) lineages
B_CONFIRMATORY = 10_000  # §6.5
B_OTHER = 2_000  # §6.5
N_SIGN_FLIPS = 10_000  # §6.3
N_C3_PERMUTATIONS = 2_000  # §6.3
FBA_CUTOFF = 0.01  # §6.1: binary call ratio_r < 0.01
L_GE06_THRESHOLD = -0.6  # §6.8: label L-ge06 = gene effect < -0.6
THRESHOLD_GRID = np.round(np.linspace(-1.0, 0.0, 21), 2)  # §6.8: Human1-paper robustness grid
BINARY_METRICS = ("mcc", "f1", "precision", "recall", "specificity", "accuracy")


class PairingError(AssertionError):
    """Methods / scores / labels are not on identical gene x line sets (§6.5 'Pairing')."""


# ----------------------------------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------------------------------
def _rng(seed) -> np.random.Generator:
    if seed is None:
        raise ValueError("an explicit seed is required (DESIGN §7: no unseeded randomness)")
    if isinstance(seed, np.random.Generator):
        return seed
    return np.random.default_rng(seed)


def _seed_sequence(seed) -> np.random.SeedSequence:
    if seed is None:
        raise ValueError("an explicit seed is required (DESIGN §7: no unseeded randomness)")
    if isinstance(seed, np.random.SeedSequence):
        return seed
    if isinstance(seed, np.random.Generator):
        return np.random.SeedSequence(int(seed.integers(2**63)))
    return np.random.SeedSequence(seed)


def harmonic(m) -> np.ndarray:
    """H_m = sum_{k=1..m} 1/k (H_0 = 0), elementwise, via digamma."""
    m = np.asarray(m, dtype=float)
    return np.where(m > 0, special.digamma(m + 1.0) + np.euler_gamma, 0.0)


def _harmonic_table(m_max: int) -> np.ndarray:
    k = np.arange(1, m_max + 1, dtype=float)
    return np.concatenate([[0.0], np.cumsum(1.0 / k)])


def _scores_labels_1d(scores, labels):
    s = np.asarray(scores, dtype=float).ravel()
    y = np.asarray(labels, dtype=float).ravel()
    if s.shape != y.shape:
        raise ValueError(f"scores {s.shape} and labels {y.shape} differ in length")
    if np.isnan(y).any():
        raise ValueError("labels contain NaN: exclude NaN pairs before scoring (§5.1)")
    if not np.isin(y, (0.0, 1.0)).all():
        raise ValueError("labels must be binary 0/1")
    if not np.isfinite(s).all():
        raise ValueError("scores must be finite on every scored pair")
    return s, y


def _sorted_blocks(s, y):
    """Sort by score descending; return sorted labels, block starts, sizes, positives per block."""
    order = np.argsort(-s, kind="stable")
    ss, ys = s[order], y[order]
    new = np.ones(ss.size, dtype=bool)
    new[1:] = ss[1:] != ss[:-1]
    starts = np.flatnonzero(new)
    sizes = np.diff(np.append(starts, ss.size))
    pos = np.add.reduceat(ys, starts) if ss.size else np.zeros(0)
    return ys, starts, sizes, pos


# ----------------------------------------------------------------------------------------------------
# §6.1 tie-aware expected AP and its null
# ----------------------------------------------------------------------------------------------------
def tap(scores, labels) -> float:
    """Tie-aware expected AP, O(n_g) form of DESIGN §6.1 (primary implementation).

    A block of n_g tied items with p_g positives, preceded by N items of which P are positive,
    contributes sum_{r=1..n_g} (p_g/n_g) [P + 1 + (r-1)(p_g-1)/(n_g-1)] / (N + r)
    (bracketed fraction = 0 when n_g = 1). tAP = sum over blocks / total positives. NaN if no positive.
    """
    s, y = _scores_labels_1d(scores, labels)
    n_pos = y.sum()
    if s.size == 0 or n_pos == 0:
        return float("nan")
    _, starts, sizes, pos = _sorted_blocks(s, y)
    p_before = np.cumsum(pos) - pos
    blk = np.repeat(np.arange(starts.size), sizes)
    i = np.arange(s.size)  # 0-based overall rank; N + r = i + 1
    r = i - starts[blk] + 1
    ng, pg, pb = sizes[blk].astype(float), pos[blk], p_before[blk]
    frac = np.where(ng > 1, (r - 1) * (pg - 1) / np.maximum(ng - 1, 1), 0.0)
    term = (pg / ng) * (pb + 1 + frac) / (i + 1)
    return float(term.sum() / n_pos)


def _log_comb(n, k):
    n, k = np.asarray(n, float), np.asarray(k, float)
    ok = (k >= 0) & (k <= n) & (n >= 0)
    out = special.gammaln(np.where(ok, n, 0) + 1) - special.gammaln(np.where(ok, k, 0) + 1) - special.gammaln(
        np.where(ok, n - k, 0) + 1
    )
    return np.where(ok, out, -np.inf)


def tap_hypergeometric(scores, labels) -> float:
    """Tie-aware expected AP, defining hypergeometric form of DESIGN §6.1 (cross-check implementation).

    The j-th positive of a block falls at in-block rank r with probability
    C(r-1, j-1) C(n_g-r, p_g-j) / C(n_g, p_g); its precision there is (P + j)/(N + r).
    Log-binomials avoid overflow. O(sum_g n_g p_g): use for verification, not for production loops.
    """
    s, y = _scores_labels_1d(scores, labels)
    n_pos = y.sum()
    if s.size == 0 or n_pos == 0:
        return float("nan")
    _, starts, sizes, pos = _sorted_blocks(s, y)
    total, n_before, p_before = 0.0, 0, 0
    for ng, pg in zip(sizes.tolist(), pos.astype(int).tolist()):
        if pg > 0:
            r = np.arange(1, ng + 1)[None, :]
            j = np.arange(1, pg + 1)[:, None]
            logp = _log_comb(r - 1, j - 1) + _log_comb(ng - r, pg - j) - _log_comb(ng, pg)
            total += float((np.exp(logp) * (p_before + j) / (n_before + r)).sum())
        n_before += ng
        p_before += pg
    return total / n_pos


def expected_ap_random(n, n_pos):
    """E[AP_random](n, n_pos): expected AP of a uniformly random ranking = tAP of a constant scorer.

    Closed form of the single-block O(n_g) sum divided by p: [H_n + (p-1)/(n-1) (n - H_n)] / n; NaN if p = 0.
    Vectorised over array inputs.
    """
    n = np.asarray(n, dtype=float)
    p = np.asarray(n_pos, dtype=float)
    if np.any(p > n) or np.any(p < 0):
        raise ValueError("need 0 <= n_pos <= n")
    hn = harmonic(n)
    c = np.where(n > 1, (p - 1) / np.where(n > 1, n - 1, 1), 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(p > 0, (hn + c * (n - hn)) / np.where(n > 0, n, 1), np.nan)
    return float(out) if out.ndim == 0 else out


def ap_step_sklearn(scores, labels) -> float:
    """sklearn step-wise AP: **sensitivity row only** (§6.1). NaN (not sklearn's 0.0) when no positive."""
    from sklearn.metrics import average_precision_score

    s, y = _scores_labels_1d(scores, labels)
    if y.sum() == 0:
        return float("nan")
    return float(average_precision_score(y, s))


def auroc(scores, labels) -> float:
    """AUROC (secondary, §6.2); sklearn's ROC is tie-aware (mid-rank). NaN if only one class."""
    from sklearn.metrics import roc_auc_score

    s, y = _scores_labels_1d(scores, labels)
    if y.sum() == 0 or y.sum() == y.size:
        return float("nan")
    return float(roc_auc_score(y, s))


def expected_tp_at_k(scores, labels, k=None) -> float:
    """Expected number of positives among the top k under random tie-breaking (§6.4 translation).

    Default k = number of positives. A block straddling k contributes p_g (k - N)/n_g.
    """
    s, y = _scores_labels_1d(scores, labels)
    k = int(y.sum()) if k is None else int(k)
    _, starts, sizes, pos = _sorted_blocks(s, y)
    take = np.clip(k - starts, 0, sizes)
    return float((pos * take / sizes).sum())


# ----------------------------------------------------------------------------------------------------
# Vectorised tAP over the columns of a matrix
# ----------------------------------------------------------------------------------------------------
def tap_columns(scores, labels):
    """tAP of every column of ``scores`` (n_items x n_cols) against ``labels`` (NaN = excluded item).

    O(n_g) per-item form, vectorised. Returns (tap, n, n_pos) arrays of length n_cols.
    """
    S = np.asarray(scores, dtype=float)
    Y = np.asarray(labels, dtype=float)
    if S.shape != Y.shape or S.ndim != 2:
        raise ValueError("scores and labels must be 2-D arrays of equal shape")
    mask = ~np.isnan(Y)
    if not np.isfinite(S[mask]).all():
        raise PairingError("non-finite score on a scored (non-NaN label) pair")
    if not np.isin(Y[mask], (0.0, 1.0)).all():
        raise ValueError("labels must be binary 0/1 or NaN")
    n_rows, n_cols = S.shape
    n = mask.sum(0)
    Y0 = np.where(mask, Y, 0.0)
    n_pos = Y0.sum(0)
    if n_rows == 0 or n_cols == 0:
        return np.full(n_cols, np.nan), n, n_pos
    S2 = np.where(mask, S, -np.inf)  # excluded items form a trailing block with no positive
    order = np.argsort(-S2, axis=0, kind="stable")
    ss = np.take_along_axis(S2, order, 0).T.ravel()
    ys = np.take_along_axis(Y0, order, 0).T.ravel()
    new = np.ones(ss.size, dtype=bool)
    new[1:] = ss[1:] != ss[:-1]
    new[::n_rows] = True
    starts = np.flatnonzero(new)
    sizes = np.diff(np.append(starts, ss.size))
    pos = np.add.reduceat(ys, starts)
    col_of_block = starts // n_rows
    cum_excl = np.cumsum(pos) - pos
    first_block = np.searchsorted(col_of_block, np.arange(n_cols))
    p_before = cum_excl - cum_excl[first_block[col_of_block]]
    blk = np.repeat(np.arange(starts.size), sizes)
    idx = np.arange(ss.size)
    i_in_col = idx % n_rows
    r = idx - starts[blk] + 1
    ng, pg, pb = sizes[blk].astype(float), pos[blk], p_before[blk]
    frac = np.where(ng > 1, (r - 1) * (pg - 1) / np.maximum(ng - 1, 1), 0.0)
    term = (pg / ng) * (pb + 1 + frac) / (i_in_col + 1)
    sums = term.reshape(n_cols, n_rows).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(n_pos > 0, sums / np.where(n_pos > 0, n_pos, 1), np.nan)
    return out, n, n_pos


class _SortedColumns:
    """Fixed per-column ranking (scores fixed) re-scored under item weights and/or new labels.

    Integer weights c_i are equivalent to c_i tied copies of item i (bootstrap resampling); labels may
    be replaced (label permutation). Each block's O(n_g) sum is evaluated in closed form:
      sum_{r=1..n} (p/n)[P + 1 + (r-1)c]/(N + r) = (p/n)[(P+1) dH + c (n - (N+1) dH)],
      dH = H_{N+n} - H_N, c = (p-1)/(n-1) (0 if n <= 1).
    """

    def __init__(self, scores, mask):
        S = np.asarray(scores, dtype=float)
        M = np.asarray(mask, dtype=bool)
        if not np.isfinite(S[M]).all():
            raise PairingError("non-finite score on a scored pair")
        self.n_rows, self.n_cols = S.shape
        S2 = np.where(M, S, -np.inf)
        self.order = np.argsort(-S2, axis=0, kind="stable")
        ss = np.take_along_axis(S2, self.order, 0).T.ravel()
        self.mask_flat = np.take_along_axis(M, self.order, 0).T.ravel().astype(float)
        self.rows_flat = self.order.T.ravel()
        new = np.ones(ss.size, dtype=bool)
        new[1:] = ss[1:] != ss[:-1]
        new[:: self.n_rows] = True
        self.starts = np.flatnonzero(new)
        self.col_of_block = self.starts // self.n_rows
        self.first_block = np.searchsorted(self.col_of_block, np.arange(self.n_cols))

    def sorted_labels(self, labels):
        Y0 = np.nan_to_num(np.asarray(labels, dtype=float), nan=0.0)
        return np.take_along_axis(Y0, self.order, 0).T.ravel() * self.mask_flat

    def tap(self, weights, ys_flat, h_table=None):
        """weights: (b, n_rows) non-negative integers. Returns (tap, n, n_pos), each (b, n_cols)."""
        W = np.asarray(weights, dtype=float)
        if W.ndim == 1:
            W = W[None, :]
        Wf = W[:, self.rows_flat] * self.mask_flat
        nb = np.add.reduceat(Wf, self.starts, axis=1)
        pb = np.add.reduceat(Wf * ys_flat, self.starts, axis=1)
        cn = np.cumsum(nb, axis=1) - nb
        cp = np.cumsum(pb, axis=1) - pb
        fb = self.first_block[self.col_of_block]
        N = cn - cn[:, fb]
        P = cp - cp[:, fb]
        if h_table is None:
            h_table = _harmonic_table(int(np.ceil(W.sum(1).max())) + 1)
        Ni = np.rint(N).astype(np.int64)
        ni = np.rint(nb).astype(np.int64)
        dH = h_table[Ni + ni] - h_table[Ni]
        with np.errstate(invalid="ignore", divide="ignore"):
            c = np.where(nb > 1, (pb - 1) / np.where(nb > 1, nb - 1, 1), 0.0)
            contrib = np.where(pb > 0, (pb / np.where(nb > 0, nb, 1)) * ((P + 1) * dH + c * (nb - (N + 1) * dH)), 0.0)
        tot = np.add.reduceat(contrib, self.first_block, axis=1)
        n = np.add.reduceat(nb, self.first_block, axis=1)
        n_pos = np.add.reduceat(pb, self.first_block, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = np.where(n_pos > 0, tot / np.where(n_pos > 0, n_pos, 1), np.nan)
        return out, n, n_pos


def tap_weighted(scores, labels, weights) -> float:
    """tAP with integer item multiplicities (= duplicated items); harmonic closed form (cross-checked)."""
    s, y = _scores_labels_1d(scores, labels)
    w = np.asarray(weights, dtype=float).ravel()
    if w.shape != s.shape or (w < 0).any() or not np.allclose(w, np.rint(w)):
        raise ValueError("weights must be non-negative integers, one per item")
    sc = _SortedColumns(s[:, None], np.ones((s.size, 1), bool))
    out, _, _ = sc.tap(w[None, :], sc.sorted_labels(y[:, None]))
    return float(out[0, 0])


# ----------------------------------------------------------------------------------------------------
# Pairing and alignment (§6.5)
# ----------------------------------------------------------------------------------------------------
def assert_same_pairs(*frames: pd.DataFrame) -> None:
    """Raise PairingError unless every frame has the same line (index) and gene (column) sets."""
    ref = frames[0]
    for k, f in enumerate(frames[1:], start=1):
        if f.index.has_duplicates or f.columns.has_duplicates:
            raise PairingError(f"frame {k} has duplicated line or gene labels")
        if set(f.index) != set(ref.index):
            raise PairingError(
                f"line sets differ (frame 0 vs {k}): {len(set(ref.index) ^ set(f.index))} lines not shared"
            )
        if set(f.columns) != set(ref.columns):
            raise PairingError(
                f"gene sets differ (frame 0 vs {k}): {len(set(ref.columns) ^ set(f.columns))} genes not shared"
            )
    if ref.index.has_duplicates or ref.columns.has_duplicates:
        raise PairingError("frame 0 has duplicated line or gene labels")


def _align(scores: pd.DataFrame, labels: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    assert_same_pairs(labels, scores)
    S = scores.loc[labels.index, labels.columns].to_numpy(dtype=float)
    Y = labels.to_numpy(dtype=float)
    if not np.isfinite(S[~np.isnan(Y)]).all():
        raise PairingError("scores are missing (non-finite) on scored pairs; methods would not be paired")
    return S, Y


def _lineage_of(index, lineage: pd.Series) -> pd.Series:
    missing = pd.Index(index).difference(lineage.index)
    if len(missing):
        raise PairingError(f"{len(missing)} lines have no lineage")
    lin = lineage.loc[index]
    if lin.isna().any():
        raise PairingError("lineage is NaN for some lines")
    return lin


# ----------------------------------------------------------------------------------------------------
# §6.2 units and §5.2 estimability
# ----------------------------------------------------------------------------------------------------
def _unit_table(tap_values, n, n_pos, index, estimable) -> pd.DataFrame:
    n = np.asarray(n, float)
    n_pos = np.asarray(n_pos, float)
    e = np.full(n.size, np.nan)
    has = n_pos > 0
    e[has] = expected_ap_random(n[has], n_pos[has])
    t = np.where(estimable, tap_values, np.nan)
    e = np.where(estimable, e, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        df = pd.DataFrame(
            {
                "n": n.astype(int),
                "n_pos": n_pos.astype(int),
                "prevalence": np.where(n > 0, n_pos / np.where(n > 0, n, 1), np.nan),
                "estimable": estimable,
                "tap": t,
                "e_random": e,
                "lift": t / e,
                "delta_null": t - e,
            },
            index=index,
        )
    return df


def per_line(scores: pd.DataFrame, labels: pd.DataFrame, min_pos: int = MIN_POS_LINE) -> pd.DataFrame:
    """tAP_line for every line (over genes). NA (tap NaN, estimable False) if < ``min_pos`` positives (§5.2).

    Columns: n, n_pos, prevalence, estimable, tap, e_random, lift (= tap/e_random), delta_null (= tap - e_random).
    """
    S, Y = _align(scores, labels)
    t, n, n_pos = tap_columns(S.T, Y.T)
    return _unit_table(t, n, n_pos, labels.index, n_pos >= min_pos)


def per_gene(
    scores: pd.DataFrame, labels: pd.DataFrame, min_pos: int = MIN_POS_GENE, min_neg: int = MIN_NEG_GENE
) -> pd.DataFrame:
    """tAP_gene for every gene (over lines). Estimable if >= ``min_pos`` positive and >= ``min_neg`` negative lines."""
    S, Y = _align(scores, labels)
    t, n, n_pos = tap_columns(S, Y)
    out = _unit_table(t, n, n_pos, labels.columns, (n_pos >= min_pos) & (n - n_pos >= min_neg))
    out.insert(2, "n_neg", out["n"] - out["n_pos"])
    return out


def pooled(scores: pd.DataFrame, labels: pd.DataFrame, lines=None) -> pd.Series:
    """Pooled tAP over all non-NaN gene x line pairs (optionally restricted to ``lines``, e.g. one lineage)."""
    if lines is not None:
        scores, labels = scores.loc[list(lines)], labels.loc[list(lines)]
    S, Y = _align(scores, labels)
    m = ~np.isnan(Y)
    n, n_pos = int(m.sum()), int(Y[m].sum())
    t = tap(S[m], Y[m]) if n_pos > 0 else np.nan
    e = expected_ap_random(n, n_pos) if n_pos > 0 else np.nan
    return pd.Series(
        {"n": n, "n_pos": n_pos, "prevalence": n_pos / n if n else np.nan, "tap": t, "e_random": e,
         "lift": t / e, "delta_null": t - e}
    )


def m_delta_ap_gene(gene_table: pd.DataFrame) -> float:
    """mΔAP_gene = mean over estimable genes of tAP_gene - E[AP_random] (§6.2, C3)."""
    d = gene_table.loc[gene_table["estimable"], "delta_null"]
    return float(d.mean()) if len(d) else float("nan")


def cell_is_na(line_table: pd.DataFrame, gene_table: pd.DataFrame | None = None) -> bool:
    """§5.2 grid-cell rule: NA if < 50% of lines estimable or (if a gene table is given) < 30 genes estimable."""
    frac = line_table["estimable"].mean() if len(line_table) else 0.0
    if frac < CELL_MIN_LINE_FRAC:
        return True
    if gene_table is not None and int(gene_table["estimable"].sum()) < CELL_MIN_GENES:
        return True
    return False


def paired_line_delta(table_a: pd.DataFrame, table_b: pd.DataFrame) -> pd.Series:
    """Per-line ΔtAP = tAP_a - tAP_b on the estimable lines; asserts identical pairing (§6.5)."""
    if set(table_a.index) != set(table_b.index):
        raise PairingError("line sets differ between methods")
    b = table_b.loc[table_a.index]
    for col in ("n", "n_pos", "estimable"):
        if not (table_a[col].to_numpy() == b[col].to_numpy()).all():
            raise PairingError(f"per-line '{col}' differs between methods: not scored on identical pairs")
    d = (table_a["tap"] - b["tap"])[table_a["estimable"]]
    d.name = "delta_tap"
    return d


def no_information_scores(scores: pd.DataFrame, failed_lines, unit: str) -> pd.DataFrame:
    """No-information score for failed builds (§11.4), applied identically in every arm.

    unit='line': failed lines get all genes tied (constant 0), so tAP_line = E[AP_random].
    unit='gene': a failed line takes, per gene, the method's median score for that gene over its
    non-failed lines in the same tier (``scores`` must hold exactly that tier); it enters as a tie.
    """
    failed = pd.Index(failed_lines)
    if len(failed.difference(scores.index)):
        raise PairingError("failed lines not in the score matrix")
    out = scores.copy()
    if len(failed) == 0:
        return out
    if unit == "line":
        out.loc[failed, :] = 0.0
    elif unit == "gene":
        ok = scores.drop(index=failed)
        if ok.empty:
            raise ValueError("every line failed: no median available")
        med = ok.median(axis=0)
        out.loc[failed, :] = np.broadcast_to(med.to_numpy(), (len(failed), scores.shape[1]))
    else:
        raise ValueError("unit must be 'line' or 'gene'")
    return out


# ----------------------------------------------------------------------------------------------------
# Lineage-equal and line-weighted summaries
# ----------------------------------------------------------------------------------------------------
def _clean(values: pd.Series, lineage: pd.Series) -> tuple[pd.Series, pd.Series]:
    v = pd.Series(values, dtype=float).dropna()
    return v, _lineage_of(v.index, lineage)


def lineage_weights(lineage: pd.Series) -> pd.Series:
    """w_l = 1 / (K n_k): lineage-equal-weighted mean = sum_l w_l v_l."""
    counts = lineage.map(lineage.value_counts())
    return 1.0 / (lineage.nunique() * counts.astype(float))


def lineage_means(values: pd.Series, lineage: pd.Series) -> pd.Series:
    v, lin = _clean(values, lineage)
    return v.groupby(lin).mean()


def lineage_equal_mean(values: pd.Series, lineage: pd.Series) -> float:
    """Δ̄_lin: mean over lineages of the mean over that lineage's (non-NaN) lines (§6.2 primary summary)."""
    m = lineage_means(values, lineage)
    return float(m.mean()) if len(m) else float("nan")


def line_weighted_mean(values: pd.Series) -> float:
    """Δ̄_line: mean over (non-NaN) lines (secondary)."""
    v = pd.Series(values, dtype=float).dropna()
    return float(v.mean()) if len(v) else float("nan")


def per_lineage_summary(
    line_table: pd.DataFrame, lineage: pd.Series, reference: pd.DataFrame | None = None
) -> pd.DataFrame:
    """Per-lineage breakdown (§6.2): headline = mean tAP_line, with E[AP_random], π, lift and Δ vs reference.

    lift = mean tAP_line / mean E[AP_random] over the lineage's estimable lines (ratio of means).
    ``delta_vs_ref`` = mean paired line Δ against ``reference`` (e.g. G*). Never rank lineages by raw AP.
    """
    lin = _lineage_of(line_table.index, lineage)
    est = line_table[line_table["estimable"]]
    g = est.groupby(lin.loc[est.index])
    out = pd.DataFrame(
        {
            "n_lines": lin.value_counts(),
            "n_estimable": g.size(),
            "mean_tap": g["tap"].mean(),
            "mean_e_random": g["e_random"].mean(),
            "mean_prevalence": g["prevalence"].mean(),
        }
    )
    out["n_estimable"] = out["n_estimable"].fillna(0).astype(int)
    out["lift"] = out["mean_tap"] / out["mean_e_random"]
    out["delta_null"] = out["mean_tap"] - out["mean_e_random"]
    if reference is not None:
        d = paired_line_delta(line_table, reference)
        out["delta_vs_ref"] = d.groupby(lin.loc[d.index]).mean()
    out.index.name = "lineage"
    return out


# ----------------------------------------------------------------------------------------------------
# §6.5 bootstrap (lines resampled within lineage)
# ----------------------------------------------------------------------------------------------------
def stratified_bootstrap_counts(lineage: pd.Series, B: int, seed) -> np.ndarray:
    """(B, n_lines) resampling multiplicities: lines drawn with replacement within lineage, n_k per lineage.

    Columns follow ``lineage.index``; lineages are processed in sorted order for reproducibility.
    """
    rng = _rng(seed)
    lin = lineage.to_numpy()
    C = np.zeros((B, lin.size), dtype=np.int32)
    for k in sorted(pd.unique(lin), key=str):
        idx = np.flatnonzero(lin == k)
        C[:, idx] = rng.multinomial(idx.size, np.full(idx.size, 1.0 / idx.size), size=B)
    return C


@dataclass
class BootstrapResult:
    estimate: float
    ci_low: float
    ci_high: float
    method: str
    B: int
    n_lines: int
    n_lineages: int
    replicates: np.ndarray = field(repr=False)


def _weights(lin: pd.Series, weighting: str) -> np.ndarray:
    if weighting == "lineage":
        return lineage_weights(lin).to_numpy()
    if weighting == "line":
        return np.full(len(lin), 1.0 / len(lin))
    raise ValueError("weighting must be 'lineage' or 'line'")


def _jackknife_linear(v: np.ndarray, lin: pd.Series, weighting: str) -> np.ndarray:
    out = np.empty(v.size)
    labels = lin.to_numpy()
    for i in range(v.size):
        keep = np.ones(v.size, bool)
        keep[i] = False
        sub = pd.Series(labels[keep])
        out[i] = _weights(sub, weighting) @ v[keep]
    return out


def _percentile_or_bca(rep, est, conf, method, jack=None):
    a_lo, a_hi = (1 - conf) / 2, 1 - (1 - conf) / 2
    if method == "percentile":
        q = np.quantile(rep, [a_lo, a_hi])
        return float(q[0]), float(q[1])
    if method == "bca":
        prop = (np.sum(rep < est) + 0.5 * np.sum(rep == est)) / rep.size
        z0 = stats.norm.ppf(np.clip(prop, 1e-12, 1 - 1e-12))
        d = jack.mean() - jack
        den = 6.0 * (d**2).sum() ** 1.5
        acc = (d**3).sum() / den if den > 0 else 0.0
        qs = []
        for a in (a_lo, a_hi):
            z = stats.norm.ppf(a)
            qs.append(stats.norm.cdf(z0 + (z0 + z) / (1 - acc * (z0 + z))))
        q = np.quantile(rep, qs)
        return float(q[0]), float(q[1])
    raise ValueError("method must be 'percentile' or 'bca'")


def bootstrap_mean(
    values: pd.Series,
    lineage: pd.Series,
    *,
    B: int,
    seed,
    weighting: str = "lineage",
    conf: float = 0.95,
    method: str = "percentile",
) -> BootstrapResult:
    """Lineage-stratified line bootstrap CI of the lineage-equal (default) or line-weighted mean (§6.5).

    NaN values are dropped first (e.g. NA lines, NaN binary metrics). Percentile interval by default;
    'bca' (jackknife acceleration over lines) only if ``skew_triggers_bca`` held on training at M2.
    """
    v, lin = _clean(values, lineage)
    if v.empty:
        raise ValueError("no non-NaN values")
    w = _weights(lin, weighting)
    vv = v.to_numpy()
    est = float(w @ vv)
    C = stratified_bootstrap_counts(lin, B, seed)
    rep = C @ (w * vv)
    jack = _jackknife_linear(vv, lin, weighting) if method == "bca" else None
    lo, hi = _percentile_or_bca(rep, est, conf, method, jack)
    return BootstrapResult(est, lo, hi, method, B, v.size, lin.nunique(), rep)


def skew_triggers_bca(training_line_delta) -> bool:
    """§6.5: BCa is used only if the training-tier distribution of line ΔtAP has |skew| > 1."""
    v = pd.Series(training_line_delta, dtype=float).dropna().to_numpy()
    return bool(abs(stats.skew(v)) > 1)


# ----------------------------------------------------------------------------------------------------
# §6.3 tests: sign-flip (C1/C1b/C4/C5), lineage-level transfer (§6.4)
# ----------------------------------------------------------------------------------------------------
def sign_flip_test(delta: pd.Series, lineage: pd.Series, *, n_flip: int = N_SIGN_FLIPS, seed, batch: int = 2_000):
    """One-sided sign-flip permutation p for the lineage-equal Δ̄_lin (§6.3, C1).

    Each line's ΔtAP sign is flipped at random (independently per line, within its lineage), the
    lineage-equal-weighted mean is recomputed, and p = (1 + #{T* >= T_obs}) / (1 + n_flip).
    Returns (observed, p).
    """
    v, lin = _clean(delta, lineage)
    w = lineage_weights(lin).to_numpy()
    wd = w * v.to_numpy()
    obs = float(wd.sum())
    tol = 1e-12 * max(1.0, abs(obs))
    rng = _rng(seed)
    count, done = 0, 0
    while done < n_flip:
        b = min(batch, n_flip - done)
        signs = rng.integers(0, 2, size=(b, wd.size), dtype=np.int8) * 2 - 1
        count += int(np.sum(signs @ wd >= obs - tol))
        done += b
    return obs, (1 + count) / (1 + n_flip)


def lineage_t_interval(lineage_means_: Sequence[float], conf: float = 0.95) -> tuple[float, float, float]:
    """Two-sided t-interval (df = K - 1) on the K lineage means (§6.4 (2)). Returns (mean, low, high)."""
    m = np.asarray(lineage_means_, dtype=float)
    k = m.size
    mean = float(m.mean())
    if k < 2:
        return mean, float("nan"), float("nan")
    half = stats.t.ppf(1 - (1 - conf) / 2, k - 1) * m.std(ddof=1) / math.sqrt(k)
    return mean, float(mean - half), float(mean + half)


def exact_sign_flip_lineage(lineage_means_: Sequence[float]) -> float:
    """Exact one-sided sign-flip p over all 2^K sign configurations of the lineage means (§6.4 (2))."""
    m = np.asarray(lineage_means_, dtype=float)
    obs = m.mean()
    tol = 1e-12 * max(1.0, abs(obs))
    signs = np.array(list(product((1.0, -1.0), repeat=m.size)))
    return float(np.mean(signs @ m / m.size >= obs - tol))


def transfer_summary(delta: pd.Series, lineage: pd.Series, conf: float = 0.95, fraction: float = TRANSFER_FRACTION):
    """§6.4 transfer: lineage t-interval, exact 2^K sign-flip p, k-of-K count; 'generalises' iff t low > 0 and k >= ceil(0.7K)."""
    m = lineage_means(delta, lineage)
    mean, lo, hi = lineage_t_interval(m.to_numpy(), conf)
    k = m.size
    k_req = math.ceil(round(fraction * k, 9))
    k_pos = int((m > 0).sum())
    return {
        "lineage_means": m,
        "K": k,
        "mean": mean,
        "t_low": lo,
        "t_high": hi,
        "exact_sign_flip_p": exact_sign_flip_lineage(m.to_numpy()),
        "k_positive": k_pos,
        "k_required": k_req,
        "generalises": bool(lo > 0 and k_pos >= k_req),
    }


# ----------------------------------------------------------------------------------------------------
# §6.4 claim language and SESOIs; §6.3 gatekeeping
# ----------------------------------------------------------------------------------------------------
GATE_STATUSES = ("confirmed", "not confirmed", "pending", "estimate only")
# §6.4 frozen SESOI table. C1b's SESOI is lift-scaled and computed once on training lines at M2 with
# ``sesoi_c1b``; it must be passed explicitly (None here). C3 has no SESOI: only 'detectable' wording.
SESOI_BY_CONTRAST = {"C1": SESOI, "C1b": None, "C4": SESOI, "C5": SESOI}


def claim_tier(passed: bool | None, estimate: float, ci_low: float, sesoi: float | None, *, gate_status: str) -> str:
    """§6.4 wording, conditional on the §6.3 gatekeeping state (``gate_status`` from ``fixed_sequence``).

    A contrast that is not confirmatory ('estimate only' after an earlier failure, or 'pending') gets no
    claim tier, whatever its own test says. For a confirmatory contrast the tiers are nested:
    'detectable improvement' (test passes) < 'improves' (estimate >= SESOI) < 'meaningfully improves'
    (CI lower bound >= SESOI as well). ``sesoi=None`` (C3) allows only 'detectable' wording.
    """
    if gate_status not in GATE_STATUSES:
        raise ValueError(f"gate_status must be one of {GATE_STATUSES}")
    if gate_status == "estimate only":
        return "estimate only (earlier gate failed)"
    if gate_status == "pending":
        return "pending"
    if bool(passed) != (gate_status == "confirmed"):
        raise ValueError(f"inconsistent: passed={passed} but gate_status={gate_status!r}")
    if not passed:
        return "no detectable improvement"
    if sesoi is None:
        return "detectable improvement"
    if estimate >= sesoi and ci_low >= sesoi:
        return "meaningfully improves"
    if estimate >= sesoi:
        return "improves"
    return "detectable improvement"


def sesoi_c1b(e_random_minus_ce, e_random_all, sesoi: float = SESOI) -> float:
    """SESOI_C1b = 0.01 x mean E[AP_random] on U minus CE_fit / mean E[AP_random] on U (training lines, M2)."""
    a = np.nanmean(np.asarray(e_random_minus_ce, dtype=float))
    b = np.nanmean(np.asarray(e_random_all, dtype=float))
    return float(sesoi * a / b)


def fixed_sequence(results) -> pd.DataFrame:
    """Fixed-sequence gatekeeping (§6.3) in the frozen order C1 -> C1b -> C3 -> C4 -> C5.

    ``results``: a sequence [(contrast, passed or None if not yet run), ...] whose names must be exactly a
    prefix of ``GATE_ORDER`` in that order (anything else raises ValueError), or a mapping
    {contrast: passed}, which is put into ``GATE_ORDER`` (its keys must form a prefix of it).
    Testing stops at the first failure; later contrasts are 'estimate only'. A contrast not yet run
    (None) is 'pending' and nothing after it is confirmatory yet.
    """
    if isinstance(results, dict):
        unknown = set(results) - set(GATE_ORDER)
        if unknown:
            raise ValueError(f"unknown contrasts {sorted(unknown)}; gate order is {GATE_ORDER}")
        results = [(c, results[c]) for c in GATE_ORDER if c in results]
    results = list(results)
    names = [name for name, _ in results]
    if names != list(GATE_ORDER[: len(names)]):
        raise ValueError(f"contrasts {names} are not a prefix of the frozen gate order {GATE_ORDER}")
    rows, alive = [], True
    for name, passed in results:
        if not alive:
            rows.append((name, passed, False, "estimate only"))
        elif passed is None:
            rows.append((name, None, False, "pending"))
            alive = False
        else:
            rows.append((name, bool(passed), True, "confirmed" if passed else "not confirmed"))
            alive = bool(passed)
    return pd.DataFrame(rows, columns=["contrast", "passed", "confirmatory", "status"])


def confirmatory_contrast(
    delta: pd.Series,
    lineage: pd.Series,
    *,
    contrast: str,
    sesoi: float,
    seed,
    alpha: float = ALPHA,
    B: int = B_CONFIRMATORY,
    n_flip: int = N_SIGN_FLIPS,
    method: str = "percentile",
) -> dict:
    """C1 / C1b / C4 / C5 (§6.3-6.4): Δ̄_lin, stratified bootstrap CI, sign-flip p, test result, transfer.

    ``contrast`` and ``sesoi`` are required (no default). C1, C4 and C5 must use the frozen 0.01;
    C1b must use its own lift-scaled SESOI from ``sesoi_c1b`` (passing 0.01 for C1b raises).
    No claim tier is returned here: the wording depends on the gatekeeping state, so use
    ``gatekept_claims`` over all contrasts. ``seed`` is spawned into [bootstrap, sign-flip] streams.
    Test passes iff p <= alpha (one-sided 0.025).
    """
    if contrast not in SESOI_BY_CONTRAST:
        raise ValueError(f"contrast must be one of {tuple(SESOI_BY_CONTRAST)} (C3 uses c3_decision)")
    frozen = SESOI_BY_CONTRAST[contrast]
    if frozen is not None and sesoi != frozen:
        raise ValueError(f"{contrast} has the frozen SESOI {frozen}, got {sesoi}")
    if frozen is None and not (0 < sesoi < SESOI):
        raise ValueError(f"C1b needs its lift-scaled SESOI from sesoi_c1b (0 < SESOI_C1b < {SESOI}), got {sesoi}")
    ss_boot, ss_flip = _seed_sequence(seed).spawn(2)
    boot = bootstrap_mean(delta, lineage, B=B, seed=ss_boot, method=method)
    obs, p = sign_flip_test(delta, lineage, n_flip=n_flip, seed=ss_flip)
    return {
        "contrast": contrast,
        "estimate": boot.estimate,
        "ci_low": boot.ci_low,
        "ci_high": boot.ci_high,
        "p_sign_flip": p,
        "passed": bool(p <= alpha),
        "sesoi": sesoi,
        "line_weighted": line_weighted_mean(delta),
        "n_lines": boot.n_lines,
        "n_lineages": boot.n_lineages,
        "transfer": transfer_summary(delta, lineage),
    }


def gatekept_claims(results: dict) -> pd.DataFrame:
    """Apply the §6.3 fixed sequence and the §6.4 wording to a set of contrast results.

    ``results`` maps contrast name -> result dict (from ``confirmatory_contrast`` or ``c3_decision``,
    with keys passed, estimate, ci_low and, except C3, sesoi) or None if not yet run. Keys must form a
    prefix of ``GATE_ORDER``. Returns one row per contrast with its gate status and claim.
    """
    seq = fixed_sequence({c: (None if r is None else bool(r["passed"])) for c, r in results.items()})
    claims = []
    for row in seq.itertuples(index=False):
        r = results[row.contrast]
        if r is None:
            claims.append((np.nan, np.nan, None, claim_tier(None, np.nan, np.nan, None, gate_status=row.status)))
            continue
        if r.get("contrast", row.contrast) != row.contrast:
            raise ValueError(f"result for {row.contrast} is labelled {r['contrast']!r}")
        sesoi = None if row.contrast == "C3" else r["sesoi"]
        claims.append(
            (r["estimate"], r["ci_low"], sesoi,
             claim_tier(r["passed"], r["estimate"], r["ci_low"], sesoi, gate_status=row.status))
        )
    extra = pd.DataFrame(claims, columns=["estimate", "ci_low", "sesoi", "claim"])
    return pd.concat([seq, extra], axis=1)


def headroom_closed(
    delta_p1: pd.Series, delta_b2: pd.Series, lineage: pd.Series, *, B: int = B_OTHER, seed, conf: float = 0.95
) -> dict:
    """Δ̄_lin(P1 - G*) / Δ̄_lin(B2 - G*) with a shared-resample bootstrap CI (§6.4).

    ``reportable`` is False unless the denominator's CI excludes 0 (then report both raw contrasts).
    """
    a = pd.Series(delta_p1, dtype=float)
    b = pd.Series(delta_b2, dtype=float)
    if set(a.dropna().index) != set(b.dropna().index):
        raise PairingError("headroom contrasts must be on the same lines")
    a = a.dropna()
    b = b.loc[a.index]
    lin = _lineage_of(a.index, lineage)
    w = lineage_weights(lin).to_numpy()
    C = stratified_bootstrap_counts(lin, B, seed)
    num, den = C @ (w * a.to_numpy()), C @ (w * b.to_numpy())
    est_den = float(w @ b.to_numpy())
    est = float(w @ a.to_numpy()) / est_den
    q = [(1 - conf) / 2, 1 - (1 - conf) / 2]
    den_lo, den_hi = np.quantile(den, q)
    with np.errstate(divide="ignore", invalid="ignore"):
        lo, hi = np.quantile(num / den, q)
    return {
        "estimate": est,
        "ci_low": float(lo),
        "ci_high": float(hi),
        "denominator": est_den,
        "denominator_ci": (float(den_lo), float(den_hi)),
        "reportable": bool(den_lo > 0 or den_hi < 0),
    }


# ----------------------------------------------------------------------------------------------------
# C3: per-gene contrast, line bootstrap and within-lineage label-permutation null
# ----------------------------------------------------------------------------------------------------
def _batch_size(n_rows: int, n_cols: int, budget: float = 2e7) -> int:
    return max(1, int(budget // max(1, n_rows * n_cols)))


def gene_bootstrap(
    scores_a: pd.DataFrame,
    scores_b: pd.DataFrame | None,
    labels: pd.DataFrame,
    lineage: pd.Series,
    *,
    B: int,
    seed,
    conf: float = 0.95,
) -> dict:
    """Line bootstrap (within lineage) of mΔAP_gene(a) [- mΔAP_gene(b)] (§6.2-6.3, C3).

    Eligible genes are fixed on the observed data (>= 5 positive and >= 5 negative lines). In each
    replicate, tAP_gene is recomputed over the resampled lines (duplicates = tied copies); genes with
    0 positives or 0 negatives in the replicate are dropped, and the mean count dropped is reported.
    """
    ga = per_gene(scores_a, labels)
    elig = ga.index[ga["estimable"]]
    est = m_delta_ap_gene(ga)
    if scores_b is not None:
        assert_same_pairs(labels, scores_b)
        est -= m_delta_ap_gene(per_gene(scores_b, labels))
    lin = _lineage_of(labels.index, lineage)
    Y = labels.loc[:, elig].to_numpy(dtype=float)
    M = ~np.isnan(Y)
    structs = [_SortedColumns(s.loc[labels.index, elig].to_numpy(float), M) for s in (scores_a, scores_b) if s is not None]
    ys = [st.sorted_labels(Y) for st in structs]
    C = stratified_bootstrap_counts(lin, B, seed)
    h = _harmonic_table(len(lin) + 1)
    reps, dropped = np.empty(B), np.empty(B)
    bs = _batch_size(*Y.shape)
    for s0 in range(0, B, bs):
        cb = C[s0 : s0 + bs]
        vals = []
        for st, y in zip(structs, ys):
            t, n, p = st.tap(cb, y, h)
            keep = (p > 0) & (p < n)
            e = np.where(keep, expected_ap_random(np.where(keep, n, 1), np.where(keep, p, 1)), np.nan)
            d = np.where(keep, t - e, np.nan)
            with np.errstate(invalid="ignore"):
                vals.append(np.nanmean(d, axis=1))
        reps[s0 : s0 + bs] = vals[0] - (vals[1] if len(vals) > 1 else 0.0)
        dropped[s0 : s0 + bs] = (~keep).sum(1)
    lo, hi = np.quantile(reps, [(1 - conf) / 2, 1 - (1 - conf) / 2])
    return {
        "estimate": float(est),
        "ci_low": float(lo),
        "ci_high": float(hi),
        "n_genes_eligible": int(len(elig)),
        "mean_genes_dropped": float(dropped.mean()),
        "replicates": reps,
    }


def permute_within_lineage(labels: np.ndarray, lineage: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Permute each column's non-NaN labels among that column's non-NaN rows, within each lineage."""
    Y = np.asarray(labels, dtype=float)
    out = Y.copy()
    for k in sorted(pd.unique(lineage), key=str):
        idx = np.flatnonzero(lineage == k)
        Yk = Y[idx]
        nanm = np.isnan(Yk)
        keys = rng.random(Yk.shape)
        keys[nanm] = np.inf
        src = np.argsort(keys, axis=0, kind="stable")  # non-NaN rows first, random order
        dst = np.argsort(nanm, axis=0, kind="stable")  # non-NaN rows first, original order
        newk = np.empty_like(Yk)
        np.put_along_axis(newk, dst, np.take_along_axis(Yk, src, 0), 0)
        out[idx] = newk
    return out


def c3_permutation_null(
    scores: pd.DataFrame,
    labels: pd.DataFrame,
    lineage: pd.Series,
    *,
    n_perm: int = N_C3_PERMUTATIONS,
    seed,
    percentile: float = C3_NULL_PERCENTILE,
) -> dict:
    """Within-lineage label-permutation null of mΔAP_gene (§6.3 C3).

    Each eligible gene's labels are permuted across lines within each lineage (NaN pairs stay put, so
    per-gene counts and eligibility are unchanged); permutation streams are spawned from ``seed``.
    Pass iff observed mΔAP_gene > the 97.5th percentile of the null (+1e-12 floating-point tolerance).
    """
    g = per_gene(scores, labels)
    elig = g.index[g["estimable"]]
    obs = m_delta_ap_gene(g)
    e = g.loc[elig, "e_random"].to_numpy()
    lin = _lineage_of(labels.index, lineage).to_numpy()
    Y = labels.loc[:, elig].to_numpy(dtype=float)
    st = _SortedColumns(scores.loc[labels.index, elig].to_numpy(float), ~np.isnan(Y))
    h = _harmonic_table(len(lin) + 1)
    ones = np.ones((1, len(lin)))
    null = np.empty(n_perm)
    for i, child in enumerate(_seed_sequence(seed).spawn(n_perm)):
        Yp = permute_within_lineage(Y, lin, np.random.default_rng(child))
        t, _, _ = st.tap(ones, st.sorted_labels(Yp), h)
        null[i] = np.mean(t[0] - e)
    thr = float(np.percentile(null, percentile))
    tol = 1e-12  # floating-point ties (e.g. a line-constant scorer: observed = every null draw = 0)
    return {
        "observed": obs,
        "threshold": thr,
        "passed": bool(obs > thr + tol),
        "p": float((1 + np.sum(null >= obs - tol)) / (1 + n_perm)),
        "n_genes_eligible": int(len(elig)),
        "null": null,
    }


def c3_decision(
    scores_p1: pd.DataFrame,
    scores_b3: pd.DataFrame,
    labels: pd.DataFrame,
    lineage: pd.Series,
    *,
    seed,
    B: int = B_CONFIRMATORY,
    n_perm: int = N_C3_PERMUTATIONS,
    conf: float = 0.95,
) -> dict:
    """C3 decision (§6.3): mΔAP_gene(P1) - mΔAP_gene(B3), all genes.

    Passes iff (a) the line-bootstrap CI lower bound of the paired P1 - B3 contrast is > 0 **and**
    (b) mΔAP_gene(P1) exceeds the 97.5th percentile of its within-lineage label-permutation null.
    P1, B3 and labels must be on identical gene x line sets. ``seed`` is spawned into
    [bootstrap, permutation] streams. The result feeds ``gatekept_claims`` (C3 has no SESOI).
    """
    assert_same_pairs(labels, scores_p1, scores_b3)
    for name, s in (("P1", scores_p1), ("B3", scores_b3)):
        _align(s, labels)  # raises PairingError on a missing score at a scored pair
    ss_boot, ss_perm = _seed_sequence(seed).spawn(2)
    boot = gene_bootstrap(scores_p1, scores_b3, labels, lineage, B=B, seed=ss_boot, conf=conf)
    null = c3_permutation_null(scores_p1, labels, lineage, n_perm=n_perm, seed=ss_perm)
    boot_ok = bool(boot["ci_low"] > 0)
    return {
        "contrast": "C3",
        "estimate": boot["estimate"],
        "ci_low": boot["ci_low"],
        "ci_high": boot["ci_high"],
        "bootstrap_passed": boot_ok,
        "m_delta_ap_gene_p1": null["observed"],
        "null_threshold": null["threshold"],
        "permutation_p": null["p"],
        "permutation_passed": null["passed"],
        "passed": bool(boot_ok and null["passed"]),
        "n_genes_eligible": boot["n_genes_eligible"],
        "mean_genes_dropped": boot["mean_genes_dropped"],
    }


def c3_claim_condition(c3_all: dict, c3_expressed: dict) -> bool:
    """§6.3: 'P1 captures line-specific dependency beyond expression' needs C3 to pass and the
    expressed-pairs version to have the same sign with CI lower bound > 0."""
    return bool(
        c3_all["passed"]
        and np.sign(c3_expressed["estimate"]) == np.sign(c3_all["estimate"])
        and c3_expressed["ci_low"] > 0
    )


# ----------------------------------------------------------------------------------------------------
# §6.8 binary metrics (Human1 / Human2 paper protocols) and §6.1 binary calls
# ----------------------------------------------------------------------------------------------------
def binary_metrics(tp, fp, tn, fn, zero_division=np.nan) -> dict:
    """MCC, F1, precision, recall (= sensitivity), specificity, accuracy from confusion counts (vectorised).

    A zero denominator gives NaN (repository-script convention, §6.8); ``zero_division=0`` gives the
    sensitivity row. F1 = 2TP / (2TP + FP + FN).
    """
    tp, fp, tn, fn = (np.asarray(x, dtype=float) for x in (tp, fp, tn, fn))

    def div(a, b):
        with np.errstate(invalid="ignore", divide="ignore"):
            out = np.where(b > 0, a / np.where(b > 0, b, 1), np.nan)
        if not (isinstance(zero_division, float) and np.isnan(zero_division)):
            out = np.where(b > 0, out, zero_division)
        return out

    res = {
        "mcc": div(tp * tn - fp * fn, np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))),
        "f1": div(2 * tp, 2 * tp + fp + fn),
        "precision": div(tp, tp + fp),
        "recall": div(tp, tp + fn),
        "specificity": div(tn, tn + fp),
        "accuracy": div(tp + tn, tp + tn + fp + fn),
    }
    return {k: (float(v) if v.ndim == 0 else v) for k, v in res.items()}


def confusion_counts(calls, truth) -> tuple[int, int, int, int]:
    """(TP, FP, TN, FN) of boolean ``calls`` against {0,1,NaN} ``truth``; NaN-truth pairs are excluded."""
    c = np.asarray(calls)
    t = np.asarray(truth, dtype=float)
    if c.shape != t.shape:
        raise PairingError("calls and truth differ in shape")
    if c.dtype != bool:
        cf = c.astype(float)
        if np.isnan(cf).any() or not np.isin(cf, (0.0, 1.0)).all():
            raise ValueError("calls must be boolean / 0-1 with no missing value")
        c = cf.astype(bool)
    m = ~np.isnan(t)
    c, t = c[m], t[m].astype(bool)
    return int((c & t).sum()), int((c & ~t).sum()), int((~c & ~t).sum()), int((~c & t).sum())


def binary_per_line(calls: pd.DataFrame, truth: pd.DataFrame, zero_division=np.nan) -> pd.DataFrame:
    """Per-line confusion matrix and §6.8 metrics. ``calls`` and ``truth`` are lines x genes on identical sets.

    The gene set (parent-model genes with a CRISPR column; parent genes absent from a line model called
    non-essential) is the caller's column selection; pairs with NaN truth are excluded per line.
    """
    assert_same_pairs(truth, calls)
    Cm = calls.loc[truth.index, truth.columns].to_numpy()
    if Cm.dtype != bool:
        Cf = Cm.astype(float)
        if np.isnan(Cf).any():
            raise ValueError("calls contain missing values")
        Cm = Cf.astype(bool)
    T = truth.to_numpy(dtype=float)
    m = ~np.isnan(T)
    Tb = np.where(m, T, 0.0).astype(bool)
    tp = (Cm & Tb & m).sum(1)
    fp = (Cm & ~Tb & m).sum(1)
    tn = (~Cm & ~Tb & m).sum(1)
    fn = (~Cm & Tb & m).sum(1)
    out = pd.DataFrame({"tp": tp, "fp": fp, "tn": tn, "fn": fn}, index=truth.index)
    for k, v in binary_metrics(tp, fp, tn, fn, zero_division).items():
        out[k] = v
    return out


def labels_from_effect(effect: pd.DataFrame, threshold: float = L_GE06_THRESHOLD) -> pd.DataFrame:
    """Binary label = gene effect < threshold (strict); NaN effect stays NaN. E.g. L-ge06 at -0.6."""
    lab = (effect < threshold).astype(float)
    return lab.where(effect.notna())


def binary_threshold_grid(
    calls: pd.DataFrame, effect: pd.DataFrame, thresholds=THRESHOLD_GRID, zero_division=np.nan
) -> pd.DataFrame:
    """Per-line binary metrics at every label threshold of the Human1-paper grid (-1.0 ... 0.0, 21 steps)."""
    parts = []
    for thr in thresholds:
        df = binary_per_line(calls, labels_from_effect(effect, float(thr)), zero_division)
        df.insert(0, "threshold", float(thr))
        parts.append(df)
    return pd.concat(parts)


def summarise_binary(per_line_table: pd.DataFrame, lineage: pd.Series, metrics=BINARY_METRICS) -> pd.DataFrame:
    """§6.8 summaries per metric: lineage-equal mean, line-weighted mean, per-line median, NaN count.

    Two rows per metric: nan_policy='exclude' (primary: NaN lines dropped and counted) and 'zero'
    (sensitivity row: NaN set to 0).
    """
    rows = []
    for metric in metrics:
        v = per_line_table[metric].astype(float)
        n_nan = int(v.isna().sum())
        for policy, vv in (("exclude", v), ("zero", v.fillna(0.0))):
            rows.append(
                {
                    "metric": metric,
                    "nan_policy": policy,
                    "lineage_equal_mean": lineage_equal_mean(vv, lineage),
                    "line_mean": line_weighted_mean(vv),
                    "median": float(vv.median()) if vv.notna().any() else float("nan"),
                    "n_lines_used": int(vv.notna().sum()),
                    "n_nan": n_nan,
                }
            )
    return pd.DataFrame(rows)


def fba_calls(ratio_r, cutoff: float = FBA_CUTOFF):
    """§6.1 fixed-cutoff binary call for FBA-type methods: essential iff ratio_r < 0.01."""
    if isinstance(ratio_r, (pd.DataFrame, pd.Series)):
        return ratio_r < cutoff
    return np.asarray(ratio_r, dtype=float) < cutoff


def rate_threshold(fit_scores, rate: float) -> float:
    """Score threshold for dichotomising a non-FBA scorer at its fitting-line positive rate (§6.1, §6.8).

    Returns t = the ceil(rate * n)-th largest fitting score; call 'essential' iff score >= t
    (ties at t make the realised rate >= ``rate``).
    """
    s = np.sort(np.asarray(fit_scores, dtype=float).ravel())[::-1]
    s = s[np.isfinite(s)]
    if not 0 < rate <= 1 or s.size == 0:
        raise ValueError("need 0 < rate <= 1 and finite fitting scores")
    k = max(1, math.ceil(round(rate * s.size, 9)))
    return float(s[k - 1])


def calls_at_threshold(scores, threshold: float):
    """Binary calls score >= threshold (pair with ``rate_threshold``)."""
    if isinstance(scores, (pd.DataFrame, pd.Series)):
        return scores >= threshold
    return np.asarray(scores, dtype=float) >= threshold
