"""Stage-A analysis library (DESIGN.md §9 W4, §13; amendments A4, A7, A10). Pure functions; no label file is opened
here: label frames are passed in by the calling script (``labels.train()`` / ``labels.fold()`` only).

Contents
--------
* label-blind gates (growth, viability tasks, repair size <= 5% of reactions), per build and per configuration;
* reaction / task / KO profiles from the Stage-A outputs and Jaccard distances;
* the method-vs-tissue variance partition (PERMANOVA-type marginal R^2 on sqrt-Jaccard distances, permutation
  p-values, lineage-stratified line bootstrap for H4.1);
* score frames (e = 1 - ratio_r over the gene universe U, Entrez columns, W1 convention);
* the panel training gate (TG: Δ̄_lin >= 0.005, paired line-bootstrap CI lower bound > 0, >= 8 of 10 panel folds
  positive, a fold "not positive" if < 2 of its lines are estimable) and nested panel LOLO selection.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from nextgem import metrics
from nextgem.config import SEEDS

REPAIR_MAX_FRAC = 0.05  # §9 W4 label-blind gate
CONFIG_MIN_PASS_FRAC = 0.95  # a configuration passes a gate if >= 95% of its builds pass (cf. §11.4 5% flag)
TG_MIN_DELTA = 0.005  # §9 TG item 1
TG_PANEL_MIN_FOLDS = 8  # §9: >= 8 of 10 panel folds
TG_PANEL_MIN_EST = 2  # a panel fold is "not positive" if < 2 of its 3 lines are estimable (r3 N-1)
RULE_ORDER = {"R1": 0, "R2": 1, "R3": 2, "lcp": 3}
GENOMIC_ORDER = {"none": 0, "loss": 1}


# ============================================================================= gates
def viability_task_indices(tasks) -> list[int]:
    """Positions of the viability tasks (GR + the 9 ER tasks) in the 57-task profile string."""
    return [i for i, t in enumerate(tasks) if t.id in ("GR", "ER")]


def build_gates(builds: pd.DataFrame, viability_idx: list[int], repair_max: float = REPAIR_MAX_FRAC) -> pd.DataFrame:
    """Per build x genomic variant: gate_growth, gate_viability, gate_repair, gate_all (label-blind, §9 W4)."""
    b = builds.copy()
    failed = b["failed"].fillna(False).astype(bool)
    b["gate_growth"] = ~failed & (b["wt_growth"].fillna(0) > 0)
    prof = b["task_profile"].fillna("")
    b["viability_passed"] = prof.map(lambda s: sum(1 for i in viability_idx if len(s) > i and s[i] == "1"))
    b["gate_viability"] = b["viability_passed"] == len(viability_idx)
    b["gate_repair"] = b["repair_frac_reactions"].fillna(0.0) <= repair_max
    b["gate_all"] = b.gate_growth & b.gate_viability & b.gate_repair
    return b


def config_gates(gated: pd.DataFrame, keys: list[str], min_frac: float = CONFIG_MIN_PASS_FRAC) -> pd.DataFrame:
    g = gated.groupby(keys, dropna=False)
    out = g.agg(n_builds=("gate_all", "size"), frac_growth=("gate_growth", "mean"),
                frac_viability=("gate_viability", "mean"), frac_repair=("gate_repair", "mean"),
                frac_all=("gate_all", "mean")).reset_index()
    out["config_gate_pass"] = out.frac_all >= min_frac
    return out


# ============================================================================= profiles and distances
def unpack_mask(hexmask: str, n: int) -> np.ndarray:
    return np.unpackbits(np.frombuffer(bytes.fromhex(hexmask), dtype=np.uint8))[:n].astype(bool)


def jaccard_distance_sq(X: np.ndarray) -> np.ndarray:
    """Squared sqrt-Jaccard distance (= Jaccard dissimilarity 1 - J) between the rows of a boolean matrix.

    sqrt(1 - J) is Euclidean-embeddable (Gower & Legendre 1986), so the PERMANOVA sums of squares below
    are non-negative and additive. Two empty rows have distance 0.
    """
    Xf = X.astype(np.float64)
    inter = Xf @ Xf.T
    s = Xf.sum(1)
    union = s[:, None] + s[None, :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        J = np.where(union > 0, inter / union, 1.0)
    return 1.0 - J


def ss_within(D2: np.ndarray, groups: np.ndarray) -> float:
    tot = 0.0
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        if idx.size:
            tot += D2[np.ix_(idx, idx)].sum() / (2.0 * idx.size)
    return tot


def r2(D2: np.ndarray, groups) -> float:
    """Marginal PERMANOVA R^2 of a grouping: SS_between / SS_total (Anderson 2001)."""
    groups = np.asarray(groups)
    sst = D2.sum() / (2.0 * D2.shape[0])
    return float((sst - ss_within(D2, groups)) / sst) if sst > 0 else float("nan")


def permanova(D2: np.ndarray, meta: pd.DataFrame, factor: str, *, perms: int = 999, seed=SEEDS[0],
              block: str | None = None, unit: str | None = None) -> dict:
    """Marginal R^2 of ``factor`` with a permutation p-value (pseudo-F is monotone in R^2 for a fixed design).

    ``unit``: permute the factor at the level of this column (e.g. lineage labels permuted across lines, each line
    keeping all its builds). ``block``: permute within blocks (e.g. method labels within each line).
    """
    rng = np.random.default_rng(seed)
    obs = r2(D2, meta[factor].to_numpy())
    null = np.empty(perms)
    for p in range(perms):
        lab = meta[factor].to_numpy().copy()
        if unit is not None:
            units = meta[unit].to_numpy()
            u_vals = pd.unique(units)
            u_lab = pd.Series(meta.groupby(unit, sort=False)[factor].first().reindex(u_vals).to_numpy())
            perm = dict(zip(u_vals, rng.permutation(u_lab.to_numpy())))
            lab = np.array([perm[u] for u in units])
        elif block is not None:
            for _, idx in meta.groupby(block).indices.items():
                lab[idx] = rng.permutation(lab[idx])
        else:
            lab = rng.permutation(lab)
        null[p] = r2(D2, lab)
    return {"factor": factor, "r2": obs, "p_perm": float((1 + (null >= obs).sum()) / (1 + perms)),
            "null_q95": float(np.quantile(null, 0.95)), "null_mean": float(null.mean()), "perms": perms}


def bootstrap_r2_diff(D2: np.ndarray, meta: pd.DataFrame, f_a: str, f_b: str, *, line: str = "line",
                      lineage: str = "lineage", B: int = metrics.B_OTHER, seed=SEEDS[0]) -> dict:
    """H4.1: R^2(f_a) - R^2(f_b) with a lineage-stratified line bootstrap (lines resampled within lineage,
    each resampled line bringing all its builds)."""
    lines = meta[[line, lineage]].drop_duplicates().reset_index(drop=True)
    idx_of = meta.groupby(line).indices
    rng = np.random.default_rng(seed)
    obs = r2(D2, meta[f_a].to_numpy()) - r2(D2, meta[f_b].to_numpy())
    reps = np.empty(B)
    by_lin = lines.groupby(lineage)[line].apply(list).to_dict()
    for b in range(B):
        pick = []
        for lin_, ls in by_lin.items():
            for li in rng.choice(len(ls), size=len(ls), replace=True):
                pick.append(ls[li])
        idx = np.concatenate([idx_of[x] for x in pick])
        sub = D2[np.ix_(idx, idx)]
        m = meta.iloc[idx]
        # resampled duplicates of a line are distinct units for the lineage/line factors
        reps[b] = r2(sub, m[f_a].to_numpy()) - r2(sub, m[f_b].to_numpy())
    lo, hi = np.quantile(reps, [0.025, 0.975])
    return {"diff": float(obs), "lo": float(lo), "hi": float(hi), "B": B}


# ============================================================================= scores
def entrez_gene_map(baseline_ko: pd.DataFrame) -> pd.Series:
    """Entrez -> Human2 Ensembl gene, exactly as W1's generic table (drop_duplicates on Entrez, gene order)."""
    ko = baseline_ko[["gene", "entrez"]].dropna(subset=["entrez"]).sort_values("gene")
    ko = ko.assign(entrez=ko.entrez.astype(int)).drop_duplicates("entrez")
    return ko.set_index("entrez").gene


def score_frame(ratio: pd.DataFrame, genes_U, emap: pd.Series) -> pd.DataFrame:
    """lines x U (Entrez) score e = 1 - ratio_r; U genes without a parent gene score 0 (§3.3)."""
    cols = list(genes_U)
    g = emap.reindex(cols)
    out = np.zeros((ratio.shape[0], len(cols)))
    have = g.notna().to_numpy()
    R = ratio.reindex(columns=g[have].to_numpy()).to_numpy(dtype=float)
    out[:, have] = np.round(1.0 - np.nan_to_num(R, nan=1.0), 6)
    return pd.DataFrame(out, index=ratio.index, columns=cols)


# ============================================================================= training gate and nested selection
def tg_panel(delta: pd.Series, lineage: pd.Series, *, B: int = metrics.B_OTHER, seed=SEEDS[0],
             min_delta: float = TG_MIN_DELTA, min_folds: int = TG_PANEL_MIN_FOLDS,
             min_est: int = TG_PANEL_MIN_EST, folds=None) -> dict:
    """Panel TG on per-line Δ (NaN = line not estimable). Folds = panel lineages."""
    lin = lineage.reindex(delta.index)
    folds = sorted(lin.unique()) if folds is None else list(folds)
    d = delta.dropna()
    out = {"n_lines": int(delta.size), "n_lines_estimable": int(d.size), "n_folds": len(folds)}
    if d.empty:
        return {**out, "delta_lin": np.nan, "lo": np.nan, "hi": np.nan, "folds_positive": 0, "tg_pass": False}
    bt = metrics.bootstrap_mean(d, lin.reindex(d.index), B=B, seed=seed)
    pos = 0
    fold_means = {}
    for f in folds:
        v = d[lin.reindex(d.index) == f]
        fold_means[f] = float(v.mean()) if len(v) else np.nan
        if len(v) >= min_est and v.mean() > 0:
            pos += 1
    out.update({"delta_lin": bt.estimate, "lo": bt.ci_low, "hi": bt.ci_high, "delta_line": float(d.mean()),
                "folds_positive": pos, "fold_means": fold_means})
    out["tg_pass"] = bool(bt.estimate >= min_delta and bt.ci_low > 0 and pos >= min_folds)
    return out


def nested_lolo_select(D: pd.DataFrame, lineage: pd.Series, order: list) -> tuple[pd.Series, dict]:
    """Nested panel LOLO. ``D``: lines x configurations of per-line Δ (or tAP; NaN = not estimable).

    For every panel lineage k, the configuration with the highest lineage-equal mean over the OTHER lineages is
    selected (ties: first in ``order``, the simplicity order), and its Δ on k's lines is the outer value. Labels of
    k never enter k's choice. Returns (outer per-line Δ, {k: selected configuration}).
    """
    lin = lineage.reindex(D.index)
    rank = {c: i for i, c in enumerate(order)}
    cols = sorted(D.columns, key=lambda c: rank.get(c, len(rank)))
    outer = pd.Series(np.nan, index=D.index)
    chosen = {}
    for k in sorted(lin.unique()):
        inner = lin != k
        best, best_v = None, -np.inf
        for c in cols:
            v = D.loc[inner, c]
            if v.notna().sum() == 0:
                continue
            m = metrics.lineage_equal_mean(v.dropna(), lin[inner].reindex(v.dropna().index))
            if m > best_v + 1e-15:
                best, best_v = c, m
        chosen[k] = best
        if best is not None:
            outer[lin == k] = D.loc[lin == k, best]
    return outer, chosen


def full_panel_rank(D: pd.DataFrame, lineage: pd.Series, order: list) -> pd.DataFrame:
    lin = lineage.reindex(D.index)
    rows = []
    for c in D.columns:
        v = D[c].dropna()
        rows.append({"config": c, "score_lin": metrics.lineage_equal_mean(v, lin.reindex(v.index)) if len(v) else np.nan,
                     "order": order.index(c) if c in order else math.inf})
    return pd.DataFrame(rows).sort_values(["score_lin", "order"], ascending=[False, True]).reset_index(drop=True)
