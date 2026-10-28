"""F2 (A15 post hoc; RESEARCH_TREE D18) shared helpers. Label-free.

Written by the F2 engineer as a non-author of M6, W7 and F1: score building, tie-aware AP, the lineage-stratified
bootstrap and the fold rule are re-implemented here from DESIGN §6.1, §6.2, §6.5 and §9 (TG) and are not
imported from nextgem.metrics, nextgem.w1, nextgem.w2_score or the F1/M6/W7 scripts. Only path and split helpers
(nextgem.paths, nextgem.labels.load_split / gene_universe, nextgem.w1.line_features) are shared.

Methods (column names):
  H1    generic Human1 on Ham's (G*; line-constant)          H2    generic Human2 (line-constant)
  P1    "M6:H1+Bloo/M2(f=1)"  (D11)                          F2    "F2:H1+artBloo/M2(f=1)" (D18; W7-03 sweeps)
  F1    curated P1 (D17), training recompute only: P1 with the F1 arm's ratio deltas (results/f1/sweeps).
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
RES = PROJECT / "results"
F2DIR = RES / "f2"
TRDIR = F2DIR / "train"
M10VAL = RES / "m10" / "val"
MANIFEST = F2DIR / "MANIFEST.tsv"
AUTHOR = "F2 engineer (Claude Opus; non-author of M6/W7/F1)"

PARENT = "human1"
M2 = "M2(f=1)"
P1_VARIANT = "B-loo:cofactor_pool_biomass"
F2_VARIANT = "W7:artefact_bloo"
P1_METHOD = "M6:H1+Bloo/M2(f=1)"
F2_METHOD = "F2:H1+artBloo/M2(f=1)"
F1_METHOD = "F1:curatedP1"
MEV_CORE = ("HMGCR", "MVK", "PMVK", "MVD")
LETHAL = 0.01  # ratio_r < 0.01 (§6.1 fixed cutoff)

W7_SWEEPS = RES / "w7" / "sweeps" / "ratio_r_human1.parquet"
W7_SWEEPS_TSV = RES / "w7" / "sweeps" / "sweeps.tsv"
W7_SWEEPS_JSON = RES / "w7" / "sweeps" / "sweeps.json"
W7_MEV = RES / "w7" / "sweeps" / "mev.tsv"
M6_SWEEPS = RES / "m6" / "sweeps" / "ratio_r_human1.parquet"
M6_SWEEP_E = RES / "m6" / "scores" / "sweep_e.parquet"
M6_ASSIGN = RES / "m6" / "scores" / "line_assignment.tsv"
M6_TRAIN_SCORES = RES / "m6" / "scores" / "scores_train.large.parquet"
M6_VAL_SCORES = RES / "m6" / "val" / "m6_val_scores.parquet"
M6_PER_LINE = RES / "m6" / "train" / "per_line.tsv"
M6_MEV = RES / "m6" / "scores" / "w3_mev.tsv"
W2_LINE_SWEEPS = RES / "w2" / "controls" / "line_sweeps.tsv"
F1_DELTA = RES / "f1" / "sweeps" / "ratio_r_delta_human1.parquet"

# Numbers to recompute (task brief; RESEARCH_TREE D18 and "F1 result")
REPORTED = {
    "F1_vs_P1": {"delta_lin": -0.00063, "ci": [-0.00071, -0.00055], "folds_positive": 0},
    "F2_vs_Gstar": {"delta_lin": 0.0304, "ci": [0.0294, 0.0313], "folds_positive": 12},
    "F2_vs_P1": {"delta_lin": -0.0230, "ci": [-0.0236, -0.0224], "folds_positive": 0},
    "P1_vs_Gstar": {"delta_lin": 0.0534, "ci": [0.0523, 0.0545], "folds_positive": 12},
}


# ----------------------------------------------------------------------------- io
def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_parquet(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_parquet(tmp, index=False, compression="zstd")
    os.replace(tmp, path)
    return path


def write_tsv(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
    os.replace(tmp, path)
    return path


def write_json(obj, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=float) + "\n")
    os.replace(tmp, path)
    return path


def provenance(*, inputs=(), seeds=(), solver="none", workers=1, allow_dirty=False):
    from nextgem.run import provenance as _prov

    F2DIR.mkdir(parents=True, exist_ok=True)
    return _prov(MANIFEST, inputs=inputs, seeds=seeds, solver=solver, workers=workers, allow_dirty=allow_dirty,
                 author=AUTHOR)


# ----------------------------------------------------------------------------- scores (label-free)
def universe() -> list[int]:
    from nextgem import labels

    return [int(x) for x in labels.gene_universe().entrez]


def gene_map() -> pd.Series:
    """Human1 model gene (ENSG) -> Entrez from the baseline KO table: NaN Entrez dropped; one gene per Entrez
    (the first in file order; the table has no duplicate Entrez, asserted)."""
    from nextgem import paths

    ko = pd.read_parquet(paths.baseline_dir() / f"{PARENT}_ko.parquet", columns=["gene", "entrez"])
    ko = ko.dropna(subset=["entrez"])
    ko = ko.assign(entrez=ko.entrez.astype(np.int64))
    if ko.entrez.duplicated().any() or ko.gene.duplicated().any():
        raise AssertionError("baseline KO table: duplicate gene or Entrez")
    return ko.set_index("gene").entrez


def symbol_map() -> pd.Series:
    from nextgem import paths

    ko = pd.read_parquet(paths.baseline_dir() / f"{PARENT}_ko.parquet", columns=["gene", "symbol"])
    return ko.set_index("gene").symbol


def e_rows(ratio: pd.DataFrame, gmap: pd.Series, U: list[int]) -> pd.DataFrame:
    """rows x U: e = 1 − ratio_r rounded to 6 decimals; a U gene absent from the model (or with no KO value)
    scores 0 (DESIGN §6.1, W1 convention)."""
    cols = [g for g in ratio.columns if g in gmap.index]
    r = ratio[cols].astype(float)
    r.columns = gmap.loc[cols].to_numpy()
    e = (1.0 - r).reindex(columns=U)
    return e.fillna(0.0).round(6)


def generic_e(U: list[int], parent: str) -> pd.Series:
    """Generic parent on Ham's: e per U gene from the baseline KO table (line-constant)."""
    from nextgem import paths

    ko = pd.read_parquet(paths.baseline_dir() / f"{parent}_ko.parquet", columns=["entrez", "ratio_r"])
    ko = ko.dropna(subset=["entrez"]).assign(entrez=lambda d: d.entrez.astype(np.int64))
    ko = ko.drop_duplicates("entrez").set_index("entrez").ratio_r
    return (1.0 - ko.reindex(U)).fillna(0.0).round(6)


def line_m2_sweep() -> pd.Series:
    """ModelID -> Human1 W2 sweep id of the line's M2(f=1) medium (results/m6/scores/line_assignment.tsv), asserted
    equal to the W2 gate-free assignment (results/w2/controls/line_sweeps.tsv)."""
    la = pd.read_csv(M6_ASSIGN, sep="\t")
    a = la[(la.parent == PARENT) & (la.step == M2)].set_index("ModelID").w2_sweep_id
    ls = pd.read_csv(W2_LINE_SWEEPS, sep="\t")
    b = ls[(ls.parent == PARENT) & (ls.ladder_step == M2)].set_index("ModelID").sweep_id
    if a.index.duplicated().any() or not a.sort_index().equals(b.reindex(a.index).sort_index()):
        raise AssertionError("line_assignment.tsv and W2 line_sweeps.tsv disagree on the Human1 M2(f=1) sweep")
    return a


def f2_sweep_e(U: list[int], gmap: pd.Series | None = None) -> pd.DataFrame:
    """W2 sweep id -> e row (U) for F2 = W7:artefact_bloo on the M2(f=1) signatures (W7-03 sweeps)."""
    gmap = gene_map() if gmap is None else gmap
    R = pd.read_parquet(W7_SWEEPS).set_index("sweep")
    if R.isna().all(axis=1).any():
        raise AssertionError("W7 artefact sweeps: a failed (all-NaN) sweep")
    E = e_rows(R, gmap, U)
    var, sid = zip(*(k.split("|", 1) for k in E.index))
    if set(var) != {F2_VARIANT}:
        raise AssertionError(f"unexpected variants {set(var)}")
    E.index = list(sid)
    return E


def p1_sweep_e(U: list[int]) -> pd.DataFrame:
    """W2 sweep id -> P1 e row (U), from M6's sweep-keyed e table (results/m6/scores/sweep_e.parquet)."""
    se = pd.read_parquet(M6_SWEEP_E)
    se = se[(se.parent == PARENT) & (se.variant == P1_VARIANT)]
    E = se.pivot(index="w2_sweep_id", columns="entrez", values="e").reindex(columns=U)
    if E.isna().any().any():
        raise AssertionError("P1 sweep_e: missing entries")
    return E


def line_matrix(sweep_e: pd.DataFrame, sweep_of: pd.Series, lines) -> pd.DataFrame:
    lines = list(lines)
    sid = sweep_of.reindex(lines)
    if sid.isna().any() or not set(sid) <= set(sweep_e.index):
        raise AssertionError("a line has no M2(f=1) sweep in this cell")
    M = sweep_e.loc[sid.to_numpy()]
    M.index = lines
    return M


def f1_sweep_e(U: list[int], p1e: pd.DataFrame, gmap: pd.Series | None = None) -> tuple[pd.DataFrame, dict]:
    """W2 sweep id -> F1 (curated P1) e row: P1's row with the F1 arm's changed ratios (results/f1/sweeps). The
    stored P1 ratio of every changed entry must reproduce P1's e (checks the composition and the gene map)."""
    gmap = gene_map() if gmap is None else gmap
    d = pd.read_parquet(F1_DELTA)
    d = d[d.arm == "F1"]
    d = d[d.gene.isin(gmap.index)].assign(entrez=lambda x: gmap.reindex(x.gene).to_numpy())
    d = d[d.entrez.isin(U)]
    E = p1e.copy()
    ri, ci = E.index.get_indexer(d.w2_sweep_id), E.columns.get_indexer(d.entrez)
    if (ri < 0).any() or (ci < 0).any():
        raise AssertionError("F1 delta: sweep id or gene not among P1 sweeps / U")
    p1_from_delta = np.round(1.0 - d.ratio_r_P1.to_numpy(float), 6)
    p1_stored = E.to_numpy()[ri, ci]
    max_dev = float(np.max(np.abs(p1_from_delta - p1_stored))) if len(d) else 0.0
    vals = np.round(1.0 - d.ratio_r.to_numpy(float), 6)
    for (sid, ent), v in zip(zip(d.w2_sweep_id, d.entrez), vals):
        E.at[sid, ent] = v
    info = {"n_delta_entries": int(len(d)), "n_genes": int(d.entrez.nunique()),
            "max_dev_P1_ratio_vs_sweep_e": max_dev}
    if max_dev > 1e-6:
        raise AssertionError(f"F1 delta file's P1 ratios disagree with P1 sweep_e (max {max_dev})")
    return E, info


# ----------------------------------------------------------------------------- metrics (own implementation)
def tap(scores: np.ndarray, y: np.ndarray) -> float:
    """Tie-aware AP (DESIGN §6.1; McSherry & Najork 2008), O(n) block form: a tied block of size n_g with p_g
    positives, preceded by N items of which P positive, contributes sum_r (p_g/n_g)(P + 1 + (r−1)(p_g−1)/(n_g−1))
    / (N + r), the fraction taken as 0 when n_g = 1; divided by the total number of positives."""
    s = np.asarray(scores, dtype=float)
    yy = np.asarray(y, dtype=bool)
    npos = int(yy.sum())
    if npos == 0:
        return float("nan")
    order = np.argsort(-s, kind="stable")
    s, yy = s[order], yy[order]
    n = len(s)
    new = np.empty(n, dtype=bool)
    new[0] = True
    new[1:] = s[1:] != s[:-1]
    bid = np.cumsum(new) - 1
    starts = np.flatnonzero(new)
    ng = np.diff(np.append(starts, n))
    pg = np.add.reduceat(yy.astype(float), starts)
    P = np.concatenate([[0.0], np.cumsum(pg)[:-1]])
    N = starts.astype(float)
    r = np.arange(n) - starts[bid] + 1.0
    ngi, pgi = ng[bid].astype(float), pg[bid]
    frac = np.where(ngi > 1, (r - 1.0) * (pgi - 1.0) / np.maximum(ngi - 1.0, 1.0), 0.0)
    terms = (pgi / ngi) * (P[bid] + 1.0 + frac) / (N[bid] + r)
    return float(terms.sum() / npos)


def per_line_tap(S: pd.DataFrame, Y: pd.DataFrame, min_pos: int = 5) -> pd.DataFrame:
    """Per line: n, n_pos, estimable (>= min_pos positives, §5.2) and tAP over the line's labelled U genes."""
    S = S.reindex(index=Y.index, columns=Y.columns)
    if S.isna().any().any():
        raise AssertionError("scores missing for a labelled line or gene")
    rows = []
    Sv, Yv = S.to_numpy(float), Y.to_numpy(float)
    for i, l in enumerate(Y.index):
        m = ~np.isnan(Yv[i])
        yy = Yv[i, m] > 0.5
        npos = int(yy.sum())
        rows.append({"ModelID": l, "n": int(m.sum()), "n_pos": npos, "estimable": npos >= min_pos,
                     "tap": tap(Sv[i, m], yy) if npos >= min_pos else np.nan})
    return pd.DataFrame(rows).set_index("ModelID")


def lineage_equal_mean(v: pd.Series, lineage: pd.Series) -> float:
    v = v.dropna()
    return float(v.groupby(lineage.reindex(v.index)).mean().mean())


def boot_lineage(v: pd.Series, lineage: pd.Series, B: int, seed: int) -> np.ndarray:
    """Lineage-stratified line bootstrap (§6.5): lines resampled with replacement within each lineage; each
    replicate is the lineage-equal mean. Lineages are drawn in sorted order from one generator."""
    v = v.dropna()
    lin = lineage.reindex(v.index)
    rng = np.random.default_rng(seed)
    reps = np.zeros(B)
    groups = sorted(lin.unique())
    for g in groups:
        x = v[lin == g].to_numpy(float)
        idx = rng.integers(0, len(x), size=(B, len(x)))
        reps += x[idx].mean(axis=1)
    return reps / len(groups)


def contrast(a: pd.DataFrame, b: pd.DataFrame, lineage: pd.Series, folds, B: int, seed: int,
             min_lines: int = 5) -> dict:
    """Paired per-line ΔtAP (a − b) on lines estimable in both; Δ̄_lin, percentile CI, line-weighted mean, and
    folds positive (§9 TG: a fold counts iff it has >= min_lines lines and a mean Δ > 0)."""
    ok = a.estimable & b.estimable.reindex(a.index).fillna(False)
    if not (a.index.equals(b.index) and (a.estimable == b.estimable).all()):
        raise AssertionError("methods are not paired on identical lines")
    d = (a.tap - b.tap)[ok]
    reps = boot_lineage(d, lineage, B, seed)
    lin = lineage.reindex(d.index)
    fm = d.groupby(lin).agg(["mean", "size"])
    pos = sum(1 for f in folds if f in fm.index and fm.loc[f, "size"] >= min_lines and fm.loc[f, "mean"] > 0)
    return {"delta_lin": lineage_equal_mean(d, lineage), "ci": [float(np.percentile(reps, 2.5)),
            float(np.percentile(reps, 97.5))], "se_boot": float(np.std(reps, ddof=1)),
            "delta_line": float(d.mean()), "folds_positive": int(pos), "n_folds": len(list(folds)),
            "n_lines": int(len(d)), "n_lineages": int(lin.nunique())}
