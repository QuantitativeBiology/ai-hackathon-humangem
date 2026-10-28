"""W7 auto-research loop library (DESIGN.md §W7, §8 B6, §9 TG, §17 A12a). Hybrids of gene priors and GEM outputs.

This module never opens a label file. Label frames reach it only as arguments: ``fit_long`` from
``labels.fold(k)`` (inside ``nextgem.cv.lolo``) or ``labels.train()``; held-out labels are used only inside the
harness's ``score`` callback (:func:`score_lines`), which returns per-line tAP tables.

Leakage rule enforced in code (:func:`assert_training_only`): every model is fitted on rows whose ModelID is a
training-tier line of the frozen split; predictions for other lines (validation look files) use their label-free
features only. Feature rows of non-fitting lines never enter a fit (tested by poisoning them in
``tests/test_w7.py``).

Feature whitelist (DESIGN §W7, PI decision 16): expression, CN, mutation, medium class / serum, GEM outputs, within-fold
B2 / B2m, neighbour and paralog features. Nothing label-derived from outside the fitting lines.

Scores are line x gene frames (columns = Entrez int of U), higher = more essential.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from nextgem import paths

LABEL = "L_prob"
LABEL_COLS = ["gene_effect", "L_prob"]  # fit_components needs gene_effect (B2c); only L_prob is a target
HGB_M2 = dict(max_iter=200, learning_rate=0.1, max_leaf_nodes=31, early_stopping=False)  # B6 as fixed at M2
P1_COL = "M6:H1+Bloo/M2(f=1)"
GEM_LINE_COLS = {  # feature name -> column of results/m6/scores/scores_<tier>.large.parquet
    "e_P1": P1_COL,
    "e_Gstar": "H1",
    "e_H1Bloo_M0": "M6:H1+Bloo/M0",
    "e_H1_M2": "M6:H1/M2(f=1)",
}
CALL = 0.99  # e >= 0.99 <=> ratio_r < 0.01 (§6.1 fixed cutoff)


def project_root() -> Path:
    return paths.PROJECT_ROOT


def w7_dir() -> Path:
    return paths.results_dir() / "w7"


# ----------------------------------------------------------------------------- ideas (model specs)
@dataclass(frozen=True)
class Idea:
    """A screened model. ``extra``: feature names added to the B6 features; ``kind``: 'clf' (HGB classifier,
    prob of L_prob) | 'resid' (HGB regressor on y - B2, score = B2 + prediction) | 'stack' (logit(B6) + beta*e_P1,
    beta chosen by nested inner LOLO from ``beta_grid``); ``mono``: monotone increasing constraint on b2 and b2m."""

    id: str
    extra: tuple = ()
    kind: str = "clf"
    hgb: tuple = tuple(sorted(HGB_M2.items()))
    mono: bool = False
    beta_grid: tuple = ()
    drop: tuple = ()  # B6 features removed
    mono_extra: tuple = ()  # extra features with a monotone increasing constraint

    def params(self) -> dict:
        return dict(self.hgb)


GEM_ALL = ("e_P1", "e_Gstar", "e_H1Bloo_M0", "e_H1_M2", "p1_lethal", "p1_effect", "p57_H1", "pgr_H1", "in_H1",
           "n_rxn_H1", "n_closed_H1")

W2MED = ("w2_medium", "w2_serum_pct", "w2_serum_free", "w2_documented")  # W2 tiered medium annotation
DEV = ("expr_dev", "cn_dev", "neigh_dev", "paralog_dev")  # line-relative omics (scripts/w7_common.deviation_features)

IDEAS: dict[str, Idea] = {
    "W7-01": Idea("W7-01"),  # H-b = B6 (ablation / reference)
    "W7-02": Idea("W7-02", extra=GEM_ALL),  # H-a
    "W7-04": Idea("W7-04", kind="stack", beta_grid=(0.0, 0.25, 0.5, 1.0, 2.0, 4.0)),  # H-d
    "W7-05": Idea("W7-05", mono=True),
    "W7-06": Idea("W7-06", kind="resid"),
    "W7-07": Idea("W7-07", kind="resid", extra=GEM_ALL),
    "W7-08": Idea("W7-08", extra=GEM_ALL + ("e_ART", "art_lethal")),
    "W7-09": Idea("W7-09", extra=DEV),
    "W7-10": Idea("W7-10", kind="resid", extra=DEV + GEM_ALL),
    "W7-11": Idea("W7-11", mono=True, extra=DEV),
    "W7-12": Idea("W7-12", mono=True, extra=GEM_ALL),
    "W7-13": Idea("W7-13", mono=True, extra=W2MED),
    "W7-14": Idea("W7-14", mono=True, extra=GEM_ALL,
                  mono_extra=("e_P1", "e_Gstar", "e_H1Bloo_M0", "e_H1_M2", "p1_lethal", "p57_H1")),
    "W7-15": Idea("W7-15", mono=True, extra=("knn_prior",), mono_extra=("knn_prior",)),
    "W7-16": Idea("W7-16", kind="knn_only"),
    "W7-17": Idea("W7-17", mono=True, extra=("knn_prior_k10", "knn_prior", "knn_prior_k200") + DEV,
                  mono_extra=("knn_prior_k10", "knn_prior", "knn_prior_k200")),
}


def register(idea: Idea) -> None:
    IDEAS[idea.id] = idea


# ----------------------------------------------------------------------------- label-free inputs
@dataclass
class Inputs:
    genes: list
    omics: dict  # expr, cn, mut_dam, mut_hot: lines x genes
    aggs: tuple  # (neighbour mean expression, max paralog expression): lines x genes
    gfeat: pd.DataFrame  # w1 gene features (n_isozyme_genes_h2 ...)
    lf: pd.DataFrame  # w1 line features (base_medium_class, serum)
    class_codes: dict
    gem_line: dict = field(default_factory=dict)  # name -> lines x genes (GEM outputs, label-free)
    gem_gene: pd.DataFrame | None = None  # genes x gene-level GEM features
    tier_of: pd.Series | None = None  # ModelID -> tier of the frozen split
    lineage_of: pd.Series | None = None

    @property
    def line_class(self) -> pd.Series:
        return self.lf.base_medium_class


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def m6_scores_file(tier: str) -> Path:
    return paths.results_dir() / "m6" / "scores" / f"scores_{tier}.large.parquet"


def check_m6_sha(tier: str) -> str:
    """The M6 per-tier score file must carry the SHA-256 pinned in results/m6/scores/SHA256SUMS.tsv."""
    f = m6_scores_file(tier)
    pins = pd.read_csv(f.parent / "SHA256SUMS.tsv", sep="\t")
    pin = pins.loc[pins.file.str.endswith(f.name), "sha256"]
    got = sha256(f)
    if len(pin) != 1 or pin.iloc[0] != got:
        raise RuntimeError(f"{f.name}: SHA-256 differs from SHA256SUMS.tsv")
    return got


def load_gem_line(tiers=("train",)) -> dict:
    """name -> lines x U (all lines of the given tiers) from the SHA-pinned M6 score files (label-free)."""
    frames = []
    for t in tiers:
        check_m6_sha(t)
        frames.append(pd.read_parquet(m6_scores_file(t), columns=["ModelID", "entrez", *GEM_LINE_COLS.values()]))
    df = pd.concat(frames, ignore_index=True)
    out = {}
    for name, col in GEM_LINE_COLS.items():
        w = df.pivot(index="ModelID", columns="entrez", values=col)
        w.columns = w.columns.astype(int)
        out[name] = w.astype("float32")
    p1 = out["e_P1"]
    out["p1_lethal"] = (p1 >= CALL).astype("float32")
    out["p1_effect"] = (p1 > 1e-6).astype("float32")
    return out


def gene_gem_features(genes, gen: pd.DataFrame, rxn: pd.DataFrame) -> pd.DataFrame:
    """genes x {p57_H1, pgr_H1, in_H1, n_rxn_H1, n_closed_H1} (label-free; ``gen`` = w1.generic_table,
    ``rxn`` = per-Entrez reaction counts of Human1)."""
    g = pd.DataFrame(index=pd.Index(list(genes), name="entrez"))
    g["p57_H1"] = gen["p57_H1"].reindex(g.index).fillna(False).astype("float32")
    g["pgr_H1"] = gen["pgr_H1"].reindex(g.index).fillna(False).astype("float32")
    g["in_H1"] = gen["in_H1"].reindex(g.index).fillna(False).astype("float32")
    g["n_rxn_H1"] = rxn["n_rxn"].reindex(g.index).fillna(0).astype("float32")
    g["n_closed_H1"] = rxn["n_closed"].reindex(g.index).fillna(0).astype("float32")
    return g


def human1_reaction_counts() -> pd.DataFrame:
    """Per Entrez (Human1): reactions with the gene in the GPR, and reactions closed by its single KO (label-free)."""
    from nextgem import w3
    from nextgem import w2_score as ws
    from nextgem.model import load_model

    m = load_model("human1")
    clo = w3.gene_closures(m)
    ko = pd.read_parquet(paths.baseline_dir() / "human1_ko.parquet", columns=["gene", "entrez", "ratio_r"])
    emap = ws.gene_to_entrez(ko)
    rows = [{"entrez": int(emap[g.id]), "n_rxn": len(g.reactions), "n_closed": len(clo.get(g.id, []))}
            for g in m.genes if g.id in emap.index]
    return pd.DataFrame(rows).drop_duplicates("entrez").set_index("entrez")


# ----------------------------------------------------------------------------- leakage guard
class LeakageError(AssertionError):
    pass


def assert_training_only(model_ids, tier_of: pd.Series) -> None:
    """Every fitting line must be a training-tier line of the frozen split."""
    ids = pd.Index(pd.unique(np.asarray(model_ids)))
    t = tier_of.reindex(ids)
    bad = ids[(t != "train").to_numpy()]
    if len(bad):
        raise LeakageError(f"{len(bad)} non-training line(s) in a fit (e.g. {list(bad[:2])})")


# ----------------------------------------------------------------------------- design matrices
def _w1_b6_module():
    spec = importlib.util.spec_from_file_location("w1_b6", project_root() / "scripts" / "w1_b6.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def feature_names(idea: Idea) -> list[str]:
    from nextgem import w1

    return [f for f in w1.B6_FEATURES if f not in idea.drop] + list(idea.extra)


def design(inp: Inputs, lines, b2_rows, b2m_rows, idea: Idea, local: dict | None = None) -> np.ndarray:
    """(n_lines * n_genes, n_features) float32, row-major over lines then genes; B6 features then ``idea.extra``.
    ``local``: fit-specific line x gene features (e.g. the within-fold kNN prior)."""
    from nextgem import w1

    lines, genes = list(lines), inp.genes
    X = w1.b6_matrix(lines, genes, b2_rows, b2m_rows, inp.omics, inp.aggs, inp.gfeat, inp.lf, inp.class_codes)
    keep = [i for i, f in enumerate(w1.B6_FEATURES) if f not in idea.drop]
    cols = [X[:, keep]]
    n, m = len(lines), len(genes)
    for f in idea.extra:
        if local is not None and f in local:
            M = local[f].reindex(index=lines, columns=genes).to_numpy(dtype="float32")
            cols.append(M.reshape(-1, 1))
        elif f in inp.gem_line:
            M = inp.gem_line[f].reindex(index=lines, columns=genes).to_numpy(dtype="float32")
            cols.append(M.reshape(-1, 1))
        elif inp.gem_gene is not None and f in inp.gem_gene.columns:
            v = inp.gem_gene[f].reindex(genes).to_numpy(dtype="float32")
            cols.append(np.tile(v, n).reshape(-1, 1))
        else:
            raise KeyError(f"feature {f!r} not available")
    out = np.concatenate(cols, axis=1)
    assert out.shape[0] == n * m
    return out


def _med_index(idea: Idea) -> int:
    return feature_names(idea).index("medium_class")


CATEGORICAL = ("medium_class", "w2_medium")  # categorical features (codes; -1 / NaN = missing)


def _clean_medium(X: np.ndarray, j: int, idea: Idea | None = None) -> np.ndarray:
    cols = [j] if idea is None else [i for i, f in enumerate(feature_names(idea)) if f in CATEGORICAL]
    for c in cols:
        X[:, c] = np.where(np.isnan(X[:, c]) | (X[:, c] < 0), np.nan, X[:, c])
    return X


KNN_K = 50
KNN_GENES = 2000
KNN_TIEBREAK = 1e-3
KNN_FEATURES = {"knn_prior": KNN_K, "knn_prior_k10": 10, "knn_prior_k200": 200}


def knn_space(inp: Inputs) -> pd.DataFrame:
    """Lines x 2,000 most variable U genes (variance over training lines) of log2(TPM+1), centred by the training-line
    gene means, then row-standardised (dot product / n = Pearson correlation). Label-free."""
    if getattr(inp, "_knn_Z", None) is not None:
        return inp._knn_Z
    X = inp.omics["expr"].reindex(columns=inp.genes).astype("float64")
    tr = [l for l in X.index if inp.tier_of is None or inp.tier_of.get(l) == "train"]
    v = X.loc[tr].var(axis=0).fillna(0.0)
    top = v.sort_values(ascending=False, kind="mergesort").index[:KNN_GENES]
    Z = X[top] - X.loc[tr, top].mean(axis=0)
    Z = Z.fillna(0.0)
    Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1, ddof=0).replace(0, 1.0), axis=0)
    inp._knn_Z = Z
    return Z


def knn_prior(inp: Inputs, Yw: pd.DataFrame, targets, exclude_lineage: pd.Series | None = None,
              fit_lineage: pd.Series | None = None, k: int = KNN_K) -> pd.DataFrame:
    """Within-fold kNN gene prior: for each target line, the mean label (NaN-ignoring) over its k most
    expression-correlated fitting lines (rows of ``Yw``, which hold fitting-line labels only). A target never uses
    itself; with ``exclude_lineage`` (fitting rows, stacked) it never uses lines of its own lineage. Targets without
    expression get NaN."""
    Z = knn_space(inp)
    F = [l for l in Yw.index if l in Z.index]
    T = [l for l in targets if l in Z.index]
    C = (Z.loc[T].to_numpy() @ Z.loc[F].to_numpy().T) / Z.shape[1]
    Fa = np.asarray(F)
    if exclude_lineage is not None:
        flin = fit_lineage.reindex(F).to_numpy()
        tlin = exclude_lineage.reindex(T).to_numpy()
        C = np.where(flin[None, :] == tlin[:, None], -np.inf, C)
    C = np.where(Fa[None, :] == np.asarray(T)[:, None], -np.inf, C)
    Yf = Yw.loc[F].to_numpy(dtype="float64")
    out = np.full((len(T), Yw.shape[1]), np.nan)
    kk = min(k, len(F))
    for i in range(len(T)):
        idx = np.argpartition(-C[i], kk - 1)[:kk]
        idx = idx[np.isfinite(C[i, idx])]
        if len(idx):
            out[i] = np.nanmean(Yf[idx], axis=0)
    res = pd.DataFrame(out, index=T, columns=Yw.columns).reindex(list(targets))
    return res


@dataclass
class Fit:
    idea: Idea
    est: object
    fitted: object  # w1.Fitted (B2/B2m of all fitting lines, used for held-out priors)
    n_rows: int
    fit_lines: list
    knn_Y: pd.DataFrame | None = None  # fitting-line labels (wide, in memory only) for the kNN prior


def fit(inp: Inputs, fit_long: pd.DataFrame, idea: Idea, seed: int, stacked=None) -> Fit:
    """Fit ``idea`` on ``fit_long`` (training labels minus held-out folds). B2/B2m features of fitting rows are
    stacked leave-own-lineage-out (as B6 at M2), so no row sees a prior containing its own label."""
    from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor

    from nextgem import w1

    if inp.tier_of is not None:
        assert_training_only(fit_long.ModelID, inp.tier_of)
    fitted = w1.fit_components(fit_long, inp.line_class, labels_=(LABEL,))
    if stacked is None:
        stacked = _w1_b6_module().stacked_priors(fit_long, inp.line_class, inp.genes)
    b2s, b2ms = stacked
    fl = list(b2s.index)
    Yw = w1.wide(fit_long, LABEL, fl, inp.genes)
    local, knn_Y = None, None
    if any(f in KNN_FEATURES for f in idea.extra) or idea.kind == "knn_only":
        knn_Y = Yw
        lin = fit_long.drop_duplicates("ModelID").set_index("ModelID").lineage
        local = {f: knn_prior(inp, knn_Y, fl, exclude_lineage=lin.reindex(fl), fit_lineage=lin, k=KNN_FEATURES[f])
                 for f in KNN_FEATURES if f in idea.extra or (f == "knn_prior" and idea.kind == "knn_only")}
    if idea.kind == "knn_only":
        return Fit(idea, None, fitted, 0, fl, knn_Y)
    X = design(inp, fl, b2s, b2ms, idea, local)
    Y = Yw.to_numpy().reshape(-1)
    ok = ~np.isnan(Y)
    X = _clean_medium(X[ok], _med_index(idea), idea)
    names = feature_names(idea)
    p = idea.params()
    kw = dict(random_state=seed, categorical_features=[i for i, f in enumerate(names) if f in CATEGORICAL])
    if idea.mono or idea.mono_extra:
        up = (("b2", "b2m") if idea.mono else ()) + tuple(idea.mono_extra)
        kw["monotonic_cst"] = [1 if f in up else 0 for f in names]
    if idea.kind in ("clf", "stack"):
        est = HistGradientBoostingClassifier(**p, **kw)
        est.fit(X, Y[ok].astype(int))
    elif idea.kind == "resid":
        b2col = names.index("b2")
        est = HistGradientBoostingRegressor(**p, **kw)
        est.fit(X, Y[ok] - np.nan_to_num(X[:, b2col]))
    else:
        raise ValueError(idea.kind)
    return Fit(idea, est, fitted, int(ok.sum()), fl, knn_Y)


def predict(inp: Inputs, f: Fit, lines) -> pd.DataFrame:
    """Scores for ``lines`` (held-out training lines, or validation lines for a look file): label-free features of
    those lines + priors from all fitting lines."""
    from nextgem import w1

    lines = list(lines)
    b2, b2m = w1.heldout_priors(f.fitted, lines, inp.line_class, inp.genes, label=LABEL)
    local = None
    if f.knn_Y is not None:
        local = {k: knn_prior(inp, f.knn_Y, lines, k=v) for k, v in KNN_FEATURES.items()
                 if k in f.idea.extra or (k == "knn_prior" and f.idea.kind == "knn_only")}
    if f.idea.kind == "knn_only":
        k = local["knn_prior"].reindex(index=lines, columns=inp.genes)
        s = (k.fillna(b2) + KNN_TIEBREAK * b2).to_numpy()
        return pd.DataFrame(np.round(s, 8), index=lines, columns=inp.genes)
    X = _clean_medium(design(inp, lines, b2, b2m, f.idea, local), _med_index(f.idea), f.idea)
    if f.idea.kind in ("clf", "stack"):
        s = f.est.predict_proba(X)[:, 1]
    else:
        s = np.nan_to_num(X[:, feature_names(f.idea).index("b2")]) + f.est.predict(X)
    return pd.DataFrame(np.round(s, 8).reshape(len(lines), len(inp.genes)), index=lines, columns=inp.genes)


def stack_scores(p_b6: pd.DataFrame, e_p1: pd.DataFrame, beta: float) -> pd.DataFrame:
    """H-d: logit(p_B6) + beta * e_P1 (p clipped to [1e-6, 1 - 1e-6])."""
    p = p_b6.clip(1e-6, 1 - 1e-6)
    e = e_p1.reindex(index=p.index, columns=p.columns)
    return np.log(p / (1 - p)) + beta * e


# ----------------------------------------------------------------------------- scoring (inside the harness)
def eval_lines(heldout_ids, inp: Inputs) -> list:
    """Held-out lines that are evaluated: those with a default expression profile (W1 primary set)."""
    expr = set(inp.omics["expr"].index)
    return [l for l in heldout_ids if l in expr]


def score_lines(preds: dict, held: pd.DataFrame, genes) -> dict:
    """``cv.lolo`` score callback: per-line tAP tables (no label values leave). preds: method -> lines x genes."""
    from nextgem import metrics, w1

    out = {}
    for name, S in preds.items():
        Y = w1.wide(held, LABEL, list(S.index), genes)
        out[name] = metrics.per_line(S.reindex(columns=genes), Y)
    return out


# ----------------------------------------------------------------------------- summaries (per-line tables in)
def tg(delta: pd.Series, lineage: pd.Series, folds, *, B: int, seed: int) -> dict:
    """§9 TG against a reference from a per-line Δ series: Δ̄_lin >= 0.005, CI lower > 0, >= 9/12 folds positive."""
    from nextgem import metrics, w2_score

    d = delta.dropna()
    lin = lineage.reindex(d.index)
    bt = metrics.bootstrap_mean(d, lin, B=B, seed=seed)
    fp = w2_score.folds_positive(d, lineage, folds)
    se = float(np.std(bt.replicates, ddof=1)) if getattr(bt, "replicates", None) is not None else np.nan
    return {"delta_lin": bt.estimate, "delta_lo": bt.ci_low, "delta_hi": bt.ci_high, "delta_se": se,
            "delta_line": float(d.mean()), "folds_positive": fp["folds_positive"], "n_folds": fp["n_folds"],
            "fold_means": fp["fold_means"],
            "tg_pass": bool(bt.estimate >= 0.005 and bt.ci_low > 0 and fp["folds_positive"] >= 9)}


def lin_mean(tap: pd.Series, lineage: pd.Series) -> float:
    from nextgem import metrics

    t = tap.dropna()
    return metrics.lineage_equal_mean(t, lineage.reindex(t.index))


def config_sha(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()
