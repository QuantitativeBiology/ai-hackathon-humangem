"""W1 baseline B6 (DESIGN.md §8): non-GEM ML, HistGradientBoosting on B2, B2m, own expression, CN, damaging and
hotspot mutation, neighbour-enzyme expression (H2 shared non-currency metabolite), max paralog expression,
isozyme count, medium class and serum flag. No lineage features. Label = L_prob (primary label; B6 is
evaluated in L_prob cells only).

Training: LOLO over the 12 eligible training lineages through ``nextgem.cv.lolo`` (the fit function sees only
``labels.fold(k)``). Inside a fold the B2/B2m features of fitting rows are *stacked* (leave-own-lineage-out
within the fitting set) so the learner never sees a prior that contains the row's own label; held-out rows
get priors from all fitting lines. Hyperparameters are fixed defaults at M2 (the grid is W7's):
max_iter 200, learning_rate 0.1, max_leaf_nodes 31, no early stopping, random_state SEEDS[0].
A final model is fitted on all training lines (for the validation look; DESIGN §4.2: fitted components from
the fitting tier only).

Usage: OMP_NUM_THREADS=8 .venv/bin/python scripts/w1_b6.py
Outputs (git-ignored): results/w1/tmp/b6_oof.parquet (training eligible lines x U, out-of-fold scores),
results/w1/tmp/b6_full.pkl (model + feature metadata). Committed: results/w1/b6_run.json.
"""

from __future__ import annotations

import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w1_common import OUT_DIR, PROJECT, TMP_DIR, atomic_write_bytes, atomic_write_parquet, current, run_main  # noqa: E402

LABEL = "L_prob"
HGB = dict(max_iter=200, learning_rate=0.1, max_leaf_nodes=31, early_stopping=False)


def stacked_priors(fit_long: pd.DataFrame, line_class: pd.Series, genes, label=LABEL):
    """Fitting rows: B2 / B2m from the fitting set minus the row's own lineage."""
    from nextgem.w1 import B2M_MIN_CLASS_LINES

    d = fit_long[["ModelID", "entrez", "lineage", label]].dropna(subset=[label]).copy()
    d["cls"] = d.ModelID.map(line_class).fillna("unknown")
    d["y"] = d[label].astype(float)
    pos = d.pivot_table(index="lineage", columns="entrez", values="y", aggfunc="sum").reindex(columns=genes).fillna(0)
    cnt = d.pivot_table(index="lineage", columns="entrez", values="y", aggfunc="count").reindex(columns=genes).fillna(0)
    posc = d.groupby(["cls", "lineage", "entrez"]).y.agg(["sum", "count"])
    lines_meta = fit_long[["ModelID", "lineage"]].drop_duplicates().set_index("ModelID").lineage
    cls_lines = pd.DataFrame({"lineage": lines_meta, "cls": lines_meta.index.map(line_class).fillna("unknown")})
    b2_rows, b2m_rows = {}, {}
    tot_p, tot_n = pos.sum(), cnt.sum()
    cache = {}
    for j in pos.index:
        p = (tot_p - pos.loc[j]) / (tot_n - cnt.loc[j]).replace(0, np.nan)
        cache[j] = p
    for l, j in lines_meta.items():
        b2_rows[l] = cache[j]
    # class-specific
    cpos = posc["sum"].unstack("entrez").reindex(columns=genes).fillna(0)
    ccnt = posc["count"].unstack("entrez").reindex(columns=genes).fillna(0)
    ccache = {}
    for l, j in lines_meta.items():
        c = cls_lines.loc[l, "cls"]
        key = (c, j)
        if key not in ccache:
            n_class_lines = int(((cls_lines.cls == c) & (cls_lines.lineage != j)).sum())
            if n_class_lines >= B2M_MIN_CLASS_LINES:
                sel_p = cpos.loc[c] if c in cpos.index.get_level_values(0) else None
                sp = sel_p.drop(index=j, errors="ignore").sum() if sel_p is not None else None
                sn = ccnt.loc[c].drop(index=j, errors="ignore").sum() if sel_p is not None else None
                v = (sp / sn.replace(0, np.nan)) if sp is not None else cache[j]
                ccache[key] = v.fillna(cache[j])
            else:
                ccache[key] = cache[j]
        b2m_rows[l] = ccache[key]
    return pd.DataFrame(b2_rows).T, pd.DataFrame(b2m_rows).T


def main() -> int:
    t0, c0 = time.time(), time.process_time()
    from sklearn.ensemble import HistGradientBoostingClassifier

    from nextgem import cv, labels, w1
    from nextgem.config import SEEDS

    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    omics = w1.load_omics()
    lf = w1.line_features()
    gf = w1.gene_features()
    nb = json.loads((TMP_DIR / "gene_neighbours.json").read_text())
    aggs = w1.expression_aggregates(omics["expr"], genes, nb["neighbours_h2"], nb["paralogs"])
    line_class = lf.base_medium_class
    class_codes = {c: i for i, c in enumerate(sorted(line_class.dropna().unique()))}
    med_idx = w1.B6_FEATURES.index("medium_class")
    cols = ["gene_effect", "dep_prob", "L_prob", "L_ge05", "L_ge06", "L_ge10"]
    log = []

    def fit(fit_long):
        fitted = w1.fit_components(fit_long, line_class, labels_=(LABEL,))
        b2s, b2ms = stacked_priors(fit_long, line_class, genes)
        fl = list(b2s.index)
        X = w1.b6_matrix(fl, genes, b2s, b2ms, omics, aggs, gf, lf, class_codes)
        Y = w1.wide(fit_long, LABEL, fl, genes).to_numpy().reshape(-1)
        ok = ~np.isnan(Y)
        clf = HistGradientBoostingClassifier(**HGB, random_state=SEEDS[0], categorical_features=[med_idx])
        Xo = X[ok]
        Xo[:, med_idx] = np.where(np.isnan(Xo[:, med_idx]), -1, Xo[:, med_idx])
        Xo[:, med_idx] = np.where(Xo[:, med_idx] < 0, np.nan, Xo[:, med_idx])
        clf.fit(Xo, Y[ok].astype(int))
        return clf, fitted, int(ok.sum())

    def predict(clf, fitted, lines):
        b2, b2m = w1.heldout_priors(fitted, lines, line_class, genes)
        X = w1.b6_matrix(lines, genes, b2, b2m, omics, aggs, gf, lf, class_codes)
        X[:, med_idx] = np.where(X[:, med_idx] < 0, np.nan, X[:, med_idx])
        p = clf.predict_proba(X)[:, 1].reshape(len(lines), len(genes))
        return pd.DataFrame(np.round(p, 8), index=lines, columns=genes)

    def fit_predict(fit_long, heldout_ids):
        tf = time.time()
        clf, fitted, n = fit(fit_long)
        pred = predict(clf, fitted, heldout_ids)
        log.append({"n_fit_pairs": n, "n_heldout_lines": len(heldout_ids), "fit_s": round(time.time() - tf, 1)})
        print(f"fold done: {log[-1]}", flush=True)
        return pred

    res = cv.lolo(fit_predict, lambda preds, held: preds, columns=cols)
    oof = pd.concat(res.values()).sort_index()
    oof.columns = oof.columns.astype(str)
    atomic_write_parquet(oof.reset_index().rename(columns={"index": "ModelID"}), TMP_DIR / "b6_oof.parquet")

    full_long = labels.train(columns=cols)
    clf, fitted, n = fit(full_long)
    atomic_write_bytes(TMP_DIR / "b6_full.pkl", pickle.dumps({"clf": clf, "class_codes": class_codes, "genes": genes,
                                                              "features": w1.B6_FEATURES, "hgb": HGB, "seed": SEEDS[0]}))
    rec = {"label": LABEL, "folds": {k: v for k, v in zip(res.keys(), log)}, "full_fit_pairs": n, "hgb": HGB,
           "features": w1.B6_FEATURES, "seed": SEEDS[0], "wall_s": round(time.time() - t0, 1)}
    (OUT_DIR / "b6_run.json").write_text(json.dumps(rec, indent=1))
    current().output(TMP_DIR / "b6_oof.parquet", TMP_DIR / "b6_full.pkl", OUT_DIR / "b6_run.json")
    print(json.dumps(rec, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(run_main(main, inputs=[TMP_DIR / f"{k}.parquet" for k in ("expr", "cn", "mut_dam", "mut_hot")] + [
        OUT_DIR / "line_features.tsv", OUT_DIR / "gene_features.tsv", TMP_DIR / "gene_neighbours.json",
        PROJECT / "results" / "splits" / "lineage_split.tsv"], seeds=[11], workers=8, reads_labels=True))
