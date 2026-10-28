"""Non-author review recompute of W7-05 (monotone B6) vs B6 and B2 (notes/reviews/M7_W7_review.md).

Independent of src/nextgem/w7.py and scripts/w7_*.py (not imported). Reused after reading (training-only feature
construction): nextgem.w1 (load_omics, line_features, gene_features, expression_aggregates, b6_matrix,
fit_components, heldout_priors, wide) and scripts/w1_b6.stacked_priors (fit-set-only, leave-own-lineage-out).
Labels: training tier only, via nextgem.cv.lolo (labels.fold(k); held-out labels only inside the score callback)
and labels.train() for the validation-score refit. No validation or test label is read.

Modes:
  lolo   : 12-fold LOLO x SEEDS[0:3] x {B6, W7-05}; per-line tAP tables -> tmp/; aggregates -> rv_w7_05_lolo.json
  valfit : fit {B6, W7-05} x SEEDS[0:3] on labels.train(), predict the validation lines of the look parquet
           (label-free features) and compare with the author's look scores -> rv_w7_05_valfit.json

Usage: OMP_NUM_THREADS=1 .venv/bin/python results/w7/review/rv_w7_05.py lolo [--workers 4] [--seeds 11,23,37]
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import concurrent.futures as cf  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import multiprocessing  # noqa: E402
import resource  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PROJECT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
TMP = HERE / "tmp"
LABEL = "L_prob"
HGB = dict(max_iter=200, learning_rate=0.1, max_leaf_nodes=31, early_stopping=False)
G = {}


def _w1_b6():
    spec = importlib.util.spec_from_file_location("w1_b6_rv", PROJECT / "scripts" / "w1_b6.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def setup():
    from nextgem import labels, w1

    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    omics = w1.load_omics()
    lf = w1.line_features()
    gf = w1.gene_features()
    nb = json.loads((PROJECT / "results" / "w1" / "tmp" / "gene_neighbours.json").read_text())
    aggs = w1.expression_aggregates(omics["expr"], genes, nb["neighbours_h2"], nb["paralogs"])
    split = labels.load_split().set_index("ModelID")
    G.update(genes=genes, omics=omics, lf=lf, gf=gf, aggs=aggs, split=split,
             codes={c: i for i, c in enumerate(sorted(lf.base_medium_class.dropna().unique()))},
             med=w1.B6_FEATURES.index("medium_class"), b2i=[w1.B6_FEATURES.index("b2"), w1.B6_FEATURES.index("b2m")],
             stacked=_w1_b6().stacked_priors)


def _X(lines, b2, b2m):
    from nextgem import w1

    X = w1.b6_matrix(lines, G["genes"], b2, b2m, G["omics"], G["aggs"], G["gf"], G["lf"], G["codes"])
    j = G["med"]
    X[:, j] = np.where(np.isnan(X[:, j]) | (X[:, j] < 0), np.nan, X[:, j])
    return X


def fit_models(fit_long, seed):
    """B6 and W7-05 on the fitting lines only (all must be training-tier)."""
    from sklearn.ensemble import HistGradientBoostingClassifier

    from nextgem import w1

    tiers = G["split"].tier.reindex(pd.unique(fit_long.ModelID))
    if not (tiers == "train").all():
        raise AssertionError("non-training line in a fit")
    fitted = w1.fit_components(fit_long, G["lf"].base_medium_class, labels_=(LABEL,))
    b2s, b2ms = G["stacked"](fit_long, G["lf"].base_medium_class, G["genes"])
    fl = list(b2s.index)
    X = _X(fl, b2s, b2ms)
    y = w1.wide(fit_long, LABEL, fl, G["genes"]).to_numpy().reshape(-1)
    ok = ~np.isnan(y)
    X, y = X[ok], y[ok].astype(int)
    out = {}
    for name, mono in (("B6", False), ("W7-05", True)):
        kw = dict(random_state=seed, categorical_features=[G["med"]])
        if mono:
            cst = [0] * X.shape[1]
            for i in G["b2i"]:
                cst[i] = 1
            kw["monotonic_cst"] = cst
        out[name] = HistGradientBoostingClassifier(**HGB, **kw).fit(X, y)
    return out, fitted


def predict(models, fitted, lines):
    from nextgem import w1

    b2, b2m = w1.heldout_priors(fitted, lines, G["lf"].base_medium_class, G["genes"], label=LABEL)
    X = _X(lines, b2, b2m)
    res = {n: pd.DataFrame(np.round(m.predict_proba(X)[:, 1], 8).reshape(len(lines), -1), index=lines,
                           columns=G["genes"]) for n, m in models.items()}
    res["B2"] = b2.fillna(0.0)
    return res


def _lolo_job(job):
    from nextgem import cv, metrics, w1

    fold, seed = job
    t0 = time.process_time()

    def fp(fit_long, ids):
        ev = [l for l in ids if l in G["omics"]["expr"].index]  # W1 primary set: lines with expression
        models, fitted = fit_models(fit_long, seed)
        return predict(models, fitted, ev)

    def score(preds, held):
        return {n: metrics.per_line(S, w1.wide(held, LABEL, list(S.index), G["genes"])) for n, S in preds.items()}

    res = cv.lolo(fp, score, columns=["gene_effect", LABEL], folds=[fold])[fold]
    out = []
    for n, t in res.items():
        t = t[["n", "n_pos", "estimable", "tap"]].copy()
        t.index.name = "ModelID"
        out.append(t.reset_index().assign(method=n, seed=seed, fold=fold))
    print(f"{fold} seed {seed}: {time.process_time() - t0:.0f} s CPU", flush=True)
    return pd.concat(out, ignore_index=True)


def summarise(df):
    from nextgem import labels, metrics
    from nextgem.config import SEEDS

    lineage = G["split"].lineage
    folds = labels.folds()
    B2 = df[df.method == "B2"].drop_duplicates("ModelID")  # seed-independent
    T = {m: d.groupby("ModelID").agg(n=("n", "first"), n_pos=("n_pos", "first"), estimable=("estimable", "first"),
                                     tap=("tap", "mean"))
         for m, d in [("B6", df[df.method == "B6"]), ("W7-05", df[df.method == "W7-05"]), ("B2", B2)]}
    out = {"n_lines": int(len(T["B6"])), "n_estimable": int(T["B6"].estimable.sum()),
           "seeds": sorted(int(s) for s in df.seed.unique())}
    for m, t in T.items():
        e = t[t.estimable]
        out[f"tap_lin_{m}"] = metrics.lineage_equal_mean(e.tap, lineage.reindex(e.index))
        if m != "B2":
            out[f"tap_lin_by_seed_{m}"] = []
            for s in out["seeds"]:
                d = df[(df.method == m) & (df.seed == s)].set_index("ModelID")
                d = d[d.estimable]
                out[f"tap_lin_by_seed_{m}"].append(metrics.lineage_equal_mean(d.tap, lineage.reindex(d.index)))
    for a, b in (("W7-05", "B6"), ("W7-05", "B2"), ("B6", "B2")):
        d = metrics.paired_line_delta(T[a], T[b])
        bt = metrics.bootstrap_mean(d, lineage.reindex(d.index), B=metrics.B_OTHER, seed=SEEDS[0])
        g = d.groupby(lineage.reindex(d.index)).agg(["mean", "size"])
        pos = sum(1 for f in folds if f in g.index and g.loc[f, "size"] >= 5 and g.loc[f, "mean"] > 0)
        out[f"d_{a}_vs_{b}"] = {"est": bt.estimate, "lo": bt.ci_low, "hi": bt.ci_high, "folds_pos": pos,
                                "n_folds": len(folds),
                                "tg_pass": bool(bt.estimate >= 0.005 and bt.ci_low > 0 and pos >= 9)}
    return out


def main_lolo(a):
    from nextgem import labels

    setup()
    seeds = [int(s) for s in a.seeds.split(",")]
    jobs = [(f, s) for s in seeds for f in labels.folds()]
    t0 = time.time()
    with cf.ProcessPoolExecutor(a.workers, mp_context=multiprocessing.get_context("fork")) as pool:
        df = pd.concat(list(pool.map(_lolo_job, jobs, chunksize=1)), ignore_index=True)
    TMP.mkdir(parents=True, exist_ok=True)
    df.to_parquet(TMP / "rv_perline_lolo.parquet", index=False)
    out = summarise(df)
    # bitwise comparison with the author's per-line tables (training per-line tAP only)
    au = pd.concat([pd.read_parquet(PROJECT / "results/w7/tmp" / f) for f in
                    ("perline_r1_W7-01.parquet", "perline_r2b_W7-05.parquet")], ignore_index=True)
    au["method"] = au.method.replace({"W7-01": "B6"})
    cmp = {}
    for m in ("B6", "W7-05", "B2"):
        x = df[df.method == m].drop_duplicates(["ModelID", "seed"] if m != "B2" else ["ModelID"])
        y = au[au.method == m].drop_duplicates(["ModelID", "seed"] if m != "B2" else ["ModelID"])
        keys = ["ModelID", "seed"] if m != "B2" else ["ModelID"]
        j = x.merge(y, on=keys, suffixes=("", "_au"))
        cmp[m] = {"n_matched": int(len(j)), "max_abs_diff_tap": float(np.nanmax(np.abs(j.tap - j.tap_au)))}
    out["vs_author_perline"] = cmp
    ru = [resource.getrusage(r) for r in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN)]
    out["cpu_h"] = sum(r.ru_utime + r.ru_stime for r in ru) / 3600
    out["wall_h"] = (time.time() - t0) / 3600
    (HERE / "rv_w7_05_lolo.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))


def _val_job(seed):
    from nextgem import labels

    fit_long = labels.train(columns=["gene_effect", LABEL])
    models, fitted = fit_models(fit_long, seed)
    lines = G["val_lines"]
    p = predict(models, fitted, lines)
    return seed, {k: v for k, v in p.items() if k != "B2"}


def main_valfit(a):
    setup()
    look = pd.read_parquet(PROJECT / "results/w7/val/w7_val_scores_W7_shortlist_1.parquet")
    lines = sorted(look.ModelID.unique())
    t = G["split"].tier.reindex(lines)
    assert (t == "validation").all(), "look lines must be validation-tier"
    G["val_lines"] = lines
    seeds = [int(s) for s in a.seeds.split(",")]
    t0 = time.time()
    with cf.ProcessPoolExecutor(min(a.workers, len(seeds)), mp_context=multiprocessing.get_context("fork")) as pool:
        res = list(pool.map(_val_job, seeds))
    out = {"n_val_lines": len(lines), "seeds": seeds, "n_rows_look": int(len(look))}
    for m, col in (("B6", "W7:B6"), ("W7-05", "W7:W7-05")):
        M = sum(r[1][m] for r in res) / len(res)
        long = M.stack().rename("rv").reset_index()
        long.columns = ["ModelID", "entrez", "rv"]
        j = look[["ModelID", "entrez", col]].merge(long, on=["ModelID", "entrez"], how="left", validate="1:1")
        out[col] = {"n_missing": int(j.rv.isna().sum()), "max_abs_diff": float(np.nanmax(np.abs(j[col] - j.rv))),
                    "bitwise_equal": bool(np.array_equal(j[col].to_numpy(), j.rv.to_numpy()))}
    ru = [resource.getrusage(r) for r in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN)]
    out["cpu_h"] = sum(r.ru_utime + r.ru_stime for r in ru) / 3600
    out["wall_h"] = (time.time() - t0) / 3600
    (HERE / "rv_w7_05_valfit.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["lolo", "valfit"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seeds", default="11,23,37")
    a = ap.parse_args()
    if a.workers > 5:
        raise SystemExit("<= 5 workers (review limit: <= 6 processes)")
    main_lolo(a) if a.mode == "lolo" else main_valfit(a)
