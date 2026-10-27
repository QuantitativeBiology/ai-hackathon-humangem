"""W5 independent review, part 2 (TRAINING LABELS via nextgem.labels.train() only; reviewer's own scoring code).

Scores (line-constant generic-model scores, so LOLO = per-lineage evaluation of a fixed score; nothing is fitted):
  e = 1 - ratio_r per U gene (Entrez; first row per Entrez in results/baseline order; absent from parent -> 0).
  H1, H2            results/baseline/<parent>_ko.parquet (released models)
  H1_TYMS           reviewer's full GLPK sweep of H1 with MAR04644 irreversible (r1_h1_tyms_full.parquet)
  H1_W5a            H1_TYMS with the genes of edited GPRs replaced by the reviewer's W5a ratios (r1)
  H2_W5a            H2 with the genes of edited GPRs replaced by the reviewer's W5a ratios (r1)
  H?_W5a_noPOLR     diagnostic: W5a without the RNA-polymerase AND->OR merges (not a candidate)
  H?_author_W5a     the author's ko_<parent>.parquet ratio_W5a column (consistency row)
Lines: training tier, non-rare stratum, with an expression profile (W1/W5 convention), 12 LOLO lineages.
Metric: reviewer's own tie-aware expected AP (checked against brute-force enumeration of within-tie orderings and
against nextgem.metrics.tap). Line estimable if >= 5 positives. Delta_lin = mean over lineages of the mean line Delta.
CI: lineage-stratified line bootstrap, B = 2000, seed = SEEDS[0]. Only aggregates are printed or written.
"""
from __future__ import annotations

import itertools
import json
import sys
import time

import numpy as np
import pandas as pd

from nextgem import labels, metrics  # metrics only for the cross-check of the reviewer's tAP
from nextgem.config import SEEDS

from nextgem import paths  # noqa: E402

OUT = paths.results_dir() / "w5" / "review"
T0, C0 = time.time(), time.process_time()


# ----------------------------------------------------------------------------- reviewer's tAP
def my_tap(s: np.ndarray, y: np.ndarray) -> float:
    """Expected AP over uniformly random orderings within tied score blocks (higher score first)."""
    P = y.sum()
    if P == 0:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    s, y = s[order], y[order]
    # block boundaries
    brk = np.flatnonzero(np.r_[True, s[1:] != s[:-1], True])
    total, N, Pb = 0.0, 0, 0
    for a, b in zip(brk[:-1], brk[1:]):
        n = b - a
        p = int(y[a:b].sum())
        if p:
            r = np.arange(1, n + 1, dtype=float)
            # P(item at in-block rank r is positive) = p/n; given that, expected # positives at ranks <= r within the
            # block = 1 + (r-1)(p-1)/(n-1)
            inner = 1 + (r - 1) * (p - 1) / (n - 1) if n > 1 else np.ones(1)
            total += float(((p / n) * (Pb + inner) / (N + r)).sum())
        N += n
        Pb += p
    return total / P


def brute_tap(s, y):
    """Average AP over all orderings that respect the score order (permute items inside tie blocks)."""
    vals = sorted(set(s), reverse=True)
    blocks = [[i for i in range(len(s)) if s[i] == v] for v in vals]
    aps = []
    for perm in itertools.product(*[list(itertools.permutations(b)) for b in blocks]):
        order = [i for blk in perm for i in blk]
        yy = np.asarray(y)[order]
        hits = np.cumsum(yy)
        aps.append(float((hits / np.arange(1, len(yy) + 1))[yy == 1].sum() / yy.sum()))
    return float(np.mean(aps))


def self_test():
    rng = np.random.default_rng(12345)
    worst_b, worst_m = 0.0, 0.0
    for _ in range(300):
        n = int(rng.integers(2, 9))
        s = rng.integers(0, 3, n).astype(float)
        y = rng.integers(0, 2, n)
        if y.sum() == 0:
            y[0] = 1
        a = my_tap(s, y)
        worst_b = max(worst_b, abs(a - brute_tap(list(s), list(y))))
        worst_m = max(worst_m, abs(a - metrics.tap(s, y.astype(float))))
    for _ in range(50):  # larger random cases vs metrics.tap only
        n = 3000
        s = np.round(rng.random(n), 1)
        y = (rng.random(n) < 0.1).astype(float)
        worst_m = max(worst_m, abs(my_tap(s, y) - metrics.tap(s, y)))
    return worst_b, worst_m


# ----------------------------------------------------------------------------- scores
def e_vector(gene_ratio: pd.Series, gmap: pd.DataFrame, U) -> pd.Series:
    """gene_ratio: Ensembl -> ratio_r.  gmap: baseline (gene, entrez) in baseline order."""
    df = gmap.assign(r=gmap.gene.map(gene_ratio)).dropna(subset=["entrez"])
    df = df.assign(entrez=df.entrez.astype(int)).drop_duplicates("entrez").set_index("entrez").r
    return (1.0 - df.reindex(U)).fillna(0.0).round(6)


def main():
    wb, wm = self_test()
    print(f"tAP self-test: max |mine - brute| = {wb:.2e}; max |mine - metrics.tap| = {wm:.2e}", flush=True)
    assert wb < 1e-12 and wm < 1e-9
    U = sorted(labels.gene_universe().entrez.astype(int))
    split = labels.load_split().set_index("ModelID")
    lf = pd.read_csv(paths.results_dir() / "w1" / "line_features.tsv", sep="\t").set_index("ModelID")
    E = sorted(l for l in split.index[(split.tier == "train") & (split.stratum != "rare")] if bool(lf.loc[l, "has_expression_matrix"]))
    lineage = split.lineage
    base = {p: pd.read_parquet(paths.baseline_dir() / f"{p}_ko.parquet", columns=["gene", "entrez", "ratio_r"]) for p in ("human1", "human2")}
    r1 = pd.read_csv(OUT / "r1_ko_affected_genes.tsv", sep="\t")
    tyms = pd.read_parquet(OUT / "r1_h1_tyms_full.parquet").set_index("gene").ratio_tyms
    au = {p: pd.read_parquet(paths.results_dir() / "w5" / f"ko_{p}.parquet") for p in ("human1", "human2")}
    sc = {}
    checks = {}
    for p, S in (("human1", "H1"), ("human2", "H2")):
        gm = base[p][["gene", "entrez"]]
        g0 = base[p].set_index("gene").ratio_r
        sc[S] = e_vector(g0, gm, U)
        aff = r1[r1.parent == p].set_index("gene")
        start = tyms if p == "human1" else g0
        if p == "human1":
            sc["H1_TYMS"] = e_vector(tyms, gm, U)
            a = au[p].set_index("gene")
            checks["H1_TYMS_full_vs_author_ratio_base_maxabs"] = float((tyms.reindex(a.index) - a.ratio_base).abs().max())
            checks["H1_TYMS_genes_changed_vs_released"] = sorted(a.symbol[(tyms.reindex(a.index) != g0.reindex(a.index))].astype(str))
        else:
            checks["H2_generic_glpk_vs_baseline_maxabs_affected"] = float((aff.ratio_generic_glpk - g0.reindex(aff.index)).abs().max())
        for arm, col in (("W5a", "ratio_W5a"), ("W5a_noPOLR", "ratio_W5a_noPOLR")):
            r = start.copy()
            r.loc[aff.index] = aff[col]
            sc[f"{S}_{arm}"] = e_vector(r, gm, U)
        a = au[p].set_index("gene")
        checks[f"{S}_W5a_vs_author_maxabs"] = float((aff.ratio_W5a - a.ratio_W5a.reindex(aff.index)).abs().max())
        sc[f"{S}_author_W5a"] = e_vector(a.ratio_W5a, gm, U)
        checks[f"{S}_W5a_vs_author_score_maxabs_U"] = float((sc[f"{S}_W5a"] - sc[f"{S}_author_W5a"]).abs().max())
    # ---- labels (training tier only)
    tl = labels.train(columns=["L_prob", "L_ge05"])
    res = {"self_test": {"max_abs_vs_brute": wb, "max_abs_vs_metrics_tap": wm}, "checks": checks, "n_lines": len(E)}
    rows = []
    for lab in ("L_prob", "L_ge05"):
        Y = tl.pivot(index="ModelID", columns="entrez", values=lab).reindex(index=E, columns=U).astype(float)
        taps = {}
        for name, e in sc.items():
            ev = e.reindex(U).to_numpy()
            t = np.full(len(E), np.nan)
            for i, l in enumerate(E):
                y = Y.loc[l].to_numpy()
                m = ~np.isnan(y)
                if y[m].sum() >= 5:
                    t[i] = my_tap(ev[m], y[m])
            taps[name] = pd.Series(t, index=E)
        lin = lineage.reindex(E)
        contrasts = [("H1_TYMS", "H1"), ("H1_W5a", "H1"), ("H2_W5a", "H2"), ("H1_W5a", "H1_TYMS"),
                     ("H1_W5a_noPOLR", "H1"), ("H2_W5a_noPOLR", "H2"), ("H1_author_W5a", "H1"), ("H2_author_W5a", "H2")]
        rng_seed = SEEDS[0]
        for a, b in contrasts:
            d = (taps[a] - taps[b]).dropna()
            ln = lin.reindex(d.index)
            per = d.groupby(ln).mean()
            nest = d.groupby(ln).size()
            est = float(per.mean())
            pos = int(((per > 0) & (nest >= 5)).sum())
            rng = np.random.default_rng(rng_seed)
            groups = [d[ln == k].to_numpy() for k in per.index]
            boots = np.empty(2000)
            for bi in range(2000):
                boots[bi] = np.mean([g[rng.integers(0, len(g), len(g))].mean() for g in groups])
            lo, hi = np.percentile(boots, [2.5, 97.5])
            rows.append({"label": lab, "method": a, "ref": b, "delta_lin": round(est, 5), "ci_lo": round(float(lo), 5),
                         "ci_hi": round(float(hi), 5), "folds_positive": pos, "n_folds": int(per.size),
                         "n_lines_est": int(d.size), "TG_pass": bool(est >= 0.005 and lo > 0 and pos >= 9)})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "r2_train_deltas.tsv", sep="\t", index=False)
    res["deltas"] = out.to_dict(orient="records")
    res["wall_s"] = round(time.time() - T0, 1)
    res["cpu_s"] = round(time.process_time() - C0, 1)
    (OUT / "r2_summary.json").write_text(json.dumps(res, indent=1, default=float))
    print(out.to_string(index=False))
    print(json.dumps(checks, indent=1, default=float))


if __name__ == "__main__":
    sys.exit(main())
