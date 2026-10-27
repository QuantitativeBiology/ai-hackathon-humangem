"""W1 §6.8 binary endpoints on the training tier (DESIGN.md §6.1, §6.8; M2): generic H1 and H2 under the
Human1-paper and Human2-paper protocols, plus our fixed-cutoff MCC. Labels only via ``nextgem.labels``.

- Call rules: P57 (KO breaks >= 1 of the 57 essential tasks; results/w1/p57_<parent>.parquet), PGR (GR task
  fails; Human1-paper biomass-only variant), FBA (ratio_r < 0.01; §6.1 fixed cutoff).
- Labels: L-ge06 (gene effect < -0.6, the papers' threshold); the Human1-paper grid (21 thresholds,
  -1.0 ... 0.0); L_prob for our fixed-cutoff MCC.
- Gene set (§6.8): the parent's genes with a CRISPR column (U ∩ parent) and a non-NaN label in the line;
  U-wide (absent genes called non-essential) is a sensitivity row.
- Lines: eligible-lineage training lines with expression (primary, as for tAP); all eligible-lineage lines
  (secondary). Summaries: lineage-equal mean, line mean, median; NaN lines excluded and counted (primary)
  or set to 0 (sensitivity); 95% CI by lineage-stratified line bootstrap (B = 2,000, SEEDS[0]).
- Context only (flagged "not a paper protocol"): B2 and B3 dichotomised at the fitting-line positive rate
  (LOLO).
- Controls: constant non-essential caller (MCC NaN; 0 in the sensitivity row), oracle caller (MCC 1).

Usage: .venv/bin/python scripts/w1_binary.py
Outputs: results/w1/train_binary.tsv, train_binary_grid.tsv, train_binary_paired.tsv
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w1_common import BASELINE, OUT_DIR, PROJECT, TMP_DIR, atomic_write_tsv, current, run_main  # noqa: E402

LABEL_COLS = ["gene_effect", "dep_prob", "L_prob", "L_ge05", "L_ge06", "L_ge10"]


def main() -> int:
    t0, c0 = time.time(), time.process_time()
    from nextgem import cv, labels, metrics, w1
    from nextgem.config import SEEDS

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    gen = w1.generic_table(U.entrez)
    expr = w1.load_omics()["expr"]
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    E_all = sorted(elig.index)
    E = [l for l in E_all if lf.loc[l, "has_expression_matrix"]]
    lineage = split.lineage
    tl = labels.train(columns=LABEL_COLS)
    GE = w1.wide(tl, "gene_effect", E_all, genes)
    L06 = w1.wide(tl, "L_ge06", E_all, genes)
    LP = w1.wide(tl, "L_prob", E_all, genes)

    rows, grid_rows, perline = [], [], {}

    def record(protocol, method, rule, label, gene_set, lineset, truth, calls_frame=None, calls=None, gset=None):
        lines = E if lineset == "expr" else E_all
        T = truth.loc[lines, gset]
        if calls_frame is None:
            summ, pl = w1.binary_endpoints(calls, T, lineage, gset)
        else:
            C = calls_frame.loc[lines, gset]
            pl = metrics.binary_per_line(C, T)
            summ = metrics.summarise_binary(pl, lineage.reindex(pl.index))
            summ["ci_lo"] = np.nan
            summ["ci_hi"] = np.nan
            for m in metrics.BINARY_METRICS:
                v = pl[m].dropna()
                if len(v):
                    bt = metrics.bootstrap_mean(v, lineage.reindex(v.index), B=metrics.B_OTHER, seed=SEEDS[0])
                    i = (summ.metric == m) & (summ.nan_policy == "exclude")
                    summ.loc[i, ["ci_lo", "ci_hi"]] = [bt.ci_low, bt.ci_high]
        tot = pl[["tp", "fp", "tn", "fn"]].sum()
        summ = summ.assign(protocol=protocol, method=method, call_rule=rule, label=label, gene_set=gene_set, lineset=lineset,
                           n_genes=len(gset), n_called=int(calls.reindex(gset).sum()) if calls is not None else np.nan,
                           TP=int(tot.tp), FP=int(tot.fp), TN=int(tot.tn), FN=int(tot.fn))
        rows.append(summ)
        perline[(protocol, method, rule, label, gene_set, lineset)] = pl
        return pl

    for s in ("H1", "H2"):
        pg = w1.parent_gene_set(gen, s)
        for lineset in ("expr", "all"):
            for gname, gset in (("parent", pg), ("U", genes)):
                if lineset == "all" and gname == "U":
                    continue
                # Human1-paper protocol: P57 and PGR x L-ge06
                record("Human1-paper", s, "P57", "L_ge06", gname, lineset, L06, calls=gen[f"p57_{s}"], gset=gset)
                record("Human1-paper", s, "PGR", "L_ge06", gname, lineset, L06, calls=gen[f"pgr_{s}"], gset=gset)
                # Human2-paper protocol: P57 x L-ge06 (identical computation; reported under its own name)
                record("Human2-paper", s, "P57", "L_ge06", gname, lineset, L06, calls=gen[f"p57_{s}"], gset=gset)
                # our fixed cutoff (§6.1): ratio_r < 0.01, primary label and L-ge06
                record("fixed-cutoff", s, "FBA<0.01", "L_prob", gname, lineset, LP, calls=gen[f"fba_call_{s}"], gset=gset)
                record("fixed-cutoff", s, "FBA<0.01", "L_ge06", gname, lineset, L06, calls=gen[f"fba_call_{s}"], gset=gset)
        # Human1-paper robustness grid (parent genes, expression lines)
        for rule in ("P57", "PGR"):
            C = pd.DataFrame(np.tile(gen[f"{rule.lower()}_{s}"].reindex(pg).to_numpy(bool), (len(E), 1)), index=E, columns=pg)
            g = metrics.binary_threshold_grid(C, GE.loc[E, pg])
            for thr, d in g.groupby("threshold"):
                sm = metrics.summarise_binary(d, lineage.reindex(d.index))
                sm = sm[sm.nan_policy == "exclude"].assign(method=s, call_rule=rule, threshold=thr)
                grid_rows.append(sm)

    # controls (H2 gene set, expression lines)
    pg2 = w1.parent_gene_set(gen, "H2")
    record("control", "constant_nonessential", "none", "L_ge06", "parent", "expr", L06,
           calls=pd.Series(False, index=gen.index), gset=pg2)
    record("control", "constant_essential", "all", "L_ge06", "parent", "expr", L06,
           calls=pd.Series(True, index=gen.index), gset=pg2)
    oracle = L06.fillna(0).astype(bool)
    record("control", "oracle", "label", "L_ge06", "parent", "expr", L06, calls_frame=oracle, gset=pg2)

    # context: non-GEM baselines dichotomised at the fitting-line positive rate (LOLO), not a paper protocol
    line_class = lf.base_medium_class

    def fit_predict(fit_long, ids):
        f = w1.fit_components(fit_long, line_class, labels_=("L_ge06",))
        ids = [i for i in ids if i in set(E)]
        p = f.prior["L_ge06"].reindex(genes).fillna(0.0)
        rate = f.rate["L_ge06"]
        thr_b2 = metrics.rate_threshold(np.tile(p.to_numpy(), len(f.fit_lines)), rate)
        fit_expr = expr.reindex(index=[l for l in f.fit_lines if l in expr.index], columns=genes).fillna(0.0)
        thr_b3 = metrics.rate_threshold(fit_expr.to_numpy(), rate)
        b2 = pd.DataFrame(np.tile(metrics.calls_at_threshold(p, thr_b2).to_numpy(), (len(ids), 1)), index=ids, columns=genes)
        b3 = metrics.calls_at_threshold(expr.reindex(index=ids, columns=genes).fillna(0.0), thr_b3)
        return b2, b3

    folds = cv.lolo(fit_predict, lambda p, held: p, columns=LABEL_COLS)
    B2c = pd.concat([v[0] for v in folds.values()]).loc[E]
    B3c = pd.concat([v[1] for v in folds.values()]).loc[E]
    for name, C in (("B2[L_ge06]", B2c), ("B3", B3c)):
        for gname, gset in (("H2-parent", pg2), ("U", genes)):
            record("context (not a paper protocol)", name, "fitting-rate threshold", "L_ge06", gname, "expr", L06,
                   calls_frame=C, gset=gset)

    out = pd.concat(rows, ignore_index=True)
    front = ["protocol", "method", "call_rule", "label", "gene_set", "lineset", "metric", "nan_policy", "lineage_equal_mean",
             "ci_lo", "ci_hi", "line_mean", "median", "n_lines_used", "n_nan"]
    out = out[front + [c for c in out.columns if c not in front]]
    atomic_write_tsv(out, OUT_DIR / "train_binary.tsv")
    atomic_write_tsv(pd.concat(grid_rows, ignore_index=True), OUT_DIR / "train_binary_grid.tsv")

    # paired per-line differences H2 - H1 (each on its own parent gene set), expression lines
    pr = []
    for (proto, rule, lab) in (("Human1-paper", "P57", "L_ge06"), ("Human1-paper", "PGR", "L_ge06"),
                               ("Human2-paper", "P57", "L_ge06"), ("fixed-cutoff", "FBA<0.01", "L_prob")):
        a = perline[(proto, "H2", rule, lab, "parent", "expr")]
        b = perline[(proto, "H1", rule, lab, "parent", "expr")]
        for m in ("mcc", "f1", "precision", "recall", "accuracy"):
            d = (a[m] - b[m]).dropna()
            bt = metrics.bootstrap_mean(d, lineage.reindex(d.index), B=metrics.B_OTHER, seed=SEEDS[0])
            pr.append({"protocol": proto, "call_rule": rule, "label": lab, "metric": m, "contrast": "H2 - H1",
                       "delta_lin": bt.estimate, "ci_lo": bt.ci_low, "ci_hi": bt.ci_high,
                       "delta_line": float(d.mean()), "n_lines": len(d)})
    atomic_write_tsv(pd.DataFrame(pr), OUT_DIR / "train_binary_paired.tsv")
    current().output(OUT_DIR / "train_binary.tsv", OUT_DIR / "train_binary_grid.tsv", OUT_DIR / "train_binary_paired.tsv")
    show = out[(out.metric == "mcc") & (out.nan_policy == "exclude") & (out.lineset == "expr")]
    print(show[["protocol", "method", "call_rule", "label", "gene_set", "lineage_equal_mean", "ci_lo", "ci_hi", "line_mean",
                "n_nan"]].to_string())
    print(f"runtime {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(run_main(main, inputs=[BASELINE / "human1_ko.parquet", BASELINE / "human2_ko.parquet",
        OUT_DIR / "p57_human1.parquet", OUT_DIR / "p57_human2.parquet", TMP_DIR / "expr.parquet", OUT_DIR / "line_features.tsv",
        PROJECT / "results" / "splits" / "lineage_split.tsv"], seeds=[11], reads_labels=True))
