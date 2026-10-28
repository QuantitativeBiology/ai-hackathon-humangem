"""F1 training scoring (post hoc / exploratory, DESIGN §17 A15; RESEARCH_TREE D17). Committed before its first
label-reading run (LESSONS).

Labels only through ``nextgem.labels.train()`` and ``labels.fold(k)`` (via ``nextgem.cv.lolo``, for fold-scoped
CE_fit). No validation or test label is read; ``--synthetic-labels`` replaces every label by Bernoulli(0.12) draws
and opens no label file (code dry run; outputs under results/f1/tmp/ only).

Lines: the W1/M6 primary training set (686 expression lines, 12 lineages). Scores: e = 1 − ratio_r of the arm's
KO sweep under the line's M2(f=1) medium signature (W2 line_sweeps.tsv), W1 Entrez convention (absent gene -> 0,
6 decimals); P1 = the M6 cell H1+Bloo/M2 (asserted bitwise against the M6 e-matrix); arms reconstructed from
results/f1/sweeps (scripts/f1_sweeps.py).

Pre-declared contrasts (fixed in this file before any label is read):
  primary (TG, §9 rule as in m6_common.tg_pass): F1 (GPR rule + TYMS) − P1 on L_prob, all U genes:
      Δ̄_lin >= 0.005, 95% CI lower bound > 0 (lineage-stratified bootstrap, B = 2,000, SEEDS[0]), >= 9/12 folds > 0;
  components (reported separately, same statistics): TYMS − P1, GPR − P1, F1 − TYMS (GPR increment on the TYMS base);
  descriptive: every arm − G* (generic H1/M0);
  label cells for every contrast: all U genes; minus CE_fit (fold-scoped, LOLO); expressed pairs (log2(TPM+1) >= 1).
  TG verdicts are computed for every contrast on the 'all' cell; only F1 − P1 decides F1.
Binary (§6.8 fixed cutoff, ratio_r < 0.01): MCC and F1-score on L_ge06 and L_prob over the parent's genes with a
CRISPR column, per arm, and the paired Δ vs P1.

Outputs (results/f1/train/): cells.tsv, contrasts.tsv, binary.tsv, binary_paired.tsv, per_line.tsv (tAP and
estimability per line; no label counts), verification.tsv, summary.json.

Usage:
  .venv/bin/python scripts/f1_score_train.py --synthetic-labels      # dry run, no label file opened
  .venv/bin/python scripts/f1_score_train.py                         # training labels (from a clean committed tree)
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f1_common as fc  # noqa: E402
import m6_common as mc  # noqa: E402

W2 = mc.PROJECT / "results" / "w2"
LABEL_COLS = ["gene_effect", "L_prob", "L_ge06"]
YCELLS = ("L_prob", "L_prob_noCE", "L_prob_expr")
YNAME = {"L_prob": "all", "L_prob_noCE": "noCE", "L_prob_expr": "expressed_pairs"}
SYNTH_RATE = 0.12
BIN_METRICS = ("mcc", "f1")
METHODS = ["H1"] + [a.key for a in fc.ARMS]  # H1 = G* (generic Human1 on Ham's)
CONTRASTS = [("F1", "P1", "primary"), ("TYMS", "P1", "component"), ("GPR", "P1", "component"),
             ("F1", "TYMS", "component"), ("P1", "H1", "descriptive"), ("TYMS", "H1", "descriptive"),
             ("GPR", "H1", "descriptive"), ("F1", "H1", "descriptive")]
G: dict = {}


def arm_e_matrices(U):
    """(arm -> e-matrix (W2 sweep id x U), Entrez map, line_sweeps, P1 ratio frame). Label-free."""
    from f1_sweeps import arm_ratios, p1_rows
    from nextgem import paths, w2_score as ws

    ls = pd.read_csv(W2 / "controls" / "line_sweeps.tsv", sep="\t")
    sigs = sorted(set(mc.line_sweep_ids(ls, fc.PARENT, mc.M2)))
    P1 = p1_rows(sigs)
    delta = pd.read_parquet(fc.SWDIR / "ratio_r_delta_human1.parquet")
    meta = pd.read_csv(fc.SWDIR / "sweeps.tsv", sep="\t")
    ko = pd.read_parquet(paths.baseline_dir() / f"{fc.PARENT}_ko.parquet", columns=["gene", "entrez", "ratio_r"])
    emap = ws.gene_to_entrez(ko)
    return {a.key: mc.e_from_ratio(arm_ratios(a.key, P1, delta, meta), emap, U) for a in fc.ARMS}, emap, ls, P1


def line_matrix(E_sig: pd.DataFrame, ls: pd.DataFrame, lines) -> pd.DataFrame:
    keys = mc.line_sweep_ids(ls, fc.PARENT, mc.M2).reindex(list(lines))
    if keys.isna().any():
        raise KeyError(f"{int(keys.isna().sum())} lines have no M2(f=1) signature")
    M = E_sig.loc[keys.to_numpy()]
    M.index = list(lines)
    return M


def load(synthetic: bool) -> None:
    from m6_scores import e_matrices
    from nextgem import cv, labels, metrics, w1
    from nextgem.config import SEEDS

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    Udf = labels.gene_universe()
    U = list(Udf.entrez.astype(int))
    gen = w1.generic_table(Udf.entrez)
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    E = sorted(l for l in elig.index if lf.loc[l, "has_expression_matrix"])
    Earm, emap, ls, P1 = arm_e_matrices(U)
    S, failed = {}, {}
    for a in fc.ARMS:
        M = line_matrix(Earm[a.key], ls, E)
        bad = list(M.index[M.isna().any(axis=1)])
        failed[a.key] = bad
        S[a.key] = metrics.no_information_scores(M.fillna(0.0), bad, "line")
    S["H1"] = pd.DataFrame(np.tile(gen["e_H1"].reindex(U).to_numpy(), (len(E), 1)), index=E, columns=U)
    failed["H1"] = []
    # verification: P1 equals the M6 cell H1+Bloo/M2 bitwise; arms differ from P1 only where the sweeps differ
    Esw6, _, _ = e_matrices(U, ls)
    M6 = mc.cell_matrix(mc.CELL["H1+Bloo/M2"], ls, Esw6["human1"], E)
    checks = [{"check": "P1 score matrix equals the M6 cell H1+Bloo/M2 on the 686 lines (bitwise)",
               "pass": bool(np.array_equal(S["P1"].to_numpy(), M6.to_numpy()))},
              {"check": "no failed sweep line in any arm", "pass": all(not v for v in failed.values())}]
    expr = w1.load_omics()["expr"]
    folds = list(labels.folds())
    if synthetic:
        rng = np.random.default_rng(SEEDS[0])
        Y = {lab: pd.DataFrame((rng.random((len(E), len(U))) < SYNTH_RATE).astype(float), index=E, columns=U)
             for lab in ("L_prob", "L_ge06")}
        line_ce = {l: set() for l in E}
    else:
        line_class = lf.base_medium_class

        def fit_predict(fit_long, ids):
            return w1.fit_components(fit_long, line_class, labels_=("L_prob",)).ce["L_prob"], ids

        res = cv.lolo(fit_predict, lambda pr, held: pr, columns=["gene_effect", "L_prob"])
        line_ce = {}
        for k, (ce, ids) in res.items():
            for l in ids:
                line_ce[l] = ce
        tl = labels.train(columns=LABEL_COLS)
        Y = {lab: w1.wide(tl, lab, E, U) for lab in ("L_prob", "L_ge06")}
        del tl
    Y["L_prob_noCE"] = Y["L_prob"].where(w1.ce_mask(E, U, line_ce))
    Y["L_prob_expr"] = w1.apply_stratum(Y["L_prob"], ("pairs", "expressed"), lf, expr, None, None)
    G.update(split=split, lf=lf, U=U, gen=gen, E=E, lineage=split.lineage, ls=ls, S=S, failed=failed, Y=Y,
             folds=folds, checks=checks)


def cell_rows(LT: dict) -> pd.DataFrame:
    from nextgem import metrics, w1

    rows = []
    for k in METHODS:
        row, _, _ = w1.method_summary(G["S"][k], G["Y"]["L_prob"], G["lineage"], gene_unit=True)
        r = {"method_key": k, "method": fc.ARM[k].method if k in fc.ARM else "H1",
             "role": fc.ARM[k].role if k in fc.ARM else "G*", "n_failed_lines": len(G["failed"][k])}
        r.update({c: row[c] for c in ("n_lines", "n_lines_estimable", "n_lineages", "tap_lin", "tap_lin_lo", "tap_lin_hi",
                                      "lift_lin", "tap_line", "m_delta_ap_gene", "n_genes_estimable")})
        for y in ("L_prob_noCE", "L_prob_expr"):
            lt = LT[(k, y)]
            r[f"tap_lin_{YNAME[y]}"] = metrics.lineage_equal_mean(lt[lt.estimable].tap, G["lineage"])
        rows.append(r)
    return pd.DataFrame(rows)


def contrast_rows(LT: dict) -> pd.DataFrame:
    from nextgem import metrics, w1, w2_score as ws
    from nextgem.config import SEEDS

    rows = []
    for a, b, role in CONTRASTS:
        for y in YCELLS:
            d = metrics.paired_line_delta(LT[(a, y)], LT[(b, y)])
            r = {"arm": a, "reference": b, "role": role, "label_cell": YNAME[y], "n_lines": int(d.notna().sum())}
            r.update(w1.delta_summary(LT[(a, y)], LT[(b, y)], G["lineage"]))
            bt = metrics.bootstrap_mean(d, G["lineage"].reindex(d.index), B=metrics.B_OTHER, seed=SEEDS[0])
            r["se_boot"] = float(np.std(bt.replicates, ddof=1))
            r.update(ws.folds_positive(d, G["lineage"], G["folds"]))
            r["tg_pass"] = mc.tg_pass(r["delta_lin"], r["delta_lo"], r["folds_positive"]) if y == "L_prob" else None
            rows.append(r)
    return pd.DataFrame(rows)


def binary_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    from nextgem import metrics, w1
    from nextgem.config import SEEDS

    lin = G["lineage"]
    gset = w1.parent_gene_set(G["gen"], "H1")
    rows, PL = [], {}
    for k in METHODS:
        C = (G["S"][k] >= mc.CALL).loc[G["E"], gset]
        for lab in ("L_ge06", "L_prob"):
            T = G["Y"][lab].loc[G["E"], gset]
            pl = metrics.binary_per_line(C, T)
            sm = metrics.summarise_binary(pl, lin.reindex(pl.index))
            sm = sm[sm.metric.isin(BIN_METRICS)].assign(method_key=k, call_rule="FBA<0.01", label=lab, n_genes=len(gset))
            rows.append(sm)
            PL[(k, lab)] = pl
    pr = []
    for a, b, role in CONTRASTS:
        for lab in ("L_ge06", "L_prob"):
            for m in BIN_METRICS:
                d = (PL[(a, lab)][m] - PL[(b, lab)][m]).dropna()
                bt = metrics.bootstrap_mean(d, lin.reindex(d.index), B=metrics.B_OTHER, seed=SEEDS[0])
                pr.append({"arm": a, "reference": b, "role": role, "call_rule": "FBA<0.01", "label": lab, "metric": m,
                           "delta_lin": bt.estimate, "ci_lo": bt.ci_low, "ci_hi": bt.ci_high, "n_lines": len(d)})
    return pd.concat(rows, ignore_index=True), pd.DataFrame(pr)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--synthetic-labels", action="store_true", help="Bernoulli labels, no label file; tmp/ output")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    from nextgem.config import SEEDS

    out_dir = fc.TMP / "score_dry" if a.synthetic_labels else fc.TRDIR
    inputs = [W2 / "controls" / "line_sweeps.tsv", fc.SWDIR / "ratio_r_delta_human1.parquet", fc.SWDIR / "sweeps.tsv",
              fc.SWDIR / "spot_checks.tsv", fc.curated_rules_path(), fc.F1 / "edits_human1.tsv",
              mc.SWDIR / "ratio_r_human1.parquet", mc.SWDIR / "catalog.tsv"]
    with fc.provenance(inputs=inputs, seeds=[SEEDS[0]], solver="none (analysis)", workers=1,
                       allow_dirty=a.allow_dirty or a.synthetic_labels) as prov:
        prov.note(synthetic_labels=a.synthetic_labels,
                  labels=("synthetic Bernoulli(0.12); no label file" if a.synthetic_labels else
                          "labels.train() + cv.lolo/labels.fold(k) (CE_fit) only; training tier"),
                  label_files_frozen_sha256=mc.frozen_label_hashes())
        outs = run(out_dir, synthetic=a.synthetic_labels)
        prov.output(*outs)
    return 0


def run(out_dir: Path, synthetic: bool) -> list[Path]:
    from nextgem import metrics

    t0 = time.time()
    load(synthetic=synthetic)
    print(f"loaded: {len(G['E'])} lines, {len(G['U'])} genes, {time.time() - t0:.0f} s", flush=True)
    LT = {(k, y): metrics.per_line(G["S"][k], G["Y"][y]) for k in METHODS for y in YCELLS}
    cells = cell_rows(LT)
    con = contrast_rows(LT)
    binary, bpaired = binary_tables()
    checks = list(G["checks"])
    if not synthetic:  # P1 reproduces the committed M6 training per-line tAP (all, noCE)
        m6pl = pd.read_csv(mc.TRDIR / "per_line.tsv", sep="\t")
        for gs, y in (("all", "L_prob"), ("noCE", "L_prob_noCE")):
            w = m6pl[(m6pl.cell == "H1+Bloo/M2") & (m6pl.label_cell == gs)].set_index("ModelID")
            lt = LT[("P1", y)]
            same = set(w.index) == set(lt.index) and bool((w.estimable.reindex(lt.index) == lt.estimable).all())
            dif = float(np.nanmax(np.abs(w.tap.reindex(lt.index).to_numpy(float) - lt.tap.to_numpy(float))))
            checks.append({"check": f"P1 tAP_line equals M6 train per_line (H1+Bloo/M2, {gs})", "pass": same and dif < 1e-6,
                           "max_abs_diff": dif})
    ver = pd.DataFrame(checks)
    out_dir.mkdir(parents=True, exist_ok=True)
    outs = [fc.atomic_write_tsv(cells, out_dir / "cells.tsv"), fc.atomic_write_tsv(con, out_dir / "contrasts.tsv"),
            fc.atomic_write_tsv(binary, out_dir / "binary.tsv"), fc.atomic_write_tsv(bpaired, out_dir / "binary_paired.tsv"),
            fc.atomic_write_tsv(ver, out_dir / "verification.tsv")]
    pl = []
    for (k, y), lt in LT.items():
        pl.append(lt[["estimable", "tap"]].assign(method_key=k, label_cell=YNAME[y]))
    pl = pd.concat(pl).rename_axis("ModelID").reset_index()
    pl["lineage"] = pl.ModelID.map(G["lineage"])
    outs.append(fc.atomic_write_tsv(pl[["method_key", "label_cell", "ModelID", "lineage", "estimable", "tap"]],
                                    out_dir / "per_line.tsv"))
    prim = con[con.label_cell == "all"].set_index(["arm", "reference"])

    def tg(a, b):
        r = prim.loc[(a, b)]
        return {"delta_lin": r.delta_lin, "ci": [r.delta_lo, r.delta_hi], "folds_positive": int(r.folds_positive),
                "tg_pass": bool(r.tg_pass)}

    nce = con[con.label_cell == "noCE"].set_index(["arm", "reference"])
    summ = {"post_hoc": "A15 post hoc / exploratory", "synthetic_labels": synthetic, "n_lines": len(G["E"]),
            "n_lineages": int(G["lineage"].reindex(G["E"]).nunique()),
            "tg_rule": "delta_lin >= 0.005, CI lower > 0, >= 9/12 folds positive (L_prob, all U genes)",
            "F1_vs_P1": tg("F1", "P1"), "TYMS_vs_P1": tg("TYMS", "P1"), "GPR_vs_P1": tg("GPR", "P1"),
            "F1_vs_TYMS": tg("F1", "TYMS"),
            "noCE": {f"{a}_vs_{b}": {"delta_lin": nce.loc[(a, b)].delta_lin,
                                     "ci": [nce.loc[(a, b)].delta_lo, nce.loc[(a, b)].delta_hi],
                                     "folds_positive": int(nce.loc[(a, b)].folds_positive)}
                     for a, b, _ in CONTRASTS if b != "H1"},
            "F1_passes_TG": bool(prim.loc[("F1", "P1")].tg_pass),
            "checks_pass": bool(ver["pass"].all()), "wall_s": round(time.time() - t0, 1)}
    outs.append(fc.atomic_write_text(json.dumps(summ, indent=1, default=float) + "\n", out_dir / "summary.json"))
    if not ver["pass"].all():
        print(ver.to_string(), file=sys.stderr)
        raise RuntimeError("verification failed")
    print(f"written to {out_dir} ({time.time() - t0:.0f} s)", flush=True)
    return outs


if __name__ == "__main__":
    raise SystemExit(main())
