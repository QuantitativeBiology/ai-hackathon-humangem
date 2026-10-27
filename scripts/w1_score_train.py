"""W1 on the training tier (DESIGN.md §6, §8, W1; M2). Labels only through ``nextgem.labels`` (train / fold API).

Evaluation lines: training lines of the 12 eligible lineages (rare lines fit, never evaluated). Primary
set = those with a default expression profile (DESIGN §4.3); all eligible-lineage lines are a secondary row
for the generic models and label-only baselines.
Every fitted component (B2, B2m, B2c, CE_fit, selective genes, B6) is LOLO: the scores of a line in lineage k
come from training minus k (``nextgem.cv.lolo``). Generic H1/H2 (and B0, B1, B3, B8) are not fitted.

Grid (§6.6): labels L_prob (primary per the kappa rule), L_ge05, L_ge10 x gene sets {all, minus CE_fit,
minus DepMap common essentials} x universes {U, U∩} x strata (§5.4) x units {line (lineage-equal and
line-weighted), gene (mΔAP_gene), pooled}; value, E[AP_random], prevalence, lift, 95% CI (lineage-stratified
line bootstrap, B = 2,000, SEEDS[0]) and paired deltas vs generic H1 and H2. NA per §5.2.

Also: per-line / per-lineage tables (primary cells), the per-gene error table (W1 outputs, ext-PI SF1),
addressable-set mAP and recall ceiling (§5.4), SESOI_C1b (frozen, §6.4), the M2 gate, the BCa skew check.

Usage: OMP_NUM_THREADS=1 .venv/bin/python scripts/w1_score_train.py [--workers 12]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w1_common import BASELINE, OUT_DIR, PROJECT, TMP_DIR, atomic_write_tsv, current, run_main  # noqa: E402

G: dict = {}  # shared with forked workers
LABEL_COLS = ["gene_effect", "dep_prob", "L_prob", "L_ge05", "L_ge06", "L_ge10"]
GENESETS = ("all", "noCE", "noDepMapCE")
UNIVERSES = ("U", "Ucap")


def build():
    from nextgem import cv, labels, w1

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    ucap = [g for g, a, b in zip(U.entrez, U.in_human1, U.in_human2) if a and b]
    gen = w1.generic_table(U.entrez)
    omics = w1.load_omics()
    line_class = lf.base_medium_class
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    E_all = sorted(elig.index)
    E = [l for l in E_all if lf.loc[l, "has_expression_matrix"]]
    lineage = split.lineage

    b6 = pd.read_parquet(TMP_DIR / "b6_oof.parquet").set_index("ModelID")
    b6.columns = b6.columns.astype(int)

    def fit_predict(fit_long, ids):
        return w1.fit_components(fit_long, line_class), ids

    folds = cv.lolo(fit_predict, lambda p, held: p, columns=LABEL_COLS)
    parts: dict[str, list] = {}
    line_ce = {lab: {} for lab in w1.GRID_LABELS + ("L_ge06",)}
    line_prior = {lab: {} for lab in w1.GRID_LABELS}
    ce_sizes = {}
    for k, (fitted, ids) in folds.items():
        ids = [i for i in ids if i in set(E_all)]
        S = w1.score_matrices(ids, genes, fitted, gen, line_class, omics["expr"], b6=b6)
        for name, df in S.items():
            parts.setdefault(name, []).append(df)
        for lab in line_ce:
            for l in ids:
                line_ce[lab][l] = fitted.ce[lab]
        for lab in line_prior:
            for l in ids:
                line_prior[lab][l] = fitted.prior[lab]
        ce_sizes[k] = {lab: len(fitted.ce[lab]) for lab in fitted.ce}
    scores = {n: pd.concat(v).loc[E_all] for n, v in parts.items()}
    train_long = labels.train(columns=LABEL_COLS)
    Y = {lab: w1.wide(train_long, lab, E_all, genes) for lab in w1.GRID_LABELS + ("L_ge06",)}
    GE = w1.wide(train_long, "gene_effect", E_all, genes)
    depmap_ce = labels.label_derived_list("depmap_inferred_common_essentials", purpose="reporting")
    panels = w1.gene_panels()
    gf = w1.gene_features()
    G.update(dict(split=split, lf=lf, genes=genes, ucap=ucap, gen=gen, omics=omics, E_all=E_all, E=E, lineage=lineage,
                  scores=scores, Y=Y, GE=GE, line_ce=line_ce, line_prior=line_prior, depmap_ce=depmap_ce,
                  strata=w1.stratum_specs(panels, gf), panels=panels, gf=gf, ce_sizes=ce_sizes, train_long=train_long,
                  U=U))


def cell_Y(label, geneset, universe, stratum, lines):
    from nextgem import w1

    Y = G["Y"][label].loc[lines]
    cols = G["genes"] if universe == "U" else G["ucap"]
    Y = Y.loc[:, cols]
    if geneset == "noCE":
        Y = Y.where(w1.ce_mask(Y.index, Y.columns, G["line_ce"][label]))
    elif geneset == "noDepMapCE":
        Y = Y.loc[:, [g for g in Y.columns if g not in G["depmap_ce"]]]
    sel = w1.selective_mask(Y.index, Y.columns, G["line_prior"][label]) if stratum == "selective_genes" else None
    return w1.apply_stratum(Y, G["strata"][stratum], G["lf"], G["omics"]["expr"], G["omics"]["cn"], sel)


def run_cell(job):
    from nextgem import w1

    label, geneset, universe, stratum, lineset = job
    lines = G["E"] if lineset == "expr" else G["E_all"]
    Y = cell_Y(label, geneset, universe, stratum, lines)
    names = w1.methods_for(label, G["scores"].keys())
    if lineset == "all":
        names = [n for n in names if n in ("H1", "H2", "H1_task", "H2_task", "B0") or n.startswith(("B2[", "B2c", "B1_"))]
    primary = (label == "L_prob" and geneset in ("all", "noCE") and universe == "U" and stratum == "none" and lineset == "expr")
    rows, lts = [], {}
    for n in names:
        S = G["scores"][n].loc[Y.index, Y.columns]
        row, lt, _ = w1.method_summary(S, Y, G["lineage"], gene_unit=True, auroc=primary)
        row.update({"label": label, "geneset": geneset, "universe": universe, "stratum": stratum, "lineset": lineset,
                    "method": n})
        rows.append(row)
        lts[n] = lt
    for ref in ("H1", "H2"):
        for r in rows:
            if r["method"] == ref:
                continue
            d = w1.delta_summary(lts[r["method"]], lts[ref], G["lineage"])
            r.update({f"{k}_vs_{ref}": v for k, v in d.items()})
    out_lines = None
    if primary:
        out_lines = pd.concat([lt.assign(method=n) for n, lt in lts.items()]).rename_axis("ModelID").reset_index()
        out_lines["geneset"] = geneset
    return rows, out_lines


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args(argv)
    t0, c0 = time.time(), time.process_time()
    from nextgem import labels, metrics, w1

    build()
    print(f"built: {len(G['E_all'])} eligible lines, {len(G['E'])} with expression, {len(G['scores'])} methods, "
          f"{time.time() - t0:.0f} s", flush=True)
    jobs = [(lab, gs, u, st, "expr") for lab in w1.GRID_LABELS for gs in GENESETS for u in UNIVERSES for st in G["strata"]]
    jobs += [("L_prob", gs, "U", "none", "all") for gs in GENESETS]
    rows, per_line = [], []
    with mp.get_context("fork").Pool(a.workers) as pool:
        for i, (r, pl) in enumerate(pool.imap_unordered(run_cell, jobs, chunksize=1), 1):
            rows += r
            if pl is not None:
                per_line.append(pl)
            if i % 25 == 0:
                print(f"{i}/{len(jobs)} cells, {time.time() - t0:.0f} s", flush=True)
    grid = pd.DataFrame(rows)
    key = ["label", "geneset", "universe", "stratum", "lineset", "method"]
    # B1: mean over SEEDS[0:5] as one row per cell
    b1 = grid[grid.method.str.startswith("B1_")]
    num = b1.select_dtypes("number").columns
    b1m = b1.groupby(key[:-1])[list(num)].mean().reset_index().assign(method="B1_mean")
    b1m["na"] = b1.groupby(key[:-1]).na.any().to_numpy()
    grid = pd.concat([grid, b1m], ignore_index=True)
    front = key + ["na", "n_lines", "n_lines_estimable", "n_lineages", "tap_lin", "tap_lin_lo", "tap_lin_hi", "e_random_lin",
                   "lift_lin", "prevalence_lin"]
    grid = grid[front + [c for c in grid.columns if c not in front]].sort_values(key)
    atomic_write_tsv(grid, OUT_DIR / "train_grid.tsv")
    pl = pd.concat(per_line, ignore_index=True)
    pl["lineage"] = pl.ModelID.map(G["lineage"])
    atomic_write_tsv(pl, OUT_DIR / "train_per_line.tsv")

    # per-lineage (primary cells)
    lin_rows = []
    for gs, sub in pl.groupby("geneset"):
        ref = sub[sub.method == "H2"].set_index("ModelID")
        for n, s in sub.groupby("method"):
            t = metrics.per_lineage_summary(s.set_index("ModelID"), G["lineage"], reference=ref)
            lin_rows.append(t.reset_index().assign(method=n, geneset=gs, delta_ref="H2"))
    atomic_write_tsv(pd.concat(lin_rows, ignore_index=True), OUT_DIR / "train_per_lineage.tsv")

    extras = derived_outputs(grid, pl)
    rec = {"n_cells": len(jobs), "n_rows": len(grid), "wall_s": round(time.time() - t0, 1), "workers": a.workers,
           "ce_fit_sizes_per_fold": G["ce_sizes"], **extras}
    (OUT_DIR / "train_summary.json").write_text(json.dumps(rec, indent=1, default=float))
    current().output(*[OUT_DIR / f for f in ("train_grid.tsv", "train_per_line.tsv", "train_per_lineage.tsv",
                                             "train_per_gene.tsv", "train_addressable.tsv", "sesoi_c1b.json",
                                             "train_summary.json")])
    print(json.dumps({k: v for k, v in rec.items() if k != "ce_fit_sizes_per_fold"}, indent=1, default=float))
    return 0


def derived_outputs(grid: pd.DataFrame, pl: pd.DataFrame) -> dict:
    from nextgem import labels, metrics, w1

    out = {}
    prim = grid[(grid.label == "L_prob") & (grid.geneset == "all") & (grid.universe == "U") & (grid.stratum == "none")
                & (grid.lineset == "expr")].set_index("method")
    # ---- M2 gate (W1 decision): H2 lift on training > 1
    out["gate"] = {"H2_lift_lin": prim.loc["H2", "lift_lin"], "H2_pooled_lift": prim.loc["H2", "pooled_lift"],
                   "H1_lift_lin": prim.loc["H1", "lift_lin"], "H1_pooled_lift": prim.loc["H1", "pooled_lift"],
                   "pass": bool(prim.loc["H2", "lift_lin"] > 1 and prim.loc["H2", "pooled_lift"] > 1)}
    # ---- BCa check (§6.5) on the training line deltas available at M2
    p = pl[pl.geneset == "all"].pivot(index="ModelID", columns="method", values="tap")
    out["bca_skew_check"] = {f"{a}-{b}": {"skew": float(pd.Series(p[a] - p[b]).dropna().skew()),
                                          "triggers_bca": metrics.skew_triggers_bca(p[a] - p[b])}
                             for a, b in (("H2", "H1"), ("B2[L_prob]", "H2"), ("B8+_H2", "H2"), ("B6[L_prob]", "H2"))}
    # ---- SESOI_C1b (§6.4): constant scorer on all training lines, CE_fit from all training lines, L_prob
    tl = G["train_long"]
    lines = sorted(tl.ModelID.unique())
    Y = w1.wide(tl, "L_prob", lines, G["genes"])
    ce = labels.ce_fit(tl, "L_prob")
    const = pd.DataFrame(0.0, index=Y.index, columns=Y.columns)
    lt_all = metrics.per_line(const, Y)
    Yn = Y.loc[:, [g for g in Y.columns if g not in ce]]
    lt_noce = metrics.per_line(const.loc[:, Yn.columns], Yn)
    s = metrics.sesoi_c1b(lt_noce.e_random, lt_all.e_random)
    sesoi = {"SESOI_C1b": s, "mean_e_random_U_minus_CE_fit": float(np.nanmean(lt_noce.e_random)),
             "mean_e_random_U": float(np.nanmean(lt_all.e_random)), "n_training_lines": len(lines),
             "n_lines_estimable_U": int(lt_all.estimable.sum()), "n_lines_estimable_noCE": int(lt_noce.estimable.sum()),
             "CE_fit_size": len(ce), "label": "L_prob", "rule": "DESIGN 6.4: 0.01 x mean E[AP_random] on U minus CE_fit / "
             "mean E[AP_random] on U; training lines (all 767), estimable lines (>= 5 positives) of each set; CE_fit from "
             "all training lines", "frozen_at": "M2 (before any validation look)"}
    (OUT_DIR / "sesoi_c1b.json").write_text(json.dumps(sesoi, indent=1))
    out["sesoi_c1b"] = sesoi
    # ---- addressable set: recall ceiling (fraction of positive pairs outside the set), E lines
    gf = G["gf"]
    rc = []
    for lab in w1.GRID_LABELS:
        Yl = G["Y"][lab].loc[G["E"]]
        pos = (Yl == 1)
        for s in ("human1", "human2"):
            addr = set(gf.index[gf[f"addressable_{s}"]])
            inside = pos.loc[:, pos.columns.isin(list(addr))].to_numpy().sum()
            tot = pos.to_numpy().sum()
            # without CE_fit (LOLO)
            keep = w1.ce_mask(Yl.index, Yl.columns, G["line_ce"][lab]).to_numpy()
            pos2 = pos.to_numpy() & keep
            inside2 = pos2[:, Yl.columns.isin(list(addr))].sum()
            rc.append({"label": lab, "parent": s, "n_addressable_U": len(addr), "positive_pairs": int(tot),
                       "frac_positives_in_addressable": inside / tot, "recall_ceiling": inside / tot,
                       "frac_outside": 1 - inside / tot, "positive_pairs_noCE": int(pos2.sum()),
                       "recall_ceiling_noCE": inside2 / pos2.sum()})
    rcd = pd.DataFrame(rc)
    atomic_write_tsv(rcd, OUT_DIR / "train_addressable.tsv")
    out["addressable"] = rcd.to_dict(orient="records")
    # ---- per-gene error table (ext-PI SF1), training only, L_prob
    per_gene_table(prim)
    return out


def per_gene_table(prim):
    from nextgem import labels, metrics

    Y = G["Y"]["L_prob"].loc[G["E"]]
    gen = G["gen"]
    frac = Y.mean()
    npos = Y.sum()
    n = Y.notna().sum()
    ce = labels.ce_fit(G["train_long"], "L_prob")  # CE_fit from all training lines
    panels = G["panels"]
    pan = panels.groupby("entrez").apply(lambda d: ";".join(sorted(set(d.section.astype(str) + ":" + d.subpanel.astype(str).str.split(":").str[0]))),
                                         include_groups=False)
    sym = G["U"].set_index("entrez")
    sy = {}
    for p in ("human1", "human2"):
        mp_ = pd.read_csv(Path(__file__).resolve().parents[1] / "results" / "baseline" / f"{p}_gene_mapping.tsv", sep="\t", dtype=str)
        sy.update({int(e): s for e, s in zip(mp_.entrez, mp_.symbol) if isinstance(e, str)})
    df = pd.DataFrame(index=pd.Index(G["genes"], name="entrez"))
    df["symbol"] = [sy.get(e, "") for e in df.index]
    for s in ("H1", "H2"):
        df[f"in_{s}"] = gen[f"in_{s}"]
        df[f"ratio_r_{s}"] = gen[f"ratio_{s}"]
        df[f"fba_call_{s}"] = gen[f"fba_call_{s}"]
        if f"p57_{s}" in gen:
            df[f"p57_{s}"] = gen[f"p57_{s}"]
    df["n_lines"] = n
    df["n_pos"] = npos
    df["frac_pos_train"] = frac
    df["ce_fit_train"] = df.index.isin(list(ce))
    for s in ("H1", "H2"):
        called = df[f"fba_call_{s}"]
        df[f"ppv_if_called_{s}"] = np.where(called, frac, np.nan)
        df[f"miss_rate_if_uncalled_{s}"] = np.where(~called, frac, np.nan)
    df["panels"] = pan.reindex(df.index).fillna("")
    for s in ("H1", "H2"):
        called = df[f"fba_call_{s}"]
        cls = np.where(called & (frac >= 0.5), "TP-majority",
                       np.where(called & (frac < 0.5), "FP-majority",
                                np.where(~called & (frac >= 0.5), "FN-majority", "TN-majority")))
        df[f"error_class_{s}"] = cls
    # per-gene tAP for the main methods (training expression lines)
    for m in ("H1", "H2", "B2[L_prob]", "B3", "B6[L_prob]"):
        gt = metrics.per_gene(G["scores"][m].loc[Y.index, Y.columns], Y)
        df[f"tap_gene_{m}"] = gt.tap
    df["e_random_gene"] = gt.e_random
    atomic_write_tsv(df.reset_index(), OUT_DIR / "train_per_gene.tsv")


if __name__ == "__main__":
    raise SystemExit(run_main(main, inputs=[BASELINE / "human1_ko.parquet", BASELINE / "human2_ko.parquet",
        OUT_DIR / "p57_human1.parquet", OUT_DIR / "p57_human2.parquet", TMP_DIR / "expr.parquet", TMP_DIR / "cn.parquet",
        TMP_DIR / "b6_oof.parquet", OUT_DIR / "line_features.tsv", OUT_DIR / "gene_features.tsv",
        PROJECT / "data" / "controls" / "gene_panels.tsv", PROJECT / "results" / "splits" / "lineage_split.tsv"],
        seeds=[11, 23, 37, 41, 53], workers=14, reads_labels=True))
