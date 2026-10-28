"""M6 integration factorial on the training tier (DESIGN §6, §6.8, §9 TG, §10.1, §17 A11-A12; RESEARCH_TREE D10).

Labels only through ``nextgem.labels.train()`` and ``labels.fold(k)`` (via ``nextgem.cv.lolo``, for fold-scoped
CE_fit). No validation or test label is read; ``--synthetic-labels`` replaces every label by Bernoulli(0.12) draws and
opens no label file (code dry run; outputs under results/m6/tmp/ only).

Lines: the W1 primary training set (eligible lineages, expression profile): 686 lines, 12 lineages. Every line is
scored (M2(f=1) on lines outside the attribution set = W2's integrated-pipeline medium).
Cells (m6_common.CELLS): G* = H1/M0 (= W1 generic H1, asserted); candidates (1) H1+B-loo/M0, (2) H1/M2(f=1),
(3) H1+B-loo/M2(f=1); H2 counterparts descriptive.
Cells of the tAP grid: L_prob x {all U genes, minus CE_fit (fold-scoped, LOLO), expressed pairs (log2(TPM+1) >= 1)}.

Reuses ``nextgem.w1`` (method_summary, delta_summary, ce_mask, apply_stratum, parent_gene_set, wide, fit_components,
generic_table), ``nextgem.metrics`` (per_line, paired_line_delta, bootstrap_mean, binary_per_line, summarise_binary,
no_information_scores) and ``nextgem.w2_score.folds_positive`` (TG fold rule). B = 2,000, SEEDS[0].

Outputs (results/m6/train/): cells.tsv, contrasts.tsv, interaction.tsv, binary.tsv, binary_paired.tsv,
per_line.tsv, verification.tsv, summary.json.

§6.8 binary endpoints (L-ge06, parent genes with a CRISPR column, expression lines; MCC, F1, precision, recall,
accuracy):
  P57 / PGR (paper protocol): each task with its own inputs, so the GR task is growth on Ham's under the cell's
    biomass; the 56 non-GR tasks close every boundary and never carry biomass flux, so their failure sets are the W3
    capacity-matrix sets for every cell. Consequence: P57/PGR do not depend on the medium (cell (2) = G*, cell (3) =
    cell (1) under these rules, by construction).
  P57|cell-medium, PGR|cell-medium (ours, sensitivity): the GR task evaluated on the line's own cell medium
    (KO growth = ratio_r x WT of the line's sweep < 1, the GR task's LB).
  FBA<0.01 (§6.1 fixed cutoff; L-ge06 and L_prob): ratio_r < 0.01 of the line's sweep.
GR calls from a sweep: KO growth < 1 (the GR task's biomass LB); on T0|M0 this equals the W3 GR task column and W1's
PGR, and P57(T0) equals W1's P57 (asserted).

Usage:
  .venv/bin/python scripts/m6_score_train.py --synthetic-labels      # dry run, no label file opened
  .venv/bin/python scripts/m6_score_train.py                         # training labels
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
import m6_common as mc  # noqa: E402
from m6_scores import W2, W3, e_matrices, sweep_ratios  # noqa: E402

LABEL_COLS = ["gene_effect", "L_prob", "L_ge06"]
YCELLS = ("L_prob", "L_prob_noCE", "L_prob_expr")
YNAME = {"L_prob": "all", "L_prob_noCE": "noCE", "L_prob_expr": "expressed_pairs"}
SYNTH_RATE = 0.12
BIN_METRICS = ("mcc", "f1", "precision", "recall", "accuracy")
PAIRWISE = [("H1+Bloo/M0", "H1"), ("H1/M2", "H1"), ("H1+Bloo/M2", "H1"),
            ("H1/M2", "H1+Bloo/M0"), ("H1+Bloo/M2", "H1+Bloo/M0"), ("H1+Bloo/M2", "H1/M2"),
            ("H2", "H1"), ("H2+Bloo/M0", "H1"), ("H2/M2", "H1"), ("H2+Bloo/M2", "H1"),
            ("H2+Bloo/M0", "H2"), ("H2/M2", "H2"), ("H2+Bloo/M2", "H2"),
            ("H2/M2", "H2+Bloo/M0"), ("H2+Bloo/M2", "H2+Bloo/M0"), ("H2+Bloo/M2", "H2/M2")]
G: dict = {}


# ============================================================================= inputs
def load(synthetic: bool) -> None:
    from nextgem import cv, labels, metrics, w1
    from nextgem.config import SEEDS

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    Udf = labels.gene_universe()
    U = list(Udf.entrez.astype(int))
    gen = w1.generic_table(Udf.entrez)
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    E = sorted(l for l in elig.index if lf.loc[l, "has_expression_matrix"])
    lineage = split.lineage
    ls = pd.read_csv(W2 / "controls" / "line_sweeps.tsv", sep="\t")
    Esw, emap, meta = e_matrices(U, ls)
    R, _ = sweep_ratios()
    S, failed = {}, {}
    for c in mc.CELLS:
        M = mc.cell_matrix(c, ls, Esw[c.parent], E)
        bad = list(M.index[M.isna().any(axis=1)])
        failed[c.key] = bad
        S[c.key] = metrics.no_information_scores(M.fillna(0.0), bad, "line")  # §11.4, identical in every arm
    checks = []
    for s in ("H1", "H2"):
        v = gen[f"e_{s}"].reindex(U).to_numpy()
        ok = bool(np.array_equal(S[s].to_numpy(), np.tile(v, (len(E), 1))))
        checks.append({"check": f"{s}/M0 scores equal the W1 generic e_{s} on the 686 lines", "pass": ok})
    if not all(c["pass"] for c in checks):
        raise RuntimeError(checks)
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
    G.update(split=split, lf=lf, U=U, gen=gen, E=E, lineage=lineage, ls=ls, Esw=Esw, emap=emap, R=R, S=S,
             failed=failed, Y=Y, folds=folds, checks=checks, meta=meta)


# ============================================================================= tAP
def per_line_tables() -> dict:
    from nextgem import metrics

    return {(c.key, y): metrics.per_line(G["S"][c.key], G["Y"][y]) for c in mc.CELLS for y in YCELLS}


def cell_rows(LT: dict) -> pd.DataFrame:
    from nextgem import metrics, w1

    rows = []
    for c in mc.CELLS:
        row, _, _ = w1.method_summary(G["S"][c.key], G["Y"]["L_prob"], G["lineage"], gene_unit=True)
        r = {"cell": c.key, "method": c.method, "parent": c.parent, "variant": c.variant, "step": c.step, "role": c.role,
             "label": mc.CANDIDATE_LABEL.get(c.key, ""), "n_failed_lines": len(G["failed"][c.key])}
        r.update({k: row[k] for k in ("n_lines", "n_lines_estimable", "n_lineages", "tap_lin", "tap_lin_lo", "tap_lin_hi",
                                      "e_random_lin", "lift_lin", "tap_line", "m_delta_ap_gene", "n_genes_estimable",
                                      "pooled_tap", "pooled_lift")})
        for y in ("L_prob_noCE", "L_prob_expr"):
            lt = LT[(c.key, y)]
            est = lt[lt.estimable]
            r[f"tap_lin_{YNAME[y]}"] = metrics.lineage_equal_mean(est.tap, G["lineage"])
        rows.append(r)
    return pd.DataFrame(rows)


def contrast_rows(LT: dict) -> pd.DataFrame:
    from nextgem import metrics, w1, w2_score as ws
    from nextgem.config import SEEDS

    rows = []
    for a, b in PAIRWISE:
        for y in YCELLS:
            d = metrics.paired_line_delta(LT[(a, y)], LT[(b, y)])
            r = {"cell": a, "reference": b, "label_cell": YNAME[y], "n_lines": int(d.notna().sum())}
            r.update(w1.delta_summary(LT[(a, y)], LT[(b, y)], G["lineage"]))
            bt = metrics.bootstrap_mean(d, G["lineage"].reindex(d.index), B=metrics.B_OTHER, seed=SEEDS[0])
            r["se_boot"] = float(np.std(bt.replicates, ddof=1))
            r.update(ws.folds_positive(d, G["lineage"], G["folds"]))
            r["tg_pass"] = mc.tg_pass(r["delta_lin"], r["delta_lo"], r["folds_positive"]) if y == "L_prob" else None
            rows.append(r)
    return pd.DataFrame(rows)


def interaction_rows(LT: dict) -> pd.DataFrame:
    from nextgem import metrics, w1, w2_score as ws
    from nextgem.config import SEEDS

    rows = []
    for s in ("H1", "H2"):
        keys = [s, f"{s}+Bloo/M0", f"{s}/M2", f"{s}+Bloo/M2"]
        for y in YCELLS:
            taps = [LT[(k, y)].tap.where(LT[(k, y)].estimable) for k in keys]
            I = mc.interaction(*taps)
            bt = metrics.bootstrap_mean(I, G["lineage"].reindex(I.index), B=metrics.B_OTHER, seed=SEEDS[0])
            dd = {k: w1.delta_summary(LT[(k, y)], LT[(s, y)], G["lineage"])["delta_lin"] for k in keys[1:]}
            r = {"parent": s, "label_cell": YNAME[y], "n_lines": int(I.size),
                 "interaction_lin": bt.estimate, "lo": bt.ci_low, "hi": bt.ci_high,
                 "interaction_line": metrics.line_weighted_mean(I),
                 "delta_bloo_M0": dd[keys[1]], "delta_M2": dd[keys[2]], "delta_joint": dd[keys[3]],
                 "additive_prediction": dd[keys[1]] + dd[keys[2]]}
            r.update(ws.folds_positive(I, G["lineage"], G["folds"]))
            r["folds_negative"] = sum(1 for v in json.loads(r["fold_means"]).values() if v is not None and v < 0)
            rows.append(r)
    return pd.DataFrame(rows)


def training_one_se(con: pd.DataFrame, LT: dict) -> dict:
    """Descriptive only (selection happens in the validation look): the paired one-SE rule applied to the training
    Δ̄_lin vs G* of the TG-passing candidates, in D10 order."""
    from nextgem import metrics
    from nextgem.config import SEEDS

    prim = con[(con.label_cell == "all") & (con.reference == "H1")].set_index("cell")
    cands = [c.key for c in mc.CANDIDATES]
    el = [c for c in cands if bool(prim.loc[c, "tg_pass"])]
    if not el:
        return {"eligible_tg": [], "pick": None}
    best = max(el, key=lambda c: prim.loc[c, "delta_lin"])
    out = {"eligible_tg": el, "best": best, "rows": []}
    pick = None
    for c in el:
        d = metrics.paired_line_delta(LT[(best, "L_prob")], LT[(c, "L_prob")])
        bt = metrics.bootstrap_mean(d, G["lineage"].reindex(d.index), B=metrics.B_OTHER, seed=SEEDS[0])
        se = float(np.std(bt.replicates, ddof=1))
        within = bool(bt.estimate <= se)
        out["rows"].append({"cell": c, "diff_best_minus_cell": bt.estimate, "se": se, "within_1se": within})
        if within and pick is None:
            pick = c
    out["pick"] = pick
    return out


# ============================================================================= §6.8 binary endpoints
def task_sets() -> dict:
    """parent -> {'nongr': set(model genes), 'gr_task': set(model genes)} from the W3 Essential capacity matrix."""
    out = {}
    for p in mc.PARENTS:
        long = pd.read_parquet(W3 / "tasks" / f"essential_{p}.parquet", columns=["task", "seed", "gene", "feasible"])
        long = long[long.seed == -1]
        meta = pd.read_csv(W3 / "tasks" / f"essential_{p}_tasks.tsv", sep="\t")
        gr = meta.task[meta.category == "GR"]
        assert len(gr) == 1
        gr = gr.iloc[0]
        nongr = set(long.gene[(~long.feasible) & (long.task != gr)])
        out[p] = {"nongr": nongr, "gr_task": set(long.gene[(~long.feasible) & (long.task == gr)]), "n_tasks": long.task.nunique()}
    return out


def to_U_bool(genes: set, emap: pd.Series, U) -> pd.Series:
    s = pd.Series(emap.index.isin(list(genes)), index=emap.to_numpy())
    return s.reindex(U, fill_value=False).astype(bool)


def binary_calls() -> tuple[dict, list]:
    """{(cell, rule): bool frame (lines x U)} and verification rows."""
    cat = pd.read_csv(mc.SWDIR / "catalog.tsv", sep="\t").set_index("key")
    ts = task_sets()
    U, E, gen = G["U"], G["E"], G["gen"]
    checks, calls = [], {}
    gr_of_key = {}
    for p in mc.PARENTS:
        R = G["R"][p]
        for k in R.index:
            gr = mc.gr_fail(R.loc[k], float(cat.loc[k, "wt_growth"]))
            gr_of_key[k] = to_U_bool(set(gr.index[gr.to_numpy()]), G["emap"][p], U)
        s = mc.SHORT[p]
        m0 = sorted(set(mc.line_sweep_ids(G["ls"], p, mc.M0)))[0]
        gr_t0 = mc.gr_fail(R.loc[mc.sweep_key(mc.T0, m0)], float(cat.loc[mc.sweep_key(mc.T0, m0), "wt_growth"]))
        gr_t0 = set(gr_t0.index[gr_t0.to_numpy()])
        checks.append({"check": f"{s}: GR from the T0|M0 sweep (KO growth < 1) equals the W3 GR task column",
                       "pass": gr_t0 == ts[p]["gr_task"]})
        nongr_U = to_U_bool(ts[p]["nongr"], G["emap"][p], U)
        p57_t0 = nongr_U | gr_of_key[mc.sweep_key(mc.T0, m0)]
        checks.append({"check": f"{s}: P57(T0) equals W1 p57_{s} on U", "pass": bool((p57_t0 == gen[f"p57_{s}"].reindex(U)).all())})
        checks.append({"check": f"{s}: PGR(T0) equals W1 pgr_{s} on U",
                       "pass": bool((gr_of_key[mc.sweep_key(mc.T0, m0)] == gen[f"pgr_{s}"].reindex(U)).all())})
        G.setdefault("nongr_U", {})[p] = nongr_U
        G.setdefault("m0_sid", {})[p] = m0
    for c in mc.CELLS:
        p = c.parent
        lit = gr_of_key[mc.sweep_key(c.variant, G["m0_sid"][p])]  # GR task = Ham's inputs under the cell's biomass
        const = lambda v: pd.DataFrame(np.tile(v.to_numpy(bool), (len(E), 1)), index=E, columns=U)  # noqa: E731
        calls[(c.key, "P57")] = const(G["nongr_U"][p] | lit)
        calls[(c.key, "PGR")] = const(lit)
        keys = mc.cell_sweep_keys(c, G["ls"]).reindex(E)
        GRm = pd.DataFrame(np.vstack([gr_of_key[k].to_numpy(bool) for k in keys]), index=E, columns=U)
        calls[(c.key, "PGR|cell-medium")] = GRm
        calls[(c.key, "P57|cell-medium")] = GRm | const(G["nongr_U"][p])
        calls[(c.key, "FBA<0.01")] = G["S"][c.key] >= mc.CALL
    return calls, checks


def binary_tables(calls: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    from nextgem import metrics, w1
    from nextgem.config import SEEDS

    lin = G["lineage"]
    specs = [("Human1/2-paper", "P57", "L_ge06"), ("Human1-paper", "PGR", "L_ge06"),
             ("ours (sensitivity)", "P57|cell-medium", "L_ge06"), ("ours (sensitivity)", "PGR|cell-medium", "L_ge06"),
             ("fixed-cutoff", "FBA<0.01", "L_ge06"), ("fixed-cutoff", "FBA<0.01", "L_prob")]
    rows, PL = [], {}
    for c in mc.CELLS:
        gset = w1.parent_gene_set(G["gen"], mc.SHORT[c.parent])
        for proto, rule, lab in specs:
            T = G["Y"][lab].loc[G["E"], gset]
            C = calls[(c.key, rule)].loc[G["E"], gset]
            pl = metrics.binary_per_line(C, T)
            sm = metrics.summarise_binary(pl, lin.reindex(pl.index))
            sm["ci_lo"], sm["ci_hi"] = np.nan, np.nan
            for m in BIN_METRICS:
                v = pl[m].dropna()
                bt = metrics.bootstrap_mean(v, lin.reindex(v.index), B=metrics.B_OTHER, seed=SEEDS[0])
                i = (sm.metric == m) & (sm.nan_policy == "exclude")
                sm.loc[i, ["ci_lo", "ci_hi"]] = [bt.ci_low, bt.ci_high]
            sm = sm[sm.metric.isin(BIN_METRICS)].assign(cell=c.key, method=c.method, protocol=proto, call_rule=rule,
                                                        label=lab, n_genes=len(gset))
            rows.append(sm)
            PL[(c.key, rule, lab)] = pl
    out = pd.concat(rows, ignore_index=True)
    pr = []
    for c in mc.CELLS:
        if c.key == "H1":
            continue
        for proto, rule, lab in specs:
            a, b = PL[(c.key, rule, lab)], PL[("H1", rule, lab)]
            for m in BIN_METRICS:
                d = (a[m] - b[m]).dropna()
                bt = metrics.bootstrap_mean(d, lin.reindex(d.index), B=metrics.B_OTHER, seed=SEEDS[0])
                pr.append({"cell": c.key, "reference": "H1", "protocol": proto, "call_rule": rule, "label": lab,
                           "metric": m, "delta_lin": bt.estimate, "ci_lo": bt.ci_low, "ci_hi": bt.ci_high,
                           "delta_line": float(d.mean()), "n_lines": len(d)})
    return out, pd.DataFrame(pr), PL


# ============================================================================= main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--synthetic-labels", action="store_true", help="Bernoulli labels, no label file; tmp/ output")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    from nextgem.config import SEEDS

    out_dir = mc.TMP / "score_dry" if a.synthetic_labels else mc.TRDIR
    inputs = [W2 / "controls" / "line_sweeps.tsv", mc.SWDIR / "catalog.tsv", mc.SWDIR / "sweeps.tsv",
              mc.SCDIR / "sweep_e.parquet", mc.SCDIR / "SHA256SUMS.tsv"]
    inputs += [mc.SWDIR / f"ratio_r_{p}.parquet" for p in mc.PARENTS]
    inputs += [W3 / "tasks" / f"essential_{p}.parquet" for p in mc.PARENTS]
    with mc.provenance(inputs=inputs, seeds=[SEEDS[0]], solver="none (analysis)", workers=1,
                       allow_dirty=a.allow_dirty or a.synthetic_labels) as prov:
        prov.note(synthetic_labels=a.synthetic_labels,
                  labels=("synthetic Bernoulli(0.12); no label file" if a.synthetic_labels else
                          "labels.train() + cv.lolo/labels.fold(k) (CE_fit) only; training tier"),
                  label_files_frozen_sha256=mc.frozen_label_hashes())
        outs = run(out_dir)
        prov.output(*outs)
    return 0


def run(out_dir: Path) -> list[Path]:
    t0 = time.time()
    load(synthetic="score_dry" in out_dir.parts)
    print(f"loaded: {len(G['E'])} lines, {len(G['U'])} genes, {time.time() - t0:.0f} s", flush=True)
    LT = per_line_tables()
    cells = cell_rows(LT)
    con = contrast_rows(LT)
    inter = interaction_rows(LT)
    one_se = training_one_se(con, LT)
    print(f"tAP done, {time.time() - t0:.0f} s", flush=True)
    calls, bchecks = binary_calls()
    binary, bpaired, _ = binary_tables(calls)
    print(f"binary done, {time.time() - t0:.0f} s", flush=True)
    checks = G["checks"] + bchecks
    # G* reproduces the W1 training values (tAP and the P57 MCC)
    w1pl = pd.read_csv(mc.PROJECT / "results" / "w1" / "train_per_line.tsv", sep="\t")
    if not ("score_dry" in out_dir.parts):
        for gs, y in (("all", "L_prob"), ("noCE", "L_prob_noCE")):
            w = w1pl[(w1pl.method == "H1") & (w1pl.geneset == gs)].set_index("ModelID")
            lt = LT[("H1", y)]
            same = set(w.index) == set(lt.index) and bool((w.estimable.reindex(lt.index) == lt.estimable).all())
            dif = float(np.nanmax(np.abs(w.tap.reindex(lt.index).to_numpy(float) - lt.tap.to_numpy(float))))
            checks.append({"check": f"G* tAP_line equals W1 train_per_line (H1, {gs})", "pass": same and dif < 1e-6,
                           "max_abs_diff": dif})
        w1b = pd.read_csv(mc.PROJECT / "results" / "w1" / "train_binary.tsv", sep="\t")
        for rule in ("P57", "PGR"):
            x = w1b[(w1b.protocol == "Human1-paper") & (w1b.method == "H1") & (w1b.call_rule == rule) & (w1b.label == "L_ge06")
                    & (w1b.gene_set == "parent") & (w1b.lineset == "expr") & (w1b.metric == "mcc") & (w1b.nan_policy == "exclude")]
            y = binary[(binary.cell == "H1") & (binary.call_rule == rule) & (binary.label == "L_ge06") & (binary.metric == "mcc")
                       & (binary.nan_policy == "exclude")]
            dif = float(abs(x.lineage_equal_mean.iloc[0] - y.lineage_equal_mean.iloc[0]))
            checks.append({"check": f"G* {rule} x L-ge06 MCC equals W1 train_binary", "pass": dif < 1e-6, "max_abs_diff": dif})  # W1 TSV: 6 significant digits
    ver = pd.DataFrame(checks)
    out_dir.mkdir(parents=True, exist_ok=True)
    outs = [mc.atomic_write_tsv(cells, out_dir / "cells.tsv"), mc.atomic_write_tsv(con, out_dir / "contrasts.tsv"),
            mc.atomic_write_tsv(inter, out_dir / "interaction.tsv"), mc.atomic_write_tsv(binary, out_dir / "binary.tsv"),
            mc.atomic_write_tsv(bpaired, out_dir / "binary_paired.tsv"), mc.atomic_write_tsv(ver, out_dir / "verification.tsv")]
    pl = []
    for (k, y), lt in LT.items():
        pl.append(lt[["n", "n_pos", "estimable", "tap", "e_random"]].assign(cell=k, label_cell=YNAME[y]))
    pl = pd.concat(pl).rename_axis("ModelID").reset_index()
    pl["lineage"] = pl.ModelID.map(G["lineage"])
    outs.append(mc.atomic_write_tsv(pl[["cell", "label_cell", "ModelID", "lineage", "n", "n_pos", "estimable", "tap",
                                        "e_random"]], out_dir / "per_line.tsv"))
    prim = con[(con.label_cell == "all") & (con.reference == "H1")].set_index("cell")
    summ = {"n_lines": len(G["E"]), "n_lineages": int(G["lineage"].reindex(G["E"]).nunique()),
            "failed_lines": {k: len(v) for k, v in G["failed"].items()},
            "tg": {c.key: {"label": mc.CANDIDATE_LABEL[c.key], "delta_lin": prim.loc[c.key, "delta_lin"],
                           "ci": [prim.loc[c.key, "delta_lo"], prim.loc[c.key, "delta_hi"]],
                           "folds_positive": int(prim.loc[c.key, "folds_positive"]), "tg_pass": bool(prim.loc[c.key, "tg_pass"])}
                   for c in mc.CANDIDATES},
            "training_one_se_descriptive": one_se, "checks_pass": bool(ver["pass"].all()), "wall_s": round(time.time() - t0, 1)}
    outs.append(mc.atomic_write_text(json.dumps(summ, indent=1, default=float) + "\n", out_dir / "summary.json"))
    if not ver["pass"].all():
        print(ver.to_string(), file=sys.stderr)
        raise RuntimeError("verification failed")
    print(f"written to {out_dir} ({time.time() - t0:.0f} s)", flush=True)
    return outs


if __name__ == "__main__":
    raise SystemExit(main())
