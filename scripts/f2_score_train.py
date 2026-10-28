"""F2 training recompute and F1 non-author recompute (A15 post hoc; RESEARCH_TREE D17, D18). Reads TRAINING labels
only, through nextgem.labels.train() and nextgem.labels.fold(k); never raw label files, validation or test labels.

Own scoring code (scripts/f2_common.py: tAP, lineage-stratified bootstrap, fold rule); it does not import the
W1/M6/W7/F1 scoring code. Set-up as in M6/W7: the 686 non-rare training lines with an expression profile, 12 LOLO
lineages, L_prob, all U genes; paired lineage-stratified line bootstrap B = 2,000 from SEEDS[0]; a fold is positive
iff it has >= 5 lines and a mean Δ > 0. "noCE" removes, for each line, CE_fit computed on labels.fold(k) of the
line's own lineage k (training minus that lineage; nested LOLO for a no-fit scorer reduces to this fold-scoped
mask). FBA scores have no fitted parameter, so the all-genes Δ is the same with or without nesting.

Contrasts (fixed here, before the first label read):
  F1_vs_P1     curated P1 − P1      (non-author recompute of the F1 result: −0.00063 [−0.00071, −0.00055], 0/12)
  F2_vs_Gstar  F2 − G* (H1/M0)      (W7 r4: +0.0304 [0.0294, 0.0313], 12/12)
  F2_vs_P1     F2 − P1              (W7 r4: −0.0230 [−0.0236, −0.0224], 0/12)
  P1_vs_Gstar  P1 − G*              (control: M6 +0.0534 [0.0523, 0.0545], 12/12)
Agreement: the point estimate equals the reported value at its reported precision, the folds count is equal, and
each CI end is within 0.0002 (the bootstrap draws differ from the authors' generator).

Verification before any contrast: P1's line x U scores built here (M6 sweep_e + line assignment) equal M6's
scores_train column bitwise; with real labels, G*'s and P1's per-line tAP equal results/m6/train/per_line.tsv
(<= 1e-6; that file is rounded to 6 significant digits).

Outputs: results/f2/train/summary.json, results/f2/train/contrasts.tsv (aggregates only; no per-line table).
Usage: .venv/bin/python scripts/f2_score_train.py [--synthetic-labels] [--allow-dirty]
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
import f2_common as fc  # noqa: E402

LABEL = "L_prob"
CONTRASTS = (("F1_vs_P1", "F1", "P1"), ("F2_vs_Gstar", "F2", "Gstar"), ("F2_vs_P1", "F2", "P1"),
             ("P1_vs_Gstar", "P1", "Gstar"))
CI_TOL = 2e-4
SYNTH_RATE, SYNTH_NAN = 0.12, 0.02


def decimals(x: float) -> int:
    s = f"{x:.10f}".rstrip("0")
    return len(s.split(".")[1]) if "." in s else 0


def agreement(name: str, got: dict) -> dict:
    rep = fc.REPORTED[name]
    k = decimals(rep["delta_lin"])
    est_ok = round(got["delta_lin"], k) == round(rep["delta_lin"], k)
    ci_ok = all(abs(a - b) <= CI_TOL for a, b in zip(got["ci"], rep["ci"]))
    return {"reported": rep, "estimate_agrees": bool(est_ok), "ci_within_tol": bool(ci_ok),
            "folds_agree": got["folds_positive"] == rep["folds_positive"],
            "agrees": bool(est_ok and ci_ok and got["folds_positive"] == rep["folds_positive"])}


def score_inputs():
    """Label-free: training lines, lineages, folds and the four lines x U score matrices."""
    from nextgem import labels, w1

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    E = sorted(l for l in split.index[(split.tier == "train") & (split.stratum != "rare")]
               if bool(lf.loc[l, "has_expression_matrix"]))
    lineage = split.lineage
    folds = list(labels.folds())
    if sorted(lineage.reindex(E).unique()) != sorted(folds):
        raise AssertionError("training lines' lineages differ from the LOLO folds")
    U, gmap = fc.universe(), fc.gene_map()
    sweep_of = fc.line_m2_sweep()
    p1e = fc.p1_sweep_e(U)
    S = {"Gstar": pd.DataFrame(np.tile(fc.generic_e(U, "human1").to_numpy(float), (len(E), 1)), index=E, columns=U),
         "P1": fc.line_matrix(p1e, sweep_of, E)}
    f1e, f1info = fc.f1_sweep_e(U, p1e, gmap)
    S["F1"] = fc.line_matrix(f1e, sweep_of, E)
    S["F2"] = fc.line_matrix(fc.f2_sweep_e(U, gmap), sweep_of, E)
    # P1 and G* against M6's training scores (bitwise)
    m6 = pd.read_parquet(fc.M6_TRAIN_SCORES, columns=["ModelID", "entrez", "H1", fc.P1_METHOD])
    m6 = m6[m6.ModelID.isin(E)]
    ver = {}
    for k, c in (("P1", fc.P1_METHOD), ("Gstar", "H1")):
        W = m6.pivot(index="ModelID", columns="entrez", values=c).reindex(index=E, columns=U)
        ver[f"{k}_scores_equal_M6_bitwise"] = bool(np.array_equal(W.to_numpy(), S[k].to_numpy()))
        if not ver[f"{k}_scores_equal_M6_bitwise"]:
            raise AssertionError(f"{k}: training scores differ from M6 scores_train")
    ver.update({f"F1_{k}": v for k, v in f1info.items()})
    ver["pairs_F1_ne_P1"] = int((S["F1"].to_numpy() != S["P1"].to_numpy()).sum())
    ver["pairs_F2_ne_P1"] = int((S["F2"].to_numpy() != S["P1"].to_numpy()).sum())
    return E, lineage, folds, U, S, ver


def real_labels(E, U, folds, lineage):
    from nextgem import labels

    tr = labels.train(columns=[LABEL])
    Y = tr.pivot(index="ModelID", columns="entrez", values=LABEL).reindex(index=E, columns=U).astype(float)
    ce = {}
    for k in folds:
        fit = labels.fold(k, columns=[LABEL])
        if set(fit.lineage) & {k}:
            raise AssertionError(f"fold {k}: held-out lineage in the fitting labels")
        frac = fit[["entrez", LABEL]].dropna().astype({LABEL: float}).groupby("entrez")[LABEL].mean()
        ce[k] = set(frac.index[frac >= 0.9].astype(int))
    return Y, ce


def synthetic_labels(E, U, folds, lineage, seed: int = 101):
    """Bernoulli(0.12) labels with 2% NaN on the training lines x U; no label file is opened."""
    rng = np.random.default_rng(seed)
    y = (rng.random((len(E), len(U))) < SYNTH_RATE).astype(float)
    y[rng.random(y.shape) < SYNTH_NAN] = np.nan
    Y = pd.DataFrame(y, index=E, columns=U)
    ce = {}
    for k in folds:
        fit = Y[lineage.reindex(E).to_numpy() != k]
        frac = fit.mean(axis=0, skipna=True)
        ce[k] = set(frac.index[frac >= 0.9].astype(int))
    return Y, ce


def no_ce(Y: pd.DataFrame, ce: dict, lineage) -> pd.DataFrame:
    lin = lineage.reindex(Y.index)
    keep = np.vstack([~np.isin(np.asarray(Y.columns), list(ce[lin[l]])) for l in Y.index])
    return Y.where(keep)


def run(a) -> dict:
    from nextgem.config import SEEDS

    t0 = time.time()
    E, lineage, folds, U, S, ver = score_inputs()
    Y, ce = (synthetic_labels if a.synthetic_labels else real_labels)(E, U, folds, lineage)
    Ys = {"all": Y, "noCE": no_ce(Y, ce, lineage)}
    LT = {(m, g): fc.per_line_tap(S[m], Ys[g]) for m in S for g in Ys}
    if not a.synthetic_labels:
        pl = pd.read_csv(fc.M6_PER_LINE, sep="\t")
        for m, cell in (("Gstar", "H1"), ("P1", "H1+Bloo/M2")):
            for g in ("all", "noCE"):
                ref = pl[(pl.cell == cell) & (pl.label_cell == g)].set_index("ModelID").reindex(E)
                mine = LT[(m, g)]
                dev = float(np.nanmax(np.abs(ref.tap.to_numpy(float) - mine.tap.to_numpy(float))))
                same = bool((ref.estimable.astype(bool).to_numpy() == mine.estimable.to_numpy()).all())
                ver[f"{m}_{g}_tap_line_max_dev_vs_M6"] = dev
                if dev > 1e-6 or not same:
                    raise AssertionError(f"{m}/{g}: per-line tAP differs from M6 per_line.tsv ({dev})")
    rows, out = [], {}
    for name, x, y in CONTRASTS:
        for g in Ys:
            r = fc.contrast(LT[(x, g)], LT[(y, g)], lineage, folds, B=2000, seed=SEEDS[0])
            r.update(contrast=name, geneset=g)
            rows.append(r)
            if g == "all":
                out[name] = {**r, **({} if a.synthetic_labels else {"check": agreement(name, r)})}
            else:
                out[name]["noCE"] = {k: r[k] for k in ("delta_lin", "ci", "folds_positive")}
    tl = {m: {g: fc.lineage_equal_mean(LT[(m, g)].tap, lineage) for g in Ys} for m in S}
    n_est = int(LT[("P1", "all")].estimable.sum())
    return {"synthetic_labels": bool(a.synthetic_labels), "label": LABEL, "n_lines": len(E),
            "n_lines_estimable": n_est, "n_lineages": int(lineage.reindex(E).nunique()), "n_genes_U": len(U),
            "B": 2000, "seed": SEEDS[0], "ci_tolerance": CI_TOL, "verification": ver, "tap_lin": tl,
            "contrasts": out, "all_agree": None if a.synthetic_labels else all(
                v["check"]["agrees"] for v in out.values()), "wall_s": round(time.time() - t0, 1),
            "_rows": rows}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--synthetic-labels", action="store_true", help="Bernoulli(0.12) labels; no label is read")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    from nextgem.config import SEEDS

    outdir = fc.TRDIR if not a.synthetic_labels else fc.F2DIR / "tmp" / "train_synthetic"
    inputs = [fc.W7_SWEEPS, fc.M6_SWEEP_E, fc.M6_ASSIGN, fc.M6_TRAIN_SCORES, fc.M6_PER_LINE, fc.F1_DELTA,
              fc.W2_LINE_SWEEPS]
    with fc.provenance(inputs=inputs, seeds=[SEEDS[0]], allow_dirty=a.allow_dirty) as prov:
        res = run(a)
        rows = pd.DataFrame(res.pop("_rows"))
        rows["ci_lo"], rows["ci_hi"] = rows.ci.str[0], rows.ci.str[1]
        rows = rows.drop(columns="ci")[["contrast", "geneset", "delta_lin", "ci_lo", "ci_hi", "se_boot", "delta_line",
                                         "folds_positive", "n_folds", "n_lines", "n_lineages"]]
        outs = [fc.write_json(res, outdir / "summary.json"), fc.write_tsv(rows, outdir / "contrasts.tsv")]
        prov.output(*outs)
        prov.note(synthetic_labels=res["synthetic_labels"], all_agree=res["all_agree"])
        print(rows.to_string(index=False))
        print(json.dumps({k: res[k] for k in ("synthetic_labels", "n_lines", "n_lines_estimable", "n_lineages",
                                              "verification", "all_agree")}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
