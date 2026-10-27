"""W4 Stage-A review (independent reviewer): own recomputation of tAP_line, Δ̄_lin vs generic Human1 (G*) and the
nested panel-LOLO selection for F2, F3 and F5, TRAINING TIER ONLY.

Labels: nested selection goes through nextgem.cv.lolo (fit_predict sees labels.fold(k) only; the held-out lineage's
labels are used only inside the harness's score callback). labels.train() is used once, for the full-panel
per-configuration values that are compared with the author's table. Never validation()/test(); nothing under
data/locked/. Writes aggregates only (per configuration and per fold), no gene x line values.

Own code: tie-aware AP (DESIGN §6.1, O(n_g) form), score construction (e = round(1 - ratio_r, 6) over U, U genes
without a parent gene -> 0), estimability (>= 5 positives), lineage-equal means, fold-positivity rule, lineage-
stratified line bootstrap. Not used from the author: w4_analysis, metrics, w1.
"""

from __future__ import annotations

import json
import sys
import time

import numpy as np
import pandas as pd

from nextgem import cv, labels, paths
from nextgem.config import SEEDS

W4 = paths.results_dir() / "w4"
ST, AN, OUT = W4 / "stage_a", W4 / "stage_a_analysis", W4 / "review"


# ------------------------------------------------------------------ own tie-aware AP
def tap(score: np.ndarray, y: np.ndarray) -> float:
    m = ~np.isnan(y)
    s, y = score[m], y[m].astype(int)
    P_tot = y.sum()
    if P_tot == 0:
        return np.nan
    order = np.argsort(-s, kind="mergesort")
    s, y = s[order], y[order]
    # blocks of equal score
    brk = np.flatnonzero(np.diff(s) != 0) + 1
    starts = np.r_[0, brk]
    ends = np.r_[brk, len(s)]
    N = P = 0
    tot = 0.0
    for a, b in zip(starts, ends):
        n_g = b - a
        p_g = int(y[a:b].sum())
        if p_g:
            r = np.arange(1, n_g + 1)
            frac = (r - 1) * (p_g - 1) / (n_g - 1) if n_g > 1 else np.zeros(n_g)
            tot += float(np.sum((p_g / n_g) * (P + 1 + frac) / (N + r)))
        N += n_g
        P += p_g
    return tot / P_tot


def per_line_tap(S: pd.DataFrame, Y: pd.DataFrame) -> pd.Series:
    out = {}
    for l in S.index:
        y = Y.loc[l].to_numpy(dtype=float)
        npos = np.nansum(y)
        out[l] = tap(S.loc[l].to_numpy(dtype=float), y) if npos >= 5 else np.nan
    return pd.Series(out)


def lin_mean(v: pd.Series, lin: pd.Series) -> float:
    v = v.dropna()
    return float(v.groupby(lin.reindex(v.index)).mean().mean())


def boot_ci(v: pd.Series, lin: pd.Series, B=2000, seed=SEEDS[0]):
    v = v.dropna()
    lv = lin.reindex(v.index)
    rng = np.random.default_rng(seed)
    groups = {k: g.to_numpy() for k, g in v.groupby(lv)}
    reps = np.empty(B)
    for b in range(B):
        reps[b] = np.mean([rng.choice(x, size=len(x), replace=True).mean() for x in groups.values()])
    return float(np.quantile(reps, 0.025)), float(np.quantile(reps, 0.975))


def main():
    t0 = time.time()
    split = labels.load_split().set_index("ModelID")
    panel = pd.read_csv(W4 / "stage_a_panel.tsv", sep="\t")
    lines = list(panel.ModelID)
    lin = split.lineage.reindex(lines)
    U = sorted(labels.gene_universe().entrez.astype(int))

    # ---- score frames (own)
    def emap(parent):
        k = pd.read_parquet(paths.baseline_dir() / f"{parent}_ko.parquet", columns=["gene", "entrez", "ratio_r"]).dropna(subset=["entrez"])
        k["entrez"] = k.entrez.astype(int)
        assert not k.entrez.duplicated().any()
        return k.set_index("entrez")
    h2, h1 = emap("human2"), emap("human1")
    gene_of = h2.gene  # entrez -> H2 Ensembl

    def generic(p):
        r = p.ratio_r.reindex(U)
        e = np.round(1.0 - r.fillna(1.0), 6)
        return pd.DataFrame(np.tile(e.to_numpy(), (len(lines), 1)), index=lines, columns=U)

    ko = {}
    for r in ("A1", "A2", "A4", "A5"):
        k = pd.read_parquet(ST / f"ko_{r}.parquet", columns=["build_id", "gene", "ratio_r"])
        k["build_id"] = k.build_id.astype(str)
        for bid, g in k.groupby("build_id"):
            ko[bid] = g.set_index("gene").ratio_r
    gU = gene_of.reindex(U)

    def frame(bids: dict) -> pd.DataFrame:
        rows = {}
        for l, bid in bids.items():
            rr = ko[bid].reindex(gU.to_numpy()).to_numpy(dtype=float)
            e = np.round(1.0 - np.nan_to_num(rr, nan=1.0), 6)
            e[gU.isna().to_numpy()] = 0.0
            rows[l] = e
        return pd.DataFrame(rows, index=U).T.loc[lines]

    cfgs = {}
    for fam in ("F2", "F3", "F4"):
        for rule in ("R1", "R2", "R3"):
            row = "A2" if rule == "R2" else "A1"
            for gen in ("none", "loss"):
                cfgs[f"{fam}|{rule}|{gen}"] = {l: f"{row}|{l}|{fam}|{rule}|rna|real|{gen}" for l in lines}
    gm = pd.read_csv(ST / "a5_grid_map.tsv", sep="\t")
    for (wp, wc, tau), g in gm.groupby(["w_p", "w_c", "tau"]):
        cfgs[f"F5|wp{wp:g}_wc{wc:g}_tau{tau:g}"] = {l: f"{b}|none" for l, b in zip(g.line, g.build_id)}
    S = {c: frame(b) for c, b in cfgs.items()}
    S["H1"], S["H2"] = generic(h1), generic(h2)
    order = {"F2": [f"F2|{r}|{g}" for r in ("R1", "R2", "R3") for g in ("none", "loss")],
             "F3": [f"F3|{r}|{g}" for r in ("R1", "R2", "R3") for g in ("none", "loss")],
             "F4": [f"F4|{r}|{g}" for r in ("R1", "R2", "R3") for g in ("none", "loss")],
             "F5": sorted([c for c in cfgs if c.startswith("F5")],
                          key=lambda c: (float(c.split("wp")[1].split("_")[0]), float(c.split("wc")[1].split("_")[0]),
                                         -float(c.split("tau")[1])))}

    # ---- (a) full-panel per-config values via labels.train() (for comparison with the author's tables)
    tl = labels.train(columns=["L_prob"])
    tl = tl[tl.ModelID.isin(lines)]
    Y = tl.pivot(index="ModelID", columns="entrez", values="L_prob").reindex(index=lines, columns=U).astype(float)
    T = {c: per_line_tap(S[c], Y) for c in S}
    full = {}
    for c in S:
        d = T[c] - T["H1"]
        full[c] = {"tap_lin": lin_mean(T[c], lin), "delta_vs_H1": lin_mean(d, lin), "delta_vs_H2": lin_mean(T[c] - T["H2"], lin)}
    auth_pl = pd.read_csv(AN / "train_per_line.tsv", sep="\t", usecols=["config", "ModelID", "tap"])
    amap = {"H1": "generic_H1", "H2": "generic_H2"}
    for c in cfgs:
        fam = c.split("|")[0]
        if fam == "F5":
            amap[c] = c.replace("F5|", "F5|wr1_")
        else:
            _, rule, gen = c.split("|")
            amap[c] = f"{fam}|{rule}|rna|real|{gen}"
    maxdiff = {}
    for c, a in amap.items():
        av = auth_pl[auth_pl.config == a].set_index("ModelID").tap.reindex(lines)
        maxdiff[c] = float((av - T[c]).abs().max()) if av.notna().any() else None
    del Y, tl  # no label frame kept beyond this point

    # ---- (b) nested panel LOLO through the cv.lolo harness (labels.fold(k) for the choice)
    pl_lins = sorted(lin.unique())

    def make_fit(fam):
        def fit_predict(fit, heldout_ids):
            fl = [l for l in lines if l in set(fit.ModelID)]
            f = fit[fit.ModelID.isin(fl)]
            Yf = f.pivot(index="ModelID", columns="entrez", values="L_prob").reindex(index=fl, columns=U).astype(float)
            ref = per_line_tap(S["H1"].loc[fl], Yf)
            best, bv = None, -np.inf
            for c in order[fam]:
                v = lin_mean(per_line_tap(S[c].loc[fl], Yf) - ref, lin)
                if v > bv + 1e-15:
                    best, bv = c, v
            hid = [l for l in heldout_ids if l in set(lines)]
            return {"chosen": best, "inner_value": bv, "ids": hid}
        return fit_predict

    def score(pred, held):
        ids = pred["ids"]
        h = held[held.ModelID.isin(ids)]
        Yh = h.pivot(index="ModelID", columns="entrez", values="L_prob").reindex(index=ids, columns=U).astype(float)
        d = per_line_tap(S[pred["chosen"]].loc[ids], Yh) - per_line_tap(S["H1"].loc[ids], Yh)
        return {"chosen": pred["chosen"], "inner_value": pred["inner_value"],
                "fold_mean": float(d.mean()) if d.notna().any() else np.nan, "n_est": int(d.notna().sum()),
                "per_line_delta": d.to_dict()}

    nested = {}
    for fam in ("F2", "F3", "F4", "F5"):
        res = cv.lolo(make_fit(fam), score, columns=["L_prob"], folds=pl_lins)
        outer = pd.Series({l: v for r in res.values() for l, v in r["per_line_delta"].items()})
        lo, hi = boot_ci(outer, lin)
        pos = sum(1 for r in res.values() if r["n_est"] >= 2 and r["fold_mean"] > 0)
        dl = lin_mean(outer, lin)
        nested[fam] = {"outer_delta_lin_vs_H1": dl, "boot_lo": lo, "boot_hi": hi, "folds_positive": pos,
                       "tg_pass": bool(dl >= 0.005 and lo > 0 and pos >= 8),
                       "chosen_per_fold": {k: r["chosen"] for k, r in res.items()},
                       "fold_means": {k: r["fold_mean"] for k, r in res.items()}}
    # author's numbers
    sel = json.loads((AN / "train_selection.json").read_text())["families"]
    cmp_ = {}
    for fam in ("F2", "F3", "F4", "F5"):
        a = sel[fam]["vs_H1"]
        own = nested[fam]
        own_ch = own["chosen_per_fold"]
        a_ch = a["chosen_per_fold"]
        norm = lambda c: c.replace("|rna|real|", "|").replace("F5|wr1_", "F5|") if c else c
        cmp_[fam] = {"author_outer": a["outer_delta_lin"], "own_outer": own["outer_delta_lin_vs_H1"],
                     "abs_diff": abs(a["outer_delta_lin"] - own["outer_delta_lin_vs_H1"]),
                     "author_ci": [a["lo"], a["hi"]], "own_ci": [own["boot_lo"], own["boot_hi"]],
                     "author_folds_pos": a["folds_positive"], "own_folds_pos": own["folds_positive"],
                     "author_tg": a["tg_pass"], "own_tg": own["tg_pass"],
                     "same_choice_every_fold": all(norm(a_ch.get(k)) == own_ch.get(k) for k in own_ch)}
    out = {"n_panel_lines": len(lines), "n_U": len(U),
           "full_panel": {c: full[c] for c in ["H1", "H2"] + order["F2"] + order["F3"] + order["F4"]},
           "F5_best_full_panel": max(order["F5"], key=lambda c: full[c]["delta_vs_H1"]),
           "per_line_tap_max_abs_diff_vs_author": maxdiff, "nested": nested, "comparison": cmp_,
           "wall_s": round(time.time() - t0)}
    for fam in nested:  # drop per-line details from the written output
        nested[fam].pop("per_line_delta", None)
    (OUT / "r13_stagea_train.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({"comparison": cmp_, "max_per_line_diff": max(v for v in maxdiff.values() if v is not None),
                      "n_configs_compared": sum(v is not None for v in maxdiff.values()),
                      "F5_best_full_panel": out["F5_best_full_panel"],
                      "H1_tap_lin": full["H1"]["tap_lin"], "H2_tap_lin": full["H2"]["tap_lin"],
                      "chosen": {f: nested[f]["chosen_per_fold"] for f in nested}}, indent=1, default=float))


if __name__ == "__main__":
    sys.exit(main())
