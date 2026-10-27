"""Independent reviewer recomputation of W1 training-tier headline numbers (M2 review).

Training labels only, via nextgem.labels.train() / labels.fold(k). Own tie-aware AP (checked against
brute-force enumeration first), own KO join, own MCC, own CE_fit / E[AP_random], own GPR-based
addressable set. Prints and writes scalars only (training-tier aggregates).

Usage: .venv/bin/python results/w1/review/recompute_w1.py
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / "recompute_w1.json"


# ----------------------------------------------------------------------------- own tie-aware AP
def tap_mine(s: np.ndarray, y: np.ndarray) -> float:
    """Expected AP over uniform random orderings inside tied score blocks (own derivation:
    precision at a positive's rank is linear in the number of block positives above it)."""
    s = np.asarray(s, float)
    y = np.asarray(y, float)
    k = y.sum()
    if k == 0:
        return float("nan")
    order = np.argsort(-s, kind="stable")
    s, y = s[order], y[order]
    _, start, counts = np.unique(-s, return_index=True, return_counts=True)
    N = P = 0.0
    total = 0.0
    for st, ng in zip(start, counts):
        pg = y[st:st + ng].sum()
        if pg > 0:
            r = np.arange(1, ng + 1, dtype=float)
            ex = (r - 1) * (pg - 1) / (ng - 1) if ng > 1 else np.zeros(1)
            total += (pg / ng) * np.sum((P + 1 + ex) / (N + r))
        N += ng
        P += pg
    return total / k


def ap_ordering(y_ordered) -> float:
    hits, acc = 0, 0.0
    for i, v in enumerate(y_ordered, 1):
        if v:
            hits += 1
            acc += hits / i
    return acc / hits


def tap_brute(s, y) -> float:
    s, y = np.asarray(s, float), np.asarray(y, int)
    vals = sorted(set(s), reverse=True)
    blocks = [list(y[s == v]) for v in vals]
    perms = [sorted(set(itertools.permutations(b))) for b in blocks]
    # weight each distinct block permutation by its multiplicity (uniform over item orderings)
    tot, wsum = 0.0, 0
    for combo in itertools.product(*[[(p, 1) for p in ps] for ps in perms]):
        seq = [v for p, _ in combo for v in p]
        tot += ap_ordering(seq)
        wsum += 1
    return tot / wsum  # distinct permutations of a 0/1 block are equally likely under uniform orderings


def e_ap_random(n: int, k: int) -> float:
    if k == 0:
        return float("nan")
    r = np.arange(1, n + 1, dtype=float)
    ex = (r - 1) * (k - 1) / (n - 1) if n > 1 else np.zeros(1)
    return float((k / n) * np.sum((1 + ex) / r) / k)


def self_check() -> dict:
    from nextgem import metrics

    rng = np.random.default_rng(2026)
    worst_mine, worst_lib, n = 0.0, 0.0, 0
    for _ in range(300):
        m = int(rng.integers(3, 9))
        s = rng.integers(0, int(rng.integers(1, 4)) + 1, m).astype(float)
        y = (rng.random(m) < 0.4).astype(int)
        if y.sum() == 0:
            y[0] = 1
        b = tap_brute(s, y)
        worst_mine = max(worst_mine, abs(tap_mine(s, y) - b))
        worst_lib = max(worst_lib, abs(metrics.tap(s, y) - b))
        n += 1
    # larger random cases: mine vs library
    worst_big = 0.0
    for _ in range(200):
        m = int(rng.integers(50, 3000))
        s = np.round(rng.random(m) * int(rng.integers(1, 40)))
        y = (rng.random(m) < rng.uniform(0.02, 0.5)).astype(int)
        if y.sum() == 0:
            continue
        worst_big = max(worst_big, abs(tap_mine(s, y) - metrics.tap(s, y)))
    # constant scorer = E[AP_random]
    worst_null = max(abs(tap_mine(np.zeros(n_), np.r_[np.ones(k_), np.zeros(n_ - k_)]) - e_ap_random(n_, k_))
                     for n_, k_ in [(10, 3), (2866, 350), (2700, 5), (100, 100)])
    return {"brute_cases": n, "max_abs_mine_vs_brute": worst_mine, "max_abs_metrics_vs_brute": worst_lib,
            "max_abs_mine_vs_metrics_large": worst_big, "max_abs_const_vs_e_random": worst_null}


# ----------------------------------------------------------------------------- helpers
def lineage_equal(v: pd.Series, lin: pd.Series) -> float:
    v = v.dropna()
    return float(v.groupby(lin.reindex(v.index)).mean().mean())


def mcc(tp, fp, tn, fn) -> float:
    d = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float("nan") if d == 0 else (tp * tn - fp * fn) / d


def main() -> int:
    from nextgem import labels

    res = {"self_check": self_check()}
    print("self-check", res["self_check"], flush=True)

    split = labels.load_split().set_index("ModelID")
    U = labels.gene_universe()
    genes = sorted(int(g) for g in U.entrez)
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    expr_lines = sorted(elig.index[elig.has_expression.astype(bool)])
    lf = pd.read_csv(ROOT / "results/w1/line_features.tsv", sep="\t").set_index("ModelID")
    expr_lines_author = sorted(l for l in elig.index if lf.loc[l, "has_expression_matrix"])
    res["n_eligible_lines"] = len(elig)
    res["n_expr_lines_split_flag"] = len(expr_lines)
    res["n_expr_lines_author_flag"] = len(expr_lines_author)
    res["expr_flag_sets_equal"] = expr_lines == expr_lines_author
    lin = split.lineage

    # own KO join (Entrez), e = 1 - clip(round(ratio_r, 6), 0, 1); absent -> 0
    score, inpar, p57 = {}, {}, {}
    for p, s in (("human1", "H1"), ("human2", "H2")):
        ko = pd.read_parquet(ROOT / f"results/baseline/{p}_ko.parquet", columns=["entrez", "ratio_r"]).dropna(subset=["entrez"])
        ko["entrez"] = ko.entrez.astype(int)
        assert not ko.entrez.duplicated().any(), p
        r = ko.set_index("entrez").ratio_r
        e = (1 - r.round(6).clip(0, 1)).reindex(genes).fillna(0.0)
        score[s] = e.to_numpy()
        inpar[s] = np.isin(genes, r.index)
        t = pd.read_parquet(ROOT / f"results/w1/p57_{p}.parquet", columns=["entrez", "p57"]).dropna(subset=["entrez"])
        t["entrez"] = t.entrez.astype(int)
        p57[s] = t.drop_duplicates("entrez").set_index("entrez").p57.astype(bool).reindex(genes, fill_value=False).to_numpy()
    res["n_U_in_H1"], res["n_U_in_H2"] = int(inpar["H1"].sum()), int(inpar["H2"].sum())

    tr = labels.train(columns=["L_prob", "L_ge06"])
    Y = tr.pivot(index="ModelID", columns="entrez", values="L_prob").reindex(columns=genes)
    Y06 = tr.pivot(index="ModelID", columns="entrez", values="L_ge06").reindex(columns=genes)

    # ---- generic tAP_line on L_prob, all U, 686 expression lines
    rows = {}
    for l in expr_lines:
        y = Y.loc[l].to_numpy()
        ok = ~np.isnan(y)
        k = y[ok].sum()
        if k < 5:
            rows[l] = (np.nan, np.nan)
            continue
        rows[l] = (tap_mine(score["H1"][ok], y[ok]), tap_mine(score["H2"][ok], y[ok]))
    t = pd.DataFrame(rows, index=["H1", "H2"]).T
    res["n_estimable"] = int(t.H1.notna().sum())
    res["tap_lin_H1"] = lineage_equal(t.H1, lin)
    res["tap_lin_H2"] = lineage_equal(t.H2, lin)
    res["delta_lin_H2_minus_H1"] = lineage_equal(t.H2 - t.H1, lin)
    res["delta_lineages_negative"] = int(((t.H2 - t.H1).groupby(lin.reindex(t.index)).mean() < 0).sum())
    print({k: res[k] for k in ("tap_lin_H1", "tap_lin_H2", "delta_lin_H2_minus_H1")}, flush=True)

    # ---- P57 x L_ge06 MCC, parent gene set, 686 expression lines
    for s in ("H1", "H2"):
        cols = np.array(genes)[inpar[s]]
        call = p57[s][inpar[s]]
        vals = {}
        for l in expr_lines:
            y = Y06.loc[l, cols].to_numpy()
            ok = ~np.isnan(y)
            yy, cc = y[ok].astype(bool), call[ok]
            vals[l] = mcc(int((yy & cc).sum()), int((~yy & cc).sum()), int((~yy & ~cc).sum()), int((yy & ~cc).sum()))
        v = pd.Series(vals)
        res[f"mcc_P57_Lge06_{s}_lineage_equal"] = lineage_equal(v, lin)
        res[f"mcc_P57_Lge06_{s}_line_mean"] = float(v.mean())
        res[f"mcc_P57_Lge06_{s}_n_nan"] = int(v.isna().sum())
        res[f"_mcc_{s}"] = v
    d = res.pop("_mcc_H2") - res.pop("_mcc_H1")
    res["mcc_P57_H2_minus_H1_lineage_equal"] = lineage_equal(d, lin)

    # ---- B2 (L_prob) under LOLO via labels.fold(k): prior from training minus k only
    b2 = {}
    for k in labels.folds():
        fit = labels.fold(k, columns=["L_prob"])
        assert not set(fit.lineage) & {k}
        prior = fit.dropna(subset=["L_prob"]).groupby("entrez").L_prob.mean().reindex(genes).fillna(0.0).to_numpy()
        for l in [x for x in expr_lines if lin[x] == k]:
            y = Y.loc[l].to_numpy()
            ok = ~np.isnan(y)
            if y[ok].sum() >= 5:
                b2[l] = tap_mine(prior[ok], y[ok])
    res["tap_lin_B2_LOLO"] = lineage_equal(pd.Series(b2), lin)

    # ---- SESOI_C1b: CE_fit from all 767 training lines; E[AP_random] on U and U minus CE_fit
    frac = tr.dropna(subset=["L_prob"]).groupby("entrez").L_prob.mean()
    ce = set(frac.index[frac >= 0.9].astype(int))
    res["CE_fit_all_training_size"] = len(ce)
    keep = ~np.isin(genes, list(ce))
    er_all, er_no = {}, {}
    for l in Y.index:
        y = Y.loc[l].to_numpy()
        ok = ~np.isnan(y)
        k_all = int(y[ok].sum())
        k_no = int(y[ok & keep].sum())
        if k_all >= 5:
            er_all[l] = e_ap_random(int(ok.sum()), k_all)
        if k_no >= 5:
            er_no[l] = e_ap_random(int((ok & keep).sum()), k_no)
    er_all, er_no = pd.Series(er_all), pd.Series(er_no)
    res["sesoi_line_mean_all767"] = 0.01 * er_no.mean() / er_all.mean()
    res["mean_e_random_U_all767"] = float(er_all.mean())
    res["mean_e_random_noCE_all767"] = float(er_no.mean())
    res["sesoi_lineage_equal_expr686"] = 0.01 * lineage_equal(er_no.reindex(expr_lines), lin) / lineage_equal(er_all.reindex(expr_lines), lin)
    res["sesoi_line_mean_eligible737"] = 0.01 * er_no.reindex(elig.index).mean() / er_all.reindex(elig.index).mean()

    # ---- addressable set (own GPR evaluation): genes whose KO closes >= 1 reaction; recall ceiling
    import cobra  # noqa: F401

    from nextgem.model import load_model

    for p, s in (("human1", "H1"), ("human2", "H2")):
        m = load_model(p)
        gmap = pd.read_csv(ROOT / f"results/baseline/{p}_gene_mapping.tsv", sep="\t", dtype=str)
        id2e = {r.gene: int(r.entrez) for r in gmap.itertuples() if isinstance(r.entrez, str) and r.entrez.isdigit()}
        addr = set()
        for g in m.genes:
            if any(not rx.gpr.eval({g.id}) for rx in g.reactions):
                if g.id in id2e:
                    addr.add(id2e[g.id])
        addr &= set(genes)
        pos = (Y.loc[expr_lines] == 1).to_numpy()
        inside = pos[:, np.isin(genes, list(addr))].sum()
        res[f"addressable_U_{s}"] = len(addr)
        res[f"positive_pairs_expr686"] = int(pos.sum())
        res[f"frac_pos_outside_addressable_{s}"] = float(1 - inside / pos.sum())
        del m

    OUT.write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
