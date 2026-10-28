"""M4 W3 review (non-author): independent training-tier rescoring of T0 and the two shortlisted B-loo variants, plus a
POST-HOC, label-free random-removal control (NOT a gate).

Own code path: no nextgem.w1 / nextgem.w3 scoring function is imported. Uses nextgem.labels (train() only),
nextgem.metrics (per_line tAP, bootstrap_mean) and nextgem.config.SEEDS. Gene scores come from
results/w3/review/rv_ratio_r_<parent>.parquet (written by rv_ko.py, own KO sweeps), e = 1 - ratio_r, absent gene -> 0.

Lines: training tier, stratum != rare, with an expression profile (results/w1/line_features.tsv). Label: L_prob.
Δ̄_lin: lineage-equal mean of per-line ΔtAP; CI: lineage-stratified bootstrap, B = 2,000, SEEDS[0].
Folds positive: lineages with >= 5 estimable lines and mean Δ > 0 (denominator 12).

Random-removal control (post hoc, label-free): per line and per seed in SEEDS[0:5], remove (e -> 0) the same number of
T0 lethal calls (e >= 0.99) as the variant removes, chosen uniformly at random among the parent's T0-called genes in U
(rng = default_rng([seed, line index])). Variant vs control Δ̄_lin is reported per seed and against the seed-mean tAP.

Outputs (aggregates only; nothing gene x line): results/w3/review/rv_scores.tsv, rv_random_removal.tsv,
rv_changed_calls.tsv (label-free gene lists).
Usage: .venv/bin/python results/w3/review/rv_score.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "w3" / "review"
CALL = 0.99
VARS = ("B-loo:cofactor_pool_biomass", "B-loo:metabolite_pool_biomass")


def gene_scores(parent: str, U: list[int]) -> tuple[dict, pd.Series]:
    ko = pd.read_parquet(ROOT / "results/baseline" / f"{parent}_ko.parquet", columns=["gene", "entrez", "symbol"])
    ko = ko.dropna(subset=["entrez"]).assign(entrez=lambda d: d.entrez.astype(int)).drop_duplicates("entrez")
    R = pd.read_parquet(OUT / f"rv_ratio_r_{parent}.parquet").set_index("variant")
    out = {}
    for v in R.index:
        r = R.loc[v].reindex(ko.gene).to_numpy()
        s = pd.Series(r, index=ko.entrez.to_numpy()).reindex(U)
        out[v] = (1.0 - s).fillna(0.0).round(6)
    in_p = pd.Series(np.isin(U, ko.entrez.to_numpy()), index=U)
    sym = pd.Series(ko.symbol.to_numpy(), index=ko.entrez.to_numpy())
    return out, in_p, sym


def main():
    from nextgem import labels, metrics
    from nextgem.config import SEEDS

    split = labels.load_split().set_index("ModelID")
    lf = pd.read_csv(ROOT / "results/w1/line_features.tsv", sep="\t").set_index("ModelID")
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    E = sorted(l for l in elig.index if bool(lf.loc[l, "has_expression_matrix"]))
    lineage = split.lineage
    U = [int(x) for x in labels.gene_universe().entrez]
    tr = labels.train(columns=["L_prob"])
    Y = tr.pivot(index="ModelID", columns="entrez", values="L_prob").reindex(index=E, columns=U).astype(float)
    del tr
    folds = sorted(elig.lineage.unique())

    def per_line(e: pd.Series | pd.DataFrame) -> pd.DataFrame:
        S = e if isinstance(e, pd.DataFrame) else pd.DataFrame(np.tile(e.to_numpy(), (len(E), 1)), index=E, columns=U)
        return metrics.per_line(S, Y)

    def delta(a: pd.DataFrame, b: pd.DataFrame) -> dict:
        assert (a.estimable == b.estimable).all() and (a.n_pos == b.n_pos).all()
        d = (a.tap - b.tap)[a.estimable]
        lin = lineage.reindex(d.index)
        bt = metrics.bootstrap_mean(d, lin, B=2000, seed=SEEDS[0])
        g = d.groupby(lin).agg(["mean", "size"])
        npos = sum(1 for f in folds if f in g.index and g.loc[f, "size"] >= 5 and g.loc[f, "mean"] > 0)
        return {"delta_lin": bt.estimate, "lo": bt.ci_low, "hi": bt.ci_high, "folds_pos": npos, "n_folds": len(folds)}

    rows, rr, changed = [], [], []
    for p in ("human1", "human2"):
        sc, in_p, sym = gene_scores(p, U)
        lt0 = per_line(sc["T0"])
        est = lt0[lt0.estimable]
        bt = metrics.bootstrap_mean(est.tap, lineage.reindex(est.index), B=2000, seed=SEEDS[0])
        rows.append({"parent": p, "mode": "T0", "tap_lin": bt.estimate, "lo": bt.ci_low, "hi": bt.ci_high,
                     "n_lines": len(E), "n_estimable": int(lt0.estimable.sum())})
        c0 = (sc["T0"] >= CALL) & in_p
        for v in VARS:
            c1 = (sc[v] >= CALL) & in_p
            rem, add = list(c0[c0 & ~c1].index), list(c1[c1 & ~c0].index)
            changed += [{"parent": p, "variant": v, "change": "removed", "entrez": g, "symbol": sym.get(g),
                         "e_T0": sc["T0"][g], "e_variant": sc[v][g]} for g in rem]
            changed += [{"parent": p, "variant": v, "change": "added", "entrez": g, "symbol": sym.get(g),
                         "e_T0": sc["T0"][g], "e_variant": sc[v][g]} for g in add]
            ltv = per_line(sc[v])
            est = ltv[ltv.estimable]
            bt = metrics.bootstrap_mean(est.tap, lineage.reindex(est.index), B=2000, seed=SEEDS[0])
            d = delta(ltv, lt0)
            n_partial = int(((sc[v] != sc["T0"]) & ~(c0 ^ c1)).sum())
            rows.append({"parent": p, "mode": v, "tap_lin": bt.estimate, "lo": bt.ci_low, "hi": bt.ci_high,
                         "n_lines": len(E), "n_estimable": int(ltv.estimable.sum()), "n_removed_calls": len(rem),
                         "n_added_calls": len(add), "n_other_score_changes": n_partial,
                         **{f"vsT0_{k}": x for k, x in d.items()}})
            # post-hoc label-free random-removal control
            pool = np.array(c0[c0].index)
            k = len(rem)
            base = sc["T0"].to_numpy()
            col = {g: i for i, g in enumerate(U)}
            pool_idx = np.array([col[g] for g in pool])
            seed_lts = []
            for s in SEEDS[0:5]:
                M = np.tile(base, (len(E), 1))
                for li in range(len(E)):
                    rng = np.random.default_rng([s, li])
                    M[li, rng.choice(pool_idx, size=k, replace=False)] = 0.0
                ltc = per_line(pd.DataFrame(M, index=E, columns=U))
                seed_lts.append(ltc)
                dc = delta(ltc, lt0)
                dv = delta(ltv, ltc)
                rr.append({"parent": p, "variant": v, "seed": s, "k_removed": k,
                           "control_vs_T0": dc["delta_lin"], "control_vs_T0_lo": dc["lo"], "control_vs_T0_hi": dc["hi"],
                           "variant_vs_control": dv["delta_lin"], "lo": dv["lo"], "hi": dv["hi"], "folds_pos": dv["folds_pos"]})
            mean_lt = seed_lts[0].copy()
            mean_lt["tap"] = np.mean([x.tap.to_numpy() for x in seed_lts], axis=0)
            dv = delta(ltv, mean_lt)
            dc = delta(mean_lt, lt0)
            rr.append({"parent": p, "variant": v, "seed": "mean(SEEDS[0:5])", "k_removed": k,
                       "control_vs_T0": dc["delta_lin"], "control_vs_T0_lo": dc["lo"], "control_vs_T0_hi": dc["hi"],
                       "variant_vs_control": dv["delta_lin"], "lo": dv["lo"], "hi": dv["hi"], "folds_pos": dv["folds_pos"]})
    pd.DataFrame(rows).to_csv(OUT / "rv_scores.tsv", sep="\t", index=False, float_format="%.6f")
    pd.DataFrame(rr).to_csv(OUT / "rv_random_removal.tsv", sep="\t", index=False, float_format="%.6f")
    pd.DataFrame(changed).to_csv(OUT / "rv_changed_calls.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).round(4).to_string())
    print(pd.DataFrame(rr).round(4).to_string())


if __name__ == "__main__":
    main()
