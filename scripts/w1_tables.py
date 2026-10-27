"""Print the W1 README tables (markdown) from the result files, so every number in results/w1/README.md is
copied from an output, not typed. Usage: .venv/bin/python scripts/w1_tables.py"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

D = Path(__file__).resolve().parents[1] / "results" / "w1"


def f(x, n=3):
    return "NA" if pd.isna(x) else f"{x:.{n}f}"


def ci(e, lo, hi, n=3):
    return f"{f(e, n)} [{f(lo, n)}, {f(hi, n)}]"


def main():
    g = pd.read_csv(D / "train_grid.tsv", sep="\t")
    order = ["H1", "H2", "H1_task", "H2_task", "B0", "B1_mean", "B2[L_prob]", "B2m[L_prob]", "B2c", "B3", "B6[L_prob]",
             "B8+_H1", "B8-_H1", "B8'_H1", "B8+_H2", "B8-_H2", "B8'_H2"]
    for gs in ("all", "noCE", "noDepMapCE"):
        p = g[(g.label == "L_prob") & (g.geneset == gs) & (g.universe == "U") & (g.stratum == "none") & (g.lineset == "expr")]
        p = p.set_index("method").reindex(order)
        print(f"\n**L_prob, {gs}, U, training expression lines (n = {int(p.loc['H2', 'n_lines'])})**\n")
        print("| method | tAP_line (lineage-equal) [95% CI] | E[AP_random] | lift | Δ vs H1 [CI] | Δ vs H2 [CI] | mΔAP_gene | pooled tAP (lift) |"
              + (" AUROC | sklearn AP |" if gs != "noDepMapCE" else ""))
        print("|---|---|---|---|---|---|---|---|" + ("---|---|" if gs != "noDepMapCE" else ""))
        for m, r in p.iterrows():
            extra = f" {f(r.get('auroc_lin'))} | {f(r.get('ap_sklearn_lin'))} |" if gs != "noDepMapCE" else ""
            print(f"| {m} | {ci(r.tap_lin, r.tap_lin_lo, r.tap_lin_hi)} | {f(r.e_random_lin)} | {f(r.lift_lin, 2)} | "
                  f"{ci(r.get('delta_lin_vs_H1'), r.get('delta_lo_vs_H1'), r.get('delta_hi_vs_H1'), 4)} | "
                  f"{ci(r.get('delta_lin_vs_H2'), r.get('delta_lo_vs_H2'), r.get('delta_hi_vs_H2'), 4)} | "
                  f"{f(r.m_delta_ap_gene, 4)} | {f(r.pooled_tap)} ({f(r.pooled_lift, 2)}) |" + extra)
    # labels x gene sets for H1/H2/B2
    print("\n**Generic models and B2 across labels and gene sets (U, expression lines; tAP_line lineage-equal, lift, NA flag)**\n")
    print("| label | gene set | H1 | H2 | H2 − H1 [CI] | B2 | B8+ (H2) |")
    print("|---|---|---|---|---|---|---|")
    for lab in ("L_prob", "L_ge05", "L_ge10"):
        for gs in ("all", "noCE", "noDepMapCE"):
            p = g[(g.label == lab) & (g.geneset == gs) & (g.universe == "U") & (g.stratum == "none") & (g.lineset == "expr")].set_index("method")
            b2 = p.loc[f"B2[{lab}]"]
            h1, h2 = p.loc["H1"], p.loc["H2"]
            na = " (NA)" if h2.na else ""
            print(f"| {lab} | {gs}{na} | {f(h1.tap_lin)} ({f(h1.lift_lin, 2)}) | {f(h2.tap_lin)} ({f(h2.lift_lin, 2)}) | "
                  f"{ci(h2.delta_lin_vs_H1, h2.delta_lo_vs_H1, h2.delta_hi_vs_H1, 4)} | {f(b2.tap_lin)} ({f(b2.lift_lin, 2)}) | "
                  f"{f(p.loc['B8+_H2'].tap_lin)} |")
    # strata
    print("\n**Strata (L_prob, all genes, U): H1, H2 tAP_line and H2 − H1**\n")
    print("| stratum | lines estimable | NA | H1 | H2 | H2 − H1 [CI] | B2 | E[AP_random] |")
    print("|---|---|---|---|---|---|---|---|")
    for st in g.stratum.unique():
        p = g[(g.label == "L_prob") & (g.geneset == "all") & (g.universe == "U") & (g.stratum == st) & (g.lineset == "expr")].set_index("method")
        h1, h2 = p.loc["H1"], p.loc["H2"]
        print(f"| {st} | {int(h2.n_lines_estimable)}/{int(h2.n_lines)} | {h2.na} | {f(h1.tap_lin)} | {f(h2.tap_lin)} | "
              f"{ci(h2.delta_lin_vs_H1, h2.delta_lo_vs_H1, h2.delta_hi_vs_H1, 4)} | {f(p.loc['B2[L_prob]'].tap_lin)} | {f(h2.e_random_lin)} |")
    # U-cap and all lines
    p = g[(g.label == "L_prob") & (g.geneset == "all") & (g.universe == "Ucap") & (g.stratum == "none") & (g.lineset == "expr")].set_index("method")
    print(f"\nU∩: H1 {f(p.loc['H1'].tap_lin)}, H2 {f(p.loc['H2'].tap_lin)}, H2 − H1 "
          f"{ci(p.loc['H2'].delta_lin_vs_H1, p.loc['H2'].delta_lo_vs_H1, p.loc['H2'].delta_hi_vs_H1, 4)}")
    p = g[(g.label == "L_prob") & (g.geneset == "all") & (g.universe == "U") & (g.stratum == "none") & (g.lineset == "all")].set_index("method")
    print(f"All eligible-lineage lines (incl. no expression, n = {int(p.loc['H2'].n_lines)}): H1 {f(p.loc['H1'].tap_lin)}, "
          f"H2 {f(p.loc['H2'].tap_lin)}, H2 − H1 {ci(p.loc['H2'].delta_lin_vs_H1, p.loc['H2'].delta_lo_vs_H1, p.loc['H2'].delta_hi_vs_H1, 4)}, "
          f"B2 {f(p.loc['B2[L_prob]'].tap_lin)}")
    # per lineage
    pl = pd.read_csv(D / "train_per_lineage.tsv", sep="\t")
    pl = pl[(pl.geneset == "all") & pl.method.isin(["H1", "H2", "B2[L_prob]", "B8+_H2"])]
    w = pl.pivot(index="lineage", columns="method", values="mean_tap")
    w["n_lines"] = pl[pl.method == "H2"].set_index("lineage").n_estimable
    w["H2 − H1"] = w["H2"] - w["H1"]
    w["e_random"] = pl[pl.method == "H2"].set_index("lineage").mean_e_random
    print("\n**Per training lineage (L_prob, all genes, U): mean tAP_line**\n")
    print("| lineage | n | E[AP_random] | H1 | H2 | H2 − H1 | B2 | B8+ (H2) |")
    print("|---|---|---|---|---|---|---|---|")
    for lin, r in w.iterrows():
        print(f"| {lin} | {int(r.n_lines)} | {f(r.e_random)} | {f(r.H1)} | {f(r.H2)} | {f(r['H2 − H1'], 4)} | {f(r['B2[L_prob]'])} | {f(r['B8+_H2'])} |")
    # binary
    b = pd.read_csv(D / "train_binary.tsv", sep="\t")
    b = b[(b.lineset == "expr") & (b.nan_policy == "exclude")]
    print("\n**§6.8 binary endpoints, training expression lines, parent gene set (lineage-equal mean [95% CI]; line mean)**\n")
    print("| protocol | model | call | label | MCC | F1 | precision | recall | specificity | accuracy | called genes |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for (proto, m, rule, lab), d in b[b.gene_set == "parent"].groupby(["protocol", "method", "call_rule", "label"], sort=False):
        d = d.set_index("metric")
        cells = [f"{ci(d.loc[k].lineage_equal_mean, d.loc[k].ci_lo, d.loc[k].ci_hi)}; {f(d.loc[k].line_mean)}" if k == "mcc"
                 else f(d.loc[k].lineage_equal_mean) for k in ("mcc", "f1", "precision", "recall", "specificity", "accuracy")]
        print(f"| {proto} | {m} | {rule} | {lab} | " + " | ".join(cells) + f" | {d.iloc[0].n_called if 'n_called' in d else ''} |")
    pr = pd.read_csv(D / "train_binary_paired.tsv", sep="\t")
    print("\n**Paired H2 − H1 (per line, own parent gene sets)**\n")
    print("| protocol | call | label | metric | Δ (lineage-equal) [95% CI] |")
    print("|---|---|---|---|---|")
    for _, r in pr.iterrows():
        print(f"| {r.protocol} | {r.call_rule} | {r.label} | {r.metric} | {ci(r.delta_lin, r.ci_lo, r.ci_hi, 4)} |")
    s = json.loads((D / "train_summary.json").read_text())
    print("\n" + json.dumps({k: s[k] for k in ("gate", "sesoi_c1b", "bca_skew_check")}, indent=1, default=float))
    a = pd.read_csv(D / "train_addressable.tsv", sep="\t")
    print("\n**Addressable set (KO closes >= 1 reaction): share of training positive pairs inside (= recall ceiling)**\n")
    print("| label | parent | addressable U genes | positive pairs | recall ceiling | recall ceiling, no CE_fit |")
    print("|---|---|---|---|---|---|")
    for _, r in a.iterrows():
        print(f"| {r.label} | {r.parent} | {r.n_addressable_U} | {r.positive_pairs} | {f(r.recall_ceiling)} | {f(r.recall_ceiling_noCE)} |")
    bg = pd.read_csv(D / "background_groups.tsv", sep="\t")
    print("\n**Background groups (screened lines; members / evaluable)**\n")
    print("| group | train | validation | test | all screened |")
    print("|---|---|---|---|---|")
    for grp, d in bg.groupby("group", sort=False):
        d = d.set_index("tier")
        cell = lambda t: f"{int(d.loc[t].n_members)} / {'yes' if d.loc[t].evaluable else 'no'} ({int(d.loc[t].n_eligible_lineages)} lin.)"
        print(f"| {grp} | {cell('train')} | {cell('validation')} | {cell('test')} | {int(d.loc['all'].n_members)} |")


if __name__ == "__main__":
    main()
