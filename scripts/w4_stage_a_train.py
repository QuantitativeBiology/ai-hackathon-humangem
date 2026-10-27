"""Stage-A training-tier scoring and selection (DESIGN.md §9 W4 "Metric and decision", §9 TG, §6; A4, A7, A8, A10).

Labels only through ``nextgem.labels.train()`` (panel lines, training tier) and ``nextgem.cv.lolo`` (CE_fit fitted on
training minus the held-out lineage). Never calls ``labels.validation()`` or ``labels.test()``.

Scores: e = 1 - ratio_r over U (W1 convention; a U gene without a Human2 gene scores 0). Primary cell: L_prob, all
genes, U, panel lines (30 training lines, 10 lineages). Secondary: minus CE_fit (LOLO), expressed pairs (§5.4),
mΔAP_gene (gene unit, vs B3), §6.8 binary endpoints at the fixed cutoff (ratio_r < 0.01 x L_prob; and the
Human1-paper biomass-only protocol PGR x L_ge06, PGR = the FBA lethal set on Ham's).

References: G* = Human1 (W1 look W1_M2_baseline) and the Stage-A parent, generic Human2 on M0. Comparators (frozen
W1 per-line values, same lines and gene set): B2[L_prob], B2m, B3, B8+/B8-/B8' for H1 and H2.

Selection (frozen rule, §9 W4): label-blind gates first (stage_a_analysis/config_summary.tsv), then nested panel
LOLO per family over its gate-passing RNA configurations {R1, R2, R3} x {no genomic, genomic}: for every held-out
panel lineage the configuration is chosen on the other 9 lineages (Δ̄_lin vs G*), and its Δ on the held-out lines
is the outer value. Panel TG on the outer values: Δ̄_lin >= 0.005, paired line-bootstrap CI lower bound > 0, >= 8 of
10 folds positive (a fold with < 2 estimable lines is not positive). <= 2 thresholds kept per family (full-panel
rank among gate-passing configurations) and the modality (protein arms compared on the 10 ProCan panel lines).
F5 (A5): nested LOLO over the 24-point grid, tagged "label-calibrated". F3b CORDA (A6'): descriptive only, 10 lines.

Outputs: results/w4/stage_a_analysis/train_*.tsv / .json (see README), MANIFEST.tsv (provenance).
Usage: .venv/bin/python scripts/w4_stage_a_train.py [--boot 2000] [--allow-dirty]
"""

from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np
import pandas as pd

from nextgem import metrics, paths
from nextgem import w4_analysis as wa
from nextgem.config import SEEDS

W4 = paths.results_dir() / "w4"
STAGE = W4 / "stage_a"
OUT = W4 / "stage_a_analysis"
LP = ["F0", "F2", "F3", "F4"]
G_STAR = "H1"  # W1 look W1_M2_baseline (results/w1/val_look_W1_M2_baseline_gstar.json)
CORDA_TAG = "descriptive, not eligible for P1 (amendment A7)"
F5_TAG = "label-calibrated (F5)"
COMPARATORS = ["B2[L_prob]", "B2m[L_prob]", "B3", "B8+_H1", "B8-_H1", "B8'_H1", "B8+_H2", "B8-_H2", "B8'_H2"]


def cfg_name(fam, rule, mod, var, gen):
    return f"{fam}|{rule}|{mod}|{var}|{gen}"


def load(a):
    from nextgem import cv, labels, w1

    split = labels.load_split().set_index("ModelID")
    panel = pd.read_csv(W4 / "stage_a_panel.tsv", sep="\t")
    lines = list(panel.ModelID)
    lineage = split.lineage.reindex(lines)
    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    gen = w1.generic_table(U.entrez)
    emap = wa.entrez_gene_map(pd.read_parquet(paths.baseline_dir() / "human2_ko.parquet", columns=["gene", "entrez"]))
    tl = labels.train(columns=["L_prob", "L_ge06"])
    tl = tl[tl.ModelID.isin(lines)]
    Y = w1.wide(tl, "L_prob", lines, genes)
    Y06 = w1.wide(tl, "L_ge06", lines, genes)
    panel_lins = sorted(lineage.unique())
    folds = cv.lolo(lambda fit, ids: (labels.ce_fit(fit, "L_prob"), ids), lambda p, held: p,
                    columns=["L_prob"], folds=panel_lins)
    line_ce = {l: ce for k, (ce, ids) in folds.items() for l in ids if l in set(lines)}
    expr = w1.load_omics()["expr"].reindex(index=lines, columns=genes)
    b = pd.read_csv(OUT / "build_gates.tsv", sep="\t", low_memory=False)
    cs = pd.read_csv(OUT / "config_summary.tsv", sep="\t")
    w1pl = pd.read_csv(paths.results_dir() / "w1" / "train_per_line.tsv", sep="\t")
    w1pl = w1pl[w1pl.ModelID.isin(lines)]
    gmap = pd.read_csv(STAGE / "a5_grid_map.tsv", sep="\t")
    return dict(lines=lines, lineage=lineage, genes=genes, gen=gen, emap=emap, Y=Y, Y06=Y06, line_ce=line_ce,
                expr=expr, builds=b, cs=cs, w1pl=w1pl, gmap=gmap, panel=panel, split=split)


def ratio_tables():
    """build_id -> Series(gene -> ratio_r) for every Stage-A build x genomic variant."""
    out = {}
    for r in ("A1", "A2", "A3", "A4", "A5", "A6p", "A7"):
        k = pd.read_parquet(STAGE / f"ko_{r}.parquet", columns=["build_id", "gene", "ratio_r"])
        k["build_id"] = k.build_id.astype(str)
        for bid, g in k.groupby("build_id"):
            out[bid] = g.set_index("gene").ratio_r
    return out


class Scorer:
    def __init__(self, C, ratios):
        self.C, self.R = C, ratios
        lines, genes = C["lines"], C["genes"]
        self.Y = C["Y"]
        self.Ync = self.Y.where(pd.DataFrame(np.vstack([~np.isin(genes, list(C["line_ce"][l])) for l in lines]),
                                             index=lines, columns=genes))
        self.Yex = self.Y.where(C["expr"] >= 1.0)
        self.h2_genes = [g for g in genes if C["gen"].loc[g, "in_H2"]]
        self.h1_genes = [g for g in genes if C["gen"].loc[g, "in_H1"]]

    def frame(self, build_ids: dict) -> pd.DataFrame:
        """build_ids: line -> build_id. Score frame lines x U."""
        lines = list(build_ids)
        ratio = pd.DataFrame({l: self.R[build_ids[l]] for l in lines}).T
        return wa.score_frame(ratio, self.C["genes"], self.C["emap"])

    def generic(self, s: str, lines) -> pd.DataFrame:
        e = self.C["gen"][f"e_{s}"]
        return pd.DataFrame(np.tile(e.reindex(self.C["genes"]).fillna(0).to_numpy(), (len(lines), 1)),
                            index=lines, columns=self.C["genes"])

    def per_line(self, S: pd.DataFrame) -> pd.DataFrame:
        idx = S.index
        a = metrics.per_line(S, self.Y.loc[idx])
        b = metrics.per_line(S, self.Ync.loc[idx])
        c = metrics.per_line(S, self.Yex.loc[idx])
        out = a[["n_pos", "estimable", "tap", "e_random", "lift"]].copy()
        out["tap_noCE"], out["est_noCE"] = b.tap, b.estimable
        out["tap_expr"], out["est_expr"] = c.tap, c.estimable
        return out

    def binary(self, S: pd.DataFrame, parent_genes) -> pd.DataFrame:
        idx = S.index
        calls = (S.loc[:, parent_genes] > 1 - metrics.FBA_CUTOFF + 1e-12)  # e > 0.99  <=>  ratio_r < 0.01
        out = {}
        for lab, Yt in (("L_prob", self.Y), ("L_ge06", self.C["Y06"])):
            pl = metrics.binary_per_line(calls, Yt.loc[idx, parent_genes])
            out[f"mcc_{lab}"] = pl.mcc
            out[f"f1_{lab}"] = pl.f1
        return pd.DataFrame(out)


def delta_stats(a: pd.Series, b: pd.Series, lineage, B, seed=SEEDS[0]) -> dict:
    d = (a - b).dropna()
    if d.empty:
        return {"delta_lin": np.nan, "lo": np.nan, "hi": np.nan}
    bt = metrics.bootstrap_mean(d, lineage.reindex(d.index), B=B, seed=seed)
    return {"delta_lin": bt.estimate, "lo": bt.ci_low, "hi": bt.ci_high}


def main(argv=None) -> int:
    from nextgem.run import provenance

    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    inputs = [STAGE / f for f in ("builds.tsv", "ko_A1.parquet", "ko_A2.parquet", "ko_A3.parquet", "ko_A4.parquet",
                                   "ko_A5.parquet", "ko_A6p.parquet", "ko_A7.parquet", "a5_grid_map.tsv", "MANIFEST.tsv")]
    inputs += [OUT / "build_gates.tsv", OUT / "config_summary.tsv", W4 / "stage_a_panel.tsv",
               paths.results_dir() / "w1" / "train_per_line.tsv", paths.baseline_dir() / "human2_ko.parquet",
               paths.baseline_dir() / "human1_ko.parquet"]
    with provenance(OUT / "MANIFEST.tsv", inputs=inputs, seeds=[SEEDS[0]], solver="none (analysis)", workers=1,
                    allow_dirty=a.allow_dirty, author="W4 engineer agent") as prov:
        prov.note(script_part="training scoring", labels="labels.train() + cv.lolo (CE_fit) only; no validation/test",
                  g_star=G_STAR, boot=a.boot)
        outs = run(a)
        for f in outs:
            prov.output(OUT / f)
    return 0


def run(a):
    t0 = time.time()
    C = load(a)
    R = ratio_tables()
    sc = Scorer(C, R)
    lines, lineage, B = C["lines"], C["lineage"], a.boot
    b = C["builds"]
    b["bid"] = b.build_id.astype(str)

    # ------------------------------------------------------------------ references, comparators, verification
    ref = {s: sc.per_line(sc.generic(s, lines)) for s in ("H1", "H2")}
    ver = []
    for s in ("H1", "H2"):
        for gs, col in (("all", "tap"), ("noCE", "tap_noCE")):
            w = C["w1pl"][(C["w1pl"].method == s) & (C["w1pl"].geneset == gs)].set_index("ModelID").tap.reindex(lines)
            diff = float((ref[s][col] - w).abs().max())
            # W1's TSV stores 6 decimals, so equality is checked to 1e-6 (rounding <= 5e-7)
            ver.append({"check": f"generic {s} tAP_line ({gs}) equals W1 train_per_line (6-decimal TSV)",
                        "max_abs_diff": diff, "pass": diff <= 1e-6})
    comp = {}
    for m in COMPARATORS:
        for gs, col in (("all", "tap"), ("noCE", "tap_noCE")):
            w = C["w1pl"][(C["w1pl"].method == m) & (C["w1pl"].geneset == gs)].set_index("ModelID").tap.reindex(lines)
            comp.setdefault(m, pd.DataFrame(index=lines))[col] = w
    bin_ref = {"H1": sc.binary(sc.generic("H1", lines), sc.h1_genes), "H2": sc.binary(sc.generic("H2", lines), sc.h2_genes)}

    # ------------------------------------------------------------------ per-line tables for every configuration
    per_line, meta = {}, {}

    def add(name, bmap, row, fam, rule, mod, var, gen, extra=None):
        S = sc.frame(bmap)
        pl = sc.per_line(S)
        bn = sc.binary(S, sc.h2_genes)
        per_line[name] = pl.join(bn)
        meta[name] = {"config": name, "row": row, "family": fam, "rule": rule, "modality": mod, "variant": var,
                      "genomic": gen, "n_lines": len(bmap), **(extra or {})}
        if name.startswith(("F5", "A5")) or fam == "F5":
            meta[name]["status_tag"] = F5_TAG
        if fam == "F3b":
            meta[name]["status_tag"] = CORDA_TAG
        mg = metrics.per_gene(S, sc.Y.loc[S.index])
        meta[name]["m_delta_ap_gene"] = metrics.m_delta_ap_gene(mg)

    sel = b[b.row.isin(["A1", "A2", "A3", "A6p", "A7"])]
    for key, g in sel.groupby(["row", "family", "rule", "modality", "variant", "genomic"]):
        row, fam, rule, mod, var, gen = key
        name = cfg_name(fam, rule, mod, var, gen) if row != "A7" else f"A7|{fam}|{rule}|{var}"
        add(name, dict(zip(g.line, g.bid)), row, fam, rule, mod, var, gen)
    # A5: 24 grid points, each line's build from the grid map (deduplicated builds)
    gm = C["gmap"]
    for (wr, wp, wc, tau), g in gm.groupby(["w_r", "w_p", "w_c", "tau"]):
        name = f"F5|wr{wr:g}_wp{wp:g}_wc{wc:g}_tau{tau:g}"
        add(name, {l: f"{bid}|none" for l, bid in zip(g.line, g.build_id)}, "A5", "F5", "lcp", "lcp-grid",
            f"wp{wp:g}_wc{wc:g}_tau{tau:g}", "none", {"w_r": wr, "w_p": wp, "w_c": wc, "tau": tau})
    # gene-unit reference for B3 (expression) on the panel
    mg_b3 = metrics.m_delta_ap_gene(metrics.per_gene(C["expr"].fillna(0.0), sc.Y))
    mg_ref = {s: metrics.m_delta_ap_gene(metrics.per_gene(sc.generic(s, lines), sc.Y)) for s in ("H1", "H2")}

    # long per-line table
    rows = []
    for n, pl in per_line.items():
        t = pl.copy()
        t.insert(0, "config", n)
        t["lineage"] = lineage.reindex(t.index)
        rows.append(t.rename_axis("ModelID").reset_index())
    for s in ("H1", "H2"):
        t = ref[s].join(bin_ref[s]).copy()
        t.insert(0, "config", f"generic_{s}")
        t["lineage"] = lineage.reindex(t.index)
        rows.append(t.rename_axis("ModelID").reset_index())
    pd.concat(rows).to_csv(OUT / "train_per_line.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ per-configuration scores (non-nested)
    cs = C["cs"]
    gate = {cfg_name(r.family, r.rule, r.modality, r.variant, r.genomic): bool(r.config_gate_pass)
            for r in cs.itertuples() if r.row in ("A1", "A2", "A3", "A6p")}
    crow = []
    for n, pl in per_line.items():
        m = meta[n]
        lin = lineage.reindex(pl.index)
        est = pl[pl.estimable]
        r_ = {**m, "gate_pass": gate.get(n, True if m["family"] == "F5" else None),
              "tap_lin": metrics.lineage_equal_mean(est.tap, lin.reindex(est.index)),
              "lift_lin": metrics.lineage_equal_mean(est.lift, lin.reindex(est.index)),
              "n_lines_estimable": int(pl.estimable.sum()),
              "m_delta_ap_gene_minus_B3": m["m_delta_ap_gene"] - mg_b3}
        for s in ("H1", "H2"):
            idx = pl.index
            tg = wa.tg_panel(pl.tap - ref[s].tap.reindex(idx), lineage, B=B)
            r_.update({f"delta_vs_{s}": tg["delta_lin"], f"lo_vs_{s}": tg["lo"], f"hi_vs_{s}": tg["hi"],
                       f"folds_pos_vs_{s}": tg["folds_positive"], f"tg_vs_{s}": tg["tg_pass"]})
            ds = delta_stats(pl.tap_noCE, ref[s].tap_noCE.reindex(idx), lineage, B)
            r_.update({f"delta_noCE_vs_{s}": ds["delta_lin"], f"lo_noCE_vs_{s}": ds["lo"]})
            de = delta_stats(pl.tap_expr, ref[s].tap_expr.reindex(idx), lineage, B)
            r_.update({f"delta_expr_vs_{s}": de["delta_lin"], f"lo_expr_vs_{s}": de["lo"]})
            dm = delta_stats(pl.mcc_L_prob, bin_ref[s].mcc_L_prob.reindex(idx), lineage, B)
            r_.update({f"delta_mcc_vs_{s}": dm["delta_lin"], f"lo_mcc_vs_{s}": dm["lo"]})
        for cm in ("B2[L_prob]", "B8+_H1", "B8-_H1", "B8'_H1", "B8+_H2", "B8-_H2", "B8'_H2", "B3"):
            ds = delta_stats(pl.tap, comp[cm].tap.reindex(pl.index), lineage, B)
            r_.update({f"delta_vs_{cm}": ds["delta_lin"], f"lo_vs_{cm}": ds["lo"]})
        r_["mcc_L_prob_lin"] = metrics.lineage_equal_mean(pl.mcc_L_prob.dropna(), lin.reindex(pl.mcc_L_prob.dropna().index))
        r_["mcc_L_ge06_lin"] = metrics.lineage_equal_mean(pl.mcc_L_ge06.dropna(), lin.reindex(pl.mcc_L_ge06.dropna().index))
        crow.append(r_)
    for s in ("H1", "H2"):
        pl = ref[s].join(bin_ref[s])
        crow.append({"config": f"generic_{s}", "family": "generic", "tap_lin": metrics.lineage_equal_mean(pl.tap, lineage),
                     "lift_lin": metrics.lineage_equal_mean(pl.lift, lineage), "m_delta_ap_gene": mg_ref[s],
                     "m_delta_ap_gene_minus_B3": mg_ref[s] - mg_b3,
                     "mcc_L_prob_lin": metrics.lineage_equal_mean(pl.mcc_L_prob.dropna(), lineage.reindex(pl.mcc_L_prob.dropna().index)),
                     "mcc_L_ge06_lin": metrics.lineage_equal_mean(pl.mcc_L_ge06.dropna(), lineage.reindex(pl.mcc_L_ge06.dropna().index))})
    for cm in COMPARATORS:
        pl = comp[cm]
        crow.append({"config": cm, "family": "comparator (W1, frozen)",
                     "tap_lin": metrics.lineage_equal_mean(pl.tap.dropna(), lineage.reindex(pl.tap.dropna().index)),
                     **{f"delta_vs_{s}": delta_stats(pl.tap, ref[s].tap, lineage, B)["delta_lin"] for s in ("H1", "H2")}})
    crow.append({"config": "B3 gene unit", "family": "comparator", "m_delta_ap_gene": mg_b3})
    ctab = pd.DataFrame(crow)
    ctab.to_csv(OUT / "train_config_scores.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ nested panel LOLO selection per family
    order_rna = [(r, g) for r in ("R1", "R2", "R3") for g in ("none", "loss")]
    fam_out, sel_rows = {}, []
    for fam in LP + ["F5"]:
        if fam == "F5":
            cands = [n for n in per_line if meta[n]["family"] == "F5"]
            cands = sorted(cands, key=lambda n: (meta[n]["w_p"], meta[n]["w_c"], -meta[n]["tau"]))  # simplest first
        else:
            cands = [cfg_name(fam, r, "rna", "real", g) for r, g in order_rna]
            cands = [n for n in cands if n in per_line]
        gated = [n for n in cands if gate.get(n, True)]
        res = {"family": fam, "candidates": len(cands), "gate_passing": len(gated),
               "status_tag": F5_TAG if fam == "F5" else ""}
        use = gated if gated else cands  # gate-failing families are scored descriptively only
        res["descriptive_only_gate_failed"] = not gated
        for s in ("H1", "H2"):
            D = pd.DataFrame({n: per_line[n].tap - ref[s].tap.reindex(per_line[n].index) for n in use})
            outer, chosen = wa.nested_lolo_select(D, lineage, order=use)
            tg = wa.tg_panel(outer, lineage, B=B)
            rank = wa.full_panel_rank(D, lineage, order=use)
            res[f"vs_{s}"] = {"outer_delta_lin": tg["delta_lin"], "lo": tg["lo"], "hi": tg["hi"],
                              "folds_positive": tg["folds_positive"], "tg_pass": tg["tg_pass"] and bool(gated),
                              "chosen_per_fold": chosen, "full_panel_best": rank.config.iloc[0],
                              "full_panel_rank": rank.head(6).to_dict("records")}
            for k, c in chosen.items():
                sel_rows.append({"family": fam, "reference": s, "heldout_lineage": k, "chosen": c,
                                 "outer_delta": tg["fold_means"].get(k)})
        # <= 2 thresholds kept (full-panel rank vs G*, gate-passing, distinct thresholds)
        if fam != "F5":
            kept, rules = [], set()
            for c in res["vs_H1"]["full_panel_rank"]:
                r_ = meta[c["config"]]["rule"]
                if r_ not in rules:
                    kept.append(c["config"])
                    rules.add(r_)
                if len(kept) == 2:
                    break
            res["kept_thresholds"] = kept
        fam_out[fam] = res
    pd.DataFrame(sel_rows).to_csv(OUT / "train_nested_selection.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ modality (A3) on the 10 ProCan panel lines
    mrows = []
    for fam in LP:
        for rule in ("R1", "R2", "R3"):
            for gen in ("none", "loss"):
                base = cfg_name(fam, rule, "rna", "real", gen)
                if base not in per_line:
                    continue
                for mod in ("protein", "rna_protein"):
                    n = cfg_name(fam, rule, mod, "real", gen)
                    if n not in per_line:
                        continue
                    pl = per_line[n]
                    ds = delta_stats(pl.tap, per_line[base].tap.reindex(pl.index), lineage, B)
                    dg = delta_stats(pl.tap, ref[G_STAR].tap.reindex(pl.index), lineage, B)
                    mrows.append({"family": fam, "rule": rule, "genomic": gen, "modality": mod, "n_lines": len(pl),
                                  "n_lineages": int(lineage.reindex(pl.index).nunique()), "gate_pass": gate.get(n),
                                  "delta_vs_rna": ds["delta_lin"], "lo_vs_rna": ds["lo"], "hi_vs_rna": ds["hi"],
                                  "delta_vs_H1": dg["delta_lin"]})
    mtab = pd.DataFrame(mrows)
    mtab.to_csv(OUT / "train_modality.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ CORDA (A6'), A7 stability, permuted control
    a6 = pd.read_csv(STAGE / "a6prime_lines.tsv", sep="\t").ModelID.tolist()
    crows = []
    for rule in ("R1", "R3"):
        cn = cfg_name("F3b", rule, "rna", "real", "none")
        for fam in LP + ["F3b"]:
            n = cfg_name(fam, rule, "rna", "real", "none")
            pl = per_line[n].reindex(a6)
            crows.append({"rule": rule, "family": fam, "n_lines": int(pl.tap.notna().sum()),
                          "tap_line_mean": float(pl.tap.mean()),
                          "delta_vs_H1": float((pl.tap - ref["H1"].tap.reindex(a6)).mean()),
                          "delta_vs_H2": float((pl.tap - ref["H2"].tap.reindex(a6)).mean()),
                          "delta_corda_minus_family": float((per_line[cn].tap.reindex(a6) - pl.tap).mean()) if fam != "F3b" else 0.0,
                          "mcc_L_prob": float(pl.mcc_L_prob.mean()),
                          "status_tag": CORDA_TAG if fam == "F3b" else ""})
    pd.DataFrame(crows).to_csv(OUT / "train_corda_a6prime.tsv", sep="\t", index=False)
    arows = []
    for fam in ("F2", "F3"):
        base = per_line[cfg_name(fam, "R1", "rna", "real", "none")]
        seeds = [n for n in per_line if n.startswith(f"A7|{fam}|")]
        T = pd.DataFrame({n: per_line[n].tap for n in seeds})
        idx = T.index
        arows.append({"family": fam, "n_lines": len(idx), "n_seeds": len(seeds),
                      "tap_sd_across_seeds_mean": float(T.std(axis=1, ddof=1).mean()),
                      "tap_range_across_seeds_max": float((T.max(1) - T.min(1)).max()),
                      "mean_delta_seed_minus_unperturbed": float((T.mean(1) - base.tap.reindex(idx)).mean()),
                      "unperturbed_tap_line_mean": float(base.tap.reindex(idx).mean())})
    pd.DataFrame(arows).to_csv(OUT / "train_a7_stability.tsv", sep="\t", index=False)
    prows = []
    for fam in LP:
        for rule in ("R1", "R3"):
            r_, p_ = per_line[cfg_name(fam, rule, "rna", "real", "none")], per_line[cfg_name(fam, rule, "rna", "perm", "none")]
            ds = delta_stats(r_.tap, p_.tap.reindex(r_.index), lineage, B)
            dsp = delta_stats(p_.tap, ref[G_STAR].tap.reindex(p_.index), lineage, B)
            prows.append({"family": fam, "rule": rule, "delta_real_minus_perm": ds["delta_lin"], "lo": ds["lo"],
                          "hi": ds["hi"], "perm_delta_vs_H1": dsp["delta_lin"]})
    pd.DataFrame(prows).to_csv(OUT / "train_permuted_control.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ shortlist
    elig = [f for f in LP + ["F5"] if fam_out[f]["vs_H1"]["tg_pass"]]
    ranked = sorted(elig, key=lambda f: -fam_out[f]["vs_H1"]["outer_delta_lin"])
    all_ranked = sorted([f for f in LP + ["F5"] if not fam_out[f]["descriptive_only_gate_failed"]],
                        key=lambda f: -fam_out[f]["vs_H1"]["outer_delta_lin"])
    shortlist = {
        "rule": "TG-passing families (nested panel LOLO vs G* = H1, >= 8 of 10, gates first); <= 2 configurations, "
                "ranked by nested outer Δ̄_lin; the best configuration is carried even if it fails TG (flagged)",
        "tg_passing_families": ranked,
        "w4_shortlist": [{"family": f, "config": fam_out[f]["vs_H1"]["full_panel_best"],
                          "outer_delta_vs_H1": fam_out[f]["vs_H1"]["outer_delta_lin"],
                          "tag": fam_out[f]["status_tag"]} for f in ranked[:2]],
        "best_if_none_pass": None if ranked else {"family": all_ranked[0], "config": fam_out[all_ranked[0]]["vs_H1"]["full_panel_best"],
                                                  "flag": "fails TG; carried to M6 per §9 W4"},
        "stage_b_configs": [fam_out[f]["vs_H1"]["full_panel_best"] for f in (ranked or all_ranked[:1])][:3],
        "kept_thresholds": {f: fam_out[f].get("kept_thresholds") for f in LP},
        "families": fam_out,
        "verification": ver, "m_delta_ap_gene_B3": mg_b3, "m_delta_ap_gene_generic": mg_ref,
        "wall_s": round(time.time() - t0),
    }
    (OUT / "train_selection.json").write_text(json.dumps(shortlist, indent=1, default=float))
    pd.DataFrame(ver).to_csv(OUT / "train_verification.tsv", sep="\t", index=False)
    print(json.dumps({k: shortlist[k] for k in ("tg_passing_families", "w4_shortlist", "best_if_none_pass",
                                                "stage_b_configs", "kept_thresholds")}, indent=1, default=float))
    for f, r in fam_out.items():
        v = r["vs_H1"]
        print(f, "gate-passing", r["gate_passing"], "| nested vs H1:", round(v["outer_delta_lin"], 4),
              [round(v["lo"], 4), round(v["hi"], 4)], "folds+", v["folds_positive"], "TG", v["tg_pass"],
              "| vs H2:", round(r["vs_H2"]["outer_delta_lin"], 4), "TG", r["vs_H2"]["tg_pass"], "| best", v["full_panel_best"])
    print("verification:", ver)
    return ["train_per_line.tsv", "train_config_scores.tsv", "train_nested_selection.tsv", "train_modality.tsv",
            "train_corda_a6prime.tsv", "train_a7_stability.tsv", "train_permuted_control.tsv",
            "train_selection.json", "train_verification.tsv"]


if __name__ == "__main__":
    sys.exit(main())
