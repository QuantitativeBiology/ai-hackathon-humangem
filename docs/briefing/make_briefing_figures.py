"""Figures for the nextgem research briefing (BRIEFING.md).

Reads only files committed to the nextgem repository (results/, autoresearch/, docs/explainer/figures/
figure_data.json) and the repository's git history (commit times). It never opens a DepMap label file or
anything under data/.

Outputs, for each figure: <out>/<name>.svg (text kept as text) and <out>/<name>.png (150 dpi), plus
<out>/briefing_figure_data.json with every plotted number and the file it came from.

Figures are drawn on an opaque white card so they read the same in light and dark viewers.
Colours: the validated reference categorical palette also used by the repository explainer
(blue = P1 / training, orange = Human2 / validation / safer variant, green = test, grey = reference).

Usage:
  /data/tmpA/nextgem/.venv/bin/python make_briefing_figures.py --repo /data/tmpA/nextgem --out figures
"""
from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

INK, MUTED, GRID = "#1d1d1b", "#5f5e58", "#e3e2dc"
S1, S2, S3, NEUTRAL = "#2a78d6", "#eb6834", "#1baf7a", "#8d8c85"
BAND, STALL, CARD = "#d7e6f8", "#efe3d6", "#ffffff"

plt.rcParams.update({
    "svg.fonttype": "none", "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 10,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.8, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.facecolor": CARD, "figure.facecolor": CARD, "savefig.facecolor": CARD, "legend.frameon": False,
    "svg.hashsalt": "nextgem-briefing",
})

DATA: dict = {}


def save(fig, out: Path, name: str):
    fig.savefig(out / f"{name}.svg", bbox_inches="tight", pad_inches=0.15, metadata={"Date": None})
    fig.savefig(out / f"{name}.png", dpi=150, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def tsv(repo, rel, **kw):
    return pd.read_csv(repo / rel, sep="\t", comment=kw.pop("comment", None), **kw)


def js(repo, rel):
    return json.loads((repo / rel).read_text())


# ------------------------------------------------------------------------------------------ fig 1
def fig_design(repo, out):
    s = tsv(repo, "results/splits/split_summary.tsv")
    meta = js(repo, "results/splits/split_meta.json")
    vl = tsv(repo, "results/validation_ledger.tsv")
    tl = tsv(repo, "results/test_ledger.tsv")
    vlooks = vl[vl.event == "completed"].look_id.tolist()
    tlooks = tl[tl.event == "completed"].look_id.tolist()

    def tier(t):
        x = s[s.tier == t]
        main = x[x.stratum != "rare"].sort_values("n_lines", ascending=False)
        rare = x[x.stratum == "rare"]
        return main, rare, int(x.n_lines.sum()), int(x.n_expression.sum())

    fig, ax = plt.subplots(figsize=(10.4, 5.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(-7, 60)
    ax.axis("off")

    def box(x, y, w, h, fc, ec, lw=1.2):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", fc=fc, ec=ec, lw=lw))

    # source box
    box(0.5, 21, 15, 21, "#f4f3ef", MUTED)
    ax.text(8, 40.6, f"{meta['n_lines']:,} CRISPR-\nscreened lines", ha="center", va="top", fontsize=10, fontweight="bold")
    ax.text(8, 34.0, "DepMap 26Q3,\nsplit once by\nlineage (tissue);\n" f"seed {meta['seed']}\nrevealed after\nthe design freeze",
            ha="center", va="top", fontsize=7.8, color=MUTED)
    cols = [("train", "TRAINING", S1, 19.5), ("validation", "VALIDATION (sealed)", S2, 46.5), ("test", "TEST (sealed)", S3, 73.5)]
    info = {}
    for t, title, col, x0 in cols:
        main, rare, n, ne = tier(t)
        info[t] = {"lineages": main.lineage.tolist(), "n_lines": n, "n_expression": ne, "rare_lines": int(rare.n_lines.sum())}
        box(x0, 12, 25, 45, CARD, col, lw=1.8)
        ax.text(x0 + 12.5, 55.3, title, ha="center", va="top", fontsize=10.5, fontweight="bold", color=col)
        ax.text(x0 + 12.5, 52.0, f"{len(main)} lineages, {n} lines" + (f"\n(incl. {int(rare.n_lines.sum())} lines of rare lineages)" if len(rare) else ""),
                ha="center", va="top", fontsize=8.6)
        lines = [f"{r.lineage} ({r.n_lines})" for r in main.itertuples()]
        ax.text(x0 + 1.2, 47.6 if len(rare) else 49.0, "\n".join(lines), ha="left", va="top", fontsize=7.4, color=MUTED, linespacing=1.15)
        ax.annotate("", xy=(x0 - 0.3, 31.5), xytext=(16.2, 31.5), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.0,
                    connectionstyle="arc3,rad=0.0")) if t == "train" else None
    # what each tier is used for
    ax.text(32, 23.6, "all development:\nlabels via a fold-scoped API;\n12 leave-one-lineage-out folds;\n"
            "training gate: gain ≥ 0.005,\nCI > 0, ≥ 9/12 folds positive", ha="center", va="top", fontsize=7.8,
            bbox=dict(boxstyle="round,pad=0.4", fc="#eef4fc", ec="none"))
    ax.text(59, 23.6, "≤ 3 candidates per look;\nscalar-only results,\nhash-chained ledger;\n"
            f"cap 40 looks, {len(vlooks)} used", ha="center", va="top", fontsize=7.8,
            bbox=dict(boxstyle="round,pad=0.4", fc="#fdf0ea", ec="none"))
    ax.text(86, 23.6, "frozen harness, run by the PI;\nT1: headline chain\nT2: secondary contrasts\n"
            "then unseal (tables re-derived\nwith matching SHA-256)", ha="center", va="top", fontsize=7.8,
            bbox=dict(boxstyle="round,pad=0.4", fc="#e9f7f1", ec="none"))
    # look strip
    ax.plot([19.5, 98.5], [5, 5], color=GRID, lw=1)
    ax.text(19.5, 7.3, "Held-out looks, in order:", fontsize=8.5, fontweight="bold", va="bottom")
    names = {"W1_M2_baseline": "1 generic\nbaseline", "W3_M4_shortlist": "2 biomass\nfix", "M6_factorial": "3 factorial\n→ P1",
             "W7_shortlist_1": "4 ML\nshortlist", "W7_P2_select": "5 ML\nrunner-up", "REHEARSAL_P1_M6": "6 test-harness\nrehearsal",
             "M10_F2": "7 safer variant\n(post hoc)", "T1": "T1\nheadline", "T2": "T2\nsecondary"}
    order = vlooks[:6] + tlooks + vlooks[6:]
    xs = np.linspace(22, 96, len(order))
    for x, lk in zip(xs, order):
        test = lk in tlooks
        ax.plot([x], [5], marker="D" if test else "s", ms=8, color=S3 if test else S2, mfc=S3 if test else CARD, mew=1.5)
        ax.text(x, 2.9, names.get(lk, lk), ha="center", va="top", fontsize=7.3, color=INK)
    ax.text(0.5, -6.0, "Before any score: design pre-registered (frozen commit 61b0beb) after 3 review rounds; a hook blocks "
            "direct reads of label files.", fontsize=7.8, color=MUTED)
    DATA["fig1"] = {"tiers": info, "validation_looks": vlooks, "test_looks": tlooks, "seed": meta["seed"],
                    "sources": ["results/splits/split_summary.tsv", "results/splits/split_meta.json",
                                "results/validation_ledger.tsv", "results/test_ledger.tsv"]}
    save(fig, out, "fig1_benchmark_design")


# ------------------------------------------------------------------------------------------ fig 2
def tier_tables(repo):
    R = "results/"
    rows = {}
    c = tsv(repo, R + "m6/train/cells.tsv").set_index("cell")
    g = tsv(repo, R + "w1/train_grid.tsv")
    g = g[(g.label == "L_prob") & (g.geneset == "all") & (g.universe == "U") & (g.stratum == "none")
          & (g.lineset == "expr")].set_index("method")
    tr = {"H2": c.loc["H2"], "H1": c.loc["H1"], "M2": c.loc["H1/M2"], "BLOO": c.loc["H1+Bloo/M0"],
          "P1": c.loc["H1+Bloo/M2"], "B3": g.loc["B3"], "B8+": g.loc["B8+_H1"], "B2": g.loc["B2[L_prob]"]}
    rows["train"] = {k: (float(v.tap_lin), float(v.tap_lin_lo), float(v.tap_lin_hi)) for k, v in tr.items()}
    v = tsv(repo, R + "m6/val/val_look_M6_factorial.tsv", comment="#")
    v = v[(v.section == "tAP_grid") & (v.label == "L_prob") & (v.geneset == "all") & (v.stratum == "none")].set_index("method")
    w = tsv(repo, R + "w1/val_look_W1_M2_baseline.tsv", comment="#")
    w = w[(w.section == "tAP_grid") & (w.label == "L_prob") & (w.geneset == "all") & (w.stratum == "none")
          & (w.universe == "U")].set_index("method")
    va = {"H2": v.loc["H2"], "H1": v.loc["H1"], "M2": v.loc["M6:H1/M2(f=1)"], "BLOO": v.loc["M6:H1+Bloo/M0"],
          "P1": v.loc["M6:H1+Bloo/M2(f=1)"], "B3": w.loc["B3"], "B8+": w.loc["B8+_H1"], "B2": w.loc["B2[L_prob]"]}
    rows["val"] = {k: (float(x.tap_lin), float(x.tap_lin_lo), float(x.tap_lin_hi)) for k, x in va.items()}
    b = tsv(repo, R + "m6/test/unsealed/T1/w1_benchmark.tsv", comment="#")
    b = b[(b.label == "L_prob") & (b.geneset == "all") & (b.stratum == "none") & (b.universe == "U")]
    b = b.drop_duplicates("method").set_index("method")
    te = {"H2": "H2", "H1": "H1", "M2": "M2[H1]", "P1": "P1[H1]", "B3": "B3", "B8+": "B8+_H1", "B2": "B2[L_prob]"}
    rows["test"] = {k: (float(b.loc[m].tap_lin), float(b.loc[m].tap_lin_lo), float(b.loc[m].tap_lin_hi)) for k, m in te.items()}
    rows["test"]["B0"] = (float(b.loc["B0"].tap_lin), None, None)
    return rows


TIERS = {"train": ("Training (12 lineages)", S1, "o", -0.22), "val": ("Validation (4)", S2, "s", 0.0),
         "test": ("Test (7)", S3, "D", 0.22)}


def fig_tiers(repo, out):
    rows = tier_tables(repo)
    items = [("B2", "Gene prior: fraction of training lines\nin which the gene is essential (no model)"),
             ("B8+", "Generic Human1, ties broken by expression"),
             ("P1", "P1: both fixes (final method)"), ("BLOO", "Biomass fix only"), ("B3", "Expression only (no model)"),
             ("M2", "Medium fix only"), ("H1", "Generic Human1 on Ham's medium (reference)"), ("H2", "Generic Human2 on Ham's medium")]
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ys = np.arange(len(items))[::-1]
    for y, (k, lab) in zip(ys, items):
        for t, (tl, col, mk, dy) in TIERS.items():
            if k not in rows[t]:
                ax.text(0.15, y + dy, "not scored on test", fontsize=7.5, color=MUTED, va="center")
                continue
            est, lo, hi = rows[t][k]
            ax.plot([lo, hi], [y + dy] * 2, color=col, lw=1.6)
            ax.plot([est], [y + dy], marker=mk, ms=6, color=col, lw=0)
        if k == "P1":
            ax.axhspan(y - 0.45, y + 0.45, xmax=0.42, color=BAND, zorder=-5, lw=0)
    ax.axvline(rows["test"]["B0"][0], color=MUTED, lw=1, ls=(0, (2, 3)))
    ax.text(rows["test"]["B0"][0] + 0.008, len(items) - 0.35, "random ranking\n(≈ prevalence, 0.12)", fontsize=7.5, color=MUTED, va="top")
    ax.set_yticks(ys)
    ax.set_yticklabels([lab for _, lab in items], fontsize=8.8)
    ax.set_xlim(0.1, 0.95)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("tie-aware AUPRC (tAP), lineage-equal mean, 95% CI (line bootstrap)")
    h = [plt.Line2D([], [], color=c, marker=m, lw=1.6, ms=6, label=t) for t, c, m, _ in TIERS.values()]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.66, 1.0), fontsize=8.8)
    # inset zoom on GEM arms
    ins = ax.inset_axes([0.50, 0.08, 0.28, 0.52])
    zoom = [("P1", "P1"), ("BLOO", "biomass"), ("M2", "medium"), ("H1", "Human1"), ("H2", "Human2")]
    zy = np.arange(len(zoom))[::-1]
    for y, (k, lab) in zip(zy, zoom):
        for t, (tl, col, mk, dy) in TIERS.items():
            if k in rows[t]:
                est, lo, hi = rows[t][k]
                ins.plot([lo, hi], [y + dy] * 2, color=col, lw=1.3)
                ins.plot([est], [y + dy], marker=mk, ms=4.5, color=col, lw=0)
    ins.set_yticks(zy)
    ins.set_yticklabels([l for _, l in zoom], fontsize=7.5)
    ins.set_xlim(0.28, 0.36)
    ins.tick_params(labelsize=7, length=2)
    ins.grid(axis="x")
    ins.set_title("zoom: GEM arms", fontsize=7.8, loc="left")
    for sp in ins.spines.values():
        sp.set_visible(True)
        sp.set_color(GRID)
    DATA["fig2"] = {"values": rows, "sources": ["results/m6/train/cells.tsv", "results/w1/train_grid.tsv",
                    "results/m6/val/val_look_M6_factorial.tsv", "results/w1/val_look_W1_M2_baseline.tsv",
                    "results/m6/test/unsealed/T1/w1_benchmark.tsv"]}
    save(fig, out, "fig2_tap_tiers")


# ------------------------------------------------------------------------------------------ fig 3
def fig_forest(repo, out):
    R = "results/"
    rows = []  # (group, label, est, lo, hi, val_est, val_lo, val_hi, kind, note)

    def add(group, label, est, lo=None, hi=None, val=None, kind="train", note=""):
        rows.append(dict(group=group, label=label, est=float(est), lo=None if lo is None else float(lo),
                         hi=None if hi is None else float(hi), val=val, kind=kind, note=note))

    g = tsv(repo, R + "w1/train_grid.tsv")
    g = g[(g.label == "L_prob") & (g.geneset == "all") & (g.universe == "U") & (g.stratum == "none") & (g.lineset == "expr")].set_index("method")
    vg = js(repo, R + "w1/val_look_W1_M2_baseline_gstar.json")
    G0 = "Base model (vs generic Human1)"
    add(G0, "Human2 instead of Human1", g.loc["H2"].delta_lin_vs_H1, g.loc["H2"].delta_lo_vs_H1, g.loc["H2"].delta_hi_vs_H1,
        val=(vg["delta_lin_H2_vs_H1"], *vg["delta_ci"]))

    w4 = js(repo, R + "w4/stage_a_analysis/train_selection.json")["families"]
    G1 = "Context-specific extraction from expression (30-line panel, vs Human1)"
    for k, lab in [("F0", "Thresholding + repair"), ("F2", "GIMME"), ("F3", "FASTCORE"), ("F4", "E-Flux"),
                   ("F5", "Label-calibrated presence")]:
        d = w4[k]["vs_H1"]
        add(G1, lab, d["outer_delta_lin"], d["lo"], d["hi"])
    cor = tsv(repo, R + "w4/stage_a_analysis/train_corda_a6prime.tsv")
    cor = cor[(cor.rule == "R1") & (cor.family == "F3b")].iloc[0]
    add(G1, "CORDA (10 lines, no CI)", cor.delta_vs_H1)

    w5 = {r["method"]: r for r in js(repo, R + "w5/train_summary.json")["primary"]}
    G2 = "GPR curation (vs Human1)"
    for k, lab in [("H1_W5a", "Evidence-based edits as pre-registered"), ("H1_W5a_pid30", "Paralogs restricted to ≥ 30% identity"),
                   ("H1_TYMS", "Thymidylate synthase irreversible")]:
        d = w5[k]
        add(G2, lab, d["delta_lin"], d["delta_lo"], d["delta_hi"])

    st = tsv(repo, R + "w2/train/step_scores.tsv")
    st = st[(st.stratum == "attribution") & (st.parent == "human1")].set_index("step")
    sh = tsv(repo, R + "w2/train/shuffled_control.tsv")
    sh = sh[(sh.scope == "pooled") & (sh.parent == "human1")].set_index("step")
    m6v = {c["method"]: c for c in js(repo, R + "m6/val/val_look_M6_factorial_decision.json")["candidates"]}
    G3 = "Medium (466 lines with a documented medium, vs Ham's)"
    for k, lab, val in [("M1", "Documented base medium", None),
                        ("M2(f=1)", "Documented medium + serum", m6v["M6:H1/M2(f=1)"]),
                        ("M-HPLM:M2(f=1)", "Physiological: HPLM + serum", None),
                        ("M-Plasmax:M2(f=1)", "Physiological: Plasmax + serum", None),
                        ("M-open", "Control: every uptake open", None)]:
        r = st.loc[k]
        add(G3, lab, r.delta_lin_vs_M0, r.delta_lo_vs_M0, r.delta_hi_vs_M0,
            val=None if val is None else (val["delta_lin"], *val["delta_ci"]), kind="control" if k == "M-open" else "train")
    add(G3, "Control: medium shuffled between lines (no CI)", sh.loc["M2(f=1)"].shuffled_minus_m0_mean, kind="control")

    m = tsv(repo, R + "w3/train/mode_scores.tsv")
    m = m[m.parent == "human1"].set_index("mode")
    w3v = js(repo, R + "w3/val/val_look_W3_M4_shortlist_decision.json")["candidates"][0]
    G4 = "Objective: biomass, tasks, alternative objectives (vs Human1)"
    for k, lab, val in [("B-loo:cofactor_pool_biomass", "Remove cofactor pool from biomass", w3v),
                        ("B-loo:metabolite_pool_biomass", "Remove metabolite pool from biomass", None),
                        ("T3", "Essential-task calls (57-task family)", None),
                        ("T5", "CellFie-inspired task calls", None),
                        ("O-MOMA", "MOMA objective", None)]:
        r = m.loc[k]
        add(G4, lab, r.delta_lin, r.delta_lo, r.delta_hi, val=None if val is None else (val["delta_lin"], *val["delta_ci"]))
    strip = m.delta_lin.astype(float).tolist()

    tr = js(repo, R + "m6/train/summary.json")["tg"]
    f2c = tsv(repo, R + "f2/train/contrasts.tsv").set_index("contrast")
    f2v = js(repo, R + "m10/val/val_look_M10_F2_decision.json")["candidates"][0]
    G5 = "Combination (vs generic Human1)"
    add(G5, "P1 = biomass fix + medium", tr["H1+Bloo/M2"]["delta_lin"], *tr["H1+Bloo/M2"]["ci"],
        val=(m6v["M6:H1+Bloo/M2(f=1)"]["delta_lin"], *m6v["M6:H1+Bloo/M2(f=1)"]["delta_ci"]))
    r = f2c[f2c.geneset == "all"].loc["F2_vs_Gstar"]
    add(G5, "Safer biomass fix + medium (post hoc)", r.delta_lin, r.ci_lo, r.ci_hi, val=(f2v["delta_lin"], *f2v["delta_ci"]))

    m5 = js(repo, R + "m5/train/train_decision.json")["families"]
    f1 = js(repo, R + "f1/train/summary.json")["F1_vs_P1"]
    f4a = js(repo, R + "f4/train/summary.json")["F4a_vs_P1"]
    f4b = js(repo, R + "f4/eflux_panel/train/decision.json")
    G6 = "Added on top of the biomass fix or of P1 (vs that pipeline)"
    add(G6, "GIMME on biomass fix (panel)", m5["F2"]["delta_lin_vs_H1_Bloo"], *m5["F2"]["ci"])
    add(G6, "E-Flux on biomass fix (panel)", m5["F4"]["delta_lin_vs_H1_Bloo"], *m5["F4"]["ci"])
    add(G6, "Stricter GPR curation on P1 (post hoc)", f1["delta_lin"], *f1["ci"])
    add(G6, "Genomic-loss constraint on P1 (post hoc)", f4a["delta_lin"], *f4a["ci"])
    add(G6, "E-Flux on P1 (panel, post hoc)", f4b["delta_lin_vs_P1"], *f4b["ci"])

    ar = tsv(repo, "autoresearch/results.tsv").set_index("id")
    G7 = "Machine learning on gene priors (vs best gene-prior model)"
    for k, lab in [("W7-02", "Add GEM knockout features"), ("W7-04", "Stack P1 with the gene-prior model"),
                   ("W7-14", "Add sign-constrained GEM features")]:
        x = ar.loc[k]
        add(G7, lab, float(x.d_vs_comparator), float(x.ci_lo), float(x.ci_hi))

    groups = list(dict.fromkeys(r["group"] for r in rows))
    ypos, y, labels, heads = [], 0, [], []
    for gname in groups:
        heads.append((gname, y))
        y -= 1
        for r in [r for r in rows if r["group"] == gname]:
            r["y"] = y
            labels.append((y, r["label"]))
            y -= 1
        if gname == G4:
            strip_y = y
            y -= 1
        y -= 0.4
    fig, ax = plt.subplots(figsize=(9.8, 11.2))
    ax.axvspan(-0.04, 0, color="#f6f5f1", zorder=-6, lw=0)
    ax.axvline(0, color=MUTED, lw=1)
    ax.axvline(0.005, color=MUTED, lw=0.9, ls=(0, (4, 3)))
    ax.text(0.0055, 0.6, "training gate\n(+0.005)", fontsize=7.5, color=MUTED, va="bottom")
    for r in rows:
        col = NEUTRAL if r["kind"] == "control" else (S1 if r["est"] >= 0.005 and (r["lo"] is None or r["lo"] > 0) else INK)
        if r["lo"] is not None:
            ax.plot([r["lo"], r["hi"]], [r["y"] + 0.12] * 2, color=col, lw=1.6)
        ax.plot([r["est"]], [r["y"] + 0.12], marker="o", ms=6.5, color=col, mfc=col if r["lo"] is not None else CARD, mew=1.4, lw=0)
        if r["val"] is not None:
            ve, vlo, vhi = r["val"]
            ax.plot([vlo, vhi], [r["y"] - 0.2] * 2, color=S2, lw=1.4)
            ax.plot([ve], [r["y"] - 0.2], marker="s", ms=5.5, color=S2, lw=0)
        ax.text(0.0675, r["y"], f"{r['est']:+.4f}", fontsize=7.6, va="center", ha="right", color=MUTED)
    ax.plot(strip, [strip_y] * len(strip), "|", color=NEUTRAL, ms=9, mew=1)
    labels.append((strip_y, f"All {len(strip)} Human1 objective modes (one tick each)"))
    for gname, gy in heads:
        ax.text(-0.0415, gy, gname, fontsize=8.8, fontweight="bold", va="center", ha="left")
    ax.set_yticks([l[0] for l in labels])
    ax.set_yticklabels([l[1] for l in labels], fontsize=8.4)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(y + 0.6, 1.4)
    ax.set_xlim(-0.0415, 0.068)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.set_xlabel("change in tie-aware AUPRC versus the comparator named in each group (95% CI)")
    h = [plt.Line2D([], [], color=S1, marker="o", lw=1.6, label="training: gain ≥ 0.005, CI > 0"),
         plt.Line2D([], [], color=INK, marker="o", lw=1.6, label="training: below that bar"),
         plt.Line2D([], [], color=NEUTRAL, marker="o", lw=1.6, label="negative control"),
         plt.Line2D([], [], color=S2, marker="s", lw=1.4, label="held-out validation look")]
    ax.legend(handles=h, loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=4, fontsize=8.4)
    DATA["fig3"] = {"rows": rows, "w3_human1_M0_mode_deltas": strip, "sources": [
        "results/w1/train_grid.tsv", "results/w1/val_look_W1_M2_baseline_gstar.json",
        "results/w4/stage_a_analysis/train_selection.json", "results/w4/stage_a_analysis/train_corda_a6prime.tsv",
        "results/w5/train_summary.json", "results/w2/train/step_scores.tsv", "results/w2/train/shuffled_control.tsv",
        "results/w3/train/mode_scores.tsv", "results/w3/val/val_look_W3_M4_shortlist_decision.json",
        "results/m6/train/summary.json", "results/m6/val/val_look_M6_factorial_decision.json",
        "results/f2/train/contrasts.tsv", "results/m10/val/val_look_M10_F2_decision.json",
        "results/m5/train/train_decision.json", "results/f1/train/summary.json", "results/f4/train/summary.json",
        "results/f4/eflux_panel/train/decision.json", "autoresearch/results.tsv"]}
    save(fig, out, "fig3_levers_forest")


# ------------------------------------------------------------------------------------------ fig 4
def fig_factorial(repo, out):
    R = "results/"
    tr = js(repo, R + "m6/train/summary.json")["tg"]
    vd = {c["method"]: (c["delta_lin"], *c["delta_ci"]) for c in js(repo, R + "m6/val/val_look_M6_factorial_decision.json")["candidates"]}
    at = tsv(repo, R + "m6/test/unsealed/T1/attribution.tsv").set_index("contrast")
    it = tsv(repo, R + "m6/train/interaction.tsv")
    it = it[(it.parent == "H1") & (it.label_cell == "all")].iloc[0]
    cells = {
        "Training (686 lines)": {"bio": (tr["H1+Bloo/M0"]["delta_lin"], *tr["H1+Bloo/M0"]["ci"]),
                                 "med": (tr["H1/M2"]["delta_lin"], *tr["H1/M2"]["ci"]),
                                 "both": (tr["H1+Bloo/M2"]["delta_lin"], *tr["H1+Bloo/M2"]["ci"])},
        "Validation (198 lines)": {"bio": vd["M6:H1+Bloo/M0"], "med": vd["M6:H1/M2(f=1)"], "both": vd["M6:H1+Bloo/M2(f=1)"]},
        "Test (380 lines)": {"bio": None, "med": tuple(at.loc["ii", ["estimate", "ci_lo", "ci_hi"]].astype(float)),
                             "both": tuple(at.loc["i", ["estimate", "ci_lo", "ci_hi"]].astype(float))},
    }
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.6), sharey=True)
    for ax, (tier, c) in zip(axes, cells.items()):
        ax.plot([0, 1], [0, c["med"][0]], color=NEUTRAL, lw=2, marker="o", ms=7, label="biomass as released")
        if c["bio"] is not None:
            ax.plot([0, 1], [c["bio"][0], c["both"][0]], color=S1, lw=2, marker="o", ms=7, label="cofactor pool removed")
            ax.plot([1], [c["bio"][0] + c["med"][0]], marker="o", ms=8, mfc="none", mec=MUTED, mew=1.2, lw=0, label="additive prediction")
            ax.annotate(f"{c['bio'][0]:+.3f}", (0, c["bio"][0]), xytext=(-8, 0), textcoords="offset points", ha="right", va="center", fontsize=8.5)
        else:
            ax.plot([1], [c["both"][0]], color=S1, marker="o", ms=7, lw=0)
            ax.text(0.03, 0.038, "biomass fix alone\nnot scored on test", fontsize=8, color=MUTED)
        for key, col in (("med", NEUTRAL), ("both", S1)):
            est, lo, hi = c[key]
            ax.plot([1, 1], [lo, hi], color=col, lw=1.4)
            ax.annotate(f"{est:+.3f}", (1, est), xytext=(8, 0), textcoords="offset points", va="center", fontsize=8.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Ham's", "documented\n+ serum"])
        ax.set_xlim(-0.5, 1.6)
        ax.set_title(tier, fontsize=9.5, loc="left")
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    axes[0].set_ylabel("change in tAP vs generic Human1")
    axes[0].set_ylim(-0.004, 0.06)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper left", bbox_to_anchor=(0.07, 1.08), ncol=3, fontsize=8.8)
    DATA["fig4"] = {"cells": cells, "training_interaction": {k: float(it[k]) for k in ("interaction_lin", "lo", "hi", "additive_prediction")},
                    "sources": ["results/m6/train/summary.json", "results/m6/train/interaction.tsv",
                                "results/m6/val/val_look_M6_factorial_decision.json", "results/m6/test/unsealed/T1/attribution.tsv"]}
    save(fig, out, "fig4_factorial")


# ------------------------------------------------------------------------------------------ fig 5
def fig_lineage(repo, out):
    p = tsv(repo, "results/m6/test/unsealed/T1/per_lineage.tsv")
    gates = js(repo, "results/m6/test/test_look_T1_gates.json")["transfer_C1"]
    lins = sorted(p.lineage.unique(), reverse=True)
    n = p[p.method == "P1[H1]"].set_index("lineage").n_lines
    fig, ax = plt.subplots(figsize=(7.6, 3.8))
    ax.axvspan(gates["t_low"], gates["t_high"], color=BAND, lw=0, zorder=-5)
    ax.axvline(0, color=MUTED, lw=1)
    vals = {}
    for mth, lab, col, mk, hol in [("P1[H1]", "P1 (final, on Human1)", S1, "o", False), ("P1[H2]", "Same fixes on Human2", S2, "s", False),
                                   ("H2", "Generic Human2", NEUTRAL, "D", True)]:
        d = p[p.method == mth].set_index("lineage").delta_vs_ref
        vals[mth] = d.to_dict()
        ax.plot(d.reindex(lins).to_numpy(), np.arange(len(lins)), lw=0, marker=mk, ms=7, color=col, mfc="none" if hol else col, mew=1.4, label=lab)
    ax.set_yticks(np.arange(len(lins)))
    ax.set_yticklabels([f"{l} ({int(n[l])} lines)" for l in lins])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.015, 0.062)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.set_xlabel("change in tAP vs generic Human1 on Ham's medium, test tier")
    ax.text(0.004, (len(lins) - 1) / 2, f"P1 > 0 in {gates['k_positive']}/7 tissues\nband: lineage t-interval\n"
            f"[{gates['t_low']:.3f}, {gates['t_high']:.3f}]\nexact sign-flip p = 1/128", fontsize=8, color=MUTED, va="center")
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=3, fontsize=8.6)
    DATA["fig5"] = {"delta_vs_generic_human1": vals, "transfer": gates,
                    "sources": ["results/m6/test/unsealed/T1/per_lineage.tsv", "results/m6/test/test_look_T1_gates.json"]}
    save(fig, out, "fig5_lineage_gain")


# ------------------------------------------------------------------------------------------ fig 6
def fig_mcc(repo, out):
    t = tsv(repo, "results/m6/test/test_look_T1.tsv", comment="#")
    b = t[(t.section == "binary_6.8") & (t.nan_policy == "exclude") & (t.gene_set == "parent")]
    prots = [("Human1-paper", "P57", "L_ge06", "Human1-paper protocol:\nfails any of 57 tasks,\ngene effect < −0.6"),
             ("Human1-paper", "PGR", "L_ge06", "Same, growth-only calls\n(biomass < 1% of normal)"),
             ("fixed-cutoff", "FBA<0.01", "L_prob", "Growth-only calls vs\nthis study's label\n(dependency prob. > 0.5)")]
    meths = [("H1", "Generic Human1", NEUTRAL, "o"), ("H2", "Generic Human2", S2, "s"), ("P1[H1]", "P1 (Human1)", S1, "D"),
             ("P1[H2]", "P1 on Human2", S3, "^")]
    fig, (a, b2) = plt.subplots(1, 2, figsize=(10.2, 3.9), gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.55})
    vals = {}
    for i, (pr, cr, lab, _) in enumerate(prots):
        y = len(prots) - 1 - i
        for j, (m, ml, col, mk) in enumerate(meths):
            r = b[(b.protocol == pr) & (b.call_rule == cr) & (b.label == lab) & (b.method == m) & (b.metric == "mcc")].iloc[0]
            vals[f"{pr}|{cr}|{lab}|{m}"] = [float(r.lineage_equal_mean), float(r.ci_lo), float(r.ci_hi)]
            yy = y + 0.3 - j * 0.2
            a.plot([r.ci_lo, r.ci_hi], [yy, yy], color=col, lw=1.5)
            a.plot([r.lineage_equal_mean], [yy], marker=mk, ms=6, color=col, lw=0, label=ml if i == 0 else None)
    a.axvspan(0.25, 0.45, color="#f6f5f1", zorder=-6, lw=0)
    a.text(0.252, 2.62, "range reported in published GEM studies\n(≈ 0.25–0.45; best-of-settings, no tissue hold-out)", fontsize=7.4, color=MUTED, va="bottom")
    a.set_yticks(range(len(prots))[::-1])
    a.set_yticklabels([p[3] for p in prots], fontsize=8.2)
    a.tick_params(axis="y", length=0)
    a.set_xlim(0.2, 0.5)
    a.set_ylim(-0.5, 3.15)
    a.grid(axis="x")
    a.set_axisbelow(True)
    a.set_xlabel("Matthews correlation (MCC), test tier, 95% CI")
    a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.04), ncol=4, fontsize=8, handletextpad=0.3, columnspacing=0.9)
    mets = [("precision", "precision"), ("recall", "recall"), ("f1", "F1"), ("mcc", "MCC")]
    pv = {}
    for k, (m, ml, col, mk) in enumerate([meths[0], meths[2]]):
        for i, (met, lab) in enumerate(mets):
            r = b[(b.protocol == "Human1-paper") & (b.call_rule == "P57") & (b.label == "L_ge06") & (b.method == m) & (b.metric == met)].iloc[0]
            pv[f"{m}|{met}"] = float(r.lineage_equal_mean)
            xb = i + (k - 0.5) * 0.36
            b2.bar(xb, r.lineage_equal_mean, width=0.34, color=col, label=ml if i == 0 else None)
            b2.plot([xb, xb], [r.ci_lo, r.ci_hi], color=INK, lw=1.1)
            pv[f"{m}|{met}|ci"] = [float(r.ci_lo), float(r.ci_hi)]
            b2.text(xb, r.ci_hi + 0.008, f"{r.lineage_equal_mean:.3f}", ha="center", va="bottom", rotation=90, fontsize=7.2, color=INK)
    b2.set_xticks(range(len(mets)))
    b2.set_xticklabels([m[1] for m in mets])
    b2.set_ylim(0, 0.6)
    b2.grid(axis="y")
    b2.set_axisbelow(True)
    b2.set_title("Human1-paper protocol (57 tasks), 95% CI", fontsize=9, loc="left")
    b2.legend(loc="upper right", fontsize=8)
    DATA["fig6"] = {"mcc": vals, "p57_components": pv, "sources": ["results/m6/test/test_look_T1.tsv (section binary_6.8, nan_policy exclude, gene_set parent)"]}
    save(fig, out, "fig6_mcc_protocols")


# ------------------------------------------------------------------------------------------ fig 7
def fig_decomp(repo, out):
    d = tsv(repo, "results/f2/decomposition.tsv")
    groups = [("cofactor_pool", "Cofactor-pool genes (35)"), ("sterol", "Sterol synthesis (12)"),
              ("medium_flip", "Medium-dependent genes (13)"), ("other", "All other genes (2,802)"),
              ("mev_core", "Mevalonate core (4)"), ("residual_total_minus_sum", "Interaction (total − sum)")]
    fig, ax = plt.subplots(figsize=(8.0, 3.9))
    y = np.arange(len(groups))[::-1]
    vals = {}
    for off, meth, col, lab in ((0.19, "P1", S1, "P1 (whole pool removed)"), (-0.19, "F2", S2, "Safer variant: artefacts only (post hoc)")):
        sub = d[d.method == meth].set_index("group")
        v = [float(sub.loc[gk, "contribution"]) for gk, _ in groups]
        vals[meth] = dict(zip([gk for gk, _ in groups], v))
        ax.barh(y + off, v, height=0.34, color=col, label=lab)
        for yy, vv, gk in zip(y + off, v, [gk for gk, _ in groups]):
            lo, hi = sub.loc[gk, "ci_lo"], sub.loc[gk, "ci_hi"]
            if pd.notna(lo):
                ax.plot([float(lo), float(hi)], [yy, yy], color=INK, lw=1)
            xa = max(vv, 0, float(hi) if pd.notna(hi) else vv)
            ax.annotate(f"{vv:+.4f}", (xa, yy), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8)
    tot = d[d.group == "total"].set_index("method").contribution
    ax.axvline(0, color=MUTED, lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, lab in groups])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.004, 0.04)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.set_xlabel(f"contribution to the training gain vs generic Human1 (totals: P1 {tot['P1']:+.4f}, safer variant {tot['F2']:+.4f})", fontsize=8.8)
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=2, fontsize=8.6)
    DATA["fig7"] = {"contribution": vals, "totals": tot.to_dict(), "sources": ["results/f2/decomposition.tsv"]}
    save(fig, out, "fig7_decomposition")


# ------------------------------------------------------------------------------------------ fig 8
def fig_timeline(repo, out):
    fd = js(repo, "docs/explainer/figures/figure_data.json")["fig5"]["hours_since_day1_00UTC"]
    D0 = datetime(2026, 10, 26, tzinfo=timezone.utc)
    log = subprocess.run(["git", "-C", str(repo), "log", "--format=%ct"], capture_output=True, text=True, check=True).stdout.split()
    hrs = np.array([(datetime.fromtimestamp(int(t), tz=timezone.utc) - D0).total_seconds() / 3600 for t in log])
    hrs = hrs[(hrs >= 12) & (hrs <= 60.5)]
    fig, (ax0, ax) = plt.subplots(2, 1, figsize=(10.4, 4.9), sharex=True, gridspec_kw={"height_ratios": [1, 2.4], "hspace": 0.08})
    bins = np.arange(12, 61, 1)
    ax0.hist(hrs, bins=bins, color=NEUTRAL, edgecolor=CARD, lw=0.6)
    ax0.set_ylabel("commits\nper hour", fontsize=8.5)
    ax0.grid(axis="y")
    ax0.set_axisbelow(True)
    for a in (ax0, ax):
        for s0, s1, lab in fd["stalls"]:
            a.axvspan(s0, s1, color=STALL, lw=0, zorder=-5)
    lanes = ["Human input", "Design & evaluation", "Levers scored", "Held-out looks"]
    yl = {k: i for i, k in enumerate(lanes[::-1])}
    ax.plot(fd["human"], [yl["Human input"]] * len(fd["human"]), lw=0, marker="v", ms=8, color=S2)
    ax.plot(fd["design"], [yl["Design & evaluation"]] * len(fd["design"]), lw=0, marker="o", ms=7, color=NEUTRAL)
    f2_t = 2 * 24 + 10 + 26 / 60  # safer biomass variant (F2) result, day 3 10:26 UTC: validated, not adopted
    for t, ok in fd["levers"]:
        if abs(t - f2_t) < 1e-6:
            ax.plot([t], [yl["Levers scored"]], marker="o", ms=8, color=S1, mfc=BAND, mew=1.5, lw=0)
            continue
        ax.plot([t], [yl["Levers scored"]], marker="o", ms=8, color=S1, mfc=S1 if ok else CARD, mew=1.5, lw=0)
    looks = []
    for rel, kind in (("results/validation_ledger.tsv", "v"), ("results/test_ledger.tsv", "t")):
        lg = tsv(repo, rel)
        for r in lg[lg.event == "completed"].itertuples():
            ts = datetime.fromisoformat(r.timestamp.replace("+0000", "+00:00"))
            looks.append(((ts - D0).total_seconds() / 3600, kind, r.look_id))
    for t, k, _ in looks:
        ax.plot([t], [yl["Held-out looks"]], marker="s" if k == "v" else "D", ms=7, color=S3 if k == "t" else INK,
                mfc=S3 if k == "t" else CARD, mew=1.3, lw=0)
    for s0, s1, lab in fd["stalls"]:
        if lab:
            ax.text((s0 + s1) / 2, -0.8, lab, ha="center", va="center", fontsize=7.3, color=MUTED)
    notes = [(fd["human"][2], 3.38, "budget cut;\nftINIT dropped"), (fd["levers"][0][0] - 0.6, 1.38, "context models,\nGPR curation fail"),
             (fd["levers"][2][0], 1.38, "biomass fix"), ([t for t, k, l in looks if l == "T1"][0], 0.4, "T1"), (fd["human"][-1], 3.38, "follow-ups\nrequested")]
    for x, yv, txt in notes:
        ax.text(x, yv, txt, fontsize=7.2, color=INK, ha="center", va="bottom")
    ax.set_yticks(list(yl.values()))
    ax.set_yticklabels(list(yl.keys()))
    ax.tick_params(axis="y", length=0)
    ticks = [12, 18, 24, 30, 36, 42, 48, 54, 60]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["day 1\n12:00", "18:00", "day 2\n00:00", "06:00", "12:00", "18:00", "day 3\n00:00", "06:00", "12:00"])
    ax.set_xlim(12, 60.6)
    ax.set_ylim(-1.15, 4.0)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    h = [plt.Line2D([], [], marker="v", color=S2, lw=0, ms=8, label="user message that changed scope"),
         plt.Line2D([], [], marker="o", color=S1, lw=0, ms=8, label="lever kept"),
         plt.Line2D([], [], marker="o", color=S1, mfc=CARD, mew=1.5, lw=0, ms=8, label="lever failed"),
         plt.Line2D([], [], marker="o", color=S1, mfc=BAND, mew=1.5, lw=0, ms=8, label="validated, not adopted"),
         plt.Line2D([], [], marker="s", color=INK, mfc=CARD, lw=0, ms=7, label="validation look"),
         plt.Line2D([], [], marker="D", color=S3, lw=0, ms=7, label="test look"),
         plt.Rectangle((0, 0), 1, 1, color=STALL, label="stall")]
    ax0.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=7, fontsize=7.4, handlelength=1.1, columnspacing=0.8)
    DATA["fig8"] = {"events_hours_since_day1_00UTC": fd, "looks_from_ledgers": looks, "commits_in_window": int(len(hrs)),
                    "sources": ["docs/explainer/figures/figure_data.json (fig5; stall times there come from the session transcript)",
                                "git log commit times", "results/validation_ledger.tsv", "results/test_ledger.tsv"]}
    save(fig, out, "fig8_timeline")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="/data/tmpA/nextgem")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    repo, out = Path(a.repo), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    fig_design(repo, out)
    fig_tiers(repo, out)
    fig_forest(repo, out)
    fig_factorial(repo, out)
    fig_lineage(repo, out)
    fig_mcc(repo, out)
    fig_decomp(repo, out)
    fig_timeline(repo, out)
    (out / "briefing_figure_data.json").write_text(json.dumps(DATA, indent=1, default=float))
    print("ok", sorted(p.name for p in out.glob("*.svg")))


if __name__ == "__main__":
    main()
