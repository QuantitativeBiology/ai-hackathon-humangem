"""Explainer figures for P1, built only from committed result files in /data/tmpA/nextgem (no label file).

Outputs (SVG, text kept as text): <out>/fig1_tap_tiers.svg, fig2_lineage_gain.svg, fig3_factorial.svg,
fig4_decomposition.svg, fig5_timeline.svg, plus a figure_data.json with every plotted number and its source.

Colours are written as CSS variables with fallbacks, so the SVGs follow light/dark themes:
  * standalone files carry a <style> block (class .fx) with light values and a prefers-color-scheme dark override;
  * the HTML page strips that block (between the FX-STYLE markers) and defines the same variables on :root.

Usage: /data/tmpA/nextgem/.venv/bin/python make_figures.py --repo /data/tmpA/nextgem --out <dir>
"""
from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse
import io
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# sentinel colours -> CSS variables (validated reference palette, slots 1-3: all-pairs safe in both modes)
SENT = {
    "#010101": ("--fx-ink", "#1d1d1b", "#ecebe6"),
    "#020202": ("--fx-muted", "#6b6a64", "#a3a29a"),
    "#030303": ("--fx-grid", "#e3e2dc", "#3a3a37"),
    "#0a0a01": ("--fx-s1", "#2a78d6", "#3987e5"),  # P1 / training
    "#0a0a02": ("--fx-s2", "#eb6834", "#d95926"),  # Human2 / validation / F2
    "#0a0a03": ("--fx-s3", "#1baf7a", "#199e70"),  # test
    "#0a0a04": ("--fx-neutral", "#8d8c85", "#8a8983"),  # reference / generic
    "#0a0a05": ("--fx-band", "#d7e6f8", "#1f3550"),  # interval band
    "#0a0a06": ("--fx-stall", "#efe3d6", "#3d3127"),  # stalls
}
PREVIEW = os.environ.get("FX_PREVIEW") in ("light", "dark")
if PREVIEW:  # layout check only: concrete colours, PNG output
    _i = 1 if os.environ["FX_PREVIEW"] == "light" else 2
    INK, MUTED, GRID, S1, S2, S3, NEUTRAL, BAND, STALL = [v[_i] for v in SENT.values()]
else:
    INK, MUTED, GRID, S1, S2, S3, NEUTRAL, BAND, STALL = list(SENT)

plt.rcParams.update({
    "svg.fonttype": "none", "font.family": "sans-serif", "font.size": 10.5,
    "font.sans-serif": ["DejaVu Sans"], "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8, "axes.grid": False,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.facecolor": "none", "figure.facecolor": "none", "legend.frameon": False, "svg.hashsalt": "nextgem-explainer",
})

DATA: dict = {}


def style_block() -> str:
    light = ";".join(f"{v}:{l}" for v, l, _ in SENT.values())
    dark = ";".join(f"{v}:{d}" for v, _, d in SENT.values())
    return ("<!--FX-STYLE--><style>.fx{" + light + "}"
            "@media (prefers-color-scheme: dark){.fx{" + dark + "}}"
            ".fx text{font-family:'IBM Plex Sans',system-ui,-apple-system,'Segoe UI',Roboto,Arial,sans-serif}"
            "</style><!--/FX-STYLE-->")


def save(fig, out: Path, name: str, title: str, desc: str):
    if PREVIEW:
        bg = "#ffffff" if os.environ["FX_PREVIEW"] == "light" else "#161615"
        fig.savefig(out / f"{name}.{os.environ['FX_PREVIEW']}.png", dpi=110, bbox_inches="tight", facecolor=bg)
        plt.close(fig)
        return
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", transparent=True, metadata={"Date": None})
    plt.close(fig)
    s = buf.getvalue()
    s = s[s.index("<svg"):]
    for hexv, (var, light, _) in SENT.items():
        s = s.replace(hexv, f"var({var}, {light})")
    s = re.sub(r"font: (bold )?([0-9.]+px) 'DejaVu Sans'(, sans-serif)?",
               lambda m: f"font-size: {m.group(2)}" + ("; font-weight: 600" if m.group(1) else ""), s)
    s = re.sub(r"font-family: 'DejaVu Sans'(, sans-serif)?;? ?", "", s)
    s = s.replace("*{stroke-linejoin", ".fx *{stroke-linejoin")
    # unique ids per figure, so several figures can be inlined in one HTML page
    s = re.sub(r'id="([^"]+)"', lambda m: f'id="{name}-{m.group(1)}"' if not m.group(1).startswith(name) else m.group(0), s)
    s = re.sub(r'href="#([^"]+)"', rf'href="#{name}-\1"', s)
    s = re.sub(r"url\(#([^)]+)\)", rf"url(#{name}-\1)", s)
    s = s.replace("<svg ", f'<svg class="fx" role="img" aria-labelledby="{name}-t {name}-d" ', 1)
    head_end = s.index(">", s.index("<svg")) + 1
    s = (s[:head_end] + f'<title id="{name}-t">{title}</title><desc id="{name}-d">{desc}</desc>'
         + style_block() + s[head_end:])
    s = re.sub(r"<metadata>.*?</metadata>", "", s, flags=re.S)
    (out / f"{name}.svg").write_text(s)


# --------------------------------------------------------------------------------------------- data
def tier_tables(repo: Path) -> dict:
    R = repo / "results"
    rows = {}
    c = pd.read_csv(R / "m6/train/cells.tsv", sep="\t").set_index("cell")
    g = pd.read_csv(R / "w1/train_grid.tsv", sep="\t")
    g = g[(g.label == "L_prob") & (g.geneset == "all") & (g.universe == "U") & (g.stratum == "none")
          & (g.lineset == "expr")].set_index("method")
    f2t = json.loads((R / "f2/train/summary.json").read_text())["tap_lin"]["F2"]["all"]
    tr = {"H2": c.loc["H2"], "H1": c.loc["H1"], "M2": c.loc["H1/M2"], "BLOO": c.loc["H1+Bloo/M0"],
          "P1": c.loc["H1+Bloo/M2"], "P1H2": c.loc["H2+Bloo/M2"], "B3": g.loc["B3"], "B8+": g.loc["B8+_H1"],
          "B2": g.loc["B2[L_prob]"]}
    rows["train"] = {k: (float(v.tap_lin), float(v.tap_lin_lo), float(v.tap_lin_hi)) for k, v in tr.items()}
    rows["train"]["F2"] = (float(f2t), None, None)
    v = pd.read_csv(R / "m6/val/val_look_M6_factorial.tsv", sep="\t", comment="#")
    v = v[(v.section == "tAP_grid") & (v.label == "L_prob") & (v.geneset == "all") & (v.stratum == "none")].set_index("method")
    w = pd.read_csv(R / "w1/val_look_W1_M2_baseline.tsv", sep="\t", comment="#")
    w = w[(w.section == "tAP_grid") & (w.label == "L_prob") & (w.geneset == "all") & (w.stratum == "none")
          & (w.universe == "U")].set_index("method")
    m10 = pd.read_csv(R / "m10/val/val_look_M10_F2.tsv", sep="\t", comment="#")
    m10 = m10[(m10.section == "tAP_grid") & (m10.label == "L_prob") & (m10.geneset == "all") & (m10.stratum == "none")].set_index("method")
    va = {"H2": v.loc["H2"], "H1": v.loc["H1"], "M2": v.loc["M6:H1/M2(f=1)"], "BLOO": v.loc["M6:H1+Bloo/M0"],
          "P1": v.loc["M6:H1+Bloo/M2(f=1)"], "P1H2": v.loc["M6:H2+Bloo/M2(f=1)"], "B3": w.loc["B3"],
          "B8+": w.loc["B8+_H1"], "B2": w.loc["B2[L_prob]"], "F2": m10.loc["F2:H1+artBloo/M2(f=1)"]}
    rows["val"] = {k: (float(x.tap_lin), float(x.tap_lin_lo), float(x.tap_lin_hi)) for k, x in va.items()}
    b = pd.read_csv(R / "m6/test/unsealed/T1/w1_benchmark.tsv", sep="\t", comment="#")
    b = b[(b.label == "L_prob") & (b.geneset == "all") & (b.stratum == "none") & (b.universe == "U")]
    b = b.drop_duplicates("method").set_index("method")
    te = {"H2": "H2", "H1": "H1", "M2": "M2[H1]", "P1": "P1[H1]", "P1H2": "P1[H2]", "B3": "B3", "B8+": "B8+_H1",
          "B2": "B2[L_prob]"}
    rows["test"] = {k: (float(b.loc[m].tap_lin), float(b.loc[m].tap_lin_lo), float(b.loc[m].tap_lin_hi)) for k, m in te.items()}
    DATA["fig1"] = {"values": rows, "sources": [
        "results/m6/train/cells.tsv", "results/w1/train_grid.tsv (lineset expr)", "results/f2/train/summary.json",
        "results/m6/val/val_look_M6_factorial.tsv", "results/w1/val_look_W1_M2_baseline.tsv",
        "results/m10/val/val_look_M10_F2.tsv", "results/m6/test/unsealed/T1/w1_benchmark.tsv"],
        "filter": "L_prob, all genes, universe U, no stratum; lineage-equal tAP with 95% line-bootstrap CI"}
    return rows


TIER_STYLE = {"train": ("Training (12 lineages)", S1, "o", -0.22),
              "val": ("Validation (4 lineages)", S2, "s", 0.0),
              "test": ("Test (7 lineages)", S3, "D", 0.22)}


def fig1(rows, out):
    gem = [("H2", "Generic Human2"), ("H1", "Generic Human1 (reference)"), ("M2", "Medium fix only"),
           ("BLOO", "Biomass fix only"), ("P1", "P1: both fixes"), ("P1H2", "P1 built on Human2"),
           ("F2", "Safer biomass fix + medium (post hoc)")]
    base = [("P1", "P1: both fixes"), ("B3", "Expression only"), ("B8+", "Reference, ties broken by expression"),
            ("B2", "Gene prior (fraction of lines essential)")]
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 4.6), gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.08})
    for ax, items, xl in ((a, gem, (0.28, 0.365)), (b, base, (0.0, 1.0))):
        ys = np.arange(len(items))[::-1]
        for y, (k, lab) in zip(ys, items):
            for tier, (tl, col, mk, dy) in TIER_STYLE.items():
                if k not in rows[tier]:
                    if k == "BLOO" and tier == "test":
                        ax.text(xl[0] + 0.001 if ax is a else 0.01, y + dy, "not scored on test", fontsize=8.5,
                                color=MUTED, va="center")
                    continue
                est, lo, hi = rows[tier][k]
                if lo is not None:
                    ax.plot([lo, hi], [y + dy] * 2, color=col, lw=1.6, solid_capstyle="round")
                hollow = k == "F2"
                ax.plot([est], [y + dy], marker=mk, ms=6.5, color=col, mfc="none" if hollow else col, mew=1.4, lw=0)
            if k == "P1":
                ax.axhspan(y - 0.42, y + 0.42, color=BAND, zorder=-5, lw=0)
        ax.set_yticks(ys)
        ax.set_yticklabels([lab for _, lab in items] if ax is a else [lab for _, lab in items])
        ax.set_xlim(*xl)
        ax.grid(axis="x")
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)
    b.yaxis.tick_right()
    b.spines["right"].set_visible(False)
    b.spines["left"].set_visible(False)
    a.set_xlabel("tie-aware AUPRC (tAP), GEM arms")
    b.set_xlabel("tAP, same lines, full 0–1 scale")
    b.axvline(rows["test"]["H1"][0], color=NEUTRAL, lw=1, ls=(0, (3, 3)))
    b.text(rows["test"]["H1"][0] + 0.015, -0.55, "generic Human1", fontsize=8.5, color=MUTED, va="center")
    handles = [plt.Line2D([], [], color=c, marker=m, lw=1.6, ms=6.5, label=t) for t, c, m, _ in TIER_STYLE.values()]
    a.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, 1.13), ncol=3, fontsize=9.5, handlelength=1.6)
    save(fig, out, "fig1_tap_tiers", "Accuracy by pipeline and tier",
         "Lineage-equal tie-aware AUPRC with 95% CIs on training, validation and test lineages. P1 gains about 0.05 "
         "over generic Human1 on every tier, but sits far below the gene prior and below the expression tie-break.")


def fig2(repo, out):
    R = repo / "results"
    p = pd.read_csv(R / "m6/test/unsealed/T1/per_lineage.tsv", sep="\t")
    gates = json.loads((R / "m6/test/test_look_T1_gates.json").read_text())["transfer_C1"]
    lins = sorted(p.lineage.unique(), reverse=True)
    n = p[p.method == "P1[H1]"].set_index("lineage").n_lines
    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    ax.axvspan(gates["t_low"], gates["t_high"], color=BAND, lw=0, zorder=-5)
    ax.axvline(0, color=MUTED, lw=1)
    series = [("P1[H1]", "P1 (Human1)", S1, "o", False), ("P1[H2]", "P1 built on Human2", S2, "s", False),
              ("H2", "Generic Human2", NEUTRAL, "D", True)]
    vals = {}
    for m, lab, col, mk, hol in series:
        d = p[p.method == m].set_index("lineage").delta_vs_ref
        vals[m] = d.to_dict()
        ax.plot(d.reindex(lins).to_numpy(), np.arange(len(lins)), lw=0, marker=mk, ms=7, color=col,
                mfc="none" if hol else col, mew=1.4, label=lab)
    ax.set_yticks(np.arange(len(lins)))
    ax.set_yticklabels([f"{l}  ({int(n[l])} lines)" for l in lins])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.015, 0.062)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.set_xlabel("change in tAP vs generic Human1 on Ham's medium (same lines)")
    ax.text(0.006, (len(lins) - 1) / 2, f"P1 > 0 in 7/7 lineages\nshaded band: lineage t-interval\n[{gates['t_low']:.3f}, {gates['t_high']:.3f}], exact sign-flip p = 1/128",
            fontsize=8.5, color=MUTED, va="center")
    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.16), ncol=3, fontsize=9.5)
    DATA["fig2"] = {"delta_vs_generic_human1": vals, "transfer": gates,
                    "sources": ["results/m6/test/unsealed/T1/per_lineage.tsv", "results/m6/test/test_look_T1_gates.json"]}
    save(fig, out, "fig2_lineage_gain", "Per-lineage test gain",
         "Change in tAP versus generic Human1 for each of the 7 held-out test lineages. P1 gains between 0.046 and "
         "0.054 in every lineage; generic Human2 is below Human1 in every lineage.")


def fig3(repo, out):
    R = repo / "results"
    tr = json.loads((R / "m6/train/summary.json").read_text())["tg"]
    it = pd.read_csv(R / "m6/train/interaction.tsv", sep="\t")
    it = it[(it.parent == "H1") & (it.label_cell == "all")].iloc[0]
    vd = json.loads((R / "m6/val/val_look_M6_factorial_decision.json").read_text())["candidates"]
    vd = {c["method"]: (c["delta_lin"], *c["delta_ci"]) for c in vd}
    at = pd.read_csv(R / "m6/test/unsealed/T1/attribution.tsv", sep="\t").set_index("contrast")
    cells = {
        "Training": {"bio": (tr["H1+Bloo/M0"]["delta_lin"], *tr["H1+Bloo/M0"]["ci"]),
                     "med": (tr["H1/M2"]["delta_lin"], *tr["H1/M2"]["ci"]),
                     "both": (tr["H1+Bloo/M2"]["delta_lin"], *tr["H1+Bloo/M2"]["ci"])},
        "Validation": {"bio": vd["M6:H1+Bloo/M0"], "med": vd["M6:H1/M2(f=1)"], "both": vd["M6:H1+Bloo/M2(f=1)"]},
        "Test": {"bio": None, "med": tuple(at.loc["ii", ["estimate", "ci_lo", "ci_hi"]].astype(float)),
                 "both": tuple(at.loc["i", ["estimate", "ci_lo", "ci_hi"]].astype(float))},
    }
    fig, axes = plt.subplots(1, 3, figsize=(8.6, 3.5), sharey=True)
    for ax, (tier, c) in zip(axes, cells.items()):
        x = [0, 1]
        ax.plot(x, [0, c["med"][0]], color=NEUTRAL, lw=2, marker="o", ms=7, label="biomass as released")
        if c["bio"] is not None:
            ax.plot(x, [c["bio"][0], c["both"][0]], color=S1, lw=2, marker="o", ms=7, label="cofactor pool removed")
            add = c["bio"][0] + c["med"][0]
            ax.plot([1], [add], marker="o", ms=8, mfc="none", mec=MUTED, mew=1.3, lw=0, label="additive prediction")
            ax.annotate(f"{c['bio'][0]:+.3f}", (0, c["bio"][0]), xytext=(-8, 0), textcoords="offset points",
                        ha="right", va="center", fontsize=9)
        else:
            ax.plot([1], [c["both"][0]], color=S1, marker="o", ms=7, lw=0, label="cofactor pool removed")
            ax.plot([0], [0.0], lw=0)
            ax.text(0.02, 0.040, "biomass fix alone\nnot scored on test", fontsize=8.5, color=MUTED)
        for key, xx, col in (("med", 1, NEUTRAL), ("both", 1, S1)):
            est, lo, hi = c[key]
            ax.plot([xx, xx], [lo, hi], color=col, lw=1.4)
            ax.annotate(f"{est:+.3f}", (xx, est), xytext=(8, 0), textcoords="offset points", va="center", fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(["Ham's", "documented\n+ serum"])
        ax.set_xlim(-0.45, 1.6)
        ax.set_title(tier, fontsize=10.5, loc="left")
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    axes[0].set_ylabel("change in tAP vs generic Human1")
    axes[0].set_ylim(-0.004, 0.06)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper left", bbox_to_anchor=(0.06, 1.07), ncol=3, fontsize=9.5)
    DATA["fig3"] = {"cells": cells, "training_interaction": {k: float(it[k]) for k in
                    ("interaction_lin", "lo", "hi", "additive_prediction", "delta_joint")},
                    "sources": ["results/m6/train/summary.json", "results/m6/train/interaction.tsv",
                                "results/m6/val/val_look_M6_factorial_decision.json",
                                "results/m6/test/unsealed/T1/attribution.tsv"]}
    save(fig, out, "fig3_factorial", "Biomass by medium factorial",
         "Each fix helps on its own and the two add up (slightly more than additive on training: +0.003). The pattern "
         "replicates on validation; on test the medium fix alone gives +0.013 and both give +0.050.")


def fig4(repo, out):
    d = pd.read_csv(repo / "results/f2/decomposition.tsv", sep="\t")
    groups = [("cofactor_pool", "Cofactor-pool genes (35)"), ("sterol", "Sterol synthesis (12)"),
              ("medium_flip", "Medium-dependent genes (13)"), ("other", "All other genes (2,802)"),
              ("mev_core", "Mevalonate core (4)"), ("residual_total_minus_sum", "Interaction (total − sum)")]
    fig, ax = plt.subplots(figsize=(7.6, 3.9))
    y = np.arange(len(groups))[::-1]
    vals = {}
    for off, meth, col in ((0.19, "P1", S1), (-0.19, "F2", S2)):
        sub = d[d.method == meth].set_index("group")
        v = [float(sub.loc[g, "contribution"]) for g, _ in groups]
        lo = [sub.loc[g, "ci_lo"] for g, _ in groups]
        hi = [sub.loc[g, "ci_hi"] for g, _ in groups]
        vals[meth] = dict(zip([g for g, _ in groups], v))
        ax.barh(y + off, v, height=0.34, color=col, label={"P1": "P1 (whole pool removed)",
                                                          "F2": "Safer variant: artefacts only removed (post hoc)"}[meth])
        for yy, vv, l, h in zip(y + off, v, lo, hi):
            if pd.notna(l):
                ax.plot([float(l), float(h)], [yy, yy], color=INK, lw=1)
            ax.annotate(f"{vv:+.4f}", (max(vv, 0), yy), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8.5)
    tot = d[d.group == "total"].set_index("method").contribution
    ax.axvline(0, color=MUTED, lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, lab in groups])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.004, 0.039)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.set_xlabel(f"contribution to training gain vs generic Human1 (totals: P1 {tot['P1']:+.4f}, safer variant {tot['F2']:+.4f})")
    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.15), ncol=2, fontsize=9.5)
    DATA["fig4"] = {"contribution": vals, "totals": tot.to_dict(), "sources": ["results/f2/decomposition.tsv"],
                    "note": "post hoc (A15), training tier; group contribution = tAP change when only that group's scores differ from the reference"}
    save(fig, out, "fig4_decomposition", "Where the gain comes from",
         "Training-tier decomposition of the gain by gene group. Un-calling the 35 cofactor-pool genes gives 62% of "
         "P1's gain; the mevalonate artefact costs 0.0014; the safer variant loses most of the cofactor-pool gain.")


def fig5(out):
    D0 = datetime(2026, 10, 26, tzinfo=timezone.utc)

    def h(day, hh, mm=0):
        return (day - 1) * 24 + hh + mm / 60

    lanes = ["Human input", "Design & evaluation", "Levers scored", "Held-out looks"]
    yl = {k: i for i, k in enumerate(lanes[::-1])}
    fig, ax = plt.subplots(figsize=(10.2, 3.6))
    stalls = [(h(1, 22, 48), h(2, 0, 20), "limit"), (h(2, 3, 15), h(2, 7, 41), "waiting for user"),
              (h(2, 12, 22), h(2, 19, 0), "weekly limit"), (h(2, 23, 52), h(3, 0, 0), ""), (h(3, 3, 54), h(3, 5, 0), "limit")]
    for a, b, lab in stalls:
        ax.axvspan(a, b, color=STALL, lw=0, zorder=-5)
        if lab:
            ax.text((a + b) / 2, -0.85, lab, ha="center", va="center", fontsize=7.5, color=MUTED)
    human = [h(1, 13, 41), h(1, 15, 32), h(1, 16, 30), h(2, 7, 41), h(2, 7, 46), h(2, 8, 3), h(3, 9, 31)]
    design = [h(1, 14, 46), h(1, 15, 21), h(1, 17, 5), h(1, 20, 5), h(2, 0, 50)]
    levers = [(h(2, 19, 46), False), (h(2, 20, 42), False), (h(3, 0, 16), True), (h(3, 1, 31), False),
              (h(3, 1, 46), True), (h(3, 3, 13), False), (h(3, 5, 10), False), (h(3, 10, 12), False),
              (h(3, 10, 26), True), (h(3, 11, 36), False)]
    looks = [(h(2, 0, 54), "v"), (h(3, 0, 47), "v"), (h(3, 2, 26), "v"), (h(3, 5, 28), "v"), (h(3, 5, 31), "v"),
             (h(3, 7, 2), "t"), (h(3, 7, 16), "t"), (h(3, 10, 30), "v")]
    ax.plot(human, [yl["Human input"]] * len(human), lw=0, marker="v", ms=8, color=S2)
    ax.plot(design, [yl["Design & evaluation"]] * len(design), lw=0, marker="o", ms=7, color=NEUTRAL)
    for t, ok in levers:
        ax.plot([t], [yl["Levers scored"]], marker="o", ms=8, color=S1, mfc=S1 if ok else "none", mew=1.5, lw=0)
    for t, k in looks:
        ax.plot([t], [yl["Held-out looks"]], marker="s" if k == "v" else "D", ms=7, color=S3 if k == "t" else INK,
                mfc=S3 if k == "t" else "none", mew=1.3, lw=0)
    tps = [("1", h(1, 15, 55)), ("2", h(2, 0, 52)), ("3", h(2, 19, 46)), ("4", h(2, 20, 42)), ("5", h(3, 0, 30)),
           ("6", h(3, 1, 46)), ("7", h(3, 4, 45)), ("8", h(3, 4, 20)), ("E", h(3, 10, 30))]
    tops = {"7": 4.25, "8": 3.75}
    for lab, t in tps:
        yy = tops.get(lab, 4.0)
        ax.plot([t, t], [-0.45, yy - 0.18], color=GRID, lw=0.8, zorder=-4)
        ax.text(t, yy, lab, ha="center", va="center", fontsize=8.5, fontweight="bold", color=INK,
                bbox=dict(boxstyle="circle,pad=0.25", fc="none", ec=INK, lw=0.9))
    ax.set_yticks(list(yl.values()))
    ax.set_yticklabels(list(yl.keys()))
    ax.tick_params(axis="y", length=0)
    ticks = [h(1, 12), h(1, 18), h(2, 0), h(2, 6), h(2, 12), h(2, 18), h(3, 0), h(3, 6), h(3, 12)]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["day 1\n12:00", "18:00", "day 2\n00:00", "06:00", "12:00", "18:00", "day 3\n00:00", "06:00", "12:00"])
    ax.set_xlim(h(1, 12), h(3, 12, 30))
    ax.set_ylim(-1.15, 4.6)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    handles = [plt.Line2D([], [], marker="o", color=S1, lw=0, ms=8, label="lever kept"),
               plt.Line2D([], [], marker="o", color=S1, mfc="none", mew=1.5, lw=0, ms=8, label="lever failed"),
               plt.Line2D([], [], marker="s", color=INK, mfc="none", lw=0, ms=7, label="validation look"),
               plt.Line2D([], [], marker="D", color=S3, lw=0, ms=7, label="test look"),
               plt.Rectangle((0, 0), 1, 1, color=STALL, label="stall (session transcript)")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, 1.2), ncol=5, fontsize=8.8, handlelength=1.2)
    DATA["fig5"] = {"hours_since_day1_00UTC": {"human": human, "design": design, "levers": levers, "looks": looks,
                                               "stalls": stalls, "turning_points": tps},
                    "sources": ["git log (commit times)", "results/validation_ledger.tsv", "results/test_ledger.tsv",
                                "session transcript (stalls, user messages; not in the repository)"],
                    "day1": D0.isoformat()}
    save(fig, out, "fig5_timeline", "Project timeline",
         "Two and a half days of work. Data-integration levers failed first (day 2 evening); the biomass fix and the "
         "medium fix were found on day 3 after midnight; usage-limit stalls and a wait for the user are shaded.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="/data/tmpA/nextgem")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    repo, out = Path(a.repo), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = tier_tables(repo)
    fig1(rows, out)
    fig2(repo, out)
    fig3(repo, out)
    fig4(repo, out)
    fig5(out)
    if PREVIEW:
        return
    (out / "figure_data.json").write_text(json.dumps(DATA, indent=1, default=float))
    print("ok", sorted(p.name for p in out.glob("*.svg")))


if __name__ == "__main__":
    main()
