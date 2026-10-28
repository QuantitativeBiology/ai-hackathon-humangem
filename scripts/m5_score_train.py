"""M5 / amendment A12(b): training-tier scoring and decision of the context-model panel pre-check (DESIGN.md §17 A12 b).

Scores the builds of scripts/m5_build.py on the 30 Stage-A panel lines (10 training lineages) with the Stage-A scoring
machinery (scripts/w4_stage_a_train.py: ``load``, ``Scorer``, ``delta_stats``; nextgem.w4_analysis: ``score_frame``,
``nested_lolo_select``, ``tg_panel``). The only adaptation is the parent: scores map Human1 genes to U with the Human1
baseline Entrez map (the W1 generic-table convention), and the binary endpoints use the Human1 gene set.

Arms (L_prob, all genes, U; e = 1 - ratio_r; failed builds take the §11.4 no-information score in both arms):
  generic_H1          W1 generic Human1 (G*)
  generic_H1+Bloo     (1) the H1 + B-loo:cofactor_pool_biomass parent (validated W3 variant)
  F2_R1+Bloo          (2) GIMME R1 on that parent
  F4_R1+Bloo          (3) E-Flux R1 on that parent
Contrasts: (2) vs (1), (3) vs (1) [the A12(b) decision], (1) vs generic_H1 [sanity: W3 gain ≈ +0.03], (2)/(3) vs
generic_H1 [context]. Nested panel LOLO as in Stage A: one configuration per family is fixed (Stage A's), so the inner
choice is that configuration in every fold and the outer per-line Δ equals the plain Δ (asserted). Δ̄_lin is the
lineage-equal mean with a paired line bootstrap stratified by lineage (B = 2,000, SEEDS[0]); folds = panel lineages,
a fold is positive if >= 2 lines are estimable and its mean Δ > 0.

Decision (pre-registered, A12 b): a family enters Stage B only if its combined Δ̄_lin vs H1 + B-loo is >= +0.005.

Labels: only ``labels.train()`` (panel lines) and ``cv.lolo`` / ``labels.fold()`` (CE_fit per held-out panel lineage),
through ``w4_stage_a_train.load``. ``--synthetic-labels`` replaces every label by Bernoulli(0.12) draws (no label file is
opened; a code dry run, written only under a tmp/ directory).

Outputs (default results/m5/train/): train_per_line.tsv, train_config_scores.tsv, train_contrasts.tsv,
train_decision.json, train_verification.tsv; a provenance row in results/m5/MANIFEST.tsv.

Usage:
  .venv/bin/python scripts/m5_score_train.py --builds results/m5/tmp/pilot --out results/m5/tmp/score_dry --synthetic-labels
  .venv/bin/python scripts/m5_score_train.py                      # production (after scripts/m5_build.py)
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
import w4_stage_a_train as w4t  # noqa: E402  (Stage-A scoring machinery, reused unchanged)

from nextgem import metrics, paths  # noqa: E402
from nextgem import w4_analysis as wa  # noqa: E402
from nextgem.config import SEEDS  # noqa: E402

M5 = paths.results_dir() / "m5"
BUILDS = M5 / "builds"
OUT = M5 / "train"
MANIFEST = M5 / "MANIFEST.tsv"
A12B_MIN_DELTA = 0.005  # A12(b): combined Δ̄ vs H1 + B-loo >= +0.005 on the panel
W3_MODE = "B-loo:cofactor_pool_biomass"
BASE = "generic_H1+Bloo"
FAMILY_ARMS = {"F2": "F2_R1+Bloo", "F4": "F4_R1+Bloo"}
CONTRASTS = [  # (arm, reference, role)
    ("F2_R1+Bloo", BASE, "decision"),
    ("F4_R1+Bloo", BASE, "decision"),
    (BASE, "generic_H1", "sanity (W3 B-loo gain)"),
    ("F2_R1+Bloo", "generic_H1", "context"),
    ("F4_R1+Bloo", "generic_H1", "context"),
]
SYNTH_RATE = 0.12


# ----------------------------------------------------------------------------- pure helpers (tested)
def stage_b_verdict(delta_lin: float, min_delta: float = A12B_MIN_DELTA) -> bool:
    """A12(b): the family enters Stage B iff its Δ̄_lin vs H1 + B-loo is >= +0.005 (NaN never enters)."""
    return bool(np.isfinite(delta_lin) and delta_lin >= min_delta - 1e-12)


def h1_entrez_map(baseline_dir: Path | None = None) -> pd.Series:
    """Entrez -> Human1 gene (the W1 generic-table convention, as wa.entrez_gene_map does for Human2)."""
    d = baseline_dir or paths.baseline_dir()
    return wa.entrez_gene_map(pd.read_parquet(d / "human1_ko.parquet", columns=["gene", "entrez"]))


def ratio_tables(ko: pd.DataFrame) -> dict:
    ko = ko.assign(build_id=ko.build_id.astype(str))
    return {bid: g.set_index("gene").ratio_r for bid, g in ko.groupby("build_id")}


def arm_build_ids(lines, family: str) -> dict:
    if family == "generic":
        return {l: "generic|none" for l in lines}
    return {l: f"M5|{l}|{family}|R1|rna|real|none" for l in lines}


def single_config_contrast(delta: pd.Series, lineage: pd.Series, name: str, B: int, seed=SEEDS[0]) -> dict:
    """Stage A's nested panel LOLO with one fixed configuration, then the panel TG statistics on the outer Δ."""
    D = pd.DataFrame({name: delta})
    outer, chosen = wa.nested_lolo_select(D, lineage, order=[name])
    assert set(chosen.values()) <= {name, None}
    pd.testing.assert_series_equal(outer.rename(None), delta.rename(None).astype(float), check_names=False)
    tg = wa.tg_panel(outer, lineage, B=B, seed=seed)
    return {"delta_lin": tg["delta_lin"], "lo": tg["lo"], "hi": tg["hi"], "delta_line": tg.get("delta_line"),
            "folds_positive": tg["folds_positive"], "n_folds": tg["n_folds"], "n_lines": tg["n_lines"],
            "n_lines_estimable": tg["n_lines_estimable"], "fold_means": tg.get("fold_means", {}),
            "tg_rule_pass": tg["tg_pass"]}


# ----------------------------------------------------------------------------- contexts
def synthetic_context() -> dict:
    """Label-free stand-in for ``w4_stage_a_train.load`` (same keys; Bernoulli labels; empty CE_fit)."""
    from nextgem import labels, w1

    split = labels.load_split().set_index("ModelID")
    panel = pd.read_csv(w4t.W4 / "stage_a_panel.tsv", sep="\t")
    lines = list(panel.ModelID)
    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    rng = np.random.default_rng(SEEDS[0])
    Y = pd.DataFrame((rng.random((len(lines), len(genes))) < SYNTH_RATE).astype(float), index=lines, columns=genes)
    Y06 = pd.DataFrame((rng.random((len(lines), len(genes))) < SYNTH_RATE).astype(float), index=lines, columns=genes)
    return dict(lines=lines, lineage=split.lineage.reindex(lines), genes=genes, gen=w1.generic_table(U.entrez), Y=Y,
                Y06=Y06, line_ce={l: set() for l in lines},
                expr=w1.load_omics()["expr"].reindex(index=lines, columns=genes), panel=panel, split=split)


# ----------------------------------------------------------------------------- main
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--builds", help="m5_build output directory (default results/m5/builds)")
    ap.add_argument("--out", help="output directory (default results/m5/train)")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--synthetic-labels", action="store_true", help="dry run on Bernoulli labels (tmp/ output only)")
    ap.add_argument("--allow-dirty", action="store_true")
    return ap.parse_args(argv)


def _abs(p: str | None, default: Path) -> Path:
    if not p:
        return default
    return Path(p) if p.startswith("/") else paths.PROJECT_ROOT / p


def main(argv=None) -> int:
    from nextgem.run import provenance

    a = parse_args(argv)
    builds, out = _abs(a.builds, BUILDS), _abs(a.out, OUT)
    in_tmp = lambda p: "tmp" in p.resolve().relative_to(paths.PROJECT_ROOT.resolve()).parts  # noqa: E731
    if a.synthetic_labels and not in_tmp(out):
        sys.exit("--synthetic-labels writes only under a tmp/ directory")
    if not a.synthetic_labels and in_tmp(builds):
        sys.exit("real training labels are scored only on production builds (not on a pilot under tmp/)")
    out.mkdir(parents=True, exist_ok=True)
    manifest = MANIFEST if not a.synthetic_labels else out / "MANIFEST.tsv"
    inputs = [builds / f for f in ("builds.tsv", "gates.tsv", "ko.parquet", "build_summary.json")]
    inputs += [w4t.W4 / "stage_a_panel.tsv", paths.baseline_dir() / "human1_ko.parquet",
               paths.baseline_dir() / "human2_ko.parquet", paths.results_dir() / "w1" / "train_per_line.tsv",
               paths.results_dir() / "w3" / "train" / "per_line.tsv"]
    with provenance(manifest, inputs=inputs, seeds=[SEEDS[0]], solver="none (analysis)", workers=1,
                    allow_dirty=a.allow_dirty or a.synthetic_labels, author="M5 engineer agent") as prov:
        prov.note(step="m5_score_train (A12 b)", boot=a.boot, synthetic_labels=a.synthetic_labels,
                  labels="synthetic Bernoulli(0.12); no label file" if a.synthetic_labels else
                  "labels.train() + cv.lolo (CE_fit) via w4_stage_a_train.load; no validation/test")
        files = run(a, builds, out)
        for f in files:
            prov.output(out / f)
    return 0


def run(a, builds: Path, out: Path) -> list[str]:
    t0 = time.time()
    C = synthetic_context() if a.synthetic_labels else w4t.load(a)
    C["emap"] = h1_entrez_map()  # Human1 parent: Human1 genes -> U
    lines, lineage, B, genes = C["lines"], C["lineage"], a.boot, C["genes"]
    ver = []

    # the Human1 map reproduces W1's generic e_H1 on U exactly (score_frame convention == generic table)
    h1 = pd.read_parquet(paths.baseline_dir() / "human1_ko.parquet", columns=["gene", "ratio_r"]).set_index("gene").ratio_r
    S1 = wa.score_frame(pd.DataFrame([h1.to_numpy()], index=["x"], columns=h1.index), genes, C["emap"])
    d = float(np.abs(S1.iloc[0].to_numpy() - C["gen"]["e_H1"].reindex(genes).fillna(0).to_numpy()).max())
    ver.append({"check": "H1 Entrez map: score_frame(baseline H1) equals W1 generic e_H1 on U", "max_abs_diff": d,
                "pass": d == 0.0})

    ko = pd.read_parquet(builds / "ko.parquet", columns=["build_id", "gene", "ratio_r"])
    R = ratio_tables(ko)
    bt = pd.read_csv(builds / "builds.tsv", sep="\t", low_memory=False)
    gates = pd.read_csv(builds / "gates.tsv", sep="\t")
    bsum = json.loads((builds / "build_summary.json").read_text())
    blines = sorted(set(bt.line) - {"generic"})
    if set(blines) != set(lines):  # dry run on the 2-line pilot: every panel line borrows a pilot build (code path only)
        assert a.synthetic_labels, "production scoring needs builds on all 30 panel lines"
        for i, l in enumerate(l for l in lines if l not in set(blines)):
            src = blines[i % len(blines)]
            for fam in FAMILY_ARMS:
                b = arm_build_ids([src], fam)[src]
                if b in R:
                    R[arm_build_ids([l], fam)[l]] = R[b]
    sc = w4t.Scorer(C, R)
    failed = {fam: set(bt[(bt.family == fam) & bt.failed.fillna(False).astype(bool)].line) for fam in ("F2", "F4")}
    for fam in ("F2", "F4"):  # a build without a KO table is failed (timeouts included)
        failed[fam] |= {l for l, b in arm_build_ids(lines, fam).items() if b not in R}
    assert "generic|none" in R, "the generic H1 + B-loo sweep is missing"

    S = {"generic_H1": sc.generic("H1", lines), BASE: sc.frame(arm_build_ids(lines, "generic"))}
    for fam, arm in FAMILY_ARMS.items():
        ok = {l: b for l, b in arm_build_ids(lines, fam).items() if l not in failed[fam]}
        Sf = sc.frame(ok).reindex(lines).fillna(0.0) if ok else pd.DataFrame(0.0, index=lines, columns=genes)
        S[arm] = Sf
    arm_failed = {"generic_H1": set(), BASE: set(), "F2_R1+Bloo": failed["F2"], "F4_R1+Bloo": failed["F4"]}

    def tables(arm, F):
        Sx = metrics.no_information_scores(S[arm], sorted(F), "line")
        return sc.per_line(Sx).join(sc.binary(Sx, sc.h1_genes)), Sx

    # ---- per arm (own failed lines only) and per contrast (union of failed lines, identical in both arms)
    per_line, crow = {}, []
    for arm in S:
        pl, Sx = tables(arm, arm_failed[arm])
        per_line[arm] = pl
        est = pl[pl.estimable]
        mg = metrics.m_delta_ap_gene(metrics.per_gene(Sx, sc.Y.loc[Sx.index]))
        crow.append({"config": arm, "n_lines": len(pl), "n_failed": len(arm_failed[arm]),
                     "n_lines_estimable": int(pl.estimable.sum()),
                     "tap_lin": metrics.lineage_equal_mean(est.tap, lineage.reindex(est.index)),
                     "m_delta_ap_gene": mg,
                     "mcc_L_prob_lin": metrics.lineage_equal_mean(pl.mcc_L_prob.dropna(), lineage.reindex(pl.mcc_L_prob.dropna().index)),
                     "mcc_L_ge06_lin": metrics.lineage_equal_mean(pl.mcc_L_ge06.dropna(), lineage.reindex(pl.mcc_L_ge06.dropna().index))})
    ctab = pd.DataFrame(crow)

    krow, folds_rows = [], []
    for arm, ref, role in CONTRASTS:
        F = arm_failed[arm] | arm_failed[ref]
        pa, _ = tables(arm, F)
        pb, _ = tables(ref, F)
        name = f"{arm} vs {ref}"
        r = single_config_contrast(pa.tap - pb.tap.reindex(pa.index), lineage, name, B)
        ds = w4t.delta_stats(pa.tap_noCE, pb.tap_noCE.reindex(pa.index), lineage, B)
        de = w4t.delta_stats(pa.tap_expr, pb.tap_expr.reindex(pa.index), lineage, B)
        dm = w4t.delta_stats(pa.mcc_L_prob, pb.mcc_L_prob.reindex(pa.index), lineage, B)
        for k, v in r.pop("fold_means").items():
            folds_rows.append({"contrast": name, "heldout_lineage": k, "fold_mean_delta": v})
        krow.append({"contrast": name, "arm": arm, "reference": ref, "role": role, "n_failed_lines": len(F), **r,
                     "delta_noCE": ds["delta_lin"], "lo_noCE": ds["lo"], "hi_noCE": ds["hi"],
                     "delta_expr": de["delta_lin"], "lo_expr": de["lo"], "hi_expr": de["hi"],
                     "delta_mcc_L_prob": dm["delta_lin"], "lo_mcc": dm["lo"], "hi_mcc": dm["hi"]})
    ktab = pd.DataFrame(krow)
    ftab = pd.DataFrame(folds_rows)

    # ---- verification against the frozen W1 / W3 per-line values (real labels only)
    if not a.synthetic_labels:
        w1pl = pd.read_csv(paths.results_dir() / "w1" / "train_per_line.tsv", sep="\t")
        w = w1pl[(w1pl.method == "H1") & (w1pl.geneset == "all")].set_index("ModelID").tap.reindex(lines)
        d = float((per_line["generic_H1"].tap - w).abs().max())
        ver.append({"check": "generic H1 tAP_line equals W1 train_per_line (6-decimal TSV)", "max_abs_diff": d,
                    "pass": d <= 1e-6})
        w3pl = pd.read_csv(paths.results_dir() / "w3" / "train" / "per_line.tsv", sep="\t",
                           usecols=["parent", "mode", "ModelID", "tap"])
        w = w3pl[(w3pl.parent == "human1") & (w3pl["mode"] == W3_MODE)].set_index("ModelID").tap.reindex(lines)
        d = float((per_line[BASE].tap - w).abs().max())
        ver.append({"check": f"generic H1+B-loo tAP_line equals W3 per_line ({W3_MODE}, human1, M0)",
                    "max_abs_diff": d, "pass": d <= 1e-6})
    gc = bsum.get("generic_check", {})
    ver.append({"check": "generic H1+B-loo KO sweep equals W3 sweep row M0|B-loo:cofactor_pool_biomass",
                "max_abs_diff": gc.get("max_abs_diff"), "pass": gc.get("n_mismatch_gt_1e-6") == 0})
    vtab = pd.DataFrame(ver)

    # ---- decision (A12 b)
    gfam = gates[gates.family.isin(list(FAMILY_ARMS))].groupby("family").gate_all.mean()
    decision = {"rule": "A12(b): a family enters Stage B only if its combined Δ̄_lin vs H1 + B-loo is >= +0.005 on the "
                        "30-line panel (nested panel LOLO, Stage-A configuration)",
                "synthetic_labels": bool(a.synthetic_labels), "families": {}}
    for fam, arm in FAMILY_ARMS.items():
        r = ktab[(ktab.arm == arm) & (ktab.reference == BASE)].iloc[0]
        decision["families"][fam] = {
            "config": f"{fam}|R1|rna|real|none on {bsum.get('parent')}", "delta_lin_vs_H1_Bloo": r.delta_lin,
            "ci": [r.lo, r.hi], "folds_positive": int(r.folds_positive), "n_folds": int(r.n_folds),
            "n_failed_builds": len(arm_failed[arm]), "build_gate_all_frac": float(gfam.get(fam, np.nan)),
            "config_gate_pass": bool(gfam.get(fam, 0) >= wa.CONFIG_MIN_PASS_FRAC),
            "enters_stage_b": stage_b_verdict(r.delta_lin)}
    s = ktab[ktab.reference == "generic_H1"].set_index("arm")
    decision["sanity_Bloo_vs_H1"] = {"delta_lin": s.loc[BASE, "delta_lin"], "ci": [s.loc[BASE, "lo"], s.loc[BASE, "hi"]],
                                     "folds_positive": int(s.loc[BASE, "folds_positive"])}
    decision["any_family_enters"] = any(f["enters_stage_b"] for f in decision["families"].values())
    decision["consequence"] = ("Stage B runs for: " + ", ".join(f for f, v in decision["families"].items() if v["enters_stage_b"])
                               if decision["any_family_enters"] else
                               "context lever not selected: the M6 factorial drops its context axis; Stage B and Stage C are skipped")
    decision["verification_all_pass"] = bool(vtab["pass"].all())
    decision["wall_s"] = round(time.time() - t0)

    # ---- write
    long = []
    for arm, pl in per_line.items():
        t = pl.copy()
        t.insert(0, "config", arm)
        t["lineage"] = lineage.reindex(t.index)
        long.append(t.rename_axis("ModelID").reset_index())
    pd.concat(long).to_csv(out / "train_per_line.tsv", sep="\t", index=False, float_format="%.6g")
    ctab.to_csv(out / "train_config_scores.tsv", sep="\t", index=False, float_format="%.6g")
    ktab.to_csv(out / "train_contrasts.tsv", sep="\t", index=False, float_format="%.6g")
    ftab.to_csv(out / "train_folds.tsv", sep="\t", index=False, float_format="%.6g")
    vtab.to_csv(out / "train_verification.tsv", sep="\t", index=False)
    (out / "train_decision.json").write_text(json.dumps(decision, indent=1, default=float))
    cols = ["contrast", "delta_lin", "lo", "hi", "folds_positive", "n_folds", "n_failed_lines", "delta_noCE", "delta_expr",
            "delta_mcc_L_prob"]
    print(ktab[cols].to_string(index=False, float_format=lambda x: f"{x:+.4f}"))
    print(ctab.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(vtab.to_string(index=False))
    print(json.dumps({k: decision[k] for k in ("families", "any_family_enters", "consequence")}, indent=1, default=float))
    return ["train_per_line.tsv", "train_config_scores.tsv", "train_contrasts.tsv", "train_folds.tsv",
            "train_verification.tsv", "train_decision.json"]


if __name__ == "__main__":
    sys.exit(main())
