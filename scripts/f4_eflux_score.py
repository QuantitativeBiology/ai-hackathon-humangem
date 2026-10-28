"""F4b pre-check scoring (post hoc / exploratory, DESIGN §17 A15; RESEARCH_TREE D19): E-Flux on P1 vs P1 on the 30
Stage-A panel lines (training tier). Committed before its first label-reading run; mirrors scripts/m5_score_train.py.

Machinery: Stage A's / M5's (scripts/w4_stage_a_train.py ``load``, ``Scorer``, ``delta_stats``; scripts/m5_score_train.py
``h1_entrez_map``, ``single_config_contrast``, ``synthetic_context``; nextgem.w4_analysis ``nested_lolo_select``,
``tg_panel``). Scores: e = 1 − ratio_r over U, Human1 genes mapped with the W1 generic-table Entrez convention.

Arms (L_prob, all genes, U; failed builds take the §11.4 no-information score in both arms of a contrast):
  generic_H1  W1 generic Human1 (G*)
  P1          the line's pinned P1 sweep (M6 cell H1+Bloo/M2; the line's M2(f=1) signature)
  F4b         E-Flux F4|R1|rna|real|none on the P1 parent with the line's M2(f=1) medium (scripts/f4_eflux_build.py)
Contrasts: F4b vs P1 [the D19 pre-check decision]; P1 vs generic_H1 [sanity: the P1 gain on the panel];
F4b vs generic_H1 [context]. Nested panel LOLO as in Stage A / M5: one fixed configuration, so the inner choice is that
configuration in every fold and the outer per-line Δ equals the plain Δ (asserted). Δ̄_lin = lineage-equal mean, paired
line bootstrap stratified by lineage (B = 2,000, SEEDS[0]); folds = the 10 panel lineages.

Decision (D19, pre-registered): F4b proceeds to all training lines only if its Δ̄_lin vs P1 is >= +0.005
(f4_common.panel_proceeds). Otherwise F4b stops here. Gate statistics (E-Flux GR units artefact) are reported.

Labels: only ``labels.train()`` (panel lines) and ``cv.lolo`` / ``labels.fold()`` (CE_fit per held-out panel lineage),
through ``w4_stage_a_train.load``. ``--synthetic-labels``: Bernoulli(0.12) labels, no label file, tmp/ output only.

Outputs (default results/f4/eflux_panel/train/): per_line.tsv, config_scores.tsv, contrasts.tsv, folds.tsv,
verification.tsv, decision.json.

Usage:
  .venv/bin/python scripts/f4_eflux_score.py --builds results/f4/tmp/eflux_pilot --out results/f4/tmp/eflux_score_dry --synthetic-labels
  .venv/bin/python scripts/f4_eflux_score.py                      # after scripts/f4_eflux_build.py (clean tree)
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
import f4_common as fc  # noqa: E402

fc.blas_single_thread()

import m5_score_train as m5s  # noqa: E402
import w4_stage_a_train as w4t  # noqa: E402

from nextgem import metrics, paths  # noqa: E402
from nextgem import w4_analysis as wa  # noqa: E402
from nextgem.config import SEEDS  # noqa: E402

BUILDS = fc.EFDIR / "builds"
OUT = fc.EFDIR / "train"
CONTRASTS = [("F4b", "P1", "decision"), ("P1", "generic_H1", "sanity (P1 gain)"), ("F4b", "generic_H1", "context")]


def f4b_build_ids(lines) -> dict:
    return {l: "|".join(["F4b", l, *fc.F4B_CONFIG[:4], "none"]) for l in lines}


def main(argv=None) -> int:
    from nextgem.run import provenance

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--builds")
    ap.add_argument("--out")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--synthetic-labels", action="store_true")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    builds = m5s._abs(a.builds, BUILDS)
    out = m5s._abs(a.out, OUT)
    if a.synthetic_labels and not fc.in_tmp(out):
        sys.exit("--synthetic-labels writes only under a tmp/ directory")
    if not a.synthetic_labels and fc.in_tmp(builds):
        sys.exit("real training labels are scored only on production builds (not on a pilot under tmp/)")
    out.mkdir(parents=True, exist_ok=True)
    manifest = fc.MANIFEST if not a.synthetic_labels else out / "MANIFEST.tsv"
    inputs = [builds / f for f in ("builds.tsv", "gates.tsv", "ko.parquet", "build_summary.json", "generic_check.tsv")]
    inputs += [w4t.W4 / "stage_a_panel.tsv", paths.baseline_dir() / "human1_ko.parquet",
               paths.results_dir() / "m6" / "sweeps" / "ratio_r_human1.parquet",
               paths.results_dir() / "w3" / "sweeps" / "ratio_r_human1.parquet",
               paths.results_dir() / "m6" / "scores" / "line_assignment.tsv",
               paths.results_dir() / "m6" / "train" / "per_line.tsv", paths.results_dir() / "w1" / "train_per_line.tsv"]
    with provenance(manifest, inputs=inputs, seeds=[SEEDS[0]], solver="none (analysis)", workers=1,
                    allow_dirty=a.allow_dirty or a.synthetic_labels, author=fc.AUTHOR) as prov:
        prov.note(step="f4_eflux_score (F4b panel pre-check, D19)", boot=a.boot, synthetic_labels=a.synthetic_labels,
                  labels="synthetic Bernoulli(0.12); no label file" if a.synthetic_labels else
                  "labels.train() + cv.lolo (CE_fit) via w4_stage_a_train.load; no validation/test")
        files = run(a, builds, out)
        for f in files:
            prov.output(out / f)
    return 0


def run(a, builds: Path, out: Path) -> list[str]:
    from f4_sweeps import p1_rows

    t0 = time.time()
    C = m5s.synthetic_context() if a.synthetic_labels else w4t.load(a)
    C["emap"] = m5s.h1_entrez_map()
    lines, lineage, B = C["lines"], C["lineage"], a.boot
    sig = fc.line_signatures().reindex(lines)
    ver = []
    ko = pd.read_parquet(builds / "ko.parquet", columns=["build_id", "gene", "ratio_r"])
    R = m5s.ratio_tables(ko)
    P1 = p1_rows(sorted(set(sig)))
    for s in P1.index:
        R[f"P1|{s}"] = P1.loc[s]
    bt = pd.read_csv(builds / "builds.tsv", sep="\t", low_memory=False)
    gates = pd.read_csv(builds / "gates.tsv", sep="\t")
    bsum = json.loads((builds / "build_summary.json").read_text())
    gchk = pd.read_csv(builds / "generic_check.tsv", sep="\t")
    f4ids = f4b_build_ids(lines)
    have = {l for l, b in f4ids.items() if b in R}
    if have != set(lines):  # dry run on the pilot: missing panel lines borrow a pilot build (code path only)
        assert a.synthetic_labels, "production scoring needs builds on all 30 panel lines"
        src = sorted(have)
        for i, l in enumerate(l for l in lines if l not in have):
            R[f4ids[l]] = R[f4ids[src[i % len(src)]]]
    sc = w4t.Scorer(C, R)
    failed = set(bt[(bt.family == "F4") & bt.failed.fillna(False).astype(bool)].line)
    failed |= {l for l, b in f4ids.items() if b not in R}
    S = {"generic_H1": sc.generic("H1", lines), "P1": sc.frame({l: f"P1|{sig[l]}" for l in lines}),
         "F4b": sc.frame({l: f4ids[l] for l in lines if l not in failed}).reindex(lines).fillna(0.0)}
    arm_failed = {"generic_H1": set(), "P1": set(), "F4b": failed}

    def tables(arm, F):
        Sx = metrics.no_information_scores(S[arm], sorted(F), "line")
        return sc.per_line(Sx).join(sc.binary(Sx, sc.h1_genes)), Sx

    per_line, crow = {}, []
    for arm in S:
        pl, Sx = tables(arm, arm_failed[arm])
        per_line[arm] = pl
        est = pl[pl.estimable]
        mg = metrics.m_delta_ap_gene(metrics.per_gene(Sx, sc.Y.loc[Sx.index]))
        crow.append({"config": arm, "n_lines": len(pl), "n_failed": len(arm_failed[arm]),
                     "n_lines_estimable": int(pl.estimable.sum()),
                     "tap_lin": metrics.lineage_equal_mean(est.tap, lineage.reindex(est.index)), "m_delta_ap_gene": mg,
                     "mcc_L_prob_lin": metrics.lineage_equal_mean(pl.mcc_L_prob.dropna(),
                                                                  lineage.reindex(pl.mcc_L_prob.dropna().index))})
    ctab = pd.DataFrame(crow)
    krow, frows = [], []
    for arm, ref, role in CONTRASTS:
        F = arm_failed[arm] | arm_failed[ref]
        pa, _ = tables(arm, F)
        pb, _ = tables(ref, F)
        name = f"{arm} vs {ref}"
        r = m5s.single_config_contrast(pa.tap - pb.tap.reindex(pa.index), lineage, name, B)
        ds = w4t.delta_stats(pa.tap_noCE, pb.tap_noCE.reindex(pa.index), lineage, B)
        de = w4t.delta_stats(pa.tap_expr, pb.tap_expr.reindex(pa.index), lineage, B)
        dm = w4t.delta_stats(pa.mcc_L_prob, pb.mcc_L_prob.reindex(pa.index), lineage, B)
        for k, v in r.pop("fold_means").items():
            frows.append({"contrast": name, "heldout_lineage": k, "fold_mean_delta": v})
        krow.append({"contrast": name, "arm": arm, "reference": ref, "role": role, "n_failed_lines": len(F), **r,
                     "delta_noCE": ds["delta_lin"], "lo_noCE": ds["lo"], "hi_noCE": ds["hi"],
                     "delta_expr": de["delta_lin"], "lo_expr": de["lo"], "hi_expr": de["hi"],
                     "delta_mcc_L_prob": dm["delta_lin"], "lo_mcc": dm["lo"], "hi_mcc": dm["hi"]})
    ktab, ftab = pd.DataFrame(krow), pd.DataFrame(frows)
    if not a.synthetic_labels:
        m6pl = pd.read_csv(paths.results_dir() / "m6" / "train" / "per_line.tsv", sep="\t")
        w = m6pl[(m6pl.cell == "H1+Bloo/M2") & (m6pl.label_cell == "all")].set_index("ModelID").tap.reindex(lines)
        d = float((per_line["P1"].tap - w).abs().max())
        ver.append({"check": "P1 tAP_line on the panel equals M6 train per_line (H1+Bloo/M2, all)", "max_abs_diff": d,
                    "pass": d <= 1e-6})
        w1pl = pd.read_csv(paths.results_dir() / "w1" / "train_per_line.tsv", sep="\t")
        w = w1pl[(w1pl.method == "H1") & (w1pl.geneset == "all")].set_index("ModelID").tap.reindex(lines)
        d = float((per_line["generic_H1"].tap - w).abs().max())
        ver.append({"check": "generic H1 tAP_line equals W1 train_per_line", "max_abs_diff": d, "pass": d <= 1e-6})
    ver.append({"check": "W4-engine generic P1 builds equal the pinned P1 sweeps (<= 1e-6)",
                "max_abs_diff": float(gchk.max_abs_diff.max()), "pass": bool(gchk["pass"].all())})
    vtab = pd.DataFrame(ver)
    r = ktab[(ktab.arm == "F4b") & (ktab.reference == "P1")].iloc[0]
    gf = gates[gates.family == "F4"]
    decision = {"rule": "D19: F4b proceeds to all training lines only if its Δ̄_lin vs P1 is >= +0.005 on the 30-line "
                        "panel (nested panel LOLO, Stage A's F4 configuration)",
                "post_hoc": "A15 post hoc / exploratory", "synthetic_labels": bool(a.synthetic_labels),
                "config": "|".join(fc.F4B_CONFIG), "delta_lin_vs_P1": r.delta_lin, "ci": [r.lo, r.hi],
                "folds_positive": int(r.folds_positive), "n_folds": int(r.n_folds), "n_failed_builds": len(failed),
                "delta_noCE": r.delta_noCE, "delta_expr": r.delta_expr, "delta_mcc_L_prob": r.delta_mcc_L_prob,
                "build_gate_all_frac": float(gf.gate_all.mean()) if len(gf) else None,
                "build_gate_viability_frac": float(gf.gate_viability.mean()) if len(gf) else None,
                "config_gate_pass": bsum.get("config_gate_pass"), "n_gr_task_fail": bsum.get("n_gr_task_fail"),
                "proceeds_to_all_training_lines": fc.panel_proceeds(r.delta_lin),
                "verification_all_pass": bool(vtab["pass"].all()), "wall_s": round(time.time() - t0)}
    s = ktab[ktab.reference == "generic_H1"].set_index("arm")
    decision["sanity_P1_vs_H1"] = {"delta_lin": s.loc["P1", "delta_lin"], "ci": [s.loc["P1", "lo"], s.loc["P1", "hi"]],
                                   "folds_positive": int(s.loc["P1", "folds_positive"])}
    long = []
    for arm, pl in per_line.items():
        t = pl.copy()
        t.insert(0, "config", arm)
        t["lineage"] = lineage.reindex(t.index)
        long.append(t.rename_axis("ModelID").reset_index())
    pd.concat(long).to_csv(out / "per_line.tsv", sep="\t", index=False, float_format="%.6g")
    ctab.to_csv(out / "config_scores.tsv", sep="\t", index=False, float_format="%.6g")
    ktab.to_csv(out / "contrasts.tsv", sep="\t", index=False, float_format="%.6g")
    ftab.to_csv(out / "folds.tsv", sep="\t", index=False, float_format="%.6g")
    vtab.to_csv(out / "verification.tsv", sep="\t", index=False)
    (out / "decision.json").write_text(json.dumps(decision, indent=1, default=float))
    cols = ["contrast", "delta_lin", "lo", "hi", "folds_positive", "n_folds", "n_failed_lines", "delta_noCE",
            "delta_expr", "delta_mcc_L_prob"]
    print(ktab[cols].to_string(index=False, float_format=lambda x: f"{x:+.4f}"))
    print(ctab.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(vtab.to_string(index=False))
    print(json.dumps(decision, indent=1, default=float))
    if not vtab["pass"].all():
        raise RuntimeError("verification failed")
    return ["per_line.tsv", "config_scores.tsv", "contrasts.tsv", "folds.tsv", "verification.tsv", "decision.json"]


if __name__ == "__main__":
    sys.exit(main())
