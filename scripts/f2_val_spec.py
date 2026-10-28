"""F2 validation look spec (A15 post hoc; RESEARCH_TREE D18; prepared, NOT run): results/m10/val/look_spec_M10_F2.json.

For the reviewed generic harness (scripts/evaluate_val.py, kind "shortlist", layout line_gene; A6(2), A11g, A13):
workstream "M10" = milestone "M10" (milestone-level look, total cap only), reference H1 (= G*), one candidate F2
(§8 gate vs G*), P1 descriptive. Scores: results/m10/val/f2_val_scores.parquet written by scripts/f2_scores.py
(its SHA-256 is re-checked against results/f2/scores_checks.json). Label-free: no label is read.

--pin-train adds the non-author training recompute (results/f2/train/summary.json, contrasts.tsv) to the
provenance pins; it refuses unless that summary is from real labels and every recomputed contrast agrees. Without
it the spec pins only label-free inputs and code.

Usage: .venv/bin/python scripts/f2_val_spec.py [--pin-train] [--allow-dirty]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f2_common as fc  # noqa: E402
from f2_scores import VAL_COLUMNS, VAL_SCORES  # noqa: E402

LOOK_ID = "M10_F2"
SPEC = fc.M10VAL / f"look_spec_{LOOK_ID}.json"
PROV_FILES = [
    "results/w7/sweeps/ratio_r_human1.parquet", "results/w7/sweeps/sweeps.tsv", "results/w7/sweeps/sweeps.json",
    "results/w7/sweeps/mev.tsv", "results/f2/mev_check.tsv", "results/f2/scores_checks.json",
    "results/m6/scores/sweep_e.parquet", "results/m6/scores/line_assignment.tsv", "results/m6/scores/w3_mev.tsv",
    "results/m6/scores/SHA256SUMS.tsv", "results/m6/val/m6_val_scores.parquet", "results/w2/controls/line_sweeps.tsv",
    "results/w3/structure/variants_human1.json", "src/nextgem/w3.py", "src/nextgem/w7.py", "src/nextgem/w2_engine.py",
    "scripts/w7_artefact_sweeps.py", "scripts/m6_common.py", "scripts/m6_sweeps.py", "scripts/f2_common.py",
    "scripts/f2_scores.py", "scripts/f2_score_train.py", "scripts/f2_val_spec.py",
]
TRAIN_FILES = ["results/f2/train/summary.json", "results/f2/train/contrasts.tsv"]

DECISION_RULE = (
    "A15 post hoc / exploratory follow-up F2 (RESEARCH_TREE D18). §8 gate: the single candidate F2 (H1 + artefact-only "
    "B-loo W7-03 on M2(f=1)) is eligible iff its paired validation Δ̄_lin vs G* (H1/M0, the W1 look's G-star) has CI "
    "lower bound > 0 (L_prob, all genes, U, no stratum); m = 1, so σ·√(2 ln m) = 0 and the one-SE rule is trivial. "
    f"P1 ({fc.P1_METHOD}) is descriptive and is never selected here. The F2 question is the validated cost of the "
    "safer variant, F2 − P1, read as Δ̄_lin(F2 − H1) − Δ̄_lin(P1 − H1) on the identical primary-cell lines (exact "
    "point estimate; no CI from the harness). No claim that F2 beats P1; it never replaces T1/T2 or the headline, and "
    "no test evaluation follows.")

NOTES = (
    "Post hoc / exploratory (DESIGN §17 A15; RESEARCH_TREE D18, fixed before any F2 validation scoring). "
    "Milestone-level look M10 (A13), total cap only. F2 = Human1 + artefact-only biomass leave-out (W7-03; "
    "autoresearch/ideas.tsv): retinoid, vitamin-D and tocopherol catabolites and carnitine removed from the cofactor "
    "pool MAR10065 (results/w3/structure/variants_human1.json groups B-loo:retinoid_catabolites, "
    "B-loo:vitaminD_catabolites, B-loo:tocopherol_catabolites, B-loo:carnitine); ubiquinol, heme, lipoate, cobamide "
    "and the real cofactors stay. Medium: documented M2(f=1) per line, exactly P1's (results/m6/scores/"
    "line_assignment.tsv). Scores: e = 1 − ratio_r of the W7-03 KO sweep of the line's M2(f=1) signature "
    "(results/w7/sweeps), 198 validation expression lines x U; H1, H2 and P1 equal the M6 look's columns bitwise. "
    "Purpose: the validated cost F2 − P1 of the biologically safer variant. F2 fixes W3-MEV label-free: the "
    "mevalonate core (HMGCR, MVK, PMVK, MVD) is lethal on 19/19 serum M2(f=1) signatures and on the serum-free one "
    "(1,361/1,361 lines), against 0/19 for P1 (results/f2/mev_check.tsv). Training (W7 r4, post hoc; 686 lines, "
    "12 lineages): F2 − G* +0.0304 [0.0294, 0.0313], 12/12; F2 − P1 −0.0230 [−0.0236, −0.0224], 0/12; the non-author "
    "recompute is results/f2/train/summary.json (D18 step 1). Flags carried by F2: W3-D1 (biomass leave-out, no "
    "pre-registered negative control), W2-SHUF and W2-MOPEN (medium); W3-MEV is lifted. The F2 agent did not open "
    "T1/T2 outputs, results/m6/test/ or results/f3/ (A15 blindness).")


def build_spec(scores_sha: str, provenance: dict) -> dict:
    from nextgem.config import SEEDS

    return {
        "look_id": LOOK_ID, "kind": "shortlist", "workstream": "M10", "milestone": "M10",
        "scores": {"path": f"results/m10/val/{VAL_SCORES.name}", "sha256": scores_sha, "layout": "line_gene"},
        "methods": list(VAL_COLUMNS), "candidates": [fc.F2_METHOD], "reference": "H1",
        "labels": ["L_prob"], "genesets": ["all", "noCE"], "universes": ["U"], "strata": ["none", "expressed_pairs"],
        "B": 2000, "seed": SEEDS[0], "decision_rule": DECISION_RULE, "provenance": provenance, "notes": NOTES,
    }


def check_train() -> None:
    s = json.loads((fc.TRDIR / "summary.json").read_text())
    if s.get("synthetic_labels") or not s.get("all_agree"):
        raise SystemExit("training recompute missing, synthetic or not in agreement: refusing --pin-train")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pin-train", action="store_true")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    files = PROV_FILES + (TRAIN_FILES if a.pin_train else [])
    if a.pin_train:
        check_train()
    with fc.provenance(inputs=[fc.PROJECT / f for f in files] + [VAL_SCORES], allow_dirty=a.allow_dirty) as prov:
        checks = json.loads((fc.F2DIR / "scores_checks.json").read_text())
        sha = fc.sha256(VAL_SCORES)
        if sha != checks["val_scores_sha256"]:
            raise SystemExit("f2_val_scores.parquet differs from the one f2_scores.py checked")
        pins = {f: fc.sha256(fc.PROJECT / f) for f in files}
        spec = build_spec(sha, pins)
        SPEC.parent.mkdir(parents=True, exist_ok=True)
        SPEC.write_text(json.dumps(spec, indent=1) + "\n")
        prov.output(SPEC)
        prov.note(look_id=LOOK_ID, scores_sha256=sha, pin_train=a.pin_train)
        print(json.dumps({k: spec[k] for k in ("look_id", "workstream", "methods", "candidates", "reference")},
                         indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
