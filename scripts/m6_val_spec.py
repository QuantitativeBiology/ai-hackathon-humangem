"""M6 factorial validation look spec (prepared, NOT run): results/m6/val/look_spec_M6_factorial.json.

For the reviewed generic harness (scripts/evaluate_val.py, look kind "shortlist", layout line_gene; DESIGN §17 A6(2),
A11f-g, A12c). Scores: results/m6/val/m6_val_scores.parquet (written by scripts/m6_scores.py; the 198 validation
expression lines x U, one column per method). Label-free: no label is read.

Candidates in DESIGN §9 order (simplest first: medium < tasks < combination; D10 as corrected): (1) H1 / M2(f=1),
(2) H1 + B-loo / M0, (3) H1 + B-loo / M2(f=1); reference
G* = H1 / M0; the H2 cells (generic H2 is a harness method anyway) are descriptive methods. Provenance pins the scores'
inputs, the M6 code and the training outputs.

Workstream "M6" in results/m6/val is a milestone-level look (DESIGN §17 A13).

Usage: .venv/bin/python scripts/m6_val_spec.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as mc  # noqa: E402

LOOK_ID = "M6_factorial"
PROV_FILES = [
    "results/m6/scores/SHA256SUMS.tsv", "results/m6/scores/sweep_e.parquet", "results/m6/scores/line_assignment.tsv",
    "results/m6/scores/w3_mev.tsv", "results/m6/sweeps/catalog.tsv", "results/m6/sweeps/sweeps.tsv",
    "results/m6/sweeps/spot_checks.tsv", "results/m6/sweeps/ratio_r_human1.parquet",
    "results/m6/sweeps/ratio_r_human2.parquet", "results/m6/train/summary.json", "results/m6/train/contrasts.tsv",
    "results/m6/train/cells.tsv", "results/m6/train/interaction.tsv",
    "results/w2/controls/line_sweeps.tsv", "results/w2/ko/ratio_r_human1.parquet", "results/w2/ko/ratio_r_human2.parquet",
    "results/w2/train/shortlist.json", "results/w3/sweeps/ratio_r_human1.parquet", "results/w3/sweeps/ratio_r_human2.parquet",
    "results/w3/structure/variants_human1.json", "results/w3/structure/variants_human2.json",
    "results/w3/train/shortlist.json", "results/m5/train/train_decision.json",
    "src/nextgem/w3.py", "src/nextgem/w2_engine.py", "src/nextgem/w2_score.py", "src/nextgem/media.py",
    "scripts/m6_common.py", "scripts/m6_sweeps.py", "scripts/m6_scores.py", "scripts/m6_score_train.py",
    "scripts/m6_val_spec.py",
]


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def build_spec(scores_sha: str, provenance: dict, train_summary: dict) -> dict:
    from nextgem import metrics
    from nextgem.config import SEEDS

    cands = [c.method for c in mc.CANDIDATES]
    methods = ["H1", "H2", *cands, *[c.method for c in mc.CELLS if c.parent == "human2" and c.key != "H2"]]
    tg = train_summary.get("tg", {})
    tg_txt = "; ".join(f"{v['label']} {k}: {v['delta_lin']:+.4f} [{v['ci'][0]:+.4f}, {v['ci'][1]:+.4f}], "
                       f"{v['folds_positive']}/12, TG {'pass' if v['tg_pass'] else 'fail'}" for k, v in tg.items())
    return {
        "look_id": LOOK_ID, "kind": "shortlist", "workstream": "M6", "milestone": "M6",
        "scores": {"path": "results/m6/val/m6_val_scores.parquet", "sha256": scores_sha, "layout": "line_gene"},
        "methods": methods, "candidates": cands, "reference": "H1",
        "labels": ["L_prob"], "genesets": ["all", "noCE"], "universes": ["U"], "strata": ["none", "expressed_pairs"],
        "B": metrics.B_OTHER, "seed": SEEDS[0],
        "decision_rule": (
            "§8 gate per candidate: eligible iff its paired validation Δ̄_lin vs G* (H1/M0, the W1 look's G-star) has "
            "CI lower bound > 0 (L_prob, all genes, U, no stratum). Among eligible candidates the paired-difference "
            "one-SE rule (§9) picks the first candidate in spec order (D10, simplest first: (1) H1+B-loo/M0, "
            "(2) H1/M2(f=1), (3) H1+B-loo/M2(f=1)) whose Δ̄_lin(best − candidate) <= its bootstrap SE; that cell is P1 "
            "(§10.1 steps 3-4, A12c). If none is eligible, P1 is the best point estimate labelled 'gate not passed'. "
            "σ·√(2 ln m) with m = 3 is reported. H2 cells are descriptive (never P1)."),
        "provenance": provenance,
        "notes": (
            "M6 integration factorial (A12c; cells fixed in RESEARCH_TREE D10 before any factorial scoring). Context "
            "lever not selected (M5 pre-check, A12b), so no context axis. Scores are line-specific (layout line_gene) "
            "over the 198 validation expression lines: e = 1 − ratio_r of the KO sweep of (parent, biomass variant, the "
            "line's medium); M2(f=1) = documented base medium + supplements + serum model f = 1, gate-free, as in "
            "results/w2/controls/line_sweeps.tsv (imputed / fallback / unmappable lines: W2's integrated-pipeline "
            "medium, classic RPMI + 10% FBS; onboarded lines: onboarded medium). B-loo = W3 B-loo:cofactor_pool_biomass "
            "(results/w3/structure/variants_<parent>.json). H1 and H2 columns equal the W1 generic scores. "
            "Flags carried: W3-D1 on cells (2) and (3) (no pre-registered negative control applies to biomass "
            "leave-outs; no claim may rest on 'beats controls'); W2-SHUF and W2-MOPEN on cells (1) and (3) (no "
            "line-specific-medium claim; the gain is 'realistic culture medium vs Ham's'); W3-MEV on cell (3) (under "
            "serum media the cofactor-pool B-loo makes HMGCR/MVK/PMVK/MVD non-lethal; results/m6/scores/w3_mev.tsv). "
            f"Training (686 lines, 12 lineages; results/m6/train/): {tg_txt}. "
            "This single look is also W2's validation answer (D10). Milestone-level look (A13). Candidate order "
            "follows DESIGN §9 (medium < tasks < combination), correcting D10 before any validation scoring "
            "(notes/reviews/M6_harness_and_spec_review.md F-B1)."),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    scores = mc.VALDIR / "m6_val_scores.parquet"
    with mc.provenance(inputs=[scores] + [mc.PROJECT / f for f in PROV_FILES], solver="none", workers=1,
                       allow_dirty=a.allow_dirty) as prov:
        summ = json.loads((mc.TRDIR / "summary.json").read_text())
        prov_pins = {f: _sha(mc.PROJECT / f) for f in PROV_FILES}
        spec = build_spec(_sha(scores), prov_pins, summ)
        out = mc.atomic_write_text(json.dumps(spec, indent=1) + "\n", mc.VALDIR / f"look_spec_{LOOK_ID}.json")
        prov.output(out)
        prov.note(look_id=LOOK_ID, scores_sha256=spec["scores"]["sha256"], n_methods=len(spec["methods"]))
        print(json.dumps({k: spec[k] for k in ("look_id", "workstream", "methods", "candidates")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
