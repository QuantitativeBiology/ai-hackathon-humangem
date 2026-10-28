"""F1 validation look spec (A15, post hoc; prepared, NOT run): results/m10/val/look_spec_M10_F1_gpr_curation.json.

Only if the curated P1 passed TG on training (results/f1/train/summary.json, real labels): otherwise refuses.
For the reviewed generic harness (scripts/evaluate_val.py, look kind "shortlist", layout line_gene; A6(2), A11f-g,
A13): workstream "M10" = milestone "M10" (milestone-level look under A15, total cap only), reference H1 (= G*),
one candidate = the curated P1 (F1: GPR rule + TYMS on H1 + B-loo / M2(f=1)), descriptive methods = P1 (M6 cell
H1+Bloo/M2) and the two F1 components (TYMS only, GPR only). Label-free: no label is read.

Scores: results/m10/val/f1_val_scores.parquet: the 198 validation expression lines x U, columns H1, H2 (= W1
generic, line-constant) and the arms (e = 1 − ratio_r of the arm's sweep under the line's M2(f=1) signature).
P1's column is asserted equal to the M6 validation scores (results/m6/val/m6_val_scores.parquet).

Usage: .venv/bin/python scripts/f1_val_spec.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f1_common as fc  # noqa: E402
import m6_common as mc  # noqa: E402

LOOK_ID = "M10_F1_gpr_curation"
M10VAL = mc.PROJECT / "results" / "m10" / "val"
SCORES = M10VAL / "f1_val_scores.parquet"
PROV_FILES = [
    "results/f1/edits_human1.tsv", "results/f1/curated_rules_human1.tsv", "results/f1/closure_changes_human1.tsv",
    "results/f1/edits_summary.json", "results/f1/sweeps/ratio_r_delta_human1.parquet", "results/f1/sweeps/sweeps.tsv",
    "results/f1/sweeps/spot_checks.tsv", "results/f1/train/summary.json", "results/f1/train/contrasts.tsv",
    "results/f1/train/cells.tsv", "results/m6/sweeps/ratio_r_human1.parquet", "results/m6/sweeps/catalog.tsv",
    "results/m6/val/m6_val_scores.parquet", "results/w2/controls/line_sweeps.tsv",
    "results/w3/sweeps/ratio_r_human1.parquet", "results/w3/structure/variants_human1.json",
    "src/nextgem/w3.py", "src/nextgem/w5.py", "src/nextgem/w5_evidence.py", "src/nextgem/w2_engine.py",
    "scripts/f1_common.py", "scripts/f1_edits.py", "scripts/f1_sweeps.py", "scripts/f1_score_train.py",
    "scripts/f1_val_spec.py", "scripts/m6_common.py", "scripts/m6_sweeps.py", "scripts/m6_scores.py",
]
VAL_METHODS = ["H1", "H2", fc.CANDIDATE.method, fc.P1_METHOD, fc.ARM["TYMS"].method, fc.ARM["GPR"].method]


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def val_scores() -> pd.DataFrame:
    from f1_score_train import arm_e_matrices, line_matrix
    from nextgem import labels, w1

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    U = list(labels.gene_universe().entrez.astype(int))
    gen = w1.generic_table(labels.gene_universe().entrez)
    V = sorted(l for l in split.index[split.tier == "validation"] if lf.loc[l, "has_expression_matrix"])
    Earm, _, ls, _ = arm_e_matrices(U)
    df = pd.DataFrame({"ModelID": np.repeat(V, len(U)), "entrez": np.tile(np.asarray(U, dtype=np.int64), len(V))})
    for s in ("H1", "H2"):
        df[s] = np.tile(gen[f"e_{s}"].reindex(U).to_numpy(), len(V))
    for a in fc.ARMS:
        M = line_matrix(Earm[a.key], ls, V)
        if M.isna().any().any():
            raise RuntimeError(f"{a.key}: NaN validation scores (failed sweep)")
        df[a.method] = M.to_numpy(dtype=np.float64).ravel()
    m6 = pd.read_parquet(mc.VALDIR / "m6_val_scores.parquet", columns=["ModelID", "entrez", fc.P1_METHOD, "H1", "H2"])
    m6 = m6.set_index(["ModelID", "entrez"]).reindex(pd.MultiIndex.from_frame(df[["ModelID", "entrez"]]))
    for c in (fc.P1_METHOD, "H1", "H2"):
        if not np.array_equal(m6[c].to_numpy(), df[c].to_numpy()):
            raise RuntimeError(f"{c}: validation scores differ from results/m6/val/m6_val_scores.parquet")
    return df[["ModelID", "entrez"] + VAL_METHODS]


def build_spec(scores_sha: str, provenance: dict, train: dict) -> dict:
    from nextgem import metrics
    from nextgem.config import SEEDS

    def fmt(k):
        v = train[k]
        return f"{v['delta_lin']:+.4f} [{v['ci'][0]:+.4f}, {v['ci'][1]:+.4f}], {v['folds_positive']}/12"

    return {
        "look_id": LOOK_ID, "kind": "shortlist", "workstream": "M10", "milestone": "M10",
        "scores": {"path": "results/m10/val/f1_val_scores.parquet", "sha256": scores_sha, "layout": "line_gene"},
        "methods": VAL_METHODS, "candidates": [fc.CANDIDATE.method], "reference": "H1",
        "labels": ["L_prob"], "genesets": ["all", "noCE"], "universes": ["U"], "strata": ["none", "expressed_pairs"],
        "B": metrics.B_OTHER, "seed": SEEDS[0],
        "decision_rule": (
            "A15 post hoc / exploratory follow-up F1. §8 gate: the single candidate (curated P1 = H1 + B-loo + D17 GPR "
            "rule + TYMS on M2(f=1)) is eligible iff its paired validation Δ̄_lin vs G* (H1/M0) has CI lower bound > 0 "
            "(L_prob, all genes, U, no stratum); m = 1, so σ·√(2 ln m) = 0 and the one-SE rule is trivial. The F1 "
            "question (does the curation improve on P1?) is read from the descriptive paired contrast curated P1 − P1 "
            f"({fc.P1_METHOD}) in the same look; it never replaces T1/T2 or the headline, and no test evaluation follows."),
        "provenance": provenance,
        "notes": (
            "Post hoc / exploratory (DESIGN §17 A15; RESEARCH_TREE D17). Milestone-level look M10 (A13), total cap only. "
            "D17 rule: AND→OR merges only for Ensembl paralog pairs with min %id >= 30 and equal Human-GEM genes.tsv "
            "compartment sets and no shared CORUM/Complex Portal complex, never in exchange, transport, pseudo/lumped "
            "(subsystem Artificial/Pool reactions or composition-weighted lumped stoichiometry, e.g. MAR07161) or "
            "AND>=6-gene reactions; plus the H1 TYMS (MAR04644) irreversibility fix. Edits: results/f1/edits_human1.tsv. "
            f"Training (686 lines, 12 lineages; results/f1/train/): curated P1 − P1 {fmt('F1_vs_P1')}; TYMS only − P1 "
            f"{fmt('TYMS_vs_P1')}; GPR only − P1 {fmt('GPR_vs_P1')}. Flags inherited from P1: W3-D1, W3-MEV, W2-SHUF, "
            "W2-MOPEN. The F1 agent did not open any T1/T2 output or results/m6/test/ (A15 blindness)."),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    train = json.loads((fc.TRDIR / "summary.json").read_text())
    if train.get("synthetic_labels") or not train.get("F1_passes_TG"):
        raise SystemExit("F1 did not pass TG on training (or the summary is synthetic): no validation spec is prepared")
    with fc.provenance(inputs=[mc.PROJECT / f for f in PROV_FILES], solver="none", workers=1,
                       allow_dirty=a.allow_dirty) as prov:
        df = val_scores()
        sc = fc.atomic_write_parquet(df, SCORES)
        pins = {f: _sha(mc.PROJECT / f) for f in PROV_FILES}
        spec = build_spec(_sha(sc), pins, train)
        out = fc.atomic_write_text(json.dumps(spec, indent=1) + "\n", M10VAL / f"look_spec_{LOOK_ID}.json")
        prov.output(sc, out)
        prov.note(look_id=LOOK_ID, scores_sha256=spec["scores"]["sha256"], n_val_lines=int(df.ModelID.nunique()))
        print(json.dumps({k: spec[k] for k in ("look_id", "workstream", "methods", "candidates", "reference")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
