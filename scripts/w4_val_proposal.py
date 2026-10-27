"""Proposed validation-look config entry for W4 (amendment A6: a new look reuses the reviewed evaluate_val.py
harness; it needs a non-author review of its config entry, a label-free dry run and the pre-flight check).

THIS SCRIPT READS NO LABEL and runs no look. It
  1. assembles the proposed config for look ``W4_M5_stageA_carry`` from the Stage-A selection
     (results/w4/stage_a_analysis/train_selection.json): the configuration(s) carried per §9 W4;
  2. computes the config SHA with evaluate_val.config_sha and the SHA-256 of every input that exists;
  3. lists the label-free prerequisites that do not exist yet (Stage-B builds on validation expression lines);
  4. runs evaluate_val.preflight() on the proposed required files (label-free) and reports what blocks it;
  5. writes results/w4/stage_a_analysis/val_look_proposal.json and val_look_dryrun.json.

The W4 look cannot run through the current evaluate_val.py: it has only the W1 look type and no scoring path for
per-line context models. Adding one changes a frozen harness file, which A6 (2) says needs a full independent code
review. Under the frozen §9 W4 rule, validation builds are made only for TG-passing configurations, and Stage A has
none, so this look is a proposal for the PI, not a request to run it.

Usage: .venv/bin/python scripts/w4_val_proposal.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))


def main() -> int:
    import evaluate_val as ev  # label-free helpers only: config_sha, preflight, _sha

    from nextgem import metrics, paths
    from nextgem.config import SEEDS

    out_dir = paths.results_dir() / "w4" / "stage_a_analysis"
    sel = json.loads((out_dir / "train_selection.json").read_text())
    carried = sel["w4_shortlist"] or ([sel["best_if_none_pass"]] if sel["best_if_none_pass"] else [])
    split = pd.read_csv(paths.splits_dir() / "lineage_split.tsv", sep="\t")
    val = split[(split.tier == "validation") & split.has_expression.astype(bool)]
    inputs_rel = {
        "human1_ko": "results/baseline/human1_ko.parquet", "human2_ko": "results/baseline/human2_ko.parquet",
        "stage_a_selection": "results/w4/stage_a_analysis/train_selection.json",
        "threshold_stats": "results/w4/threshold_stats.json", "threshold_gene_means": "results/w4/threshold_gene_means.parquet",
        "stage_b_val_builds": "results/w4/stage_b/ko_validation.parquet",  # does not exist: Stage B not run
        "stage_b_val_manifest": "results/w4/stage_b/MANIFEST.tsv",
        "expr": "results/w1/tmp/expr.parquet", "sesoi_c1b": "results/w1/sesoi_c1b.json",
    }
    root = paths.state_root()
    code = list(ev.CODE_FILES) + ["src/nextgem/context.py", "src/nextgem/context_inputs.py", "src/nextgem/w4_analysis.py"]
    cfg = {
        "look_id": "W4_M5_stageA_carry", "workstream": "W4", "milestone": "M5",
        "configurations": [{"family": c["family"], "config": c["config"], "tag": c.get("tag") or c.get("flag", "")}
                           for c in carried],
        "status": ("no configuration passed the Stage-A panel TG; the best gate-passing configuration is carried to M6 "
                   "flagged (§9 W4). Under the frozen rule no validation builds are made for it until it passes TG."),
        "reference": "G* = Human1 (W1_M2_baseline); also the Stage-A parent, generic Human2 on M0",
        "primary": "Δ̄_lin(config - G*), L_prob, all genes, U, validation expression lines; paired line bootstrap "
                   f"B = {metrics.B_OTHER}, seed {SEEDS[0]}",
        "secondary": ["minus CE_fit (from all training lines, §5.3)", "expressed pairs (§5.4)",
                      "mΔAP_gene vs B3", "§6.8 fixed-cutoff and PGR binary endpoints", "Δ vs B8+/B8-/B8' and B2"],
        "outputs": "scalars only (evaluate_val.assert_scalar_only)", "B": metrics.B_OTHER, "seed": SEEDS[0],
        "n_validation_expression_lines": int(len(val)), "validation_lineages": sorted(val.lineage.unique()),
        "inputs": {k: (ev._sha(root / v) if (root / v).exists() else None) for k, v in inputs_rel.items()},
        "code": {Path(c).name: ev._sha(paths.PROJECT_ROOT / c) for c in code},
        "_inputs_rel": inputs_rel,
    }
    cfg["config_sha"] = ev.config_sha(cfg)
    required = code + list(inputs_rel.values())
    problems = ev.preflight(required)
    missing = [k for k, v in cfg["inputs"].items() if v is None]
    dry = {"look_id": cfg["look_id"], "config_sha": cfg["config_sha"], "labels_read": False,
           "missing_inputs": missing, "preflight_problems": problems,
           "harness_change_needed": "evaluate_val.py has only the W1 look; a W4 look type (per-line context-model scores "
                                    "from Stage-B validation builds) must be added and fully reviewed (A6 (2))",
           "ready": not problems and not missing}
    (out_dir / "val_look_proposal.json").write_text(json.dumps({k: v for k, v in cfg.items() if not k.startswith("_")}
                                                               | {"inputs_rel": inputs_rel}, indent=1))
    (out_dir / "val_look_dryrun.json").write_text(json.dumps(dry, indent=1))
    print(json.dumps({k: dry[k] for k in ("look_id", "config_sha", "labels_read", "missing_inputs", "ready")}, indent=1))
    print(f"{len(problems)} pre-flight problems, e.g.:", problems[:5])
    return 0


if __name__ == "__main__":
    sys.exit(main())
