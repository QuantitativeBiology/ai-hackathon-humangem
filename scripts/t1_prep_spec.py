"""T1 prep (1): assemble results/m6/p1_spec.json (template results/m6/p1_spec.template.json) from the outputs of
scripts/t1_prep_tasks.py, t1_prep_build.py and t1_prep_skew.py. Label-free (reads skew.json, a scalar decision).

Every file the spec names is pinned by its current SHA-256; the spec is checked with evaluate_test.check_p1_spec()
before it is written. It is NOT committed or frozen here: the PI commits it and adds the FROZEN.tsv rows listed by
``scripts/evaluate_test.py --print-frozen-candidates``.

Usage: .venv/bin/python scripts/t1_prep_spec.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as mc  # noqa: E402
import t1_prep_build as tb  # noqa: E402
import t1_prep_common as tc  # noqa: E402

SPEC_ID = "P1_M6"
SELECTING_LOOK = "M6_factorial"
PROVENANCE_FILES = (
    "scripts/t1_prep_common.py", "scripts/t1_prep_tasks.py", "scripts/t1_prep_build.py", "scripts/t1_prep_skew.py",
    "scripts/t1_prep_spec.py", "scripts/m6_common.py", "scripts/m6_scores.py", "scripts/m6_sweeps.py",
    "results/m6/scores/SHA256SUMS.tsv", "results/m6/scores/sweep_e.parquet", "results/m6/scores/line_assignment.tsv",
    "results/m6/sweeps/catalog.tsv", "results/m6/t1/tasks_bloo_human1.parquet", "results/m6/t1/tasks_bloo_human2.parquet",
    "results/m6/t1/tasks_bloo_human1_tasks.tsv", "results/m6/t1/tasks_bloo_human2_tasks.tsv",
    "results/m6/t1/build_summary.json", "results/m6/t1/skew.json",
    "results/w3/structure/variants_human1.json", "results/w3/structure/variants_human2.json",
)


def _entry(path: Path, layout: str | None = None) -> dict:
    e = {"path": tc.rel(path), "sha256": tc.sha256(path)}
    if layout:
        e["layout"] = layout
    return e


def build_spec() -> dict:
    skew = json.loads((tc.T1 / "skew.json").read_text())
    if skew.get("synthetic_labels"):
        raise RuntimeError("skew.json comes from synthetic labels")
    summ = json.loads((tc.T1 / "build_summary.json").read_text())
    if not summ["all_checks_pass"]:
        raise RuntimeError("build_summary.json: a build check failed")
    O = tb.OUT
    notes = (
        "P1 (RESEARCH_TREE D11, look 3/40 M6_factorial): M6 cell (3). Flags carried: W3-D1, W3-MEV (mevalonate core "
        "non-lethal under serum on ~all lines), W2-SHUF, W2-MOPEN. Methods: H1/H2 generic; P1[H1] = M6:H1+Bloo/M2(f=1); "
        "P1[H2] = M6:H2+Bloo/M2(f=1) (within-parent arm); M2[H1] = M6:H1/M2(f=1), the attribution arm (ii) 'medium "
        "alone' (M0 -> M2(f=1) on the generic parent H1 = G*). Curation (W5) and context (M5) not selected (§10.1, "
        "A12c): (iii) and (v) 'lever not selected'; (iv) 'constraints alone' = P1 with generic GPRs vs G* = P1 itself "
        "vs G* (identical to (i) by construction). Calls: p57/pgr under the literal paper protocol (each task with its "
        "own inputs, GR = growth on Ham's) from the Essential-task matrix re-run with the B-loo biomass "
        "(scripts/t1_prep_tasks.py), hence line-constant per parent; fba = ratio_r < 0.01 of the line's M2(f=1) "
        "score (line-specific). No failed build (M6 sweeps 42/42 ok; no non-finite score). Test-line score and call "
        "files are git-ignored under results/m6/tmp/t1/ and regenerable by scripts/t1_prep_build.py. Gene "
        "annotations: 6 of 8 classes (gene_less_bypass and model_contradicted_loss not annotated). bootstrap_method "
        f"from the §6.5 training skew check: skew C1 {skew['C1']['skew']:.4f}, C1b {skew['C1b']['skew']:.4f} "
        f"(n = {skew['C1']['n_lines']} training lines) -> {skew['bootstrap_method']}."
    )
    spec = {
        "spec_id": SPEC_ID,
        "kind": "p1",
        "milestone": "M6",
        "pipeline": (f"Human1 + W3 B-loo:cofactor_pool_biomass on documented medium M2(f=1), gate-free, generic "
                     f"parent, no context family (selected in the M6 validation look {SELECTING_LOOK}, D11); P1[H2] = "
                     f"the same pipeline on Human2."),
        "primary": tc.PRIMARY,
        "methods": list(tc.METHODS),
        "arms": dict(tc.ARMS),
        "scores": {t: _entry(O[f"scores_{t}"], "line_gene") for t in ("validation", "test")},
        "calls": {t: _entry(O[f"calls_{t}"], "line_gene") for t in ("validation", "test")},
        "failed_builds": {"validation": None, "test": None},
        "attribution": {k: (list(v) if isinstance(v, list) else v) for k, v in tc.ATTRIBUTION.items()},
        "task_calls": _entry(O["task_calls"]),
        "gene_annotations": _entry(O["gene_annotations"]),
        "bootstrap_method": skew["bootstrap_method"],
        "fitted_on_labels": False,
        "provenance": {p: tc.sha256(tc.PROJECT / p) for p in PROVENANCE_FILES},
        "notes": notes,
    }
    return spec


def check(spec: dict) -> None:
    s = importlib.util.spec_from_file_location("evaluate_test", tc.PROJECT / "scripts" / "evaluate_test.py")
    et = importlib.util.module_from_spec(s)
    s.loader.exec_module(et)
    et.check_p1_spec(spec)


def main(argv=None) -> int:
    spec = build_spec()
    check(spec)
    with tc.provenance(inputs=[tc.PROJECT / p for p in PROVENANCE_FILES]) as prov:
        tc.atomic_write_json(spec, tc.SPEC_PATH)
        prov.output(tc.SPEC_PATH)
    print(json.dumps({"spec": tc.rel(tc.SPEC_PATH), "sha256": tc.sha256(tc.SPEC_PATH),
                      "bootstrap_method": spec["bootstrap_method"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
