"""W5 2 x 2 ablation setup (DESIGN §9 W5 "Ablation", §10.1): {generic, curated GPR} x {generic, best constraint
pipeline}, on identical gene x line sets. Curated GPRs are applied to the *already-extracted* reaction sets
(``nextgem.w5.curated_for_extracted``); re-extraction with curated scores is a P1-only sensitivity analysis.

Status: the constraint arm is a stub (``nextgem.w5.best_constraint_pipeline`` raises NotImplementedError) until
W2/W3/W4 have shortlists. This script writes the setup and checks the curated-GPR path on the generic-constraint
row: applying the curated trees to the generic parent through ``curated_for_extracted`` + a full KO sweep must
reproduce the W5a arm of results/w5/ko_<parent>.parquet exactly (label-free, GLPK).

Output: results/w5/factorial_2x2_setup.json.  Usage: .venv/bin/python scripts/w5_factorial.py [--parent human2]
"""

from __future__ import annotations

import argparse
import json
import sys

import cobra

cobra.Configuration().solver = "glpk"

import pandas as pd  # noqa: E402

from nextgem import w5  # noqa: E402


def main(prov, a) -> int:
    t = w5.Timer()
    p = a.parent
    m = w5.load_parent(p, tyms_irreversible=(p == "human1"))
    eng = w5.KOEngine(m)
    cur = pd.read_csv(w5.w5_dir() / f"curated_W5a_{p}.tsv", sep="\t")
    trees = {rx: w5.parse_rule(r.curated_rule) for r in cur[cur.edited].itertuples() for rx in r.reactions.split(";")}
    gp = w5.curated_for_extracted(m, trees)
    ko = pd.read_parquet(w5.w5_dir() / f"ko_{p}.parquet").set_index("gene")
    got = w5.full_sweep(eng, gp, list(ko.index))
    diff = (got - ko.ratio_W5a).abs()
    try:
        w5.best_constraint_pipeline()
        stub = "implemented"
    except NotImplementedError as e:
        stub = f"stub: {e}"
    setup = {
        "design": "DESIGN §9 W5 Ablation: 2 x 2 on identical gene x line sets; broken down by edit class",
        "arms": [{"gpr": g, "constraints": c} for g, c in w5.FACTORIAL_ARMS],
        "curated_gpr": {"H2": "results/w5/curated_W5a_human2.tsv", "H1": "results/w5/curated_W5a_human1.tsv (+ MAR04644 irreversible)"},
        "apply_rule": "curated trees substituted for reactions present in the extracted model; reactions removed by "
                      "extraction stay removed; no re-extraction (P1-only sensitivity: re-extract with curated scores)",
        "edit_class_breakdown": ["RMBRANCH (pan-cancer absent)", "AND2OR (paralogs, no shared complex)", "OR2AND", "TYMS (H1)"],
        "constraint_arm": stub,
        "check_generic_row": {"parent": p, "max_abs_diff_vs_W5a_arm": float(diff.max()), "n_differ": int((diff > 0).sum()),
                              "n_genes": len(diff)},
        **t.lap(),
    }
    if setup["check_generic_row"]["n_differ"]:
        raise AssertionError(setup["check_generic_row"])
    f = w5.w5_dir() / "factorial_2x2_setup.json"
    f.write_text(json.dumps(setup, indent=1))
    prov.output(f)
    print(json.dumps(setup, indent=1))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", default="human2", choices=w5.PARENTS)
    args = ap.parse_args()
    sys.exit(w5.run_main(lambda prov: main(prov, args), solver="glpk", workers=1))
