"""W2 shortlist validation look: label-free score parquet + look spec (DESIGN §10, §11.3; A6/A11f-g). Prepares the
look; it does NOT run it. No label is read (the split file and line media are label-free).

Scores (layout ``line_gene``: the documented medium makes W2 scores line-specific): for every validation line with
an expression profile (the harness's evaluation lines), e = 1 − ratio_r of the line's cached sweep
(``results/w2/controls/line_sweeps.tsv``), mapped to U as in W1/W3 (absent gene -> 0, 6 decimals):
  H1, H2                 generic parents on M0 (must equal the W1 look's generic scores; asserted here);
  W2:<step>[H1]          the W2 shortlist (``results/w2/train/shortlist.json``) on G* = Human1 (gate candidate);
  W2:<step>[H2]          the same step on Human2 (descriptive).
Validation lines outside the attribution set get the integrated-pipeline medium of DESIGN W2 (imputed, fallback and
unmappable lines: classic RPMI + 10% FBS; onboarded lines: their onboarded medium), exactly as in line_sweeps.tsv.

Outputs (results/w2/val/): w2_val_scores.parquet, look_spec_W2_M4_shortlist.json.

Usage: .venv/bin/python scripts/w2_val_scores.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
W2 = PROJECT / "results" / "w2"
VDIR = W2 / "val"
LOOK_ID = "W2_M4_shortlist"
PARENTS = ("human1", "human2")
SHORT = {"human1": "H1", "human2": "H2"}


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    from nextgem import labels, metrics, paths, w1, w2_score as ws
    from nextgem.config import SEEDS
    from nextgem.w2_runtime import Manifest

    sl = json.loads((W2 / "train" / "shortlist.json").read_text())
    if not sl["shortlist"]:
        print("W2 is negative (no shortlist): no validation look is prepared")
        return 0
    step = sl["shortlist"][0]
    VDIR.mkdir(parents=True, exist_ok=True)
    with Manifest("w2_val_scores.py", solver="none", workers=1,
                  inputs=[W2 / "ko" / f"ratio_r_{p}.parquet" for p in PARENTS]
                  + [W2 / "controls" / "line_sweeps.tsv", W2 / "train" / "shortlist.json"]) as man:
        split = labels.load_split().set_index("ModelID")
        lf = w1.line_features()
        Udf = labels.gene_universe()
        U = list(Udf.entrez.astype(int))
        val = sorted(split.index[split.tier == "validation"])
        V = [l for l in val if lf.loc[l, "has_expression_matrix"]]
        ls = pd.read_csv(W2 / "controls" / "line_sweeps.tsv", sep="\t")
        gen = w1.generic_table(Udf.entrez)
        cols = {}
        for p in PARENTS:
            ko = pd.read_parquet(paths.baseline_dir() / f"{p}_ko.parquet", columns=["gene", "entrez"])
            Esw = ws.sweep_e_matrix(pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet"), ws.gene_to_entrez(ko), U)
            for s, name in (("M0", SHORT[p]), (step, f"W2:{step}[{SHORT[p]}]")):
                sw = ls[(ls.parent == p) & (ls.ladder_step == s)].set_index("ModelID").sweep_id
                M = ws.line_matrix(sw, Esw, V)
                cols[name] = M
            g = np.tile(gen[f"e_{SHORT[p]}"].reindex(U).to_numpy(), (len(V), 1))
            if not np.array_equal(cols[SHORT[p]].to_numpy(), g):
                raise RuntimeError(f"{SHORT[p]} M0 scores differ from the W1 generic scores")
        long = pd.DataFrame({"ModelID": np.repeat(V, len(U)), "entrez": np.tile(np.asarray(U, dtype=np.int64), len(V))})
        for name, M in cols.items():
            long[name] = M.to_numpy(dtype=np.float64).ravel()
        out = VDIR / "w2_val_scores.parquet"
        tmp = out.with_suffix(".parquet.tmp")
        long.to_parquet(tmp, index=False, compression="zstd")
        tmp.replace(out)
        lm = pd.read_csv(W2 / "media_library" / "line_media.tsv", sep="\t").set_index("ModelID")
        n_attr = int((lm.stratum.reindex(V) == "attribution").sum())
        methods = ["H1", "H2", f"W2:{step}[H1]", f"W2:{step}[H2]"]
        prov_files = ["results/w2/train/shortlist.json", "results/w2/train/step_scores.tsv",
                      "results/w2/train/shuffled_control.tsv", "results/w2/controls/line_sweeps.tsv",
                      "results/w2/ko/ratio_r_human1.parquet", "results/w2/ko/ratio_r_human2.parquet",
                      "src/nextgem/w2_score.py", "scripts/w2_score_train.py", "scripts/w2_val_scores.py"]
        flags = sl.get("candidate_flags", {})
        spec = {
            "look_id": LOOK_ID, "kind": "shortlist", "workstream": "W2", "milestone": "M4",
            "scores": {"path": str(out.relative_to(PROJECT)), "sha256": _sha(out), "layout": "line_gene"},
            "methods": methods, "candidates": [f"W2:{step}[H1]"], "reference": "H1",
            "labels": ["L_prob"], "genesets": ["all", "noCE"], "universes": ["U"], "strata": ["none", "expressed_pairs"],
            "B": metrics.B_OTHER, "seed": SEEDS[0],
            "decision_rule": ("§8 gate: the W2 candidate is eligible only if its paired validation Δ̄_lin vs G* (H1, M0) "
                              "has CI lower bound > 0 (L_prob, all genes, U, no stratum); one candidate, so the §9 one-SE "
                              "rule picks it if eligible and σ·√(2 ln m) = 0 (A11g)."),
            "provenance": {f: _sha(PROJECT / f) for f in prov_files},
            "notes": (f"W2 M4 shortlist: {step} (documented base medium + supplements + serum model at f = "
                      f"{ws.parse_step(step).f:g}, gate-free) on G* = Human1; [H2] is descriptive. Scores are line-specific "
                      f"(layout line_gene) over the {len(V)} validation expression lines; {n_attr} are in the medium-"
                      "attribution set, the others get the integrated-pipeline medium (imputed/fallback/unmappable: classic "
                      "RPMI + 10% FBS; onboarded: onboarded medium), as in line_sweeps.tsv; the harness cannot restrict the "
                      "gate cell to the attribution set (strata are W1's). Training flags carried: W2-SHUF (beats_shuffled_"
                      f"medium = {flags.get('beats_shuffled_medium')}), W2-MOPEN (M-open Δ̄_lin vs M0 = "
                      f"{flags.get('m_open_delta_vs_M0', float('nan')):.4f}), W2-CONST (post hoc: one line-constant classic RPMI + "
                      "10% FBS medium for every line scores higher than the documented per-line media on training); see "
                      "results/w2/README.md §M4. "
                      "Not to be combined with W3 B-loo variants here (M6 factorial)."),
        }
        sp = VDIR / f"look_spec_{LOOK_ID}.json"
        sp.write_text(json.dumps(spec, indent=1) + "\n")
        man.add(out, sp)
        man.extra.update(step=step, n_val_expr_lines=len(V), n_attribution=n_attr, scores_sha256=spec["scores"]["sha256"])
        print(json.dumps({"step": step, "n_val_expr_lines": len(V), "n_in_attribution_set": n_attr,
                          "rows": len(long), "scores_sha256": spec["scores"]["sha256"]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
