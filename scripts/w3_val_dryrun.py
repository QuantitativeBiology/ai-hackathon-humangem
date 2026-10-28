"""W3 proposed validation look: config entry + label-free dry run (DESIGN §10, §11.3; amendment A6 / user D7).

NO validation label is read and no ledger row is written. Under A6 a new look reuses the reviewed
``scripts/evaluate_val.py`` harness and needs (1) its config entry reviewed by a non-author agent, (2) a label-free
dry run, (3) the pre-flight. This script produces (1) as a proposal and runs (2):

1. Label-free score vectors for the look's methods on the validation expression lines: generic H1 / H2 (= T0) and the
   W3 shortlist (``results/w3/train/shortlist.json``) on G* = Human1 (and the same modes on Human2 as secondary).
   Line-constant modes only need the gene vectors; T5 would use thresholds from all training-line expression.
   Written to results/w3/val/w3_val_scores.parquet (Entrez x method; no line IDs, no labels).
2. The proposed config (``val_look_W3_M4_shortlist_config.proposed.json``) with input/code SHA-256 and its own SHA,
   in the same layout as ``evaluate_val.w1_config``.
3. Dry run: the look's scalar pipeline (w1.method_summary / delta_summary vs H1 and H2; L_prob; all genes and
   minus CE_fit [from labels.train(), as §5.3 prescribes for validation looks]; U; strata none / expressed pairs)
   on SYNTHETIC labels of the validation shape (Bernoulli, SEEDS[0], prevalence 0.124 = the training null), then
   ``evaluate_val.assert_scalar_only`` on the output. It checks that every score matrix is finite and that the
   output is scalar-only. Nothing here touches labels.validation() or data/locked/.

Usage: .venv/bin/python scripts/w3_val_dryrun.py
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w3_common import OUT_DIR, PROJECT, atomic_write_parquet, atomic_write_text, current, run_main  # noqa: E402

VDIR = OUT_DIR / "val"
LOOK_ID = "W3_M4_shortlist"
CODE = ("scripts/evaluate_val.py", "src/nextgem/labels.py", "src/nextgem/w1.py", "src/nextgem/metrics.py",
        "src/nextgem/run.py", "src/nextgem/w3.py", "scripts/w3_score_train.py")


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    t0 = time.time()
    sys.path.insert(0, str(PROJECT / "scripts"))
    import evaluate_val as ev
    import w3_score_train as st

    from nextgem import labels, metrics, w1
    from nextgem.config import SEEDS

    sl = json.loads((OUT_DIR / "train" / "shortlist.json").read_text())
    cands = [(c["parent"], c["mode"]) for c in sl.get("candidates", [])]
    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    Udf = labels.gene_universe()
    U = list(Udf.entrez.astype(int))
    modes, mt, extra = st.build_modes(U, [])
    methods = {"H1": modes[("human1", "T0")], "H2": modes[("human2", "T0")]}
    for p, m in cands:
        for q in ("human1", "human2"):
            if (q, m) in modes:
                methods[f"W3:{m}[{st.SHORT[q]}]"] = modes[(q, m)]
            else:
                raise RuntimeError(f"{m} is line-specific; extend this dry run before proposing it")
    V = pd.DataFrame(methods).rename_axis("entrez").reset_index()
    sp = atomic_write_parquet(V, VDIR / "w3_val_scores.parquet")
    val_lines = sorted(split.index[split.tier == "validation"])
    vexpr = [l for l in val_lines if lf.loc[l, "has_expression_matrix"]]
    inputs = {"w3_val_scores": str(sp.relative_to(PROJECT)), "line_features": "results/w1/line_features.tsv",
              "expr": "results/w1/tmp/expr.parquet", "gene_panels": "data/controls/gene_panels.tsv",
              "w3_shortlist": "results/w3/train/shortlist.json", "w3_mode_scores": "results/w3/train/mode_scores.tsv"}
    cfg = {"look_id": LOOK_ID, "workstream": "W3", "milestone": "M4", "labels": ["L_prob"], "genesets": ["all", "noCE"],
           "universes": ["U"], "strata": ["none", "expressed_pairs"], "B": metrics.B_OTHER, "seed": SEEDS[0],
           "methods": sorted(methods), "reference": "H1 (= G*, generic Human1 = W3 T0 on M0)",
           "decision_rule": ("§8 gate: a W3 candidate is eligible only if its paired validation Δ̄_lin vs G* (H1) has CI "
                             "lower bound > 0; among eligible candidates the paired-difference one-SE rule (§9). "
                             "Context claims: same sign on validation lines of the same context (W3), reported only."),
           "scalars_only": True, "n_validation_expression_lines": len(vexpr),
           "inputs": {k: _sha(PROJECT / v) for k, v in inputs.items()},
           "code": {Path(c).name: _sha(PROJECT / c) for c in CODE},
           "harness_change_required": ("evaluate_val.py implements only W1_M2_baseline; this look needs a 'W3_M4_shortlist' "
                                       "branch that loads w3_val_scores.parquet as line-constant methods. That is a "
                                       "harness change, so A6(2) requires a full independent code review of the diff."),
           "_inputs_rel": inputs}
    csha = ev.config_sha(cfg)
    doc = {k: v for k, v in cfg.items() if not k.startswith("_")}
    cp = atomic_write_text(json.dumps({**doc, "config_sha": csha}, indent=1, sort_keys=True),
                           VDIR / f"val_look_{LOOK_ID}_config.proposed.json")
    # ------------------------------------------------------------ dry run on synthetic labels (no label file read)
    S = {n: pd.DataFrame(np.tile(v.reindex(U).to_numpy(), (len(vexpr), 1)), index=vexpr, columns=U) for n, v in methods.items()}
    for n, df in S.items():
        if not np.isfinite(df.to_numpy()).all():
            raise RuntimeError(f"non-finite validation scores for {n}")
    rng = np.random.default_rng(SEEDS[0])
    Ysyn = pd.DataFrame((rng.random((len(vexpr), len(U))) < 0.124).astype(float), index=vexpr, columns=U)
    train_long = labels.train(columns=["gene_effect", "dep_prob", "L_prob", "L_ge05", "L_ge06", "L_ge10"])
    fitted = w1.fit_components(train_long, lf.base_medium_class)
    omics = w1.load_omics()
    strata = w1.stratum_specs(w1.gene_panels(), w1.gene_features())
    lineage = split.lineage
    rows = []
    for gs in cfg["genesets"]:
        for stt in cfg["strata"]:
            Y = Ysyn if gs == "all" else Ysyn.loc[:, [g for g in U if g not in fitted.ce["L_prob"]]]
            Y = w1.apply_stratum(Y, strata[stt], lf, omics["expr"], omics["cn"], None)
            lts, cell = {}, []
            for n in sorted(methods):
                r, lt, _ = w1.method_summary(S[n].loc[Y.index, Y.columns], Y, lineage, seed=cfg["seed"], B=cfg["B"])
                lts[n] = lt
                r.update(label="L_prob", geneset=gs, universe="U", stratum=stt, method=n)
                cell.append(r)
            for ref in ("H1", "H2"):
                for r in cell:
                    if r["method"] != ref:
                        d = w1.delta_summary(lts[r["method"]], lts[ref], lineage, seed=cfg["seed"], B=cfg["B"])
                        r.update({f"{k}_vs_{ref}": v for k, v in d.items()})
            rows += cell
    out = pd.DataFrame(rows).assign(section="tAP_grid")
    out = out[["section"] + [c for c in out.columns if c != "section"]]
    for c in ("protocol", "call_rule", "gene_set", "metric", "nan_policy"):  # binary-section keys (no rows here)
        out[c] = pd.Series([None] * len(out), dtype=object)
    syms = set()
    for p in ("human1", "human2"):
        syms |= set(pd.read_csv(PROJECT / "results" / "baseline" / f"{p}_gene_mapping.tsv", sep="\t", usecols=["symbol"]).symbol.dropna())
    ctx = {"lineage": lineage, "genes": U, "gene_symbols": syms, "all_ids": list(split.index)}
    guard_cfg = {"labels": ["L_prob"], "genesets": cfg["genesets"], "universes": cfg["universes"], "strata": cfg["strata"]}
    ev.assert_scalar_only(out, guard_cfg, ctx, n_methods={"L_prob": len(methods)})
    dr = {"look_id": LOOK_ID, "config_sha": csha, "synthetic_labels": "Bernoulli(0.124), SEEDS[0]",
          "n_methods": len(methods), "n_rows": len(out), "scalar_guard": "pass", "validation_labels_read": False,
          "ledger_rows_written": 0, "wall_s": round(time.time() - t0, 1),
          "synthetic_tap_lin_H1": float(out[(out.method == "H1") & (out.geneset == "all") & (out.stratum == "none")].tap_lin.iloc[0])}
    atomic_write_text(json.dumps(dr, indent=1), VDIR / "dry_run.json")
    current().output(sp, cp, VDIR / "dry_run.json")
    print(json.dumps(dr, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(run_main(main, seeds=[11], solver="", workers=1, reads_labels=True))
