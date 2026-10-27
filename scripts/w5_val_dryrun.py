"""Proposed W5 validation-look config and its label-free dry run (amendment A6). Reads NO validation label.

The frozen ``scripts/evaluate_val.py`` implements only the W1 look; this script writes the *proposed* config entry
for a W5 look (to be reviewed by a non-author agent; a new look branch in evaluate_val.py is a harness change and
needs the full code review, A6) and dry-runs everything that does not need validation labels:
  - validation lines (split, label-free) with a default expression profile;
  - score matrices of every proposed method on those lines (line-constant W5a/TYMS from results/w5/ko_*.parquet;
    W5b refitted on ALL training lines through labels.train(), as §5.3 prescribes for validation looks);
  - pairing assertions (identical gene x line sets, no NaN), score SHA-256 per method, config SHA.
Outputs: results/w5/proposed_val_look_W5.json, results/w5/val_dryrun.json.
Usage: .venv/bin/python scripts/w5_val_dryrun.py --methods H1_W5a [H1_W5b ...]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys

import cobra

cobra.Configuration().solver = "glpk"

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from nextgem import paths, w5  # noqa: E402
from nextgem import w5_evidence as ev  # noqa: E402
from nextgem.config import SEEDS  # noqa: E402

LOOK_ID = "W5_M3_gpr_curation"


def sha_obj(o) -> str:
    return hashlib.sha256(json.dumps(o, sort_keys=True, default=str).encode()).hexdigest()


def sha_file(p) -> str:
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def w5b_full_training(parent: str, genes, variant: str = "W5b") -> tuple[pd.Series, list]:
    """W5b rule fitted on all training lines (labels.train(); validation-look fitting set, §5.3)."""
    from nextgem import labels

    sys.path.insert(0, str(paths.SCRIPTS_DIR))
    from w5_score_train import coess_ok  # same rule as the training scorer

    cands = json.loads((w5.w5_dir() / "w5b_candidates.json").read_text())[parent]
    ent = w5.gene_table(parent).entrez.dropna().astype(int).to_dict()
    need = sorted({ent[g] for v in cands.values() for c in v["candidates"] for pr in c["pairs"] for g in pr if g in ent})
    tl = labels.train(columns=["gene_effect"])
    ge = tl[tl.entrez.isin(need)].pivot(index="ModelID", columns="entrez", values="gene_effect")
    m = w5.load_parent(parent, tyms_irreversible=(parent == "human1"))
    eng, gp = w5.KOEngine(m), w5.model_gprs(m)
    k = pd.read_parquet(w5.w5_dir() / f"ko_{parent}.parquet")
    cur = pd.read_csv(w5.w5_dir() / f"curated_W5a_{parent}.tsv", sep="\t")
    new = {rx: w5.parse_rule(r.curated_rule) for r in cur[cur.edited].itertuples() for rx in r.reactions.split(";")}
    log = []
    diag = []  # training-tier diagnostic: per candidate OR node, the minimum pair r (all training lines)
    for grp, v in cands.items():
        ps = []
        for c in v["candidates"]:
            rr_ = []
            for a_, b_ in c["pairs"]:
                ea, eb = ent.get(a_), ent.get(b_)
                if ea in ge.columns and eb in ge.columns:
                    x = ge[[ea, eb]].dropna()
                    rr_.append(float(np.corrcoef(x[ea], x[eb])[0, 1]) if len(x) >= ev.COESS_MIN_N else np.nan)
                else:
                    rr_.append(np.nan)
            diag.append({"group": grp, "anchored": c["anchored"], "n_pairs": len(rr_),
                         "min_pair_r": float(np.nanmin(rr_)) if np.isfinite(rr_).any() else None,
                         "max_pair_r": float(np.nanmax(rr_)) if np.isfinite(rr_).any() else None})
            if variant == "W5b" and not c["anchored"]:
                continue
            ok, rs = coess_ok(ge, [tuple(x) for x in c["pairs"]], ent)
            if ok:
                ps.append(tuple(c["path"]))
                log.append({"group": grp, "pairs": c["pairs"], "r": rs})
        if ps:
            t = ev.convert_or_nodes(w5.parse_rule(v["curated_rule"]), ps)
            for rx in v["reactions"]:
                new[rx] = t
    rr = w5.local_sweep(eng, gp, k.set_index("gene").ratio_base, new, 1)
    return w5.e_by_entrez(k.gene, k.entrez, rr.reindex(k.gene).to_numpy(), genes), {"edits": log, "candidates": diag}


def proposed_config(methods) -> dict:
    files = {f"ko_{p}": f"results/w5/ko_{p}.parquet" for p in w5.PARENTS}
    files.update({f"curated_W5a_{p}": f"results/w5/curated_W5a_{p}.tsv" for p in w5.PARENTS})
    files.update({f"edits_W5a_{p}": f"results/w5/edits_W5a_{p}.tsv" for p in w5.PARENTS})
    files["w5b_candidates"] = "results/w5/w5b_candidates.json"
    root = paths.PROJECT_ROOT
    cfg = {
        "look_id": LOOK_ID, "workstream": "W5", "milestone": "M3 (shortlist -> validation)",
        "harness": "scripts/evaluate_val.py (reviewed W1 harness, reused per A6); needs a new --look branch that loads the "
                   "W5 score inputs below: a harness code change => full independent code review (A6 item 2)",
        "reference": "G* = H1 (generic Human1 v1.19.0 on M0; D6)",
        "status": "NOT TO BE RUN under the frozen rule: no W5 candidate passed TG on training (results/w5/README.md). "
                  "Provided so that a PI decision (or the §10.1 factorial, if P1 carries curation) can be executed "
                  "without new W5 code; the methods listed are the training-best curated arms, not a shortlist",
        "methods": list(methods), "comparators_shown": ["H1", "H2", "B2", "B2m", "B8+_H1", "B8-_H1", "B8'_H1"],
        "W5b_fit": "co-essentiality (Chronos r >= 0.4, >= 30 lines) on ALL training lines via labels.train() (§5.3)",
        "cells": {"labels": ["L_prob", "L_ge05", "L_ge10"], "genesets": ["all", "noCE"], "universes": ["U", "Ucap"],
                  "strata": ["none", "expressed_pairs"], "lines": "validation lines with a default expression profile"},
        "outputs": "scalars only (assert_scalar_only): tAP_line lineage-equal mean, lift, paired Δ̄_lin vs H1 with "
                   "lineage-stratified bootstrap CI (B = 2,000, SEEDS[0]); fixed-cutoff MCC; no per-line/per-gene tables",
        "decision": "§8 gate: eligible iff paired validation Δ̄_lin vs G* has CI lower bound > 0; among eligible, "
                    "paired-difference one-SE rule (curation ranks after medium and tasks); winner's-curse σ√(2 ln m) reported",
        "B": 2000, "seed": SEEDS[0],
        "inputs": {k: {"path": v, "sha256": sha_file(root / v)} for k, v in files.items()},
        "code": {p: sha_file(root / p) for p in ("src/nextgem/w5.py", "src/nextgem/w5_evidence.py",
                                                 "scripts/w5_score_train.py", "scripts/w5_val_dryrun.py")},
    }
    cfg["config_sha"] = sha_obj(cfg)
    return cfg


def main(prov, a) -> int:
    from nextgem import labels

    t = w5.Timer()
    split = labels.load_split()
    lf = pd.read_csv(paths.results_dir() / "w1" / "line_features.tsv", sep="\t").set_index("ModelID")
    val = sorted(split[split.tier == "validation"].ModelID)
    V = [l for l in val if lf.loc[l, "has_expression_matrix"]]
    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    mats, w5b_logs = {}, {}
    for name in ["H1", "H2"] + list(a.methods):
        s, _, arm = name.partition("_")
        p = "human1" if s == "H1" else "human2"
        k = pd.read_parquet(w5.w5_dir() / f"ko_{p}.parquet")
        if arm in ("", "W5a", "W5a_pid30", "W5a_gpr_only", "TYMS") or arm.startswith("neg_"):
            col = {"": "ratio_generic", "TYMS": "ratio_base"}.get(arm, f"ratio_{arm}")
            e = w5.e_by_entrez(k.gene, k.entrez, k[col], genes)
        elif arm in ("W5b", "W5b_loose"):
            e, w5b_logs[name] = w5b_full_training(p, genes, arm)
        else:
            raise ValueError(name)
        mats[name] = pd.DataFrame(np.tile(e.to_numpy(), (len(V), 1)), index=V, columns=genes)
    shapes = {n: list(m.shape) for n, m in mats.items()}
    ref = mats["H1"]
    for n, m in mats.items():
        assert m.index.equals(ref.index) and m.columns.equals(ref.columns), n
        assert not m.isna().any().any(), n
    cfg = proposed_config(a.methods)
    out = w5.w5_dir()
    (out / "proposed_val_look_W5.json").write_text(json.dumps(cfg, indent=1))
    dry = {"look_id": LOOK_ID, "config_sha": cfg["config_sha"], "validation_labels_read": False,
           "n_validation_lines_expr": len(V), "n_validation_lines": len(val), "n_genes_U": len(genes),
           "shapes": shapes, "score_sha256": {n: sha_obj(m.round(6).to_numpy().tolist()) for n, m in mats.items()},
           "n_genes_differing_from_H1": {n: int((m.iloc[0] != ref.iloc[0]).sum()) for n, m in mats.items()},
           "w5b_full_training_edits": w5b_logs, "pairing_asserts": "passed", **t.lap()}
    (out / "val_dryrun.json").write_text(json.dumps(dry, indent=1, default=str))
    prov.output(out / "proposed_val_look_W5.json", out / "val_dryrun.json")
    print(json.dumps(dry, indent=1, default=str))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--methods", nargs="+", required=True)
    args = ap.parse_args()
    sys.exit(w5.run_main(lambda prov: main(prov, args), solver="glpk", workers=1, reads_labels=True))
