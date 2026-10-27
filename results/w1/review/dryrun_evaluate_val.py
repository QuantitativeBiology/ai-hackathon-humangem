"""Reviewer dry run of scripts/evaluate_val.py's *real-data* code path without any validation label.

(1) Label-free: real_context() and the score matrices (incl. the B6 full model) on the real validation
    expression lines -- nothing label-bearing for validation is touched.
(2) Pseudo-look: w1_scalars() on 4 *training* lineages treated as "validation", with components fitted on
    the remaining training lines (labels.train() only). No ledger row, no output file; prints only shape /
    pass-fail / runtime (no metric values).

Usage: .venv/bin/python results/w1/review/dryrun_evaluate_val.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def main() -> int:
    from nextgem import labels, w1

    spec = importlib.util.spec_from_file_location("evaluate_val", ROOT / "scripts" / "evaluate_val.py")
    ev = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ev)
    rep = {}
    t0 = time.time()
    cfg = ev.w1_config()
    rep["config_sha"] = ev.config_sha(cfg)
    ctx = ev.real_context(cfg)
    rep["n_val_lines"] = len(ctx["val_lines"])
    rep["n_val_expr_lines"] = len(ctx["val_expr_lines"])
    rep["t_context_s"] = round(time.time() - t0, 1)

    # (1) label-free path on real validation lines (fit on all training = the look's prepare step)
    train_long = labels.train(columns=ev.LABEL_COLS)
    t1 = time.time()
    prep = ev.prepare_w1(train_long, ctx)
    b6 = ctx["b6"](ctx["val_expr_lines"], prep["fitted"])
    S = w1.score_matrices(ctx["val_expr_lines"], ctx["genes"], prep["fitted"], ctx["gen"], ctx["line_class"],
                          ctx["omics"]["expr"], b6=b6, random_seeds=cfg["random_seeds"])
    rep["val_score_methods"] = len(S)
    rep["val_scores_all_finite"] = bool(all(np.isfinite(v.to_numpy()).all() for v in S.values()))
    rep["val_unknown_medium_class_lines"] = int(sum(ctx["line_class"].get(l, None) not in (ev.pickle.loads(
        Path(cfg["_paths"]["b6_full"]).read_bytes())["class_codes"]) for l in ctx["val_expr_lines"]))
    rep["t_label_free_scores_s"] = round(time.time() - t1, 1)

    # (2) pseudo-look on 4 training lineages
    folds = labels.folds()
    pseudo = folds[:4]
    fit_long = train_long[~train_long.lineage.isin(pseudo)].reset_index(drop=True)
    pv_long = train_long[train_long.lineage.isin(pseudo)].drop(columns="lineage").reset_index(drop=True)
    lf = ctx["lf"]
    pv_lines = sorted(l for l in pv_long.ModelID.unique() if lf.loc[l, "has_expression_matrix"])
    ctx2 = dict(ctx, val_lines=pv_lines, val_expr_lines=pv_lines)
    t2 = time.time()
    prep2 = ev.prepare_w1(fit_long, ctx2)
    out, gstar = ev.w1_scalars(pv_long, prep2, ctx2, cfg)
    rep["pseudo_lineages"] = len(pseudo)
    rep["pseudo_lines"] = len(pv_lines)
    rep["pseudo_out_rows"] = int(len(out))
    rep["pseudo_out_cols"] = int(out.shape[1])
    rep["pseudo_sections"] = sorted(out.section.unique())
    rep["pseudo_gstar_keys"] = sorted(gstar)
    rep["pseudo_scalar_guard_passed"] = True  # w1_scalars calls assert_scalar_only itself
    rep["t_pseudo_look_s"] = round(time.time() - t2, 1)
    rep["t_total_s"] = round(time.time() - t0, 1)
    (Path(__file__).resolve().parent / "dryrun_evaluate_val.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
