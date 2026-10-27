"""Reviewer re-check dry run of the revised scripts/evaluate_val.py (no validation label is read).

(1) Label-free on the real validation lines: real_context() + prepare_look() (fit on labels.train(), all
    validation-line score matrices) -- exactly the pre-'started' work of the look.
(2) Pseudo-look on 4 *training* lineages standing in for validation (fit on the remaining training lines):
    w1_scalars() incl. the hardened assert_scalar_only() and decide_gstar(). No ledger, no output table;
    prints shapes / pass-fail / runtimes only (no metric values).

Usage: .venv/bin/python results/w1/review/dryrun_evaluate_val_v2.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def main() -> int:
    from nextgem import labels

    spec = importlib.util.spec_from_file_location("evaluate_val", ROOT / "scripts" / "evaluate_val.py")
    ev = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ev)
    rep = {}
    t0 = time.time()
    cfg = ev.w1_config()
    rep["config_sha"] = ev.config_sha(cfg)
    rep["include_b6"] = cfg["include_b6"]
    ctx = ev.real_context(cfg)
    rep["n_val_lines"] = len(ctx["val_lines"])
    rep["n_val_expr_lines"] = len(ctx["val_expr_lines"])
    train_long = labels.train(columns=ev.LABEL_COLS)
    t1 = time.time()
    prep = ev.prepare_look(cfg, ctx, train_long)  # raises on non-finite scores
    rep["val_methods"] = sorted(prep["S"])
    rep["t_prepare_real_val_s"] = round(time.time() - t1, 1)

    pseudo = labels.folds()[:4]
    fit_long = train_long[~train_long.lineage.isin(pseudo)].reset_index(drop=True)
    pv_long = train_long[train_long.lineage.isin(pseudo)].drop(columns="lineage").reset_index(drop=True)
    pv_lines = sorted(l for l in pv_long.ModelID.unique() if ctx["lf"].loc[l, "has_expression_matrix"])
    ctx2 = dict(ctx, val_lines=pv_lines, val_expr_lines=pv_lines)
    t2 = time.time()
    prep2 = ev.prepare_look(cfg, ctx2, fit_long)
    out, gstar = ev.w1_scalars(pv_long, prep2, ctx2, cfg)  # calls the hardened assert_scalar_only
    rep["pseudo_lines"] = len(pv_lines)
    rep["pseudo_out_rows"] = int(len(out))
    rep["pseudo_out_rows_by_section"] = out.section.value_counts().to_dict()
    rep["pseudo_out_cols"] = int(out.shape[1])
    rep["pseudo_guard_passed"] = True
    rep["pseudo_gstar_keys"] = sorted(gstar)
    rep["pseudo_gstar_consistent_with_sign"] = gstar["G_star"] == ("Human1" if gstar["delta_lin_H2_vs_H1"] < 0 else "Human2")
    rep["t_pseudo_look_s"] = round(time.time() - t2, 1)
    rep["t_total_s"] = round(time.time() - t0, 1)
    (Path(__file__).resolve().parent / "dryrun_evaluate_val_v2.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
