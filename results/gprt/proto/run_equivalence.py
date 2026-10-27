"""Equivalence test of the GPR transform on generic Human1/Human2 (Ham's medium M0, biomass), label-free.

For each model and isozyme encoding: GPR statistics, model size before/after, transform time, WT solve cost,
and single-gene KO of every gene (u__g bounded to 0; every gene solved, no shortcut) compared with
results/baseline/human{1,2}_ko.parquet (Gurobi, Boolean GPR KO).

Usage: .venv/bin/python results/gprt/proto/run_equivalence.py --model human2 --split arm --workers 8
Writes results/gprt/equivalence_<model>_<split>.parquet and results/gprt/sizes_<model>_<split>.json.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
ROOT = HERE.parents[2]
_TM = None


def _init(blob: bytes) -> None:
    global _TM
    import cobra

    cobra.Configuration().solver = "glpk"
    _TM = pickle.loads(blob)


def _ko_chunk(genes):
    sys.path.insert(0, str(HERE))
    import gpr_transform as gt
    from optlang.interface import INFEASIBLE, OPTIMAL

    rows = []
    for g in genes:
        c0 = time.process_time()
        with _TM:
            gt.knock_out(_TM, g)
            _TM.slim_optimize()
            st = _TM.solver.status
            val = _TM.solver.objective.value if st == OPTIMAL else (0.0 if st == INFEASIBLE else math.nan)
        rows.append((g, val, st, time.process_time() - c0))
    return rows, os.getpid()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=["human1", "human2"], required=True)
    ap.add_argument("--split", choices=["arm", "paper"], default="arm")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--no-sweep", action="store_true")
    a = ap.parse_args()

    import cobra

    cobra.Configuration().solver = "glpk"  # before any model is loaded
    import pandas as pd

    sys.path.insert(0, str(HERE))
    import gpr_transform as gt
    from nextgem.knockouts import configure, rounded_ratio
    from nextgem.model import load_model

    t0 = time.time()
    m = load_model(a.model)
    configure(m)  # Ham's medium (GR task of the model's release) + biomass objective
    c0 = time.process_time()
    wt0 = m.slim_optimize()
    wt0_cpu = time.process_time() - c0
    stats = gt.gpr_stats(m)
    tm, info = gt.transform(m, split=a.split)
    c0 = time.process_time()
    wt1 = tm.slim_optimize()
    wt1_cpu_cold = time.process_time() - c0
    c0 = time.process_time()
    tm.slim_optimize()
    wt1_cpu_warm = time.process_time() - c0
    # cold re-solve of the original for a like-for-like comparison
    sizes = {
        "model": a.model, "split": a.split, "solver": "glpk", **stats,
        **{k: v for k, v in info.items() if k not in ("map", "terms", "genes")},
        "wt_growth_original_glpk": wt0, "wt_growth_transformed": wt1, "wt_abs_diff": abs(wt1 - wt0),
        "wt_cpu_s_original_cold": wt0_cpu, "wt_cpu_s_transformed_cold": wt1_cpu_cold,
        "wt_cpu_s_transformed_warm": wt1_cpu_warm,
    }
    print(json.dumps(sizes, indent=1))
    if not a.no_sweep:
        base = pd.read_parquet(ROOT / "results" / "baseline" / f"{a.model}_ko.parquet")
        genes = list(base.gene)
        chunks = [genes[i:i + 40] for i in range(0, len(genes), 40)]
        blob = pickle.dumps(tm)
        w0 = time.time()
        with ProcessPoolExecutor(a.workers, mp_context=get_context("spawn"), initializer=_init, initargs=(blob,)) as pool:
            res = list(pool.map(_ko_chunk, chunks))
        rows = [r for rr, _ in res for r in rr]
        ko = pd.DataFrame(rows, columns=["gene", "ko_growth_t", "status_t", "cpu_s_t"])
        ko["wt_growth_t"] = wt1
        ko["ratio_t"] = ko.ko_growth_t / wt1
        ko["ratio_r_t"] = rounded_ratio(ko.ratio_t)
        df = base[["gene", "symbol", "wt_growth", "ko_growth", "ratio", "ratio_r", "status"]].merge(ko, on="gene")
        df["abs_diff_ratio"] = (df.ratio_t - df.ratio).abs()
        df["abs_diff_ratio_r"] = (df.ratio_r_t - df.ratio_r).abs()
        df["lethal_base"] = df.ratio_r < 0.01
        df["lethal_t"] = df.ratio_r_t < 0.01
        df.to_parquet(OUT / f"equivalence_{a.model}_{a.split}.parquet", index=False)
        mism = df[(df.abs_diff_ratio > 1e-6) | (df.lethal_base != df.lethal_t)]
        sizes.update({
            "sweep_genes": len(df), "sweep_wall_s": time.time() - w0, "sweep_workers": a.workers,
            "sweep_cpu_s_sum": float(df.cpu_s_t.sum()), "sweep_cpu_s_median": float(df.cpu_s_t.median()),
            "status_counts": df.status_t.value_counts().to_dict(),
            "n_mismatch_ratio_gt_1e-6": int((df.abs_diff_ratio > 1e-6).sum()),
            "n_mismatch_ratio_r": int((df.abs_diff_ratio_r > 0).sum()),
            "max_abs_diff_ratio": float(df.abs_diff_ratio.max()),
            "n_lethal_base": int(df.lethal_base.sum()), "n_lethal_t": int(df.lethal_t.sum()),
            "lethal_set_equal": bool((df.lethal_base == df.lethal_t).all()),
            "mismatch_genes": mism.symbol.fillna(mism.gene).tolist()[:50],
        })
        print({k: sizes[k] for k in list(sizes)[-12:]})
    sizes["total_wall_s"] = time.time() - t0
    (OUT / f"sizes_{a.model}_{a.split}.json").write_text(json.dumps(sizes, indent=1, default=str))


if __name__ == "__main__":
    main()
