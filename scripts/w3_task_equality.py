"""W3 §12 equality test: the P57 calls of the W3 capacity matrix equal Human-GEM ``taskEssentialGenes`` on the
generic parents (DESIGN §2, §6.8, §12 "the P57 call equals taskEssentialGenes.py on the generic parents").

For each parent: (a) Human-GEM ``find_task_essential_categories`` (reused, not re-implemented) is run on GLPK
(``cobra.Configuration().solver = "glpk"`` before loading; no Gurobi session); (b) our calls = genes failing >= 1
task in ``results/w3/tasks/essential_<parent>.parquet``. Checks: same gene set; same per-gene category sets;
equality with the reviewed W1 Gurobi P57 calls (``results/w1/p57_<parent>.parquet``); GR failures = the FBA
lethal set (ratio_r < 0.01). Label-free.

Usage: .venv/bin/python scripts/w3_task_equality.py --workers 7
Output: results/w3/tasks/p57_equality.tsv
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w3_common import OUT_DIR, PROJECT, atomic_write_tsv, current, run_main  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--models", nargs="*", default=["human1", "human2"])
    a = ap.parse_args(argv)
    import cobra

    cobra.Configuration().solver = "glpk"
    sys.path.insert(0, str(PROJECT / "data" / "Human-GEM-main" / "code" / "test"))
    from taskEssentialGenes import find_task_essential_categories  # Human-GEM code (reuse, DESIGN §2)

    from nextgem.model import GEMS, load_model

    rows = []
    for name in a.models:
        t0 = time.time()
        m = load_model(name)
        ref = find_task_essential_categories(m, GEMS[name].essential_tasks, processes=a.workers)
        del m
        cap = pd.read_parquet(OUT_DIR / "tasks" / f"essential_{name}.parquet")
        meta = pd.read_csv(OUT_DIR / "tasks" / f"essential_{name}_tasks.tsv", sep="\t")
        cat = dict(zip(meta.task, meta.category))
        fail = cap[~cap.feasible]
        mine = fail.groupby("gene").task.apply(lambda s: {cat[t] for t in s}).to_dict()
        a_, b_ = set(mine), set(ref)
        cat_eq = sum(mine[g] == ref[g] for g in a_ & b_)
        w1 = pd.read_parquet(PROJECT / "results" / "w1" / f"p57_{name}.parquet").set_index("gene")
        w1set = set(w1.index[w1.p57])
        gr_key = meta.task[meta.category == "GR"].iloc[0]
        gr_fail = set(cap.gene[(cap.task == gr_key) & ~cap.feasible])
        base = pd.read_parquet(PROJECT / "results" / "baseline" / f"{name}_ko.parquet", columns=["gene", "ratio_r"])
        lethal = set(base.gene[base.ratio_r < 0.01])
        rows += [
            {"model": name, "check": "p57_equals_taskEssentialGenes_glpk", "pass": a_ == b_,
             "detail": f"ours {len(a_)} vs Human-GEM {len(b_)}; only ours {sorted(a_ - b_)[:5]}; only theirs {sorted(b_ - a_)[:5]}"},
            {"model": name, "check": "category_sets_equal", "pass": cat_eq == len(a_ & b_),
             "detail": f"{cat_eq}/{len(a_ & b_)}"},
            {"model": name, "check": "p57_equals_w1_gurobi", "pass": a_ == w1set,
             "detail": f"ours {len(a_)} vs W1 {len(w1set)}; sym diff {sorted(a_ ^ w1set)[:5]}"},
            {"model": name, "check": "GR_failures_equal_fba_lethal", "pass": gr_fail == lethal,
             "detail": f"GR fail {len(gr_fail)} vs ratio_r<0.01 {len(lethal)}; sym diff {sorted(gr_fail ^ lethal)[:5]}"},
        ]
        current().input(GEMS[name].sbml, GEMS[name].essential_tasks, OUT_DIR / "tasks" / f"essential_{name}.parquet")
        print(f"[{name}] {time.time() - t0:.0f} s", rows[-4:], flush=True)
    out = atomic_write_tsv(pd.DataFrame(rows), OUT_DIR / "tasks" / "p57_equality.tsv")
    current().output(out)
    return 0 if all(r["pass"] for r in rows) else 1


if __name__ == "__main__":
    w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 7
    raise SystemExit(run_main(main, solver="glpk", workers=w))
