"""W4 review, part 4: consistency of the author's FASTCORE (F3) output through cobra fastcc (GLPK).

For each line: the FASTCORE support (before repair) and the final open set (after repair) are checked: the
reactions outside the set are closed on the M0 parent and cobra fastcc (flux_threshold 1, the author's eps)
returns the consistent subset; blocked = set minus consistent.
Usage: .venv/bin/python results/w4/review/r4_fastcc.py --line ACH-000957 --which support
"""

import argparse
import json
import time

import cobra

cobra.Configuration().solver = "glpk"
from cobra.flux_analysis import fastcc  # noqa: E402

from nextgem import context as cx  # noqa: E402
from nextgem import paths  # noqa: E402
from nextgem.model import load_model  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--line", required=True)
ap.add_argument("--which", choices=["support", "final_open"], required=True)
ap.add_argument("--solver-tag", default="gurobi")
a = ap.parse_args()
OUT = paths.results_dir() / "w4" / "review"
sets = json.loads((OUT / "tmp" / f"f3_sets_{a.line}_{a.solver_tag}.json").read_text())
S = set(sets[a.which])
t0, c0 = time.time(), time.process_time()
m = load_model("human2")
parent = cx.make_parent(m, "human2")
m.objective = cx.BIOMASS_REACTION  # cobra-level objective (fastcc replaces it anyway)
with m:
    for r in m.reactions:
        if r.id not in S:
            r.bounds = (0.0, 0.0)
    cm = fastcc(m, flux_threshold=1.0)
cons = {r.id for r in cm.reactions}
blocked = sorted(S - cons)
res = {"line": a.line, "which": a.which, "n_set": len(S), "n_consistent": len(cons & S), "n_blocked": len(blocked),
       "blocked_sample": blocked[:25], "core_blocked": len(set(sets["core"]) & set(blocked)),
       "wall_s": round(time.time() - t0), "cpu_s": round(time.process_time() - c0)}
(OUT / f"r4_fastcc_{a.line}_{a.which}.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res))
