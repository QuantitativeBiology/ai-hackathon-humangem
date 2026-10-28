"""W8: is a random-KO gate miss LP noise? Re-solve one model's sampled KOs with GLPK (default and tight tolerances),
GLPK exact (rational) for the mismatching genes, and Gurobi (tight tolerances). Label-free. Diagnostic only: the
gate and its tolerances (REL_TOL 1e-6, ABS_FLOOR 1e-9; RESEARCH_TREE 9.1) are NOT changed by this script.

The two models are the ones the gate compared: ``source`` = the in-memory P1 model rebuilt exactly as the build job
does (parent copy + M2(f=1) bounds + B-loo variant) and ``reloaded`` = the released SBML file read back. The gene
sample is the gate's own (``resource.random_ko_compare``: 20 genes, SEEDS[0]).

Usage: .venv/bin/python scripts/w8_ko_noise_check.py --parent human2 --w2-sweep-id human2:20cc933a4b852b6a
"""

from __future__ import annotations

import argparse
import json
import math
import os
import signal
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as mc  # noqa: E402
import w8_build_p1 as wb  # noqa: E402

OUT = wb.W8 / "ko_noise"
EXACT_CPU_S = 900


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def solve(model, gene: str | None) -> tuple[float, str, float]:
    """Objective after the KO (no rounding, no retry) and the biomass bound violation of the solution."""
    with model:
        if gene:
            model.genes.get_by_id(gene).knock_out()
        model.slim_optimize()
        st = model.solver.status
        v = model.solver.objective.value if st == "optimal" else math.nan
        bio = model.reactions.get_by_id("MAR13082")
        viol = max(0.0, bio.lower_bound - bio.flux) if st == "optimal" else math.nan
    return v, st, viol


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", required=True)
    ap.add_argument("--w2-sweep-id", required=True)
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--no-gurobi", action="store_true")
    a = ap.parse_args(argv)
    import cobra

    cobra.Configuration().solver = "glpk"
    import m6_sweeps

    from nextgem import resource as rs, w3
    from nextgem.model import GEMS, load_model
    from nextgem.run import provenance

    models = pd.read_csv(wb.PROJECT / "resource" / "models.tsv", sep="\t", dtype=str, keep_default_na=False)
    row = models[(models.parent == a.parent) & (models.w2_sweep_id == a.w2_sweep_id)].iloc[0]
    path = wb.PROJECT / "resource" / row.file
    inputs = [path, wb.PROJECT / "resource" / "models.tsv", GEMS[a.parent].sbml,
              wb.PROJECT / "results" / "w3" / "structure" / f"variants_{a.parent}.json"]
    with provenance(wb.MANIFEST, inputs=inputs, seeds=[], solver="glpk;glpk_exact;gurobi", workers=1,
                    allow_dirty=a.allow_dirty, author=wb.AUTHOR) as prov:
        t0 = time.time()
        bounds, _ = m6_sweeps.bounds_by_id({a.parent: {mc.M2: [a.w2_sweep_id]}})
        parent = load_model(a.parent)
        src = rs.build_medium_variant_model(parent, bounds[a.w2_sweep_id],
                                            w3.Variant.from_json(m6_sweeps.variant_json(a.parent, mc.BLOO)))
        rel = rs.load_checked(path)
        genes = rs.random_ko_compare(src, src, n=rs.N_RANDOM_KO)["genes"]  # the gate's sample (SEEDS[0])
        rows = []

        def run_all(tag: str, m, configure=None, genes_=None):
            if configure:
                configure(m)
            for g in [None] + list(genes_ if genes_ is not None else genes):
                c0 = time.process_time()
                v, st, viol = solve(m, g)
                rows.append(dict(solver=tag, model=m.notes.get("_w8_tag", ""), gene=g or "WT", value=v, status=st,
                                 biomass_lb_violation=viol, cpu_s=round(time.process_time() - c0, 2)))

        def glpk_tight(m):
            c = m.solver.configuration
            c.tolerances.feasibility = 1e-10
            c.presolve = False

        for tag, m in (("source", src), ("reloaded", rel)):
            m.notes["_w8_tag"] = tag
            run_all("glpk_default", m)
            run_all("glpk_feastol1e-10", m, glpk_tight)
        df = pd.DataFrame(rows)
        # gate rule on the default solves, per gene: which genes mismatch?
        d = df[df.solver == "glpk_default"].pivot(index="gene", columns="model", values="value")
        bad = [g for g in d.index if g != "WT" and not rs._rel_equal(d.loc[g, "source"], d.loc[g, "reloaded"])]
        recorded = []  # the gate's own mismatches in the build (results/w8/p1/gates.jsonl or unscreened/gates.jsonl)
        for gl in (wb.W8 / "p1" / "gates.jsonl", wb.W8 / "unscreened" / "gates.jsonl"):
            if gl.exists():
                for ln in gl.read_text().splitlines():
                    r_ = json.loads(ln)
                    if r_["parent"] == a.parent and r_["w2_sweep_id"] == a.w2_sweep_id and r_["gates"]:
                        recorded += [m_[0] for m_ in r_["gates"]["random_ko_match"][1]["mismatches"]]
        rerun_bad = list(bad)
        bad = sorted(set(bad) | set(recorded))
        # exact rational GLPK for the mismatching genes (and WT), both models
        for tag, m in (("source", src), ("reloaded", rel)):
            mm = m.copy()
            mm.notes["_w8_tag"] = tag
            try:
                signal.signal(signal.SIGPROF, _alarm)
                signal.setitimer(signal.ITIMER_PROF, EXACT_CPU_S)
                mm.solver = "glpk_exact"
                run_all("glpk_exact", mm, genes_=bad)
            except _Timeout:
                rows.append(dict(solver="glpk_exact", model=tag, gene="(timeout)", value=math.nan, status="timeout",
                                 biomass_lb_violation=math.nan, cpu_s=EXACT_CPU_S))
            finally:
                signal.setitimer(signal.ITIMER_PROF, 0.0)
        if not a.no_gurobi:
            from nextgem.w2_runtime import ensure_gurobi_slot

            ensure_gurobi_slot()

            def gurobi_tight(m):
                p = m.solver.problem.Params
                p.Threads, p.FeasibilityTol, p.OptimalityTol = 1, 1e-9, 1e-9

            for tag, m in (("source", src), ("reloaded", rel)):
                mm = m.copy()
                mm.notes["_w8_tag"] = tag
                mm.solver.configuration.tolerances.feasibility = 1e-9  # Gurobi's minimum; carried over on conversion
                mm.solver = "gurobi"
                run_all("gurobi_tol1e-9", mm, gurobi_tight)
        df = pd.DataFrame(rows)
        OUT.mkdir(parents=True, exist_ok=True)
        tsv = mc.atomic_write_tsv(df, OUT / f"{a.w2_sweep_id.replace(':', '_')}_resolves.tsv", float_format="%.12g")
        wt = df[df.gene == "WT"].set_index(["solver", "model"]).value
        summ = {"model": row.file, "genes_sampled": genes, "gate_mismatch_genes_build": recorded,
                "gate_mismatch_genes_rerun_glpk_default": rerun_bad,
                "per_solver": {}, "wall_s": round(time.time() - t0, 1)}
        for s_, sub in df[df.gene.isin(bad)].groupby("solver"):
            summ["per_solver"][s_] = {f"{r.model}:{r.gene}": {"value": r.value, "status": r.status,
                                                               "biomass_lb_violation": r.biomass_lb_violation,
                                                               "ratio": r.value / wt.get((s_, r.model), math.nan)}
                                      for r in sub.itertuples()}
        # every sampled gene: max |source - reloaded| per solver
        for s_, sub in df[df.gene != "WT"].groupby("solver"):
            p_ = sub.pivot_table(index="gene", columns="model", values="value")
            if {"source", "reloaded"} <= set(p_.columns):
                summ["per_solver"].setdefault(s_, {})["max_abs_diff_source_reloaded_all_sampled"] = float(
                    np.nanmax(np.abs(p_.source - p_.reloaded)))
        js = mc.atomic_write_text(json.dumps(summ, indent=1, default=str) + "\n",
                                  OUT / f"{a.w2_sweep_id.replace(':', '_')}_summary.json")
        prov.output(tsv, js)
        print(json.dumps(summ, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
