"""W4 re-check (B2): the per-LP time limit is installed on both solvers and a GLPK time-limit hit is
retried once and then returned as non-optimal (the build then counts as failed / the KO solve as NaN)."""
import json, sys, time
import cobra
solver = sys.argv[1]
slot = None
if solver == "gurobi":
    from nextgem import gurobi_pool
    slot = gurobi_pool.hold_slot(timeout=3600, poll=10)
cobra.Configuration().solver = solver
from nextgem import context as cx, paths
from nextgem.model import load_model
m = load_model("human2"); p = cx.make_parent(m, "human2")
res = {"solver": cx.solver_name(m)}
if solver == "gurobi":
    res["TimeLimit"] = m.solver.problem.Params.TimeLimit; res["Threads"] = m.solver.problem.Params.Threads
else:
    res["optlang_timeout"] = m.solver.configuration.timeout
    with m:
        m.reactions.get_by_id(cx.BIOMASS_REACTION).lower_bound = 0.5 * p.wt_generic
        cx._l1_objective(m, {rid: 1.0 for rid in p.rxn_ids})
        m.solver.configuration.timeout = 1  # force a hit on the full-model L1 LP
        t = time.time(); st = cx._status_retry(m); res["forced_1s_status_after_retry"] = st; res["wall_s"] = round(time.time() - t, 1)
        m.solver.configuration.timeout = 300
        cx._biomass_objective(m)
    res["restored_wt_ok"] = abs(m.slim_optimize() - p.wt_generic) < 1e-6
(paths.results_dir() / "w4" / "review" / f"r11_lp_limit_{solver}.json").write_text(json.dumps(res, indent=1)); print(res)
