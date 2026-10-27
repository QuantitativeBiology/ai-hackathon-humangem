"""W4 review, part 6: why are some consistent core reactions missing from the FASTCORE support? FVA (GLPK) of the
missing core reactions on the M0 parent: max |v| reachable vs the FASTCORE threshold 0.99 * eps (eps = 1)."""
import json, sys
import cobra
cobra.Configuration().solver = "glpk"
from cobra.flux_analysis import flux_variability_analysis
from nextgem import context as cx, paths
from nextgem.model import load_model
OUT = paths.results_dir() / "w4" / "review"
m = load_model("human2"); parent = cx.make_parent(m, "human2"); m.objective = cx.BIOMASS_REACTION
res = {}
for line in sys.argv[1:]:
    s = json.loads((OUT / "tmp" / f"f3_sets_{line}_gurobi.json").read_text())
    miss = sorted(set(s["core"]) - set(s["support"]))
    fva = flux_variability_analysis(m, miss, fraction_of_optimum=0.0, processes=1)
    mx = fva.abs().max(axis=1)
    res[line] = {"n_missing": len(miss), "n_max_abs_flux_below_0.99": int((mx < 0.99).sum()),
                 "max_abs_flux": {k: round(float(v), 6) for k, v in mx.items()},
                 "subsystems": sorted({m.reactions.get_by_id(r).subsystem for r in miss})}
(OUT / "r6_core_fva.json").write_text(json.dumps(res, indent=1)); print(json.dumps({k: {kk: v[kk] for kk in ("n_missing", "n_max_abs_flux_below_0.99", "subsystems")} for k, v in res.items()}, indent=1))
