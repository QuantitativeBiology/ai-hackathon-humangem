"""W4 review, part 7: is the E-Flux (F4) KO profile line-specific? Generic Human2/M0 KO (author's sweep, GLPK)
vs the rebuilt F2/F3/F4 KO tables of 3 panel lines (label-free)."""
import json
import cobra
cobra.Configuration().solver = "glpk"
import pandas as pd
from nextgem import context as cx, paths
from nextgem.model import load_model
T = paths.results_dir() / "w4" / "review" / "tmp"
m = load_model("human2"); p = cx.make_parent(m, "human2")
gen = cx.ko_sweep(p, cx.LineModel(line="generic", family="F0", closed=set())).set_index("gene").ratio_r
ge = set(gen.index[gen < 0.01]); res = {"generic_n_lethal": len(ge)}
for f in ("F2", "F3", "F4"):
    for line in ("ACH-000957", "ACH-002024", "ACH-003383"):
        sol = "gurobi" if f == "F3" else "glpk"
        k = pd.read_parquet(T / f"ko_{line}_{f}_{sol}.parquet").set_index("gene").ratio_r
        le = set(k.index[k < 0.01])
        res[f"{f}_{line}"] = {"n_lethal": len(le), "equal_generic_lethal": le == ge, "only_line": len(le - ge), "only_generic": len(ge - le),
                              "n_partial_0.01_0.99": int(((k >= 0.01) & (k < 0.99)).sum()), "n_partial_generic": int(((gen >= 0.01) & (gen < 0.99)).sum())}
(T.parent / "r7_generic_vs_f4.json").write_text(json.dumps(res, indent=1, default=bool)); print(json.dumps(res, indent=1, default=bool))
