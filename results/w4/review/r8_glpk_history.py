"""W4 review, part 8: GLPK history dependence of the F0 repair (ACH-002024, R1, real RNA). One fresh process:
extract once, then finish + sweep three times in a row; compare re-opened sets, open masks and KO tables."""
import json, sys
import cobra
cobra.Configuration().solver = "glpk"
import pandas as pd
from nextgem import context as cx, context_inputs as ci, paths
from nextgem.model import load_model
T = paths.results_dir() / "w4" / "review"
line = sys.argv[1] if len(sys.argv) > 1 else "ACH-002024"
m = load_model("human2"); p = cx.make_parent(m, "human2"); g = ci.gene_table(m)
expr = pd.read_parquet(T / "tmp" / "panel_expr.parquet")
s = ci.signed_scores(line, "R1", "rna", g, expr, ci.ThresholdStats("rna", 0, 0, pd.Series(dtype=float), 0, 0)).reindex(p.genes)
base = cx.extract(p, "F0", s)
runs = []
for k in range(2):
    lm = cx.LineModel(line=line, family="F0", closed=set(base.closed), off=set(base.off), info=dict(base.info))
    import time; t = time.process_time()
    cx.finish(p, lm); t1 = time.process_time(); print(f"run {k}: finish {t1 - t:.0f}s reopened {len(lm.reopened)} tasks {lm.info['repair_task_lps']}", flush=True)
    ko = cx.ko_sweep(p, lm).set_index("gene").ratio_r; print(f"run {k}: sweep {time.process_time() - t1:.0f}s", flush=True)
    runs.append((lm.reopened, ko, lm.info["wt_growth"]))
res = {"line": line, "n_reopened": [len(r[0]) for r in runs], "wt": [r[2] for r in runs],
       "reopened_equal_1_2": runs[0][0] == runs[1][0],
       "symdiff_1_2": len(runs[0][0] ^ runs[1][0]),
       "ko_n_differ_1_2": int((runs[0][1] - runs[1][1]).abs().gt(1e-9).sum()),
       "lethal_differ_1_2": int(((runs[0][1] < 0.01) != (runs[1][1] < 0.01)).sum()),
       "ko_max_abs_diff_1_2": float((runs[0][1] - runs[1][1]).abs().max())}
(T / f"r8_glpk_history_{line}.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res))
