"""Label-free pilot: gene-level E-Flux (candidate F6-GE) on 3 Stage-A panel lines, Human2, M0, GLPK.

Formulation (see notes/methods/GPRTransform_assessment.md, variant B):
  * parent: Human2 on Ham's medium M0 (biomass objective), GPR-transformed with the exact "arm" encoding;
  * gene linear value lin_g = 2**s_g if s_g >= 0 else 0, with s_g the R1 RNA signed score
    (log2(TPM+1) - 1), i.e. exactly the F4 inputs; genes without expression: no bound;
  * gene-usage bound u_g <= max(FLOOR, SCALE * lin_g / max_g lin_g) with FLOOR = 10, SCALE = 1000
    (F4's floor and scale, moved from reactions to genes); reaction, exchange and arm bounds unchanged;
  * KO of gene g: u_g = 0; ratio = KO growth / line WT growth; ratio_r = rounded to 6 dp, clipped [0, 1];
  * exact zero-usage shortcut: genes with u_g = 0 in the line's minimal-total-usage WT solution
    (gene-pFBA at WT growth) keep ratio 1 without a solve.
Compared with the reaction-level F4 builds of the same lines (results/w4/tmp/pilot_ko_gurobi_lp.parquet,
R1/rna/real/no genomic constraint). No label is read.

Usage: .venv/bin/python results/gprt/proto/pilot_gene_eflux.py [--lines ACH-... ...]
Writes results/gprt/tmp/pilot_gene_eflux_ko.parquet (per-gene, git-ignored) and
results/gprt/pilot_gene_eflux_summary.json (aggregates only).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
ROOT = HERE.parents[2]
FLOOR, SCALE = 10.0, 1000.0
LINES = ["ACH-000957", "ACH-002024", "ACH-003383"]  # first 3 panel lines = the lines of the W4 F4 pilot


def line_job(args):
    line, s_g = args
    import cobra

    cobra.Configuration().solver = "glpk"
    sys.path.insert(0, str(HERE))
    import gpr_transform as gt
    from optlang.interface import INFEASIBLE, OPTIMAL

    from nextgem.knockouts import configure
    from nextgem.model import BIOMASS_REACTION, load_model

    c_load = time.process_time()
    m = load_model("human2")
    configure(m)
    tm, info = gt.transform(m, split="arm")
    del m
    tm.solver.configuration.timeout = 300
    load_cpu = time.process_time() - c_load

    c0 = time.process_time()
    lin = {g: (2.0 ** v if v >= 0 else 0.0) for g, v in s_g.items() if not math.isnan(v)}
    cap_max = max(lin.values())
    bounds = {g: max(FLOOR, SCALE * v / cap_max) for g, v in lin.items()}
    rows = []
    with tm:
        for g, b in bounds.items():
            tm.reactions.get_by_id(gt.usage(g)).bounds = (0.0, b)
        tm.slim_optimize()
        assert tm.solver.status == OPTIMAL, tm.solver.status
        wt = tm.solver.objective.value
        wt_cpu = time.process_time() - c0
        # gene-pFBA at WT growth: minimal total usage (unweighted), for the exact zero-usage shortcut
        c1 = time.process_time()
        with tm:
            bm = tm.reactions.get_by_id(BIOMASS_REACTION)
            bm.lower_bound = wt * (1 - 1e-9)
            tm.objective = {tm.reactions.get_by_id(gt.usage(g)): 1.0 for g in info["genes"]}
            tm.objective_direction = "min"
            tm.slim_optimize()
            assert tm.solver.status == OPTIMAL, tm.solver.status
            u = gt.usage_fluxes(tm, info)
            total_usage = tm.solver.objective.value
        pfba_cpu = time.process_time() - c1
        tm.slim_optimize()  # restore a WT basis for warm starts
        c2 = time.process_time()
        n_solved = 0
        for g in info["genes"]:
            if u[g] <= 1e-9:
                rows.append((g, wt, "zero_usage", u[g], bounds.get(g, math.inf)))
                continue
            with tm:
                gt.knock_out(tm, g)
                tm.slim_optimize()
                st = tm.solver.status
                v = tm.solver.objective.value if st == OPTIMAL else (0.0 if st == INFEASIBLE else math.nan)
            n_solved += 1
            rows.append((g, v, st, u[g], bounds.get(g, math.inf)))
        sweep_cpu = time.process_time() - c2
    n_bound_active = sum(1 for g, b in bounds.items() if u[g] >= b * (1 - 1e-6))
    meta = {"line": line, "wt_growth": wt, "cap_max": cap_max, "n_genes_bounded": len(bounds),
            "n_genes_floor": sum(1 for b in bounds.values() if b == FLOOR), "total_usage_wt": total_usage,
            "n_genes_used_wt": sum(1 for x in u.values() if x > 1e-9), "n_usage_at_bound_wt": n_bound_active,
            "n_ko_solved": n_solved, "load_transform_cpu_s": load_cpu, "wt_cpu_s": wt_cpu, "pfba_cpu_s": pfba_cpu,
            "sweep_cpu_s": sweep_cpu, "line_cpu_s": wt_cpu + pfba_cpu + sweep_cpu,
            "isozyme_genes": sorted({g for t in info["terms"].values() if len(t) > 1 for tt in t for g in tt}),
            "complex_genes": sorted({g for t in info["terms"].values() for tt in t if len(tt) > 1 for g in tt})}
    return line, rows, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lines", nargs="+", default=LINES)
    a = ap.parse_args()
    import cobra

    cobra.Configuration().solver = "glpk"
    import numpy as np
    import pandas as pd

    from nextgem import context_inputs as ci
    from nextgem.knockouts import rounded_ratio
    from nextgem.model import load_model

    panel = pd.read_csv(ROOT / "results" / "w4" / "stage_a_panel.tsv", sep="\t")
    assert set(a.lines) <= set(panel.ModelID), "pilot lines must be Stage-A (training) panel lines"
    t0 = time.time()
    genes = ci.gene_table(load_model("human2"))
    rna = ci.read_expression(genes, a.lines)
    st = json.loads((ROOT / "results" / "w4" / "threshold_stats.json").read_text())
    gm = pd.read_parquet(ROOT / "results" / "w4" / "threshold_gene_means.parquet").set_index("gene")
    stats = ci.ThresholdStats("rna", st["rna"]["p25"], st["rna"]["p75"], gm.rna_mean.dropna(), st["rna"]["floor_value"],
                              st["rna"]["n_lines"])
    scores = {ln: ci.signed_scores(ln, "R1", "rna", genes, rna, stats) for ln in a.lines}
    inputs_wall = time.time() - t0
    with ProcessPoolExecutor(len(a.lines), mp_context=get_context("spawn")) as pool:
        res = list(pool.map(line_job, [(ln, scores[ln]) for ln in a.lines]))
    ref = pd.read_parquet(ROOT / "results" / "w4" / "tmp" / "pilot_ko_gurobi_lp.parquet")
    base = pd.read_parquet(ROOT / "results" / "baseline" / "human2_ko.parquet")[["gene", "symbol", "ratio_r"]]
    base = base.rename(columns={"ratio_r": "ratio_r_generic"})
    kos, summ = [], {"formulation": {"floor": FLOOR, "scale": SCALE, "lin": "2**s_g (s_g>=0) else 0; R1 RNA",
                                     "encoding": "arm", "solver": "glpk", "parent": "human2", "medium": "M0"},
                     "inputs_wall_s": inputs_wall, "wall_s": time.time() - t0, "lines": {}}
    for line, rows, meta in res:
        ko = pd.DataFrame(rows, columns=["gene", "ko_growth", "status", "u_wt", "u_bound"])
        ko["line"] = line
        ko["ratio_r"] = rounded_ratio(ko.ko_growth / meta["wt_growth"])
        f4 = ref[ref.build_id == f"{line}|F4|R1|rna|real|none"][["gene", "ratio_r"]].rename(columns={"ratio_r": "ratio_r_F4"})
        ko = ko.merge(f4, on="gene", how="left").merge(base, on="gene", how="left")
        iso, cplx = set(meta.pop("isozyme_genes")), set(meta.pop("complex_genes"))
        ko["isozyme_gene"] = ko.gene.isin(iso)
        ko["complex_gene"] = ko.gene.isin(cplx)
        kos.append(ko)
        r, r4 = ko.ratio_r, ko.ratio_r_F4
        diff = (r - r4).abs() > 1e-6
        def dist(x):
            return {"n": int(x.notna().sum()), "distinct": int(x.nunique()), "lethal_lt_0.01": int((x < 0.01).sum()),
                    "partial_0.01_1": int(((x >= 0.01) & (x < 1)).sum()), "eq_1": int((x == 1).sum()),
                    "quantiles_partial": [round(float(q), 4) for q in np.nanquantile(x[(x >= 0.01) & (x < 1)], [0, .25, .5, .75, 1])]
                    if ((x >= 0.01) & (x < 1)).any() else [],
                    "largest_tie_group": int(x.value_counts().iloc[0])}
        iso_changed = ko[ko.isozyme_gene & diff]
        iso_graded = ko[ko.isozyme_gene & (r < 1) & (r >= 0.01) & (r4 == 1)]
        meta.update({
            "gene_eflux": dist(r), "reaction_F4": dist(r4), "generic": dist(ko.ratio_r_generic),
            "n_genes_score_differs_vs_F4": int(diff.sum()),
            "lethal_only_gene_level": int(((r < 0.01) & (r4 >= 0.01)).sum()),
            "lethal_only_F4": int(((r >= 0.01) & (r4 < 0.01)).sum()),
            "spearman_vs_F4": float(r.corr(r4, method="spearman")),
            "isozyme_genes_in_model": len(iso), "isozyme_genes_score_differs_vs_F4": len(iso_changed),
            "isozyme_genes_partial_here_but_1_in_F4": len(iso_graded),
            "examples_isozyme_partial_here_but_1_in_F4": iso_graded.sort_values("ratio_r").symbol.head(15).tolist(),
            "isozyme_genes_lower_here": int((ko.isozyme_gene & (r < r4 - 1e-6)).sum()),
            "isozyme_genes_higher_here": int((ko.isozyme_gene & (r > r4 + 1e-6)).sum()),
            "status_counts": ko.status.value_counts().to_dict(),
        })
        summ["lines"][line] = meta
    allko = pd.concat(kos)
    (OUT / "tmp").mkdir(exist_ok=True)
    allko.to_parquet(OUT / "tmp" / "pilot_gene_eflux_ko.parquet", index=False)
    summ["wall_s"] = time.time() - t0
    summ["cpu_s_total_lines"] = sum(v["line_cpu_s"] + v["load_transform_cpu_s"] for v in summ["lines"].values())
    (OUT / "pilot_gene_eflux_summary.json").write_text(json.dumps(summ, indent=1, default=str))
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
