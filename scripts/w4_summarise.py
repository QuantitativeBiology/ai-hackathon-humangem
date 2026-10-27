"""Summarise the W4 / M3 pilots into results/w4/pilot_summary.{json,md} (every number in README comes from here).

Inputs (results/w4/): repair_probe12.tsv, loss_rule_panel.tsv, pilot_builds_gurobi_lp.tsv,
pilot_builds_gurobi_corda.tsv, glpk_vs_gurobi.tsv, pilot_builds_glpk_<line>.tsv, fastcc_human2_M0.json,
stage_a_panel.tsv. Label-free.
"""

import argparse
import json

import numpy as np
import pandas as pd
import psutil
from scipy.stats import spearmanr

from nextgem import paths

W4 = paths.results_dir() / "w4"
LP = ["F0", "F2", "F3", "F4"]
GB_PER_SLOT_RESERVE = 20.0  # §11.4: keep >= 20 GB free
CORDA_BANDS = (0.25, 1.0)  # CPU-h per model (DESIGN row A6)


def r(x, k=2):
    return None if x is None or pd.isna(x) else round(float(x), k)


def _bits(hexmask: str) -> np.ndarray:
    return np.unpackbits(np.frombuffer(bytes.fromhex(hexmask), dtype=np.uint8)).astype(bool)


def a6prime_comparison(stage_dir) -> pd.DataFrame:
    """Amendment A7, label-free: every CORDA (A6') build against F0/F2/F3/F4 built in A1 on the same line and
    threshold (real RNA, no genomic constraint): reaction-set Jaccard, essential-set (ratio_r < 0.01) Jaccard and
    sizes, Spearman of ratio_r over all genes, WT growth ratio, repair size and build CPU. Descriptive only."""
    from nextgem.context_inputs import A6P_TAG

    b = pd.read_csv(stage_dir / "builds.tsv", sep="\t")
    b = b[b.genomic == "none"].set_index("build_id")
    ko = pd.concat([pd.read_parquet(stage_dir / f"ko_{r}.parquet") for r in ("A1", "A6p")
                    if (stage_dir / f"ko_{r}.parquet").exists()])
    ko["build_id"] = ko.build_id.astype(str)
    out = []
    for bid, c in b[b.family == "F3b"].iterrows():
        if c.get("failed", False) or pd.isna(c.get("open_mask")):
            out.append({"line": c.line, "rule": c.rule, "family": "F3b", "corda_failed": True, "status": A6P_TAG})
            continue
        kc = ko[ko.build_id == bid].set_index("gene").ratio_r
        for f in LP:
            oid = f"A1|{c.line}|{f}|{c.rule}|rna|real|none"
            if oid not in b.index:
                continue
            o = b.loc[oid]
            ma, mb = _bits(c.open_mask), _bits(o.open_mask)
            ko_o = ko[ko.build_id == oid].set_index("gene").ratio_r.reindex(kc.index)
            ea, eb = set(kc.index[kc < 0.01]), set(ko_o.index[ko_o < 0.01])
            out.append({
                "line": c.line, "rule": c.rule, "vs_family": f,
                "reaction_jaccard": r((ma & mb).sum() / max(1, (ma | mb).sum()), 4),
                "essential_jaccard": r(len(ea & eb) / max(1, len(ea | eb)), 4), "n_essential_corda": len(ea),
                "n_essential_other": len(eb), "ratio_r_spearman": r(spearmanr(kc, ko_o).statistic, 4),
                "wt_ratio_corda": r(c.wt_ratio_generic, 3), "wt_ratio_other": r(o.wt_ratio_generic, 3),
                "n_reopened_corda": int(c.n_reopened), "n_reopened_other": int(o.n_reopened),
                "build_cpu_min_corda": r((c.extract_cpu_s + c.variant_cpu_s) / 60, 1),
                "build_cpu_min_other": r((o.extract_cpu_s + o.variant_cpu_s) / 60, 2),
                "n_open_corda": int(c.n_open), "n_open_other": int(o.n_open), "status": A6P_TAG})
    df = pd.DataFrame(out)
    df.to_csv(stage_dir / "a6prime_comparison.tsv", sep="\t", index=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage-a", default=str(W4 / "stage_a"), help="Stage-A output dir for the A6' comparison")
    a = ap.parse_args()
    from pathlib import Path

    stage_dir = Path(a.stage_a)
    S = {}
    rp = pd.read_csv(W4 / "repair_probe12.tsv", sep="\t")
    S["repair12"] = {
        "n_lines": len(rp), "n_wt_zero_before": int((rp.wt_before_repair <= 1e-9).sum()),
        "n_restored": int((~rp.failed).sum()), "wt_ratio_after_min": r(rp.wt_ratio_generic.min(), 3),
        "wt_ratio_after_max": r(rp.wt_ratio_generic.max(), 3), "reopened_min": int(rp.n_reopened.min()),
        "reopened_max": int(rp.n_reopened.max()), "reopen_frac_closed_min": r(rp.reopen_frac_of_closed.min(), 4),
        "reopen_frac_closed_max": r(rp.reopen_frac_of_closed.max(), 4), "repair_cpu_s_median": r(rp.repair_cpu_s.median(), 1),
        "tiers": rp.tier.value_counts().to_dict()}
    lr = pd.read_csv(W4 / "loss_rule_panel.tsv", sep="\t")
    S["loss_rule"] = {
        "n_lines": len(lr), "all_order_independent": bool(lr.order_independent.all()), "orders_per_line": int(lr.n_orders.iloc[0]),
        "lost_model_genes_median": r(lr.n_lost_model_genes.median(), 1), "lost_model_genes_max": int(lr.n_lost_model_genes.max()),
        "n_lines_contradicted": int((lr.n_contradicted > 0).sum()), "contradicted_genes": sorted(set(";".join(lr.contradicted.fillna("")).split(";")) - {""}),
        "n_lines_joint_drop": int((lr.n_dropped_joint > 0).sum()), "cpu_s_per_line_max": r(lr.cpu_s_per_run.max(), 2)}
    fc = json.loads((W4 / "fastcc_human2_M0.json").read_text())
    S["fastcc"] = {"consistent": len(fc["consistent"]), "reactions": fc["n_reactions"], "cpu_min": r(fc["cpu_s"] / 60, 1)}

    d = pd.read_csv(W4 / "pilot_builds_gurobi_lp.tsv", sep="\t")
    d["build_cpu_s"] = d.extract_cpu_s + d.variant_cpu_s
    fam = {}
    for f in LP:
        x = d[d.family == f]
        b, g = x[x.genomic == "none"], x[x.genomic == "loss"]
        fam[f] = {
            "build_cpu_min": r(b.build_cpu_s.mean() / 60), "build_cpu_min_max": r(b.build_cpu_s.max() / 60),
            "resweep_cpu_min": r(g.variant_cpu_s.mean() / 60), "extract_cpu_s": r(b.extract_cpu_s.mean(), 1),
            "repair_cpu_s": r(x.repair_cpu_s.mean(), 1), "sweep_cpu_s": r(x.sweep_cpu_s.mean(), 1),
            "tasks_cpu_s": r(x.tasks_cpu_s.mean(), 1), "sweep_n_solved": r(x.sweep_n_solved.mean(), 0),
            "n_reopened_range": [int(x.n_reopened.min()), int(x.n_reopened.max())],
            "repair_frac_max": r(x.repair_frac_reactions.max(), 4), "open_reactions_range": [int(x.n_open.min()), int(x.n_open.max())],
            "wt_ratio_generic_range": [r(x.wt_ratio_generic.min(), 3), r(x.wt_ratio_generic.max(), 3)],
            "lethal_range": [int(x.sweep_n_lethal.min()), int(x.sweep_n_lethal.max())],
            "tasks_passed_range": [int(x.tasks_passed.min()), int(x.tasks_passed.max())],
            "failed": int(x.failed.sum()), "n_builds": len(x), "peak_rss_gb_max": r(x.peak_rss_gb.max(), 2),
            "perm_built": int((x.variant == "perm").sum()), "perm_failed": int(x[x.variant == "perm"].failed.sum()),
            "n_lines": x.line.nunique()}
    S["families"] = fam

    cp = W4 / "pilot_builds_gurobi_corda.tsv"
    if cp.exists():
        c = pd.read_csv(cp, sep="\t")
        c["model_cpu_h"] = (c.extract_cpu_s + c.variant_cpu_s) / 3600
        h = float(c.model_cpu_h.max())
        band = "full" if h <= CORDA_BANDS[0] else ("fallback" if h <= CORDA_BANDS[1] else "drop")
        S["corda"] = {"n_models": len(c), "extract_cpu_h": [r(v, 3) for v in c.extract_cpu_s / 3600],
                      "model_cpu_h": [r(v, 3) for v in c.model_cpu_h], "band": band,
                      "peak_rss_gb_max": r(c.peak_rss_gb.max(), 2), "failed": int(c.failed.sum()),
                      "n_open": [int(v) for v in c.n_open], "n_reopened": [int(v) for v in c.n_reopened],
                      "wt_ratio_generic": [r(v, 3) for v in c.wt_ratio_generic]}
    gp = W4 / "glpk_vs_gurobi.tsv"
    if gp.exists():
        g = pd.read_csv(gp, sep="\t")
        gl = pd.concat(pd.read_csv(p, sep="\t") for p in W4.glob("pilot_builds_glpk_*.tsv"))
        out = {}
        for f, x in g.groupby("family"):
            full = x[x["mode"] == "full_pipeline"].iloc[0]
            same = x[x["mode"] == "same_model_sweep"]
            out[f] = {"full_jaccard": r(full.reaction_set_jaccard, 4), "full_n_differ": int(full.n_ratio_r_differ),
                      "full_max_abs_diff": r(full.max_abs_diff, 6),
                      "same_model_n_differ": int(same.n_ratio_r_differ.iloc[0]) if len(same) else None,
                      "same_model_max_abs_diff": r(same.max_abs_diff.iloc[0], 6) if len(same) else None,
                      "glpk_over_gurobi_sweep": r(full.glpk_sweep_cpu_s / full.gurobi_sweep_cpu_s, 2)}
            # like-for-like build cost (extraction + repair + sweep, no task profile; same line, genomic none)
            lg = gl[(gl.family == f) & (gl.genomic == "none")].iloc[0]
            rg = d[(d.family == f) & (d.line == full.line) & (d.genomic == "none") & (d.variant == "real")].iloc[0]
            out[f]["glpk_build_cpu_s"] = r(lg.extract_cpu_s + lg.variant_cpu_s, 1)
            out[f]["gurobi_build_cpu_s"] = r(rg.extract_cpu_s + rg.variant_cpu_s - rg.tasks_cpu_s, 1)
            out[f]["glpk_over_gurobi_cpu"] = r(out[f]["glpk_build_cpu_s"] / out[f]["gurobi_build_cpu_s"], 2)
            out[f]["equal"] = (out[f]["full_n_differ"] == 0) and (out[f]["same_model_n_differ"] in (0, None))
            out[f]["glpk_allowed"] = out[f]["equal"] and (out[f]["glpk_over_gurobi_cpu"] or 99) <= 2.0
        S["glpk"] = out
        gb = list(W4.glob("pilot_builds_glpk_*.tsv"))
        if gb:
            S["glpk_worker_peak_rss_gb_max"] = r(pd.concat(pd.read_csv(p, sep="\t") for p in gb).peak_rss_gb.max(), 2)

    # Stage-A projection (frozen table; n_p from the panel)
    panel = pd.read_csv(W4 / "stage_a_panel.tsv", sep="\t")
    n, n_p = len(panel), int(panel.has_procan.sum())
    builds = {"A1": 4 * n, "A2": n, "A3": 6 * n_p}  # per LP family
    resweeps = 2 * n + n + 6 * n_p  # A4 per LP family: real builds of A1-A3
    cpu = {row: sum(k * fam[f]["build_cpu_min"] for f in LP) / 60 for row, k in builds.items()}
    cpu["A4"] = sum(resweeps * fam[f]["resweep_cpu_min"] for f in LP) / 60
    cpu["A7"] = 3 * 20 * (fam["F2"]["build_cpu_min"] + fam["F3"]["build_cpu_min"]) / 60
    if "corda" in S:
        per = max(S["corda"]["model_cpu_h"])
        cpu["A6"] = {"full": 120 * per, "fallback": 30 * per, "drop": 0.0}[S["corda"]["band"]]
    S["stage_a_projection_cpu_h"] = {k: r(v, 1) for k, v in cpu.items()}
    S["stage_a_projection_cpu_h"]["total"] = r(sum(cpu.values()), 1)
    S["stage_a_builds"] = {"A1": 4 * builds["A1"], "A2": 4 * builds["A2"], "A3": 4 * builds["A3"], "A4": 4 * resweeps,
                           "A7": 120, "n_panel": n, "n_p": n_p}

    # solver lanes (solver fixed per family): CPU-h per family over A1-A4 + A7 (+ A6 CORDA on Gurobi)
    per_fam = {}
    for f in LP:
        b = (builds["A1"] + builds["A2"] + builds["A3"]) * fam[f]["build_cpu_min"] + resweeps * fam[f]["resweep_cpu_min"]
        if f in ("F2", "F3"):
            b += 60 * fam[f]["build_cpu_min"]  # A7
        per_fam[f] = b / 60
    glpk_ratio = {f: (S.get("glpk", {}).get(f, {}).get("glpk_over_gurobi_cpu") or 1.0) for f in LP}
    a6 = S["stage_a_projection_cpu_h"].get("A6", 0.0) or 0.0

    def lanes(glpk_fams, n_gur=4, n_glpk=20):
        gur = sum(per_fam[f] for f in LP if f not in glpk_fams) + a6
        glp = sum(per_fam[f] * glpk_ratio[f] for f in glpk_fams)
        return {"glpk_families": list(glpk_fams), "gurobi_lane_cpu_h": r(gur, 1), "glpk_lane_cpu_h": r(glp, 1),
                "gurobi_workers": n_gur, "glpk_workers": n_glpk,
                "wall_h": r(max(gur / n_gur, glp / n_glpk if n_glpk else 0), 1)}
    S["stage_a_lanes"] = {"strict (F0, F2 on GLPK)": lanes(("F0", "F2")),
                          "if F4 tolerance 1e-6 accepted": lanes(("F0", "F2", "F4")),
                          "all Gurobi, 4 sessions": lanes((), 4, 0)}
    S["per_family_stage_a_cpu_h"] = {f: r(v, 1) for f, v in per_fam.items()}

    # A5 (amendment A4 d): distinct builds per line after deduplication, cost from the A5 GLPK smoke test
    sm = W4 / "stage_a_smoke" / "builds.tsv"  # post-fix smoke of every row (A5 rows on 2 lines)
    if sm.exists():
        b5 = pd.read_csv(sm, sep="\t")
        b5 = b5[b5.row == "A5"]
        a5_min = float(((b5.extract_cpu_s + b5.variant_cpu_s) / 60).mean())
        acc = json.loads((W4 / "tmp" / "accepted_losses_panel.json").read_text())
        n5 = sum(24 if p else (8 if acc.get(l) else 4) for l, p in zip(panel.ModelID, panel.has_procan.astype(bool)))
        a5 = n5 * a5_min / 60
        S["a5"] = {"grid_points": 24 * n, "distinct_builds_projected": int(n5), "glpk_build_cpu_min": r(a5_min),
                   "cpu_h": r(a5, 1), "smoke_builds": len(b5), "smoke_failed": int(b5.failed.sum()),
                   "smoke_repair_frac_max": r(b5.repair_frac_reactions.max(), 4)}
        # declared solvers (A4 c): GLPK F0, F2, F4, F5; Gurobi F3 (A6 dropped)
        glp = sum(per_fam[f] * glpk_ratio[f] for f in ("F0", "F2", "F4")) + a5
        gur = per_fam["F3"]
        S["stage_a_declared"] = {"gurobi_lane_cpu_h": r(gur, 1), "glpk_lane_cpu_h": r(glp, 1),
                                 "total_cpu_h": r(gur + glp, 1),
                                 "wall_h_4_gurobi_20_glpk": r(max(gur / 4, glp / 20), 1)}

    # Projection after review fixes + A7 (A6' CORDA, 20 builds) with one build per process and Gurobi pool = 8.
    # Per-job process overhead (spawn + model load + parent setup) from the post-fix smoke manifest.
    sm_dir = W4 / "stage_a_smoke"
    if (sm_dir / "MANIFEST.tsv").exists() and (sm_dir / "builds.tsv").exists() and "a5" in S:
        man = pd.read_csv(sm_dir / "MANIFEST.tsv", sep="\t").iloc[-1]
        bs = pd.read_csv(sm_dir / "builds.tsv", sep="\t")
        extra = json.loads(man.extra)
        n_jobs = int(bs.drop_duplicates("build_id").build_id.str.rsplit("|", n=1).str[0].nunique()) + extra.get("n_lines", 2)
        build_cpu = float(extra["build_cpu_h"]) * 3600
        over_s = max(0.0, (float(man.cpu_s) - build_cpu) / n_jobs)
        corda = S.get("corda", {}).get("model_cpu_h", [1.63])
        corda_h = float(np.mean(corda))
        n_gur = 270 + 20  # F3 jobs (A1 120, A2 30, A3 60, A7 60) + A6' 20
        n_glp = 3 * 210 + 60 + S["a5"]["distinct_builds_projected"] + 30  # F0/F2/F4 A1-A3, F2 A7, A5, loss rule
        gur_h = per_fam["F3"] + 20 * corda_h + n_gur * over_s / 3600
        glp_h = (sum(per_fam[f] * glpk_ratio[f] for f in ("F0", "F2", "F4")) + S["a5"]["cpu_h"]
                 + n_glp * over_s / 3600)
        S["stage_a_final"] = {
            "per_job_overhead_cpu_s": r(over_s, 1), "n_gurobi_jobs": n_gur, "n_glpk_jobs": n_glp,
            "a6prime_cpu_h": r(20 * corda_h, 1), "gurobi_lane_cpu_h": r(gur_h, 1), "glpk_lane_cpu_h": r(glp_h, 1),
            "total_cpu_h": r(gur_h + glp_h, 1),
            "wall_h_8_gurobi_20_glpk": r(max(gur_h / 8, glp_h / 20, corda_h * np.ceil(20 / 8)), 1),
            "memory_gb_8_plus_20": r(8 * 1.6 + 20 * 1.56 + 3.5, 1),
            "smoke_wall_s": float(man.wall_s), "smoke_cpu_h": r(float(man.cpu_s) / 3600, 2)}

    # worker counts from RSS
    vm = psutil.virtual_memory()
    rss_g = max(fam[f]["peak_rss_gb_max"] for f in LP)
    rss_l = S.get("glpk_worker_peak_rss_gb_max") or rss_g
    rss_c = S.get("corda", {}).get("peak_rss_gb_max") or rss_g
    avail = vm.available / 1e9
    S["memory"] = {"host_total_gb": r(vm.total / 1e9, 1), "available_gb_now": r(avail, 1),
                   "gurobi_worker_rss_gb": rss_g, "glpk_worker_rss_gb": rss_l, "corda_worker_rss_gb": rss_c,
                   "budget_gb": r(avail - GB_PER_SLOT_RESERVE, 1),
                   "max_workers_by_memory": int((avail - GB_PER_SLOT_RESERVE) // max(rss_g, rss_l, rss_c)),
                   "plan_24_plus_4_gb": r(24 * rss_l + 4 * max(rss_g, rss_c), 1),
                   "plan_28_gb": r(28 * max(rss_g, rss_l, rss_c), 1)}
    if (stage_dir / "builds.tsv").exists() and "F3b" in set(pd.read_csv(stage_dir / "builds.tsv", sep="\t").family):
        cmp_ = a6prime_comparison(stage_dir)
        S["a6prime_comparison"] = {"dir": str(stage_dir), "rows": len(cmp_),
                                   "file": str(stage_dir / "a6prime_comparison.tsv")}
    (W4 / "pilot_summary.json").write_text(json.dumps(S, indent=1))
    print(json.dumps(S, indent=1))


if __name__ == "__main__":
    main()
