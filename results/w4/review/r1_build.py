"""W4 review, part 1: rebuild panel lines x families at R1 (real RNA, no genomic constraint) with the author's
pipeline on the declared solver, and compare with the reviewer's own GIMME / E-Flux (own.py) and an own
brute-force KO sweep under the A4 (a) convention. FASTCORE (F3): core coverage, and the support is saved for
a cobra-fastcc consistency check (r3_fastcc.py).

Usage: .venv/bin/python results/w4/review/r1_build.py --solver glpk --families F2 F4 --lines ACH-000957 ...
Gurobi runs hold one slot of nextgem.gurobi_pool for the life of the process. Label-free.
"""

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import cobra

ap = argparse.ArgumentParser()
ap.add_argument("--solver", required=True, choices=["glpk", "gurobi"])
ap.add_argument("--families", nargs="+", required=True)
ap.add_argument("--lines", nargs="+", required=True)
ap.add_argument("--tag", default="")
ap.add_argument("--no-own-ko", action="store_true")
a = ap.parse_args()

slot = None
if a.solver == "gurobi":
    from nextgem import gurobi_pool

    slot = gurobi_pool.hold_slot(timeout=3600, poll=10)
cobra.Configuration().solver = a.solver  # before any model exists: GLPK runs open no Gurobi session

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import own  # noqa: E402

from nextgem import context as cx  # noqa: E402
from nextgem import context_inputs as ci  # noqa: E402
from nextgem import paths  # noqa: E402
from nextgem.knockouts import rounded_ratio  # noqa: E402
from nextgem.model import load_model  # noqa: E402

R = paths.results_dir() / "w4"
OUT = R / "review"
TMP = OUT / "tmp"
t0 = time.time()
m = load_model("human2")
assert cx.solver_name(m) == a.solver
parent = cx.make_parent(m, "human2")
if "F3" in a.families:
    cx.consistent_reactions(parent, R / "fastcc_human2_M0.json")
genes_tab = ci.gene_table(m)
expr = pd.read_parquet(TMP / "panel_expr.parquet")
st = json.loads((R / "threshold_stats.json").read_text())
gm = pd.read_parquet(R / "threshold_gene_means.parquet").set_index("gene")
rna_stats = ci.ThresholdStats("rna", st["rna"]["p25"], st["rna"]["p75"], gm.rna_mean.dropna(), st["rna"]["floor_value"], st["rna"]["n_lines"])
pilot = pd.read_csv(R / "pilot_builds_gurobi_lp.tsv", sep="\t")
pilot_ko = pd.read_parquet(R / "tmp" / "pilot_ko_gurobi_lp.parquet")
glpk_pilot = pd.read_csv(R / "pilot_builds_glpk_ACH-000957.tsv", sep="\t")


def mask(open_set):
    return np.packbits(np.array([r in open_set for r in parent.rxn_ids])).tobytes().hex()


def jacc(a_, b_):
    return len(a_ & b_) / max(1, len(a_ | b_))


results = []
for line in a.lines:
    s_own = (expr.loc[line] - 1.0).to_dict()  # R1: log2(TPM+1) - 1, own reader
    s_auth = ci.signed_scores(line, "R1", "rna", genes_tab, expr, rna_stats).reindex(parent.genes)
    d = pd.Series(s_own).reindex(parent.genes)
    assert np.allclose(d.fillna(-999), s_auth.fillna(-999)), "own and author signed scores differ"
    off = {g for g, v in s_own.items() if v < 0}
    for fam in a.families:
        c0 = time.process_time()
        lm = cx.extract(parent, fam, s_auth)
        lm.line = line
        cx.finish(parent, lm)
        ko = cx.ko_sweep(parent, lm)
        auth_open = set(parent.rxn_ids) - lm.final_closed
        r = {"line": line, "family": fam, "solver": a.solver, "author_cpu_s": round(time.process_time() - c0, 1),
             "author_n_open": len(auth_open), "author_n_reopened": len(lm.reopened), "author_failed": bool(lm.info["failed"]),
             "author_wt": lm.info["wt_growth"], "author_n_lethal": int((ko.ratio_r < 0.01).sum()),
             "author_open_mask_sha": hashlib.sha256(mask(auth_open).encode()).hexdigest()[:16],
             "author_ko_sha": hashlib.sha256(pd.util.hash_pandas_object(ko[["gene", "ratio_r"]], index=False).values.tobytes()).hexdigest()[:16],
             "author_n_nonoptimal": int(ko.ratio.isna().sum())}
        pr = pilot[(pilot.line == line) & (pilot.family == fam) & (pilot.variant == "real") & (pilot.genomic == "none")]
        if len(pr):
            pm = pr.open_mask.iloc[0]
            bits = np.unpackbits(np.frombuffer(bytes.fromhex(pm), dtype=np.uint8))[: len(parent.rxn_ids)].astype(bool)
            p_open = {rid for rid, b in zip(parent.rxn_ids, bits) if b}
            pk = pilot_ko[pilot_ko.build_id == f"{line}|{fam}|R1|rna|real|none"].set_index("gene").ratio_r
            r.update({"vs_gurobi_pilot_open_equal": p_open == auth_open, "vs_gurobi_pilot_jaccard": round(jacc(p_open, auth_open), 5),
                      "vs_gurobi_pilot_ko_n_differ": int((ko.set_index("gene").ratio_r.reindex(pk.index) - pk).abs().gt(1e-9).sum())})
        gp = glpk_pilot[(glpk_pilot.line == line) & (glpk_pilot.family == fam) & (glpk_pilot.genomic == "none")] if "open_mask" in glpk_pilot else []
        if len(gp):
            bits = np.unpackbits(np.frombuffer(bytes.fromhex(gp.open_mask.iloc[0]), dtype=np.uint8))[: len(parent.rxn_ids)].astype(bool)
            g_open = {rid for rid, b in zip(parent.rxn_ids, bits) if b}
            r["vs_glpk_pilot_open_equal"] = g_open == auth_open
        auth_ratio = ko.set_index("gene").ratio
        auth_ess = set(ko.gene[ko.ratio_r < 0.01])
        own_res = None
        if fam == "F2":
            s_r = own.reaction_scores(m, s_own)
            closed_own, obj = own.gimme(m, s_r, parent.wt_generic, 0.9)
            own_open = set(parent.rxn_ids) - closed_own
            r.update({"own_gimme_obj": obj, "author_gimme_obj": lm.info.get("gimme_objective"),
                      "own_n_open": len(own_open), "open_jaccard_own_vs_author": round(jacc(own_open, auth_open), 6),
                      "open_symdiff": len(own_open ^ auth_open), "author_extract_closed_equal_own": lm.closed == closed_own})
            if not a.no_own_ko:
                own_res = own.ko_presence(m, closed_own, off)
        elif fam == "F4":
            lin = own.eflux_lin(s_own)
            b_own, cap_max = own.eflux_bounds(m, lin)
            diff = [rid for rid in set(b_own) | set(lm.eflux) if b_own.get(rid) != lm.eflux.get(rid)]
            r.update({"own_n_bounded": len(b_own), "author_n_bounded": len(lm.eflux), "n_bound_differ": len(diff),
                      "cap_max_own": cap_max, "cap_max_author": lm.cap_max,
                      "n_at_floor": int(sum(1 for rid, (lb, ub) in b_own.items() if max(abs(lb), abs(ub)) == 10.0)),
                      "own_n_open": len(parent.rxn_ids), "open_symdiff": len(set(parent.rxn_ids) ^ auth_open)})
            if not a.no_own_ko:
                own_res = own.ko_eflux(m, b_own, lin, cap_max)
                # literal reading of "KO zeroes the gene capacity" (no GPR closure; floor 10 applies) for the
                # genes the author calls lethal
                lit = own.ko_eflux(m, b_own, lin, cap_max, literal=True, genes=sorted(auth_ess))
                r["literal_reading_lethal_of_author_lethal"] = int(sum(1 for v in lit[1].values() if v < 0.01))
        elif fam in ("F3", "F0"):
            A = set(parent.rxn_ids) - lm.closed
            s_r = cx.reaction_scores(parent.gprs, s_auth)
            core = {rid for rid, v in s_r.items() if not math.isnan(v) and v >= 0}
            if fam == "F3":
                cons = parent.consistent
                r.update({"n_core": len(core), "n_core_consistent": len(core & cons), "n_core_cons_missing_from_support": len((core & cons) - A),
                          "author_unreachable_core": lm.info.get("unreachable_core"), "n_support": len(A),
                          "support_minus_consistent": len(A - cons), "n_open_final": len(auth_open)})
                (TMP / f"f3_sets_{line}_{a.solver}{a.tag}.json").write_text(json.dumps(
                    {"support": sorted(A), "final_open": sorted(auth_open), "core": sorted(core & cons)}))
            if not a.no_own_ko:
                own_res = own.ko_presence(m, lm.final_closed, off)
        if own_res is not None:
            wt_o, res = own_res
            own_ratio = rounded_ratio(pd.Series(res))
            own_ess = set(own_ratio.index[own_ratio < 0.01])
            cmp = pd.concat([ko.set_index("gene").ratio_r, own_ratio.rename("own")], axis=1)
            r.update({"own_wt": wt_o, "own_n_lethal": len(own_ess), "lethal_intersection": len(own_ess & auth_ess),
                      "lethal_only_author": sorted(auth_ess - own_ess)[:20], "lethal_only_own": sorted(own_ess - auth_ess)[:20],
                      "ko_ratio_n_differ_gt_1e-6": int((cmp.ratio_r - cmp.own).abs().gt(1e-6).sum()),
                      "ko_ratio_max_abs_diff": float((cmp.ratio_r - cmp.own).abs().max())})
        ko[["gene", "ratio_r", "status", "solver"]].to_parquet(TMP / f"ko_{line}_{fam}_{a.solver}{a.tag}.parquet", index=False)
        r["wall_s"] = round(time.time() - t0)
        results.append(r)
        print(json.dumps(r, default=str), flush=True)
pd.DataFrame(results).to_csv(OUT / f"r1_{a.solver}_{'_'.join(a.families)}_{'_'.join(a.lines)}{a.tag}.tsv", sep="\t", index=False)
if slot is not None:
    slot.release()
