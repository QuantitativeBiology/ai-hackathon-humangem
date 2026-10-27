"""W4 review, part 3: loss rule (own implementation vs author) and the repair rule on 3 panel lines (GLPK only).

Usage: .venv/bin/python results/w4/review/r3_loss_repair.py --lines ACH-000957 ACH-002024 ACH-003383 [--tag _b]
Label-free (CN, damaging-mutation matrix and expression inputs only).
"""

import argparse
import hashlib
import json
import math
import re
import sys
import time
from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.csv as pacsv  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from nextgem import context as cx  # noqa: E402
from nextgem import context_inputs as ci  # noqa: E402
from nextgem import paths  # noqa: E402
from nextgem.model import load_model  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--lines", nargs="+", required=True)
ap.add_argument("--tag", default="")
a = ap.parse_args()
R = paths.results_dir() / "w4"
OUT = R / "review"
TMP = OUT / "tmp"
t0 = time.time()
m = load_model("human2")
parent = cx.make_parent(m, "human2")
genes_tab = ci.gene_table(m)
expr = pd.read_parquet(TMP / "panel_expr.parquet")


def own_read(fname, lines):
    f = paths.depmap_dir() / fname
    hdr = pd.read_csv(f, nrows=0).columns
    col_of = {}
    for c in hdr:
        mm = re.match(r"^(.*) \((\d+)\)$", c)
        if mm:
            col_of.setdefault(mm.group(2), c)
    want = {g: col_of[e] for g, e in genes_tab.entrez.items() if isinstance(e, str) and e in col_of}
    meta = ["ModelID", "IsDefaultEntryForModel"]
    t = pacsv.read_csv(f, convert_options=pacsv.ConvertOptions(include_columns=meta + sorted(set(want.values())))).to_pandas()
    t = t[(t.IsDefaultEntryForModel == "Yes") & t.ModelID.isin(lines)].set_index("ModelID")
    return pd.DataFrame({g: t[c] for g, c in want.items()}, index=t.index).astype("float64")


cn = own_read(ci.CN_FILE, a.lines)
dmg = own_read(ci.DAMAGING_FILE, a.lines)
cn_auth = ci.read_cn(genes_tab, a.lines)
dmg_auth = ci.read_damaging(genes_tab, a.lines)
acc_pilot = json.loads((R / "tmp" / "accepted_losses_panel.json").read_text())
acc_smoke = {}
for d in ("stage_a_smoke", "stage_a_smoke_a5"):
    p = R / d / "accepted_losses.json"
    if p.exists():
        acc_smoke.update(json.loads(p.read_text()))

floor = 0.01 * parent.wt_generic


def growth_closed(closed):
    with m:
        for rid in closed:
            m.reactions.get_by_id(rid).bounds = (0, 0)
        v = m.slim_optimize()
        return 0.0 if m.solver.status == "infeasible" else v


def own_closure(gs):
    out = set()
    for g in gs:
        for r in m.genes.get_by_id(g).reactions:
            if not r.gpr.eval(knockouts=set(gs)):
                out.add(r.id)
    return out


rows = []
for line in a.lines:
    r = {"line": line}
    # ---- loss rule: own definition (r2-11) on own-read inputs
    d = dmg.loc[line] if line in dmg.index else pd.Series(dtype=float)
    c = cn.loc[line] if line in cn.index else pd.Series(dtype=float)
    e = expr.loc[line]
    idx = d.index.union(c.index)
    d, c, e2 = d.reindex(idx), c.reindex(idx), e.reindex(idx)
    own_lost = sorted(idx[((d >= 2) | ((d == 1) & ((c < 0.5) | (e2 < 1.0))) | (c < 0.3)).fillna(False).to_numpy()])
    auth_lost, flags = ci.lost_genes(line, dmg_auth, cn_auth, expr)
    r.update({"own_lost": len(own_lost), "author_lost": len(auth_lost), "lost_equal": own_lost == auth_lost})
    # ---- order-independent drop: author under 6 shuffled orders; own implementation
    model_lost = [g for g in auth_lost if g in set(parent.genes)]
    rng = np.random.default_rng(11)
    outs = []
    for k in range(6):
        order = list(rng.permutation(model_lost)) if k else list(model_lost)
        outs.append(cx.accept_losses(parent, order, expr.loc[line]))
    r["author_order_independent"] = all(o["accepted"] == outs[0]["accepted"] and o["contradicted"] == outs[0]["contradicted"]
                                        and o["dropped_joint"] == outs[0]["dropped_joint"] for o in outs)
    contra = sorted(g for g in model_lost if own_closure([g]) and growth_closed(own_closure([g])) < floor)
    rem = [g for g in model_lost if g not in contra]
    key = sorted(rem, key=lambda g: (-(e.get(g) if not pd.isna(e.get(g, np.nan)) else -np.inf),
                                     int(genes_tab.entrez.get(g)) if str(genes_tab.entrez.get(g)).isdigit() else 10**12))
    dropped = []
    while rem and growth_closed(own_closure(rem)) < floor:
        g = key.pop(0)
        dropped.append(g)
        rem.remove(g)
    r.update({"own_accepted": len(rem), "author_accepted": len(outs[0]["accepted"]),
              "accepted_equal_own": sorted(rem) == outs[0]["accepted"],
              "contradicted_equal_own": contra == outs[0]["contradicted"], "n_contradicted": len(contra),
              "n_dropped_joint": len(dropped),
              "accepted_equal_pilot_file": sorted(acc_pilot.get(line, [])) == outs[0]["accepted"] if line in acc_pilot else None,
              "accepted_equal_smoke": (acc_smoke[line]["accepted"] == outs[0]["accepted"]) if line in acc_smoke else None})
    accepted = set(outs[0]["accepted"])
    # ---- repair (F0), none and loss variants; determinism within process; own growth-LP check
    s = ci.signed_scores(line, "R1", "rna", genes_tab, expr,
                         ci.ThresholdStats("rna", 0, 0, pd.Series(dtype=float), 0, 0)).reindex(parent.genes)
    base = cx.extract(parent, "F0", s)
    for gvar, lost in (("none", set()), ("loss", accepted)):
        reps = []
        for k in range(2):
            lm = cx.LineModel(line=line, family="F0", closed=set(base.closed), off=set(base.off), info=dict(base.info))
            cx.finish(parent, lm, lost)
            reps.append(lm)
        lm = reps[0]
        ko = cx.ko_sweep(parent, lm)
        kh = hashlib.sha256(pd.util.hash_pandas_object(ko[["gene", "ratio_r"]], index=False).values.tobytes()).hexdigest()[:16]
        # growth before repair
        wt_before = growth_closed(lm.closed | lm.lost_closed)
        # own growth LP: min sum w|v| over repairable closed reactions s.t. biomass >= 1% generic
        s_r = lm.info["_s_r"]
        repairable = lm.closed - lm.lost_closed
        with m:
            for rid in lm.lost_closed:
                m.reactions.get_by_id(rid).bounds = (0, 0)
            m.reactions.get_by_id(cx.BIOMASS_REACTION).lower_bound = floor
            terms = []
            for rid in repairable:
                v = s_r.get(rid, math.nan)
                w = 1.0 if math.isnan(v) else 1.0 + max(0.0, -v)
                rx = m.reactions.get_by_id(rid)
                terms += [w * rx.forward_variable, w * rx.reverse_variable]
            m.objective = m.problem.Objective(sum(terms), direction="min")
            m.slim_optimize()
            own_growth_set = {rid for rid in repairable if abs(m.reactions.get_by_id(rid).flux) > 1e-6}
        prof = cx.task_profile(parent, lm)
        gen = cx.task_profile(parent, cx.LineModel(line="generic", family="F0", closed=set()))
        viab_ids = [k for k in prof if k.split("_", 1)[1] in {t.description for t in parent.tasks}]
        r.update({
            f"{gvar}_wt_before_repair": wt_before, f"{gvar}_wt_after": lm.info["wt_growth"],
            f"{gvar}_restored": lm.info["wt_growth"] >= floor, f"{gvar}_n_reopened": lm.info["n_reopened"],
            f"{gvar}_repair_frac": lm.info["repair_frac_reactions"],
            f"{gvar}_reopened_deterministic": reps[0].reopened == reps[1].reopened,
            f"{gvar}_reopen_hits_lost_closed": len(lm.reopened & lm.lost_closed),
            f"{gvar}_n_lost_closed": len(lm.lost_closed),
            f"{gvar}_own_growth_lp_n": len(own_growth_set), f"{gvar}_own_growth_subset_of_reopened": own_growth_set <= lm.reopened,
            f"{gvar}_viability_pass_line": sum(prof[k] for k in viab_ids), f"{gvar}_viability_pass_generic": sum(gen[k] for k in viab_ids),
            f"{gvar}_ko_sha": kh, f"{gvar}_open_n": lm.info["n_open"]})
    r["wall_s"] = round(time.time() - t0)
    rows.append(r)
    print(json.dumps(r, default=str), flush=True)
pd.DataFrame(rows).to_csv(OUT / f"r3_loss_repair{a.tag}.tsv", sep="\t", index=False)
