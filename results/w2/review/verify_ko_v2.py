"""Re-check after A9 (review M3-W2; label-free, GLPK only: the Gurobi pool is held by W4 Stage A).

1. M0 sweep row of results/w2/ko/ratio_r_{parent}.parquet == results/baseline ratio_r, bitwise.
2. Sweeps and labels rebuilt with W2's own sweep_context/build_sweeps/labels_frame: every parquet sweep id and every
   (parent, variant, medium_key, f, gate_open) -> sweep_id row of ko/sweep_labels.tsv is reproduced.
3. Bound assertions for A9(b) and the 2-ME fix: cysteine uptake (MAR09065) closed in the default DMEM:F12/HPLM/Plasmax
   sweeps, open in the RPMI+cys(2-ME) sweeps and in the documented-cysteine sensitivities; ATP never opened.
4. Random gene x sweep pairs among the sweeps NEW in A9 (label_changes_A9.tsv: new_sweep_existed_before_A9 == False;
   4 GLPK-run + 4 Gurobi-run per parent), re-solved WITHOUT shortcuts on GLPK (stratified as in verify_ko.py).
Writes results/w2/review/ko_verify_v2_pairs.tsv and ko_verify_v2_summary.tsv.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_ko import closures, rr, solve  # noqa: E402  (same helpers as the first review)

from nextgem import media  # noqa: E402
from nextgem.knockouts import WT_GROWTH_TOL  # noqa: E402
from nextgem.model import BIOMASS_REACTION  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
W2 = ROOT / "results" / "w2"
OUT = Path(__file__).resolve().parent
SEED = 20261027
CYS = media.EXCHANGES["cysteine"][0]


def main() -> None:
    summary = []
    labels = pd.read_csv(W2 / "ko" / "sweep_labels.tsv", sep="\t")
    sw = pd.read_csv(W2 / "ko" / "sweeps.tsv", sep="\t").set_index("sweep_id")

    # 1. M0 bitwise
    for p in ("human1", "human2"):
        R = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
        sid = labels[(labels.parent == p) & (labels.variant == "M0")].sweep_id.unique()
        base = pd.read_parquet(ROOT / "results" / "baseline" / f"{p}_ko.parquet").set_index("gene").ratio_r
        x = R.loc[sid[0]].reindex(base.index).to_numpy()
        ok = len(sid) == 1 and set(R.columns) == set(base.index) and x.tobytes() == base.to_numpy(np.float64).tobytes()
        summary.append(dict(check=f"M0_bitwise_{p}", value=ok, detail=f"{len(base)} genes; {sid[0]}"))

    # 2. rebuild sweeps and labels
    from nextgem.w2_engine import build_sweeps, labels_frame, sweep_context
    lines, sheets, doc, phys, sens, pcs, exchanges, hams = sweep_context()
    sweeps = build_sweeps(exchanges, hams, pcs, doc, phys, sens)
    by_id = {s.sweep_id: s for s in sweeps}
    for p in ("human1", "human2"):
        R = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
        mine = {k for k in by_id if k.startswith(p)}
        summary.append(dict(check=f"sweeps_rebuilt_{p}", value=set(R.index) == mine,
                            detail=f"{len(R)} in parquet; rebuilt {len(mine)}; missing {len(set(R.index) - mine)}; "
                                   f"extra {len(mine - set(R.index))}"))
    lf = labels_frame(sweeps)
    key = ["parent", "variant", "medium_key", "f", "gate_open"]
    for d in (lf, labels):
        d["medium_key"] = d.medium_key.fillna("")
        d["f"] = d.f.astype(float).fillna(-1.0)
        d["gate_open"] = d.gate_open.astype(bool)
    j = lf[key + ["sweep_id"]].merge(labels[key + ["sweep_id"]], on=key, how="outer", suffixes=("_rebuilt", "_w2"),
                                     indicator=True)
    bad = j[(j._merge != "both") | (j.sweep_id_rebuilt != j.sweep_id_w2)]
    summary.append(dict(check="labels_rebuilt", value=len(bad) == 0,
                        detail=f"{len(labels)} W2 labels, {len(lf)} rebuilt, mismatches {len(bad)}"))

    # 3. bound assertions
    def lab_bounds(p, mk, v="M1", f=None):
        r = labels[(labels.parent == p) & (labels.medium_key == mk) & (labels.variant == v) & labels.gate_open
                   & (labels.f == (-1.0 if f is None else f))]
        return [by_id[s].bounds for s in r.sweep_id.unique()]
    msgs, ok = [], True
    for p in ("human1", "human2"):
        for mk, v, f, want in (("DMEM:F12", "M1", None, False), ("RPMI+cys(2-ME)", "M1", None, True),
                               ("RPMI+cys(2-ME)", "M2", 0.01, True), ("RPMI", "M1", None, False),
                               ("M-HPLM", "M2", 0.01, False), ("M-Plasmax", "M2", 0.01, False)):
            bs = lab_bounds(p, mk, v, f)
            got = [CYS in b for b in bs]
            ok &= bool(bs) and all(g == want for g in got)
            msgs.append(f"{p}:{mk}:{v}:{f} cys_open={got} (want {want})")
        sens_doc = [k for k in sens if "cys" in k.lower() and "document" in k.lower()]
        for k in sens_doc:
            bs = lab_bounds(p, k, "M2" if k.startswith("M-") else "M1", 0.01 if k.startswith("M-") else None)
            ok &= bool(bs) and all(CYS in b for b in bs)
            msgs.append(f"{p}:{k} cys_open={[CYS in b for b in bs]} (want True)")
    atp_open = any("MAR00569" in s.bounds for s in sweeps if s.labels[0]["variant"] not in ("M-open",))
    summary.append(dict(check="cysteine_2ME_bounds", value=ok, detail="; ".join(msgs)))
    summary.append(dict(check="ATP_uptake_never_open_outside_M-open", value=not atp_open, detail=""))

    # 4. random pairs among NEW sweeps, GLPK, no shortcuts
    from nextgem.w2_runtime import load_parent
    rng = np.random.default_rng(SEED)
    # "new in A9" = sweeps that did not exist before A9 (label_changes_A9.tsv); the 196 new Gurobi sweeps carry
    # reused_from_previous_run=True because the final assembly resumed from A9 checkpoints
    lc = pd.read_csv(W2 / "ko" / "label_changes_A9.tsv", sep="\t")
    new_ids = set(lc[lc.new_sweep_existed_before_A9 == False].sweep_id_new.dropna())  # noqa: E712
    new = sw[sw.index.isin(new_ids)]
    summary.append(dict(check="new_sweeps_in_A9", value=len(new) == 216,
                        detail=f"{len(new)} new sweeps: {new.solver.value_counts().to_dict()}"))
    rows = []
    for p in ("human1", "human2"):
        m = load_parent(p, solver="glpk")
        m.objective = BIOMASS_REACTION
        clos = closures(m)
        R = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
        newp = new[new.parent == p]
        glpk_ids = sorted(newp[newp.solver == "glpk"].index)
        gur_ids = sorted(newp[newp.solver == "gurobi"].index)
        pick = list(rng.choice(glpk_ids, size=min(4, len(glpk_ids)), replace=False)) + \
            list(rng.choice(gur_ids, size=min(4, len(gur_ids)), replace=False))
        for sid in pick:
            bounds = by_id[sid].bounds
            row = R.loc[sid]
            with m:
                m.medium = bounds
                wt, _ = solve(m)
            g_all = list(row.index)
            s1 = [g for g in g_all if row[g] == 1.0 and clos[g]]
            s2 = [g for g in g_all if row[g] < 1.0]
            picks = [(g, "ratio1_closure") for g in rng.choice(s1, size=min(2, len(s1)), replace=False)]
            picks += [(g, "ratio_lt1") for g in rng.choice(s2, size=min(2, len(s2)), replace=False)]
            for g, tag in picks:
                with m:
                    m.medium = bounds
                    m.genes.get_by_id(g).knock_out()
                    kg, kst = solve(m)
                ratio = kg / wt if wt > WT_GROWTH_TOL else math.nan
                rows.append(dict(parent=p, sweep_id=sid, sweep_solver=sw.loc[sid, "solver"], gene=g, stratum=tag,
                                 wt_glpk=wt, wt_w2=sw.loc[sid, "wt_growth"], ratio_r_review=rr(ratio),
                                 ratio_r_w2=float(row[g]), abs_diff=abs(rr(ratio) - float(row[g]))))
        del m
    pairs = pd.DataFrame(rows)
    pairs.to_csv(OUT / "ko_verify_v2_pairs.tsv", sep="\t", index=False)
    summary.append(dict(check="pairs_new_sweeps_glpk_no_shortcut", value=bool((pairs.abs_diff <= 1e-5).all()),
                        detail=f"{len(pairs)} pairs over {pairs.sweep_id.nunique()} new sweeps "
                               f"({pairs[pairs.sweep_solver == 'glpk'].sweep_id.nunique()} GLPK-run); max |diff| "
                               f"{pairs.abs_diff.max():.2e}; n exact {int((pairs.abs_diff == 0).sum())}"))
    s = pd.DataFrame(summary)
    s.to_csv(OUT / "ko_verify_v2_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 400)
    print(s.to_string())
    if len(bad):
        print(bad.head(20).to_string())


if __name__ == "__main__":
    main()
