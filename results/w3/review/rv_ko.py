"""M4 W3 review (non-author): independent single-gene KO sweeps for T0 and the two shortlisted B-loo variants.

Own code path: no nextgem.w3 / w2_engine import, no flux-support shortcut. Every model gene is knocked out with
cobra's Gene.knock_out() and the biomass LP is solved (GLPK). The variant is built from
results/w3/structure/variants_<parent>.json (kind 'remove' = drop that reactant from MAR13082).
Medium M0 = Ham's (GR task inputs, bound 1000) from nextgem.model.hams_medium, asserted equal to the author's M0.
ratio_r = clip(round(KO/WT, 6), 0, 1); non-optimal KO -> 0. Label-free.

Output: results/w3/review/rv_ratio_r_<parent>.parquet (rows = variant, columns = model genes), rv_ko_meta.tsv.
Usage: .venv/bin/python results/w3/review/rv_ko.py
"""
from __future__ import annotations

import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "w3" / "review"
VARIANTS = ("T0", "B-loo:cofactor_pool_biomass", "B-loo:metabolite_pool_biomass")


def job(args):
    parent, vname = args
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem.model import hams_medium, load_model

    t0 = time.process_time()
    m = load_model(parent)
    med = hams_medium(m)
    author_m0 = json.loads((ROOT / "results/w3/structure" / f"media_{parent}.json").read_text())["M0"]
    assert med == author_m0, "M0 differs from the author's M0"
    m.medium = med
    m.objective = "MAR13082"
    bm = m.reactions.get_by_id("MAR13082")
    if vname != "T0":
        spec = next(v for v in json.loads((ROOT / "results/w3/structure" / f"variants_{parent}.json").read_text())["bloo"]
                    if v["name"] == vname)
        assert spec["kind"] == "remove"
        for mid in spec["params"]["mets"]:
            met = m.metabolites.get_by_id(mid)
            coef = bm.metabolites[met]
            bm.add_metabolites({met: -coef})  # zero the coefficient -> metabolite removed from the reaction
            assert met not in bm.metabolites
    wt = m.slim_optimize()
    assert m.solver.status == "optimal" and wt > 1e-6
    genes = sorted(g.id for g in m.genes)
    ratio = np.empty(len(genes))
    nonopt = 0
    for i, gid in enumerate(genes):
        with m:
            m.genes.get_by_id(gid).knock_out()
            v = m.slim_optimize(error_value=np.nan)
            st = m.solver.status
        if st != "optimal" or not np.isfinite(v):
            ratio[i] = 0.0
            nonopt += st != "infeasible"
        else:
            ratio[i] = v / wt
    rr = np.clip(np.round(ratio, 6), 0.0, 1.0) + 0.0
    return parent, vname, genes, rr, {"parent": parent, "variant": vname, "wt_growth": wt, "n_genes": len(genes),
                                      "n_lethal_lt_0.01": int((rr < 0.01).sum()), "n_nonoptimal_not_infeasible": nonopt,
                                      "cpu_s": round(time.process_time() - t0, 1)}


def main():
    jobs = [(p, v) for p in ("human1", "human2") for v in VARIANTS]
    with mp.get_context("spawn").Pool(6) as pool:
        res = pool.map(job, jobs, chunksize=1)
    meta = []
    for p in ("human1", "human2"):
        rows = {v: pd.Series(rr, index=g) for (pp, v, g, rr, _) in res if pp == p}
        df = pd.DataFrame(rows).T.rename_axis("variant").reset_index()
        df.to_parquet(OUT / f"rv_ratio_r_{p}.parquet", index=False)
    meta = pd.DataFrame([r[4] for r in res])
    meta.to_csv(OUT / "rv_ko_meta.tsv", sep="\t", index=False)
    print(meta.to_string())


if __name__ == "__main__":
    main()
