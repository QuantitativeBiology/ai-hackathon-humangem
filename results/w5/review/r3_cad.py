"""W5 review, part 3 (LABEL-FREE): is (S)-dihydroorotate a biomass constituent via MAR10064, and does that explain CAD
lethality under uridine? GLPK, Ham's (M0) +/- uridine / orotate exchange opened at 1000 (same as other M0 inputs).
Usage: OMP_NUM_THREADS=1 .venv/bin/python results/w5/review/r3_cad.py  -> results/w5/review/r3_cad.json
"""
import json
import time

import cobra

cobra.Configuration().solver = "glpk"
from nextgem.knockouts import configure  # noqa: E402
from nextgem.model import load_model  # noqa: E402

T0 = time.time()
out = {}


def exch_for(m, name):
    hits = []
    for r in m.boundary:
        for met in r.metabolites:
            if met.name.lower() == name.lower():
                hits.append(r.id)
    return sorted(set(hits))


def grow(m):
    v = m.slim_optimize(error_value=float("nan"))
    return 0.0 if m.solver.status == "infeasible" else round(v, 4)


for p in ("human2", "human1"):
    m = load_model(p)
    configure(m)
    o = {}
    pool = m.reactions.get_by_id("MAR10064")
    dho = [x for x in pool.metabolites if "dihydroorotate" in x.name.lower()]
    o["MAR10064_name"] = pool.name
    o["dihydroorotate_in_MAR10064"] = {x.id: [x.name, x.compartment, pool.metabolites[x]] for x in dho}
    prods = [x.id for x in pool.products]
    bm = m.reactions.get_by_id("MAR13082")
    o["MAR10064_product_in_biomass"] = {pid: bm.metabolites.get(m.metabolites.get_by_id(pid), 0) for pid in prods}
    allm = [x for x in m.metabolites if x.name.lower() in ("(s)-dihydroorotate", "dihydroorotate")]
    o["dho_reactions"] = {x.id: sorted(f"{r.id}:{r.gene_reaction_rule[:60]}:lb={r.lower_bound}" for r in x.reactions) for x in allm}
    uri = exch_for(m, "uridine")
    oro = exch_for(m, "orotate")
    o["uridine_exchange"], o["orotate_exchange"] = uri, oro
    cad = [g for g in m.genes if g.name == "CAD"] or [g for g in m.genes if "CAD" == str(getattr(g, "name", ""))]
    gid = {s: next((g.id for g in m.genes if g.name == s), None) for s in ("CAD", "DHODH", "UMPS")}
    if not all(gid.values()):
        import pandas as pd
        from nextgem.model import GEMS
        gt = pd.read_csv(GEMS[p].genes_tsv, sep="\t", dtype=str)
        gid = {s: gt.genes[gt.geneSymbols == s].iloc[0] for s in ("CAD", "DHODH", "UMPS")}
    o["genes"] = gid
    res = {}
    for cond, extra in (("hams", []), ("hams+uridine", uri), ("hams+orotate", oro)):
        with m:
            med = dict(m.medium)
            for x in extra:
                med[x] = 1000.0
            m.medium = med
            wt = grow(m)
            r = {"wt": wt}
            for s, g in gid.items():
                with m:
                    m.genes.get_by_id(g).knock_out()
                    r[f"{s}_KO_ratio"] = round(grow(m) / wt, 6) if wt > 0 else None
            # dihydroorotate removed from the biomass pool (diagnostic)
            with m:
                for x in dho:
                    pool.add_metabolites({x: -pool.metabolites[x]})
                wt2 = grow(m)
                with m:
                    m.genes.get_by_id(gid["CAD"]).knock_out()
                    r["noDHOpool_wt"] = wt2
                    r["noDHOpool_CAD_KO_ratio"] = round(grow(m) / wt2, 6) if wt2 > 0 else None
            res[cond] = r
    o["growth"] = res
    out[p] = o
out["wall_s"] = round(time.time() - T0, 1)
(open("/data/tmpA/nextgem/results/w5/review/r3_cad.json", "w")).write(json.dumps(out, indent=1, default=str))
print(json.dumps(out, indent=1, default=str))
