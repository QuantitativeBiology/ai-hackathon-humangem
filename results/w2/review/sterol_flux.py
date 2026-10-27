"""Follow-up to sterol_trace.py (review M3-W2; label-free, GLPK): under FDFT1 KO + cholesterol open (Ham's M0),
a lathosterol sink restores growth. Trace where the supplied lathosterol goes (pFBA flux downstream, 4 steps), to
identify the non-cholesterol sterol product the parents require. Writes results/w2/review/sterol_flux.tsv."""
from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"
import pandas as pd
from cobra.flux_analysis import pfba

from nextgem.model import BIOMASS_REACTION, hams_medium, load_model

ROOT = Path(__file__).resolve().parents[3]
rows = []
for parent in ("human2", "human1"):
    m = load_model(parent)
    m.objective = BIOMASS_REACTION
    gm = pd.read_csv(ROOT / "results" / "baseline" / f"{parent}_gene_mapping.tsv", sep="\t")
    gid = dict(zip(gm.symbol, gm.gene))["FDFT1"]
    m.medium = {**hams_medium(m), "MAR09285": 1000.0}
    lat = m.metabolites.get_by_id("MAM02343c")
    with m:
        m.genes.get_by_id(gid).knock_out()
        m.add_boundary(lat, type="sink", reaction_id="SK_lathosterol_c", lb=-1000, ub=0)
        sol = pfba(m)
        rows.append(dict(parent=parent, step=0, reaction="SK_lathosterol_c", flux=sol.fluxes["SK_lathosterol_c"],
                         equation="lathosterol supplied", gpr=""))
        frontier, seen = {lat}, set()
        for step in range(1, 7):
            nxt = set()
            for met in frontier:
                for r in met.reactions:
                    v = sol.fluxes[r.id]
                    if r.id in seen or abs(v) < 1e-9 or r.id.startswith("SK_"):
                        continue
                    if v * r.metabolites.get(met, 0) < 0:  # consumes met in the solution direction
                        seen.add(r.id)
                        rows.append(dict(parent=parent, step=step, reaction=r.id, flux=v,
                                         equation=r.build_reaction_string(use_metabolite_names=True)[:160],
                                         gpr=r.gene_name_reaction_rule[:80]))
                        nxt |= {x for x, c in r.metabolites.items() if c * v > 0 and x.name not in
                                ("H2O", "H+", "NADP+", "NADPH", "O2", "NAD+", "NADH", "CO2", "Pi", "PPi", "ATP", "ADP")}
            frontier = nxt
    del m
df = pd.DataFrame(rows)
df.to_csv(Path(__file__).resolve().parent / "sterol_flux.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 160); pd.set_option("display.max_rows", 200)
print(df.to_string())
