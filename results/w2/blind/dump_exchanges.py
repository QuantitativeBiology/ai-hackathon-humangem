"""Dump exchange reactions (id, metabolite id/name/formula) for Human1 and Human2 (blind extraction helper)."""
import cobra
cobra.Configuration().solver = "glpk"
import pandas as pd
from nextgem.model import load_model

rows = []
for parent in ("human1", "human2"):
    m = load_model(parent)
    for r in m.exchanges:
        for met in r.metabolites:
            rows.append(dict(parent=parent, rxn=r.id, met=met.id, name=met.name, formula=met.formula,
                             comp=met.compartment, lb=r.lower_bound, ub=r.upper_bound))
pd.DataFrame(rows).to_csv("results/w2/blind/tmp/exchanges.tsv", sep="\t", index=False)
print(len(rows))
