"""Confirmation for the sterol finding (review M3-W2; label-free, GLPK): on Ham's M0 + cholesterol (MAR09285),
do vitamin D3 (MAR09442) or calcidiol (MAR09215) uptake rescue the post-squalene KOs?
Writes results/w2/review/sterol_vitd.tsv."""
from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"
import pandas as pd

from nextgem.model import BIOMASS_REACTION, hams_medium, load_model

ROOT = Path(__file__).resolve().parents[3]
GENES = ["FDFT1", "SQLE", "LSS", "NSDHL", "HSD17B7", "EBP", "DHCR24", "DHCR7", "HMGCR"]
rows = []
for parent in ("human1", "human2"):
    m = load_model(parent)
    m.objective = BIOMASS_REACTION
    s2g = dict(pd.read_csv(ROOT / "results" / "baseline" / f"{parent}_gene_mapping.tsv", sep="\t")[["symbol", "gene"]].values)
    base = {**hams_medium(m), "MAR09285": 1000.0}
    for tag, extra in (("chol", {}), ("chol+vitD3", {"MAR09442": 1000.0}), ("chol+calcidiol", {"MAR09215": 1000.0})):
        m.medium = {**base, **extra}
        wt = m.slim_optimize()
        for gs in GENES:
            with m:
                m.genes.get_by_id(s2g[gs]).knock_out()
                rows.append(dict(parent=parent, medium=tag, gene=gs, ratio=round(m.slim_optimize(error_value=0) / wt, 6)))
df = pd.DataFrame(rows).pivot_table(index="gene", columns=["parent", "medium"], values="ratio")
df.to_csv(Path(__file__).resolve().parent / "sterol_vitd.tsv", sep="\t")
print(df.to_string())
