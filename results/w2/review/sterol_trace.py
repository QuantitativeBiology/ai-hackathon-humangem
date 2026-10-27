"""Why does serum cholesterol not rescue the post-squalene sterol genes? (review M3-W2; label-free, GLPK).

For each parent, on Ham's M0 with and without the cholesterol uptake (MAR09285) open at 1000:
  1. ratio of each sterol-branch KO (and HMGCR);
  2. for KOs still lethal with cholesterol: rescue by a sink supplying each named post-squalene intermediate
     (STEROL_NAMES, compartments c and r) -> which intermediate is required;
  3. the pool/biomass reactions that consume sterol-like metabolites (where the requirement enters biomass).
Writes results/w2/review/sterol_trace.tsv.
"""
from __future__ import annotations

import re
from pathlib import Path

import sys

import cobra

# GLPK by default; --gurobi takes one session from nextgem.gurobi_pool (faster for the sink scan)
SOLVER = "gurobi" if "--gurobi" in sys.argv else "glpk"
cobra.Configuration().solver = SOLVER

import pandas as pd

from nextgem.model import BIOMASS_REACTION, hams_medium, load_model

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / "sterol_trace.tsv"
GENES = ["HMGCR", "FDFT1", "SQLE", "LSS", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", "HSD17B7", "EBP", "SC5D", "DHCR7",
         "DHCR24"]
CHOL_EX = "MAR09285"
# named intermediates of the post-squalene pathway (Bloch and Kandutsch-Russell); compartments c and r
STEROL_NAMES = {"squalene", "(S)-2,3-epoxysqualene", "2,3-oxidosqualene", "lanosterol", "24,25-dihydrolanosterol",
                "zymosterol", "zymostenol", "lathosterol", "7-dehydrocholesterol", "desmosterol",
                "7-dehydrodesmosterol", "5alpha-cholesta-7,24-dien-3beta-ol", "cholesterol", "14-demethyllanosterol",
                "4alpha-methylzymosterol"}
SCAN_GENES = {"FDFT1", "SQLE", "LSS", "DHCR24", "EBP", "NSDHL", "HSD17B7"}


def growth(m) -> float:
    v = m.slim_optimize(error_value=0.0)
    return max(0.0, v or 0.0)


def main() -> None:
    slot = None
    if SOLVER == "gurobi":
        from nextgem.gurobi_pool import hold_slot
        slot = hold_slot(timeout=3600)
    rows = []
    for parent in ("human1", "human2"):
        m = load_model(parent)
        m.objective = BIOMASS_REACTION
        gmap = pd.read_csv(ROOT / "results" / "baseline" / f"{parent}_gene_mapping.tsv", sep="\t")
        s2g = dict(zip(gmap.symbol, gmap.gene))
        med = hams_medium(m)
        sterols = [x for x in m.metabolites if x.name in STEROL_NAMES and x.compartment in ("c", "r")]
        # where do sterol-like metabolites enter biomass? (multi-reactant pool reactions consuming them)
        for x in sterols:
            for r in x.reactions:
                if r.metabolites[x] < 0 and (r.id == BIOMASS_REACTION or len(r.reactants) > 5):
                    rows.append(dict(parent=parent, test="pool_consumer", gene="", metabolite=f"{x.name}[{x.compartment}]",
                                     value=r.metabolites[x], detail=f"{r.id} -> {', '.join(p.name for p in r.products)}"))
        for chol in (False, True):
            m.medium = {**med, CHOL_EX: 1000.0} if chol else dict(med)
            wt = growth(m)
            for gs in GENES:
                gid = s2g.get(gs)
                if gid is None or gid not in m.genes:
                    rows.append(dict(parent=parent, test="ko", gene=gs, cholesterol=chol, value=None, detail="absent"))
                    continue
                with m:
                    m.genes.get_by_id(gid).knock_out()
                    g = growth(m)
                    rows.append(dict(parent=parent, test="ko", gene=gs, cholesterol=chol, value=round(g / wt, 6)))
                    if chol and g < 0.01 * wt and gs in SCAN_GENES:
                        for x in sterols:
                            with m:
                                m.add_boundary(x, type="sink", reaction_id=f"SKrev_{x.id}", lb=-1000, ub=0)
                                gr = growth(m)
                            if gr > 0.01 * wt:
                                rows.append(dict(parent=parent, test="rescue_by_sink", gene=gs, cholesterol=chol,
                                                 metabolite=f"{x.name}[{x.compartment}] {x.id}",
                                                 value=round(gr / wt, 6)))
        del m
    if slot is not None:
        slot.release()
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 120)
    pd.set_option("display.max_rows", 500)
    print(df.to_string())


if __name__ == "__main__":
    main()
