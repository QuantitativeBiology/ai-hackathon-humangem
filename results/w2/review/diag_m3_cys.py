"""Diagnostics for the W2 review (label-free, GLPK): M3 calibration mechanism, SLC7A11/GSH, ATP supplement.

A. M3 (RPMI, M2 f=0.01, glucose 0.5, glutamine 0.1): growth and pFBA at O2 = 0.01 (O2-limited) and O2 = 1
   (physiological order); top O2-consuming reactions; top carbon-supplying uptakes (C atoms/h).
B. SLC7A11 on RPMI M1: cystine[e] reactions and GPRs; KO SLC7A11 with GSH open vs closed.
C. ATP supplement: WT growth of RPMI M1 with and without the ATP uptake opened by the RPMI+HT/adenine/ATP key.
Writes results/w2/review/diag_m3_cys.tsv.
"""
from __future__ import annotations

import re
from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"

import pandas as pd
from cobra.flux_analysis import pfba

from nextgem import media
from nextgem.model import BIOMASS_REACTION
from nextgem.w2_runtime import load_parent

ROOT = Path(__file__).resolve().parents[3]
W2 = ROOT / "results" / "w2"
OUT = Path(__file__).resolve().parent / "diag_m3_cys.tsv"


def carbons(met) -> int:
    m = re.search(r"C(\d*)(?![a-z])", met.formula or "")
    return 0 if not m else int(m.group(1) or 1)


def main() -> None:
    sheets = media.load_sheets()
    rpmi = media.build_medium(media.parse_formulation("RPMI + 10% FBS"), sheets)
    pcs = media.floor_constraints(W2 / "floor" / "floor.tsv", W2 / "floor" / "requirements.tsv",
                                  W2 / "floor" / "floor_mult.tsv")
    rows = []
    for p in ("human1", "human2"):
        m = load_parent(p, solver="glpk")
        m.objective = BIOMASS_REACTION
        pc = pcs[p]
        # ---- A. M3 mechanism
        for o2 in (0.01, 1.0):
            b = media.uptake_bounds(None, "M3", rpmi, pc, f=0.01, calibration={"glucose": 0.5, "glutamine": 0.1,
                                                                             "O2": o2})
            with m:
                m.medium = b
                sol = pfba(m)
                gr = sol.fluxes[BIOMASS_REACTION]
                rows.append(dict(parent=p, block="M3", item=f"O2={o2}", value=gr, detail="growth (pFBA)"))
                cons = []
                for o2m in (mt for mt in m.metabolites if mt.id.startswith("MAM02630") and mt.compartment != "e"):
                    cons += [(-sol.fluxes[r.id] * r.metabolites[o2m], r, o2m.compartment) for r in o2m.reactions
                             if sol.fluxes[r.id] * r.metabolites[o2m] < -1e-9 and not r.id.startswith("MAR0")
                             or (sol.fluxes[r.id] * r.metabolites[o2m] < -1e-9 and len(r.metabolites) > 2)]
                for v, r, comp in sorted(cons, key=lambda t: -t[0])[:8]:
                    rows.append(dict(parent=p, block="M3", item=f"O2={o2} O2[{comp}] consumer {r.id}", value=v,
                                     detail=f"{r.name} | {r.subsystem}"))
                upt = []
                for r in m.exchanges:
                    v = sol.fluxes[r.id]
                    if v < -1e-9:
                        met = list(r.metabolites)[0]
                        upt.append((-v * carbons(met), -v, met.name))
                tot = sum(c for c, _, _ in upt)
                for c, v, n in sorted(upt, reverse=True)[:8]:
                    rows.append(dict(parent=p, block="M3", item=f"O2={o2} uptake {n}", value=v,
                                     detail=f"C-mol/h {c:.3g} ({c / tot:.1%} of carbon uptake)"))
        # ---- B. SLC7A11
        b1 = media.uptake_bounds(None, "M1", rpmi, pc)
        cys_e = next(mt for mt in m.metabolites if mt.name == "cystine" and mt.compartment == "e")
        for r in cys_e.reactions:
            if not r.boundary:
                rows.append(dict(parent=p, block="SLC7A11", item=f"cystine[e] reaction {r.id}", value=None,
                                 detail=f"{r.reaction[:90]} | GPR: {r.gene_name_reaction_rule[:120]}"))
        g = next((x for x in m.genes if x.name == "SLC7A11"), None)
        for gsh_open in (True, False):
            bb = dict(b1)
            if not gsh_open:
                bb.pop(media.EXCHANGES["GSH"][0], None)
            with m:
                m.medium = bb
                wt = m.slim_optimize()
                ko = None
                if g is not None:
                    with m:
                        g.knock_out()
                        ko = m.slim_optimize(error_value=0.0)
                rows.append(dict(parent=p, block="SLC7A11", item=f"RPMI M1, GSH {'open' if gsh_open else 'closed'}",
                                 value=(ko / wt) if ko is not None else None, detail=f"WT {wt:.4f}; KO {ko}"))
        # ---- C. ATP supplement
        for add in ((), ("ATP",), ("adenine", "hypoxanthine", "thymidine"), ("ATP", "adenine", "hypoxanthine",
                                                                              "thymidine")):
            bb = dict(b1)
            for c in add:
                bb[media.EXCHANGES[c][0]] = 1000.0
            with m:
                m.medium = bb
                rows.append(dict(parent=p, block="ATP", item="RPMI M1 + " + ("+".join(add) or "nothing"),
                                 value=m.slim_optimize(), detail="WT growth"))
        del m
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 140)
    print(df.to_string())


if __name__ == "__main__":
    main()
