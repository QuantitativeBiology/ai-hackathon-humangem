"""Diagnostics for the W2 review (label-free, GLPK): which floor component binds at 100 x req_c in Human2, and why
CAD stays lethal in Plasmax although uridine is supplied. Writes results/w2/review/diag_floor_cad.tsv."""
from __future__ import annotations

from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"

import pandas as pd

from nextgem import media
from nextgem.model import BIOMASS_REACTION
from nextgem.w2_runtime import load_parent

ROOT = Path(__file__).resolve().parents[3]
W2 = ROOT / "results" / "w2"
OUT = Path(__file__).resolve().parent / "diag_floor_cad.tsv"


def main() -> None:
    sheets = media.load_sheets()
    pcs = media.floor_constraints(W2 / "floor" / "floor.tsv", W2 / "floor" / "requirements.tsv",
                                  W2 / "floor" / "floor_mult.tsv")
    rows = []
    m = load_parent("human2", solver="glpk")
    m.objective = BIOMASS_REACTION
    pc = pcs["human2"]
    df12 = media.build_medium(media.parse_formulation("DMEM:F12 + 10% FBS"), sheets)
    b = media.uptake_bounds(None, "M1", df12, pc)
    with m:
        m.medium = b
        g_floor = m.slim_optimize()
    with m:
        m.medium = media.uptake_bounds(None, "M1", df12, pc, floor_override=1000.0)
        g_open = m.slim_optimize()
    rows.append(dict(block="floor", item=f"H2 DMEM:F12 M1 floor at {pc.floor_mult:g}x req vs 1000", value=g_floor,
                     detail=f"open {g_open:.6f}; rel {abs(g_open - g_floor) / g_open:.3e}"))
    for c in pc.floor:
        if c in df12.present():
            continue
        bb = dict(b)
        bb[media.EXCHANGES[c][0]] = 1000.0
        with m:
            m.medium = bb
            g = m.slim_optimize()
        rows.append(dict(block="floor", item=f"only {c} raised to 1000", value=g,
                         detail=f"gain vs floor {(g - g_floor) / g_floor:.3e}; bound was {b[media.EXCHANGES[c][0]]:.4g}"))
    # CAD in Plasmax: uridine present; what is blocked?
    plas = media.physiological_medium("Plasmax", sheets)
    bp = media.uptake_bounds(None, "M2", plas, pc, f=0.01)
    cad = next(g for g in m.genes if g.name == "CAD")
    with m:
        m.medium = bp
        wt = m.slim_optimize()
        rows.append(dict(block="CAD", item="CAD reactions (H2)", value=None,
                         detail="; ".join(f"{r.id}: {r.name} [{r.gene_name_reaction_rule}]" for r in cad.reactions)))
        with m:
            cad.knock_out()
            rows.append(dict(block="CAD", item="Plasmax M2 f=0.01 CAD KO ratio", value=m.slim_optimize(error_value=0) / wt,
                             detail=f"uridine bound {bp.get(media.EXCHANGES['uridine'][0])}"))
            for met_name in ("carbamoyl-phosphate", "N-carbamoyl-L-aspartate", "(S)-dihydroorotate", "UMP",
                             "orotate"):
                mets = [x for x in m.metabolites if x.name == met_name]
                for x in mets:
                    with m:
                        sk = m.add_boundary(x, type="sink", reaction_id=f"SKrev_{x.id}", lb=-1000, ub=0)
                        rows.append(dict(block="CAD", item=f"CAD KO + supply {met_name}[{x.compartment}]",
                                         value=m.slim_optimize(error_value=0) / wt, detail=""))
        # is carbamoyl-phosphate needed outside pyrimidines? consumers in WT pFBA
        from cobra.flux_analysis import pfba
        sol = pfba(m)
        for x in (x for x in m.metabolites if x.name == "carbamoyl-phosphate"):
            for r in x.reactions:
                v = sol.fluxes[r.id] * r.metabolites[x]
                if v < -1e-9:
                    rows.append(dict(block="CAD", item=f"carbamoyl-P[{x.compartment}] consumer {r.id}", value=-v,
                                     detail=f"{r.name} | {r.subsystem} | GPR {r.gene_name_reaction_rule}"))
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 200)
    print(df.to_string())


if __name__ == "__main__":
    main()
