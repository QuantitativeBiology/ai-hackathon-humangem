"""Re-check (A9d): classic RPMI 1640 in the v2 library vs Cantor 2017 Table S1 'Classic synthetic media' (Gibco 11875),
and the effect of keeping the req_c reference medium on the on-disk ATCC sheet (label-free, GLPK).

1. canonical composition of v2 'RPMI' vs the Cantor 11875 column (uM), per metabolite, rel. diff;
2. req_c recomputed with W2's rules (scripts/w2_medium_floor.py: FVA-min at WT for floor items; pFBA >= 99% with c
   open for serum items) on reference = RPMI-ATCC (as frozen) and reference = classic RPMI; relative difference.
Writes results/w2/review/rpmi_classic_check.tsv.
"""
from __future__ import annotations

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
OUT = Path(__file__).resolve().parent / "rpmi_classic_check.tsv"
CANTOR2CANON = {"Glucose": "glucose", "Arginine": "arginine", "Asparagine": "asparagine", "Aspartate": "aspartate",
                "Cystine": "cystine", "Glutamate": "glutamate", "Glutamine": "glutamine", "Glycine": "glycine",
                "Histidine": "histidine", "Hydroxyproline": "hydroxyproline", "Isoleucine": "isoleucine",
                "Leucine": "leucine", "Lysine": "lysine", "Methionine": "methionine",
                "Phenylalanine": "phenylalanine", "Proline": "proline", "Serine": "serine", "Threonine": "threonine",
                "Tryptophan": "tryptophan", "Tyrosine": "tyrosine", "Valine": "valine", "Biotin": "biotin",
                "Choline": "choline", "Folate": "folate", "myo-Inositol": "inositol", "Niacinamide": "nicotinamide",
                "p-Aminobenzoate": "4-aminobenzoate", "Pyridoxine": "pyridoxine", "Riboflavin": "riboflavin",
                "Thiamine": "thiamin", "Vitamin B-12": "cobalamin", "MgSO4": "sulfate", "KCl": "K+",
                "NaHCO3": "HCO3-", "Na2HPO4": "Pi", "Ca(NO3)2.4H2O": "Ca2+"}


def main() -> None:
    rows = []
    comp = media.canonical_composition(media.read_compositions())
    v2 = comp[comp.medium == "RPMI"].set_index("canonical").uM
    cantor = media.parse_cantor_classic(media.media_dir() / media.HPLM_FILE)
    for name, canon in CANTOR2CANON.items():
        c = cantor[name]
        v = float(v2.get(canon, 0.0))
        if name == "Pantothenate":
            continue
        rel = abs(v - c) / c if c else (0.0 if v == 0 else float("inf"))
        rows.append(dict(check="composition_vs_cantor_11875", item=name, cantor_uM=c, v2_uM=v, rel_diff=rel))
    # pantothenate: Cantor lists the Ca-salt molarity (0.524); v2 canonical pantothenate = 2 x salt
    rows.append(dict(check="composition_vs_cantor_11875", item="Pantothenate (Ca salt x2)",
                     cantor_uM=cantor["Pantothenate"] * 2, v2_uM=float(v2.get("pantothenate", 0)),
                     rel_diff=abs(float(v2.get("pantothenate", 0)) - 2 * cantor["Pantothenate"]) / (2 * cantor["Pantothenate"])))
    extra = sorted(set(v2.index) - set(CANTOR2CANON.values()) - {"pantothenate", "chloride", "Na+", "Mg2+",
                                                                 "nitrate", "<dye>"})
    rows.append(dict(check="composition_vs_cantor_11875", item="in v2 RPMI, not in Cantor table", v2_uM=None,
                     cantor_uM=None, rel_diff=None, note=",".join(extra)))

    # req_c on the two reference media
    sheets = media.load_sheets()
    fl = pd.read_csv(W2 / "floor" / "floor.tsv", sep="\t")
    rq = pd.read_csv(W2 / "floor" / "requirements.tsv", sep="\t")
    p_rpmi = media.parse_formulation("RPMI + 10% FBS")
    for parent in ("human1", "human2"):
        m = load_parent(parent, solver="glpk")
        m.objective = BIOMASS_REACTION
        floor = tuple(sorted(fl[fl.parent == parent].component))
        frozen = rq[rq.parent == parent].set_index("component").req_c
        for ref in (media.RPMI_ATCC, "RPMI"):
            med = media.build_medium(p_rpmi, sheets, rpmi_base=ref)
            ref_b = media.uptake_bounds(None, "M1", med, media.ParentConstraints(floor=floor, req={}),
                                        floor_override=media.UPTAKE_OPEN)
            with m:
                m.medium = ref_b
                wt = m.slim_optimize()
                m.reactions.get_by_id(BIOMASS_REACTION).lower_bound = wt * (1 - 1e-9)
                for c in floor:
                    ex = m.reactions.get_by_id(media.EXCHANGES[c][0])
                    m.objective = {ex: 1}
                    m.objective_direction = "max"
                    req = max(0.0, -m.slim_optimize())
                    rows.append(dict(check=f"req_c_ref_{ref}", parent=parent, item=c, req_c=req, wt=wt,
                                     frozen=frozen[c], rel_diff=abs(req - frozen[c]) / frozen[c]))
            for c in media.SERUM_F_SCALED:
                ex = media.EXCHANGES[c][0]
                with m:
                    m.medium = {**ref_b, ex: media.UPTAKE_OPEN}
                    m.objective = BIOMASS_REACTION
                    sol = pfba(m, fraction_of_optimum=media.PFBA_FRACTION)
                    up = max(0.0, -sol.fluxes[ex])
                rows.append(dict(check=f"req_c_ref_{ref}", parent=parent, item=c, req_c=up, wt=wt,
                                 frozen=frozen[c], rel_diff=abs(up - frozen[c]) / frozen[c]))
        del m
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(df.to_string())


if __name__ == "__main__":
    main()
