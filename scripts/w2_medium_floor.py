"""W2 step 2: biomass micronutrient floor and requirement table req_c, per parent (DESIGN W2; PI decision r2-7).

DESIGN names this script scripts/medium_floor.py; it is scripts/w2_medium_floor.py here (W2 naming).

Floor (drop-one-input on each parent):
  1. E = Ham's GR inputs whose removal from M0 stops WT growth (growth <= 1e-6 x M0 growth).
  2. For each on-disk vendor sheet S (M1 availability bounds, no floor): candidates = E not supplied by S.
     If S + candidates grows, candidates are pruned greedily (sorted order) while growth persists -> needed(S).
  3. floor = union over the four sheets of needed(S); cholesterol is never in the floor (DESIGN).
Requirements on the reference medium (RPMI sheet M1 + floor open at 1000, every other input as in the sheet):
  - floor (true auxotroph) components: req_c = FVA-min uptake of c at WT growth (biomass fixed at WT x (1 - 1e-9)).
  - synthesisable serum components: req_c = uptake of c in a pFBA solution at >= 99% of WT growth with c open at
    1000 (cobra.flux_analysis.pfba, fraction_of_optimum = 0.99).
Writes results/w2/floor/{floor.tsv, requirements.tsv, drop_one_M0.tsv, floor_by_sheet.tsv}.

Usage: .venv/bin/python scripts/w2_medium_floor.py
"""

from __future__ import annotations

import math

import pandas as pd
from cobra.flux_analysis import pfba

from nextgem import media
from nextgem.model import BIOMASS_REACTION, hams_medium
from nextgem.w2_runtime import PARENTS, Manifest, atomic_write_tsv, load_parent, w2_dir

EX2CANON = {v[0]: k for k, v in media.EXCHANGES.items()}
NO_GROWTH_REL = 1e-6
FIX_REL = 1e-9


def growth(model, bounds: dict[str, float]) -> float:
    with model:
        model.medium = bounds
        v = model.slim_optimize(error_value=0.0)
    return 0.0 if (v is None or math.isnan(v)) else max(v, 0.0)


def sheet_medium(base: str, sheets) -> media.Medium:
    return media.build_medium(media.parse_formulation(f"{base} + 10% FBS"), sheets, key=base)


def main() -> None:
    out = w2_dir() / "floor"
    out.mkdir(parents=True, exist_ok=True)
    comp = media.canonical_composition(media.read_compositions())
    sheets = media.sheet_media(comp)
    no_floor = media.ParentConstraints(floor=(), req={})
    drop_rows, sheet_rows, floor_rows, req_rows = [], [], [], []
    with Manifest("scripts/w2_medium_floor.py", solver="gurobi", workers=1) as man:
        for p in PARENTS:
            m = load_parent(p)
            m.objective = BIOMASS_REACTION
            m0 = hams_medium(m)
            g0 = growth(m, m0)
            ess = []
            for ex in sorted(m0):
                g = growth(m, {k: v for k, v in m0.items() if k != ex})
                drop_rows.append(dict(parent=p, exchange=ex, canonical=EX2CANON.get(ex, ""),
                                      growth_without=g, rel=g / g0, essential=g <= NO_GROWTH_REL * g0))
                if g <= NO_GROWTH_REL * g0:
                    ess.append(ex)
            needed_all: dict[str, list[str]] = {}
            for base in sorted(media.SHEET_FILES):
                med = sheet_medium(base, sheets)
                b = media.uptake_bounds(m, "M1", med, no_floor)
                g_bare = growth(m, b)
                cands = sorted(e for e in ess if e not in b)
                b_all = {**b, **{e: media.UPTAKE_OPEN for e in cands}}
                g_all = growth(m, b_all)
                needed = list(cands)
                if g_all > NO_GROWTH_REL * g0:
                    for e in sorted(cands):
                        trial = [x for x in needed if x != e]
                        if growth(m, {**b, **{x: media.UPTAKE_OPEN for x in trial}}) > NO_GROWTH_REL * g0:
                            needed = trial
                sheet_rows.append(dict(parent=p, sheet=base, growth_M1_no_floor=g_bare,
                                       candidates=",".join(EX2CANON[e] for e in cands),
                                       growth_with_candidates=g_all,
                                       needed=",".join(EX2CANON[e] for e in needed)))
                for e in needed:
                    needed_all.setdefault(EX2CANON[e], []).append(base)
            floor = tuple(sorted(c for c in needed_all if c != "cholesterol"))

            # ---- requirements on the reference medium (RPMI + floor open at 1000)
            ref = sheet_medium("RPMI", sheets)
            ref_b = media.uptake_bounds(m, "M1", ref, media.ParentConstraints(floor=floor, req={}),
                                        floor_override=media.UPTAKE_OPEN)
            with m:
                m.medium = ref_b
                wt = m.slim_optimize()
                m.reactions.get_by_id(BIOMASS_REACTION).lower_bound = wt * (1 - FIX_REL)
                for c in floor:
                    ex = m.reactions.get_by_id(media.EXCHANGES[c][0])
                    m.objective = {ex: 1}
                    m.objective_direction = "max"
                    vmax = m.slim_optimize()
                    req = max(0.0, -vmax)
                    req_rows.append(dict(parent=p, component=c, kind="floor", rule="FVA-min uptake at WT growth",
                                         req_c=req, ref_wt_growth=wt, n_sheets_lacking=len(needed_all[c]),
                                         sheets_lacking=",".join(needed_all[c])))
                    floor_rows.append(dict(parent=p, component=c, exchange=ex.id, req_c=req,
                                           floor_bound=media.FLOOR_MULT * req, floor_mult=media.FLOOR_MULT,
                                           needed_by_sheets=",".join(needed_all[c])))
            for c in media.SERUM_F_SCALED:
                ex = media.EXCHANGES[c][0]
                with m:
                    m.medium = {**ref_b, ex: media.UPTAKE_OPEN}
                    m.objective = BIOMASS_REACTION
                    wt_c = m.slim_optimize()
                    sol = pfba(m, fraction_of_optimum=media.PFBA_FRACTION)
                    up = max(0.0, -sol.fluxes[ex])
                req_rows.append(dict(parent=p, component=c, kind="serum_synthesisable",
                                     rule=f"pFBA uptake at >= {media.PFBA_FRACTION:.0%} WT growth, c open at 1000",
                                     req_c=up, ref_wt_growth=wt_c, flagged=(up == 0.0)))
            print(p, "M0", round(g0, 4), "essential inputs", [EX2CANON.get(e, e) for e in ess])
            print(p, "floor", floor)
            del m
        fl, rq = pd.DataFrame(floor_rows), pd.DataFrame(req_rows)
        for name, df in {"drop_one_M0.tsv": pd.DataFrame(drop_rows), "floor_by_sheet.tsv": pd.DataFrame(sheet_rows),
                         "floor.tsv": fl, "requirements.tsv": rq}.items():
            man.add(atomic_write_tsv(df, out / name))
        print(rq.to_string())
        print(pd.DataFrame(sheet_rows).to_string())


if __name__ == "__main__":
    main()
