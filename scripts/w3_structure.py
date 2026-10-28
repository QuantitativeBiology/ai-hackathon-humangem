"""W3 label-free structure items (DESIGN §9 W3; M3): NGAM reaction check, B-loo components (V-13), B-add metabolite
IDs (V-14), D-BS / D-ER demand feasibility log, variant definitions for the KO sweeps, candidate FROZEN rows.

Media: M0 = Ham's GR inputs at 1000 (``set_hams_medium``). The "documented variant" is classic RPMI-1640 under
M2(f = 1), gate open (``nextgem.media.uptake_bounds``, the frozen W2 library; A9). W2 has not run its training gate,
so no "best" variant exists yet; RPMI is the modal medium of the attribution set (557 / 835 lines) and f = 1 the
full serum model. This is a label-free proxy (README "ambiguities"); its bounds signature equals W2's sweep.

Outputs (results/w3/structure/):
  ngam_check.tsv, bloo_components.tsv, badd_metabolites.tsv, demand_products.tsv, demand_feasibility_log.tsv,
  variants_<parent>.json (Variant specs per medium, read by scripts/w3_sweeps.py), media_<parent>.json (bounds),
  structure_summary.json, frozen_candidates.tsv.

Usage: .venv/bin/python scripts/w3_structure.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w3_common import OUT_DIR, atomic_write_text, atomic_write_tsv, current, run_main  # noqa: E402

SDIR = OUT_DIR / "structure"
VARIANT_MEDIUM = ("RPMI", "M2", 1.0)  # (medium_key, ladder step, f)


def media_bounds(parents=("human1", "human2")) -> tuple[dict, dict]:
    """{parent: {"M0": bounds, "RPMI_M2f1": bounds}} plus the W2 sweep ids they correspond to."""
    from nextgem import media, w2_engine
    from nextgem.w2_runtime import w2_dir

    lines, sheets, doc, phys, sens, pcs, exchanges, hams = w2_engine.sweep_context()
    labels = pd.read_csv(w2_dir() / "ko" / "sweep_labels.tsv", sep="\t")
    key, step, f = VARIANT_MEDIUM
    out, ids = {}, {}
    for p in parents:
        b = media.uptake_bounds(None, step, doc[key], pcs[p], f=f, gate_open=True)
        sig = media.bounds_signature(b)
        row = labels[(labels.parent == p) & (labels.variant == step) & (labels.medium_key == key) & (labels.f == f)
                     & (labels.gate_open)]
        if len(row) != 1 or row.signature.iloc[0] != sig:
            raise RuntimeError(f"{p}: variant bounds signature {sig} does not match the W2 sweep {row.signature.tolist()}")
        m0 = dict(hams[p])
        m0row = labels[(labels.parent == p) & (labels.variant == "M0")]
        if m0row.signature.iloc[0] != media.bounds_signature(m0):
            raise RuntimeError(f"{p}: M0 signature mismatch")
        out[p] = {"M0": m0, "RPMI_M2f1": b}
        ids[p] = {"M0": m0row.sweep_id.iloc[0], "RPMI_M2f1": row.sweep_id.iloc[0]}
    return out, ids


def max_flux(model, bounds, rid_or_rxn, biomass_free=True) -> float:
    """Max flux of a reaction under ``bounds`` with biomass free (lb 0)."""
    from nextgem import w3

    with model:
        model.medium = bounds
        model.objective = rid_or_rxn
        st = w3.solve_status(model)
        return float(model.solver.objective.value) if st == "optimal" else 0.0


def wt_growth(model, bounds, variant) -> tuple[float, str]:
    from nextgem import w3

    with model:
        model.medium = bounds
        variant.apply(model)
        model.objective = w3.BIOMASS
        return w3.solve_growth(model)


def main() -> int:
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import w3
    from nextgem.model import GEMS, load_model
    from raven_toolbox.tasks.check import task_name_maps

    t0 = time.time()
    SDIR.mkdir(parents=True, exist_ok=True)
    bounds, sweep_ids = media_bounds()
    ngam_rows, bloo_rows, badd_rows, prod_rows, dlog_rows, frozen = [], [], [], [], [], []
    summary = {"variant_medium": dict(zip(("medium_key", "step", "f"), VARIANT_MEDIUM)), "w2_sweep_ids": sweep_ids}
    for p in ("human1", "human2"):
        m = load_model(p)
        m.objective = w3.BIOMASS
        variants: dict[str, list] = {}
        # ---------------------------------------------------------------- NGAM (MAR03964)
        r = m.reactions.get_by_id(w3.NGAM_RXN)
        eq = r.build_reaction_string(use_metabolite_names=True)
        ok = (eq.replace(" ", "") == "ATP+H2O-->ADP+H++Pi" and r.compartments == {"c"} and r.bounds == (0.0, 1000.0))
        for med, b in bounds[p].items():
            mx = max_flux(m, b, r.id)
            ngam_rows.append({"parent": p, "medium": med, "reaction": r.id, "name": r.name, "equation": eq,
                              "compartments": ",".join(sorted(r.compartments)), "bounds": str(r.bounds),
                              "gpr": r.gene_reaction_rule, "is_atp_hydrolysis_c": ok, "max_flux_biomass_free": mx})
        if not ok:
            raise RuntimeError(f"{p}: {w3.NGAM_RXN} is not cytosolic ATP hydrolysis with bounds [0, 1000]: {eq}")
        # ---------------------------------------------------------------- B-loo components (V-13)
        bio = m.reactions.get_by_id(w3.BIOMASS)
        comps = []
        for met, c in sorted(bio.metabolites.items(), key=lambda x: x[0].id):
            if c < 0 and met.id not in w3.BIOMASS_EXCLUDED:
                comps.append(("direct:" + met.name, "remove", [met.id], c))
        pool = m.reactions.get_by_id(w3.COFACTOR_POOL)
        pool_sp = {x: c for x, c in pool.metabolites.items() if c < 0}
        assigned = {}
        for gname, rx in w3.POOL_GROUPS.items():
            mem = sorted(x.id for x in pool_sp if rx.search(x.name or ""))
            for x in mem:
                if x in assigned:
                    raise RuntimeError(f"{p}: {x} in two pool groups ({assigned[x]}, {gname})")
                assigned[x] = gname
            comps.append(("pool:" + gname, "pool_remove", mem, sum(pool_sp[m.metabolites.get_by_id(x)] for x in mem)))
        for cname, (pool_id, nm, cp) in w3.EXTRA_BLOO.items():
            hit = w3.met_by_name(m, nm, cp)
            if hit is None or isinstance(hit, list) or hit not in m.reactions.get_by_id(pool_id).metabolites:
                raise RuntimeError(f"{p}: {nm}[{cp}] not a single species of {pool_id}")
            comps.append((cname, "pool_remove", [hit.id], m.reactions.get_by_id(pool_id).get_coefficient(hit.id), pool_id))
        unassigned = sorted(f"{x.id} {x.name}" for x in pool_sp if x.id not in assigned)
        summary.setdefault("pool_species_not_in_any_group", {})[p] = unassigned
        for comp in comps:
            name, kind, mets, coef = comp[:4]
            pool_id = comp[4] if len(comp) > 4 else (w3.COFACTOR_POOL if kind == "pool_remove" else w3.BIOMASS)
            if not mets:
                bloo_rows.append({"parent": p, "component": name, "kind": kind, "n_species": 0, "species": "",
                                  "names": "", "coefficient_sum": 0, "status": "skipped_absent"})
                continue
            params = {"mets": mets, "pool": pool_id} if kind == "pool_remove" else {"mets": mets}
            v = w3.Variant(f"B-loo:{name.split(':', 1)[1]}", kind, params)
            g0, st0 = wt_growth(m, bounds[p]["M0"], v)
            bloo_rows.append({"parent": p, "component": name, "kind": kind, "n_species": len(mets),
                              "species": ";".join(mets),
                              "names": ";".join(m.metabolites.get_by_id(x).name for x in mets),
                              "coefficient_sum": coef, "wt_growth_M0": g0, "status": "ok" if st0 == "optimal" else st0,
                              "reaction": pool_id,
                              "origin": ("label-free structural finding (M3_W5_review §4), added before W3 scoring"
                                         if name in w3.EXTRA_BLOO else "DESIGN W3 B-loo definition")})
            variants.setdefault("bloo", []).append(v.to_json())
        # T10 biomass of the final list: cofactor pool and glycogen terms zeroed
        t10 = w3.Variant("L-T10:biomass", "remove", {"mets": ["MAM10012c", "MAM03161c"]})
        variants["t10"] = [t10.to_json()]
        # ---------------------------------------------------------------- B-add (V-14)
        in_bio = set(bio.metabolites) | set(pool.metabolites) | set(m.reactions.MAR10064.metabolites) | set(m.reactions.MAR10062.metabolites)
        in_bio_ids = {x.id for x in in_bio}
        add_ok = {}
        for comp, spec in w3.BADD_SPECS.items():
            if spec["kind"] == "sink":
                found = None
                for nm, cp in spec["names"]:
                    h = w3.met_by_name(m, nm, cp)
                    if h is not None and not isinstance(h, list):
                        found = h
                        break
                if found is None:
                    badd_rows.append({"parent": p, "component": comp, "kind": "sink", "metabolites": "", "status": "skipped_absent_in_parent"})
                    continue
                if found.id in in_bio_ids:
                    badd_rows.append({"parent": p, "component": comp, "kind": "sink", "metabolites": found.id,
                                      "names": found.name, "status": "skipped_already_in_biomass"})
                    continue
                st = {found.id: -1.0}
            else:
                st, miss = {}, []
                for nm, cp in spec["reactants"]:
                    h = w3.met_by_name(m, nm, cp)
                    (miss.append(nm) if h is None or isinstance(h, list) else st.__setitem__(h.id, st.get(h.id, 0) - 1.0))
                for nm, cp in spec["products"]:
                    h = w3.met_by_name(m, nm, cp)
                    (miss.append(nm) if h is None or isinstance(h, list) else st.__setitem__(h.id, st.get(h.id, 0) + 1.0))
                if miss:
                    badd_rows.append({"parent": p, "component": comp, "kind": "turnover", "metabolites": "",
                                      "status": f"skipped_missing:{miss}"})
                    continue
            res = {}
            for coef_name, coef in (("1e-4", w3.BADD_COEF), ("1e-3", w3.BADD_COEF_SENS)):
                v = w3.Variant(f"B-add:{comp}@{coef_name}", "add", {"stoich": {k: c * coef for k, c in st.items()}})
                g, s = wt_growth(m, bounds[p]["M0"], v)
                res[coef_name] = (g, s, v)
            ok_ = all(s == "optimal" and g > 1e-6 for g, s, _ in res.values())
            badd_rows.append({"parent": p, "component": comp, "kind": spec["kind"], "metabolites": json.dumps(st),
                              "names": ";".join(m.metabolites.get_by_id(k).name for k in st),
                              "wt_growth_M0_1e-4": res["1e-4"][0], "wt_growth_M0_1e-3": res["1e-3"][0],
                              "status": "ok" if ok_ else "skipped_not_producible"})
            if ok_:
                add_ok[comp] = st
                for coef_name, (_, _, v) in res.items():
                    variants.setdefault("badd" if coef_name == "1e-4" else "badd_sens", []).append(v.to_json())
        for coef_name, coef, key in (("1e-4", w3.BADD_COEF, "badd"), ("1e-3", w3.BADD_COEF_SENS, "badd_sens")):
            tot = {}
            for st in add_ok.values():
                for k, c in st.items():
                    tot[k] = tot.get(k, 0.0) + c * coef
            v = w3.Variant(f"B-add:ALL@{coef_name}", "add", {"stoich": tot})
            g, s = wt_growth(m, bounds[p]["M0"], v)
            badd_rows.append({"parent": p, "component": f"ALL@{coef_name}", "kind": "all", "metabolites": json.dumps(tot),
                              "wt_growth_M0_1e-4" if coef_name == "1e-4" else "wt_growth_M0_1e-3": g,
                              "status": "ok" if s == "optimal" and g > 1e-6 else s})
            variants.setdefault(key, []).append(v.to_json())
        # ---------------------------------------------------------------- D-BS / D-ER demands per medium
        tasks = w3.load_tasks(p, "essential")
        base = m.copy()
        maps = task_name_maps(base)
        prods = []
        for _, t in tasks[tasks.category == "BS"].iterrows():
            for tok, lb, ub in t.task.outputs:
                if lb > 0:
                    ids = maps[0].get(tok.upper(), [])
                    if len(ids) != 1:
                        raise RuntimeError(f"{p}: BS output {tok} resolves to {ids}")
                    prods.append((t.key, t.description, tok, ids[0]))
        seen = set()
        bs_specs = []
        for key, desc, tok, mid in prods:
            if mid in seen:
                continue
            seen.add(mid)
            bs_specs.append((f"DM_w3_{mid}", {mid: -1.0}, tok, key, desc))
        er_specs = []
        for name, (reac, prod) in w3.DER_SPECS.items():
            st = {}
            for nm, cp in reac:
                st[w3.met_by_name(m, nm, cp).id] = st.get(w3.met_by_name(m, nm, cp).id, 0) - 1.0
            for nm, cp in prod:
                st[w3.met_by_name(m, nm, cp).id] = st.get(w3.met_by_name(m, nm, cp).id, 0) + 1.0
            er_specs.append((f"DM_w3_{name}", st, name, "ER", name))
        for med, b in bounds[p].items():
            for mode, specs in (("D-BS", bs_specs), ("D-ER", er_specs)):
                kept = []
                for rid, st, tok, key, desc in specs:
                    prod_rows.append({"parent": p, "medium": med, "mode": mode, "demand": rid, "token": tok,
                                      "stoich": json.dumps(st), "task": key, "task_description": desc})
                    with m:
                        m.medium = b
                        rx = cobra.Reaction(rid, lower_bound=0.0, upper_bound=1000.0)
                        rx.add_metabolites({m.metabolites.get_by_id(k): c for k, c in st.items()})
                        m.add_reactions([rx])
                        m.objective = rid
                        stt = w3.solve_status(m)
                        mx = float(m.solver.objective.value) if stt == "optimal" else 0.0
                    lb = w3.DEMAND_FRACTION * mx
                    row = {"parent": p, "medium": med, "mode": mode, "demand": rid, "token": tok, "max_flux": mx,
                           "lb": lb}
                    if not mx > 1e-9:
                        row.update(status="dropped_not_producible", wt_growth_after=None)
                        dlog_rows.append(row)
                        continue
                    trial = kept + [(rid, st, lb)]
                    g, s = wt_growth(m, b, w3.Variant(mode, "demands", {"demands": trial}))
                    if s == "optimal" and g > 1e-6:
                        kept = trial
                        row.update(status="kept", wt_growth_after=g)
                    else:
                        row.update(status=f"dropped_joint_infeasible({s})", wt_growth_after=g)
                    dlog_rows.append(row)
                g0, _ = wt_growth(m, b, w3.Variant("T0"))
                gD, sD = wt_growth(m, b, w3.Variant(mode, "demands", {"demands": kept}))
                summary.setdefault("demands", {}).setdefault(p, {})[f"{mode}|{med}"] = {
                    "n_kept": len(kept), "n_candidates": len(specs), "wt_growth_T0": g0, "wt_growth_with_demands": gD}
                variants.setdefault(f"{mode}|{med}", []).append(w3.Variant(mode, "demands", {"demands": kept}).to_json())
            # O-NGAM_k on this medium
            mx = [x for x in ngam_rows if x["parent"] == p and x["medium"] == med][0]["max_flux_biomass_free"]
            for k in w3.NGAM_LEVELS:
                v = w3.Variant(f"O-NGAM_{k:g}", "fix", {"fluxes": {w3.NGAM_RXN: k * mx}})
                g, s = wt_growth(m, b, v)
                summary.setdefault("ngam", {}).setdefault(p, {})[f"{k:g}|{med}"] = {"fixed_flux": k * mx, "wt_growth": g, "status": s}
                variants.setdefault(f"ngam|{med}", []).append(v.to_json())
        atomic_write_text(json.dumps(variants, indent=1), SDIR / f"variants_{p}.json")
        atomic_write_text(json.dumps(bounds[p], indent=1, sort_keys=True), SDIR / f"media_{p}.json")
        current().input(GEMS[p].sbml, GEMS[p].essential_tasks)
        current().output(SDIR / f"variants_{p}.json", SDIR / f"media_{p}.json")
        print(f"[{p}] done {time.time() - t0:.0f} s", flush=True)

    tabs = {"ngam_check.tsv": ngam_rows, "bloo_components.tsv": bloo_rows, "badd_metabolites.tsv": badd_rows,
            "demand_products.tsv": prod_rows, "demand_feasibility_log.tsv": dlog_rows}
    for f, rows in tabs.items():
        atomic_write_tsv(pd.DataFrame(rows), SDIR / f)
        current().output(SDIR / f)
    # ---------------------------------------------------------------- candidate FROZEN rows (PI adds them)
    day = time.strftime("%Y-%m-%d")
    bl = pd.DataFrame(bloo_rows)
    for _, r in bl.iterrows():
        frozen.append({"item": f"w3.bloo.{r.parent}.{r.component}", "kind": "value", "path": "",
                       "sha256_or_value": r.species, "frozen_on": day, "amendment": "",
                       "note": f"DESIGN W3 B-loo (V-13); {r.kind}; status {r.status}"})
    ba = pd.DataFrame(badd_rows)
    for _, r in ba.iterrows():
        frozen.append({"item": f"w3.badd.{r.parent}.{r.component}", "kind": "value", "path": "",
                       "sha256_or_value": r.metabolites, "frozen_on": day, "amendment": "",
                       "note": f"DESIGN W3 B-add (V-14); {r.kind}; status {r.status}"})
    atomic_write_tsv(pd.DataFrame(frozen), SDIR / "frozen_candidates_structure.tsv")
    summary["wall_s"] = round(time.time() - t0, 1)
    atomic_write_text(json.dumps(summary, indent=1, default=str), SDIR / "structure_summary.json")
    current().output(SDIR / "frozen_candidates_structure.tsv", SDIR / "structure_summary.json")
    print(json.dumps({k: v for k, v in summary.items() if k != "w2_sweep_ids"}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(run_main(main, solver="glpk", workers=1))
