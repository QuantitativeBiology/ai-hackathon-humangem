"""W2 step 3 (M3 pilots): per-medium WT-growth prechecks, floor non-binding assertion, req_c != 0, M-open control,
and the M3 growth-rate calibration (DESIGN W2 "Prechecks", §12 row W2). Label-free.

Inputs: results/w2/media_library/line_media.tsv, results/w2/floor/{floor,requirements}.tsv, data/media/*.
Outputs (results/w2/prechecks/):
  sweep_labels.tsv        every ladder label (parent x variant x medium x f x gate) -> sweep signature
  precheck_wt.tsv         WT growth per distinct sweep; pass = WT > 1e-3 x M0 growth of the parent
  floor_nonbinding.tsv    WT growth with the floor at 10 x req_c vs at 1000, per medium x variant (x f), gate open;
                          parents with any binding case are re-tested at 100 x req_c (DESIGN) -> floor/floor_mult.tsv
  req_nonzero.tsv         req_c per serum/floor component and parent; fails unless > 0 or flagged
  calibration_grid.tsv    WT growth on RPMI M2(f) with glucose, glutamine and O2 bounded on a grid (mmol/gDW/h)
  calibration.tsv         calibration verdict (uses data/media/calibration.tsv if present; else 'missing')
  precheck_summary.json

Usage: .venv/bin/python scripts/w2_prechecks.py [--workers 4]
"""

from __future__ import annotations

import argparse
import itertools
import json

import pandas as pd

from nextgem import media, paths
from nextgem import w2_engine as eng
from nextgem.w2_runtime import PARENTS, Manifest, atomic_write_tsv, w2_dir

CAL_GLC = (0.05, 0.2, 0.5, 2.0)
CAL_GLN = (0.01, 0.1)
CAL_O2 = (0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0)
CAL_F = (0.01, 1.0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    out = w2_dir() / "prechecks"
    out.mkdir(parents=True, exist_ok=True)
    with Manifest("scripts/w2_prechecks.py", solver="gurobi", workers=args.workers) as man:
        lines, sheets, doc, phys, sens, pcs, exchanges, hams = eng.sweep_context(use_mult=False)
        allmeds = {**doc, **phys, **sens}
        cal_file = paths.data_root() / "media" / "calibration.tsv"

        def floor_tasks(parents, pcs_):
            tasks, meta = [], []
            for p in parents:
                for key, med in allmeds.items():
                    for v, f in [("M1", None)] + [("M2", f) for f in media.SERUM_F_LEVELS]:
                        if key in phys and v == "M1":
                            continue
                        for mode, ov in (("floor", None), ("open", media.UPTAKE_OPEN)):
                            b = media.uptake_bounds(None, v, med, pcs_[p], f=f, floor_override=ov)
                            tasks.append((p, f"floor|{p}|{key}|{v}|{f}|{mode}", b))
                        meta.append(dict(parent=p, medium_key=key, variant=v, f=f, floor_mult=pcs_[p].floor_mult))
            return tasks, meta

        def floor_table(res, meta):
            rows = []
            for mt in meta:
                base = f"floor|{mt['parent']}|{mt['medium_key']}|{mt['variant']}|{mt['f']}"
                gf, gop = res[base + "|floor"]["growth"], res[base + "|open"]["growth"]
                rel = abs(gf - gop) / max(abs(gop), 1e-12)
                rows.append(dict(**mt, growth_floor=gf, growth_floor_open=gop, rel_diff=rel,
                                 nonbinding=rel <= media.FLOOR_NONBINDING_RTOL))
            return pd.DataFrame(rows)

        with eng.make_pool(args.workers) as pool:
            # phase A: floor at FLOOR_MULT x req_c vs 1000 (gate open, every medium x M1 / M2(f))
            t, mt = floor_tasks(PARENTS, pcs)
            flt = floor_table({r["tag"]: r for r in pool.map(eng.wt_with_override_task, t, chunksize=8)}, mt)
            # DESIGN: if it binds, the bound is raised to 100 x req and logged (applied per parent, all media)
            failing = sorted(flt[~flt.nonbinding].parent.unique())
            for p in failing:
                pcs[p].floor_mult = media.FLOOR_MULT_FALLBACK
            if failing:
                t, mt = floor_tasks(failing, pcs)
                flt = pd.concat([flt, floor_table({r["tag"]: r for r in pool.map(eng.wt_with_override_task, t,
                                                                                  chunksize=8)}, mt)])
            mult = pd.DataFrame([dict(parent=p, floor_mult=pcs[p].floor_mult,
                                      raised_from=media.FLOOR_MULT if p in failing else "",
                                      still_binding_cases=int(((flt.parent == p) & (flt.floor_mult == pcs[p].floor_mult)
                                                               & ~flt.nonbinding).sum()))
                                 for p in PARENTS])
            man.add(atomic_write_tsv(mult, w2_dir() / "floor" / "floor_mult.tsv"))

            sweeps = eng.build_sweeps(exchanges, hams, pcs, doc, phys, sens)
            labels = eng.labels_frame(sweeps)
            cal_tasks = []
            for p in PARENTS:
                for f, g, q, o in itertools.product(CAL_F, CAL_GLC, CAL_GLN, CAL_O2):
                    b = media.uptake_bounds(None, "M3", doc["RPMI"], pcs[p], f=f,
                                            calibration={"glucose": g, "glutamine": q, "O2": o})
                    cal_tasks.append((p, f"cal|{p}|{f}|{g}|{q}|{o}", b))
            wt = list(pool.map(eng.wt_task, [(s.parent, s.sweep_id, s.bounds) for s in sweeps], chunksize=4))
            cal = list(pool.map(eng.wt_with_override_task, cal_tasks, chunksize=16))

        wt = pd.DataFrame(wt)
        m0 = {p: wt.set_index("sweep_id").loc[next(s.sweep_id for s in sweeps if s.parent == p and any(
            lab["variant"] == "M0" for lab in s.labels)), "wt_growth"] for p in PARENTS}
        wt["m0_growth"] = wt.parent.map(m0)
        wt["rel_to_M0"] = wt.wt_growth / wt.m0_growth
        wt["precheck_pass"] = (wt.wt_status == "optimal") & (wt.wt_growth > media.PRECHECK_REL * wt.m0_growth)
        wt["n_open"] = wt.sweep_id.map({s.sweep_id: len(s.bounds) for s in sweeps})
        lab = labels.merge(wt[["sweep_id", "wt_growth", "rel_to_M0", "precheck_pass"]], on="sweep_id")

        rq = pd.read_csv(w2_dir() / "floor" / "requirements.tsv", sep="\t")
        rq["flagged"] = rq.get("flagged", False)
        rq["flagged"] = rq.flagged.astype("boolean").fillna(False)
        rq["req_nonzero_pass"] = (rq.req_c > 0) | rq.flagged

        calg = pd.DataFrame(cal)
        calg[["_", "parent_", "f", "glucose", "glutamine", "O2"]] = calg.tag.str.split("|", expand=True)
        calg = calg.drop(columns=["_", "parent_", "tag"]).astype({"f": float, "glucose": float, "glutamine": float,
                                                                  "O2": float})
        lo, hi = media.CALIBRATION_RANGE
        calg["in_window"] = (calg.growth >= lo) & (calg.growth <= hi)
        if cal_file.exists():
            verdict = "calibration.tsv present: evaluate"  # not reached on this data snapshot
        else:
            verdict = ("M3 UNCALIBRATED: data/media/calibration.tsv (literature-anchored glucose, glutamine and O2 "
                       "rates with citations, DESIGN W2) does not exist; M3 is excluded from selection until it is "
                       "fixed. The grid shows which rate combinations would place WT growth in 0.01-0.06 h^-1.")
        calv = pd.DataFrame([dict(parent=p, f=f, n_grid=int(((calg.parent == p) & (calg.f == f)).sum()),
                                  n_in_window=int(((calg.parent == p) & (calg.f == f) & calg.in_window).sum()),
                                  min_growth=calg[(calg.parent == p) & (calg.f == f)].growth.min(),
                                  max_growth=calg[(calg.parent == p) & (calg.f == f)].growth.max(),
                                  verdict=verdict)
                             for p in PARENTS for f in CAL_F])

        files = {"sweep_labels.tsv": lab, "precheck_wt.tsv": wt, "floor_nonbinding.tsv": flt,
                 "req_nonzero.tsv": rq, "calibration_grid.tsv": calg, "calibration.tsv": calv}
        for n, df in files.items():
            man.add(atomic_write_tsv(df, out / n))
        mopen = lab[lab.variant == "M-open"].set_index("parent").wt_growth.to_dict()
        summ = dict(
            n_labels=len(lab), n_distinct_sweeps=len(wt), m0_growth=m0, m_open_growth=mopen,
            precheck_fail=lab[~lab.precheck_pass][["parent", "variant", "medium_key", "f", "gate_open",
                                                   "wt_growth"]].to_dict("records"),
            floor_mult=mult.to_dict("records"),
            floor_nonbinding_at_final_mult=bool(flt[flt.floor_mult == flt.parent.map(
                {p: pcs[p].floor_mult for p in PARENTS})].nonbinding.all()),
            floor_binding_cases=flt[~flt.nonbinding].to_dict("records"),
            req_nonzero_all=bool(rq.req_nonzero_pass.all()), calibration=calv.to_dict("records"))
        js = out / "precheck_summary.json"
        js.write_text(json.dumps(summ, indent=1, default=str))
        man.add(js)
        man.extra = dict(n_sweeps=len(wt), n_lp=len(wt) + 2 * len(flt) + len(cal_tasks))
        print(json.dumps({k: v for k, v in summ.items() if k != "calibration"}, indent=1, default=str)[:4000])
        print(calv.to_string())


if __name__ == "__main__":
    main()
