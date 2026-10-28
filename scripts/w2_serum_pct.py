"""W2 A9(f) sensitivity at the TG-chosen f only: serum components bounded at f_eff = f × serum% / 10 instead of f
(``media.uptake_bounds(..., serum_pct_scaling=True)``), for the shortlisted step M2(f). Label-free KO sweeps.

Scope: the training-tier attribution lines (the set the sensitivity is scored on in ``w2_score_train.py``). Lines
whose parsed serum% is 10 (or that are serum-free) keep their primary sweep; one new sweep per distinct
(parent, medium key, serum%) otherwise. WT precheck as DESIGN W2 (WT > 1e-3 × M0 WT); KO engine
``w2_engine.ko_task`` on GLPK workers.

Outputs (results/w2/sens_serum_pct/): ratio_r_{parent}.parquet, sweeps.tsv, line_sweeps.tsv.

Usage: OMP_NUM_THREADS=1 .venv/bin/python scripts/w2_serum_pct.py --f 1 --workers 6
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
W2 = PROJECT / "results" / "w2"
OUT = W2 / "sens_serum_pct"


def main() -> int:
    from nextgem import labels, media, w2_engine as eng, w2_score as ws
    from nextgem.w2_runtime import PARENTS, Manifest, atomic_write_tsv

    ap = argparse.ArgumentParser()
    ap.add_argument("--f", type=float, required=True, help="the TG-chosen serum level f")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--lp-timeout", type=int, default=120, help="GLPK time limit per LP (s)")
    a = ap.parse_args()
    if a.workers > 8:
        raise SystemExit("at most 8 workers")
    step = f"M2(f={a.f:g})"
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with Manifest("w2_serum_pct.py", solver="glpk", workers=a.workers,
                  inputs=[W2 / "media_library" / "line_media.tsv", W2 / "controls" / "line_sweeps.tsv"]) as man:
        lines, sheets, doc, phys, sens, pcs, exchanges, hams = eng.sweep_context()
        split = labels.load_split().set_index("ModelID")
        lm = lines.set_index("ModelID")
        A = [l for l in lm.index[lm.stratum == "attribution"] if split.tier.get(l) == "train"]
        ls = pd.read_csv(W2 / "controls" / "line_sweeps.tsv", sep="\t")
        ls = ls[ls.ladder_step == step].set_index(["parent", "ModelID"]).sweep_id
        defs, line_rows = {}, []
        for p in PARENTS:
            for l in A:
                key, pct = lm.loc[l, "medium_key"], lm.loc[l, "serum_pct"]
                med = doc[key]
                if not med.serum or not np.isfinite(pct) or pct == 10.0:
                    line_rows.append({"parent": p, "ModelID": l, "sweep_id": ls[(p, l)], "serum_pct": pct, "new": False})
                    continue
                b = media.uptake_bounds(None, "M2", dataclasses.replace(med, serum_pct=float(pct)), pcs[p], f=a.f,
                                        gate_open=True, serum_pct_scaling=True)
                sid = f"{p}:{media.bounds_signature(b)}"
                defs.setdefault(sid, {"parent": p, "medium_key": key, "serum_pct": float(pct), "f_eff": a.f * pct / 10,
                                      "bounds": b})
                line_rows.append({"parent": p, "ModelID": l, "sweep_id": sid, "serum_pct": pct, "new": True})
        print(f"{len(defs)} new sweeps for {len(A)} training attribution lines", flush=True)
        res, meta = {}, []
        with eng.make_pool(a.workers, solver="glpk") as pool:
            genes = {p: g for p, g in zip(PARENTS, pool.map(eng.genes_task, PARENTS))}
            wt0 = {p: pool.submit(eng.wt_task, (p, "M0", hams[p])).result()["wt_growth"] for p in PARENTS}
            wts = dict(zip(defs, pool.map(eng.wt_task, [(d["parent"], sid, d["bounds"]) for sid, d in defs.items()])))
            ok = [sid for sid in defs if wts[sid]["wt_growth"] > media.PRECHECK_REL * wt0[defs[sid]["parent"]]]
            if len(ok) != len(defs):
                raise RuntimeError(f"precheck failed for {sorted(set(defs) - set(ok))}")
            timed_out = []
            for r in pool.map(ws.ko_full_task, [(defs[s]["parent"], s, defs[s]["bounds"], 5, a.lp_timeout) for s in ok]):
                if r["status"] == "ok":
                    res[r["sweep_id"]] = r
                else:
                    timed_out.append(r["sweep_id"])
            ok = [s for s in ok if s in res]
            for row in line_rows:  # a timed-out sweep: the line keeps its primary sweep, flagged
                row["timeout"] = row["sweep_id"] in timed_out
                if row["timeout"]:
                    row["sweep_id"], row["new"] = ls[(row["parent"], row["ModelID"])], False
        for p in PARENTS:
            sids = [s for s in ok if defs[s]["parent"] == p]
            mat = pd.DataFrame(np.vstack([res[s]["ratio_r"] for s in sids]) if sids else np.empty((0, len(genes[p]))),
                               index=pd.Index(sids, name="sweep_id"), columns=genes[p])
            mat.to_parquet(OUT / f"ratio_r_{p}.parquet")
            man.add(OUT / f"ratio_r_{p}.parquet")
            for s in sids:
                d, r = defs[s], res[s]
                meta.append({"sweep_id": s, "parent": p, "medium_key": d["medium_key"], "serum_pct": d["serum_pct"],
                             "f_eff": d["f_eff"], "wt_growth": r["wt_growth"], "wt_rel_M0": r["wt_growth"] / wt0[p],
                             "n_lethal": int((r["ratio_r"] < 0.01).sum()), "n_lp_solved": r["n_solved"],
                             "runtime_s": round(r["runtime_s"], 1), "shortcut_verify_max_dev": r["shortcut_verify_max_dev"]})
        man.add(atomic_write_tsv(pd.DataFrame(meta), OUT / "sweeps.tsv"),
                atomic_write_tsv(pd.DataFrame(line_rows), OUT / "line_sweeps.tsv"))
        summ = {"step": step, "f": a.f, "n_new_sweeps": len(ok), "n_lines": len(A),
                "n_line_parent_changed": int(sum(r["new"] for r in line_rows)),
                "n_timeout_sweeps": len(timed_out), "wall_s": round(time.time() - t0, 1)}
        man.extra.update(**summ)
        print(json.dumps(summ, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
