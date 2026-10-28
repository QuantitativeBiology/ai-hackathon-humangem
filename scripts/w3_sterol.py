"""W3 sterol branch (M3_W2 review §2.7; A9 j): what leaving the vitamin-D catabolite group out of the biomass
cofactor pool (B-loo:vitaminD_catabolites) does to the sterol-branch genes and to the mevalonate core, in every
medium of the frozen W2 library. Label-free.

Genes: sterol branch (final a-priori list §1.2: FDFT1, SQLE, LSS, CYP51A1, TM7SF2, MSMO1, NSDHL, HSD17B7, SC5D, EBP,
DHCR7, DHCR24) and mevalonate core (§1.1H: HMGCS1, HMGCR, MVK, PMVK, MVD, IDI1, FDPS, GGPS1).
Media: every distinct gate-open bounds set of the W2 ladder with role baseline / documented / physiological
(M0, M1 and M2(f) for the 22 documented media, M-HPLM, M-Plasmax), rebuilt with ``nextgem.w2_engine.build_sweeps``
from the frozen library, plus M0 + cholesterol uptake at 1000 (the review's clean case).
Conditions: T0 (biomass as released) and B-loo:vitaminD_catabolites (variant from results/w3/structure); for the
review's diagnosis also B-loo:retinoid_catabolites and both groups together.
Each gene KO is solved explicitly (no shortcut). Serum cholesterol is open in M2(f) of serum media, in M-HPLM and
M-Plasmax (W2 serum model).

Decision rule (DESIGN W2/W3 control, line ~520): the mevalonate-core genes that are lethal on M0 (HMGCR, MVK, PMVK,
MVD) must stay lethal (ratio_r < 0.01) in every medium under every condition; a loss is a modelling error.

Outputs: results/w3/sterol/sterol_ratios.tsv (long), sterol_summary.tsv, sterol_summary.json.
Usage: .venv/bin/python scripts/w3_sterol.py --workers 4
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import multiprocessing
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w3_common import MAX_WORKERS, OUT_DIR, atomic_write_tsv, current, run_main  # noqa: E402

STEROL = ["FDFT1", "SQLE", "LSS", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", "HSD17B7", "SC5D", "EBP", "DHCR7", "DHCR24"]
MEVALONATE = ["HMGCS1", "HMGCR", "MVK", "PMVK", "MVD", "IDI1", "FDPS", "GGPS1"]
CHOLESTEROL_EX = "MAR09285"
URIDINE_EX = "MAR09439"
# --panel pyrimidine (coordinator 2026-10-27, from notes/reviews/M3_W5_review.md §4; label-free): CAD / DHODH / UMPS
# under B-loo of (S)-dihydroorotate in the small-metabolite pool MAR10064, reported next to the vitamin-D group.
PYRIMIDINE = ["CAD", "DHODH", "UMPS"]
_W: dict = {}


def _init(parent, gene_ids, variants):
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import w3
    from nextgem.model import load_model

    m = load_model(parent)
    m.objective = w3.BIOMASS
    _W.update(model=m, genes=gene_ids, variants={k: w3.Variant.from_json(v) for k, v in variants.items()})


def _job(job):
    from nextgem import w3

    sweep, bounds, cond = job
    m = _W["model"]
    t0 = time.process_time()
    out = []
    with m:
        m.medium = bounds
        _W["variants"][cond].apply(m)
        m.objective = w3.BIOMASS
        wt, st = w3.solve_growth(m)
        for sym, gid in _W["genes"].items():
            if gid is None:
                out.append((sym, np.nan, "absent"))
                continue
            if not wt > 1e-6:
                out.append((sym, np.nan, "wt_dead"))
                continue
            with m:
                m.genes.get_by_id(gid).knock_out()
                g, s = w3.solve_growth(m)
            out.append((sym, float(np.clip(np.round(g / wt, 6), 0, 1)) + 0.0, s))
    return sweep, cond, wt, st, out, time.process_time() - t0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--panel", choices=["sterol", "pyrimidine"], default="sterol")
    a = ap.parse_args(argv)
    pyr = a.panel == "pyrimidine"
    panel_genes = PYRIMIDINE if pyr else STEROL + MEVALONATE
    sfx = "_pyrimidine" if pyr else ""
    if a.workers > MAX_WORKERS:
        raise SystemExit("at most 8 worker processes")
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import paths, w2_engine, w3

    t0 = time.time()
    lines, sheets, doc, phys, sens, pcs, exchanges, hams = w2_engine.sweep_context()
    sweeps = w2_engine.build_sweeps(exchanges, hams, pcs, doc, phys, {})
    rows, summ = [], {}
    for p in ("human1", "human2"):
        gm = pd.read_parquet(paths.baseline_dir() / f"{p}_ko.parquet", columns=["gene", "symbol"])
        gid = {s: (gm.gene[gm.symbol == s].iloc[0] if (gm.symbol == s).any() else None) for s in panel_genes}
        vj = json.loads((OUT_DIR / "structure" / f"variants_{p}.json").read_text())
        byname = {v["name"]: v for v in vj["bloo"]}
        vd, rt = byname["B-loo:vitaminD_catabolites"], byname["B-loo:retinoid_catabolites"]
        both = w3.Variant("B-loo:vitaminD+retinoid", "combo", {"variants": [w3.Variant.from_json(vd), w3.Variant.from_json(rt)]})
        variants = {"T0": w3.Variant("T0").to_json(), "B-loo:vitaminD": vd, "B-loo:retinoid": rt, "B-loo:vitaminD+retinoid": both.to_json()}
        if pyr:
            variants = {"T0": w3.Variant("T0").to_json(), "B-loo:dihydroorotate": byname["B-loo:dihydroorotate"]}
        meds = {}
        for s in sweeps:
            if s.parent != p:
                continue
            labs = [l for l in s.labels if l["gate_open"] and l["role"] in ("baseline", "documented", "physiological")]
            if labs:
                meds[s.sweep_id] = (s.bounds, labs)
        m0 = dict(hams[p])
        m0c = {**m0, CHOLESTEROL_EX: 1000.0}
        meds[f"{p}:M0+cholesterol"] = (m0c, [{"variant": "M0+chol", "medium_key": "Hams_GR", "f": None, "role": "review_case"}])
        if pyr:
            meds[f"{p}:M0+uridine"] = ({**m0, URIDINE_EX: 1000.0}, [{"variant": "M0+uridine", "medium_key": "Hams_GR", "f": None, "role": "review_case"}])
        jobs = [(sid, b, c) for sid, (b, _) in meds.items() for c in variants]
        print(f"[{p}] {len(meds)} media x {len(variants)} conditions = {len(jobs)} jobs", flush=True)
        ctx = multiprocessing.get_context("spawn")
        with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init,
                                    initargs=(p, gid, variants)) as pool:
            for sid, cond, wt, st, out, cpu in pool.map(_job, jobs, chunksize=2):
                b, labs = meds[sid]
                chol = b.get(CHOLESTEROL_EX, 0.0)
                lab = "; ".join(f"{l['variant']}:{l['medium_key']}" + (f"(f={l['f']})" if l.get("f") is not None else "") for l in labs)
                for sym, r, s in out:
                    rows.append({"parent": p, "sweep_id": sid, "labels": lab, "cholesterol_uptake_ub": chol,
                                 "serum_cholesterol": chol > 0, "uridine_uptake_ub": b.get(URIDINE_EX, 0.0),
                                 "condition": cond, "wt_growth": wt, "gene": sym,
                                 "panel": "pyrimidine" if pyr else ("sterol" if sym in STEROL else "mevalonate"),
                                 "ratio_r": r, "status": s})
        print(f"[{p}] done {time.time() - t0:.0f} s", flush=True)
    df = pd.DataFrame(rows)
    out = atomic_write_tsv(df, OUT_DIR / "sterol" / f"sterol_ratios{sfx}.tsv")
    # summary: per parent x condition x gene x cholesterol state: n media, n lethal, min/max ratio
    if pyr:
        df["uridine"] = df.uridine_uptake_ub > 0
    g = df.groupby(["parent", "condition", "panel", "gene", "uridine" if pyr else "serum_cholesterol"]).ratio_r
    s = pd.DataFrame({"n_media": g.size(), "n_lethal": g.apply(lambda x: int((x < 0.01).sum())),
                      "ratio_min": g.min(), "ratio_max": g.max()}).reset_index()
    so = atomic_write_tsv(s, OUT_DIR / "sterol" / f"sterol_summary{sfx}.tsv")
    core = ["HMGCR", "MVK", "PMVK", "MVD"]
    viol = df[df.gene.isin(core) & ~(df.ratio_r < 0.01) & df.status.ne("absent")]
    if pyr:
        for p in ("human1", "human2"):
            for cond in ("T0", "B-loo:dihydroorotate"):
                for uri in (False, True):
                    e = df[(df.parent == p) & (df.condition == cond) & (df.uridine == uri)]
                    print(p, cond, "uridine" if uri else "no uridine",
                          e.groupby("gene").ratio_r.agg(["min", "max", "size"]).to_dict("index"), flush=True)
    summ = {"n_rows": len(df), "mevalonate_core_lethal_everywhere": bool(viol.empty),
            "core_violations": viol[["parent", "sweep_id", "labels", "condition", "gene", "ratio_r"]].head(20).to_dict("records"),
            "wall_s": round(time.time() - t0, 1)}
    for p in ("human1", "human2") if not pyr else ():
        for cond in ("T0", "B-loo:vitaminD"):
            d = df[(df.parent == p) & (df.condition == cond) & (df.panel == "sterol")]
            for chol in (False, True):
                e = d[d.serum_cholesterol == chol]
                summ[f"{p}|{cond}|serum_chol={chol}|sterol_genes_lethal_in_all_media"] = sorted(
                    e.groupby("gene").ratio_r.apply(lambda x: bool((x < 0.01).all())).loc[lambda x: x].index)
    (OUT_DIR / "sterol" / f"sterol_summary{sfx}.json").write_text(json.dumps(summ, indent=1, default=str))
    current().input(OUT_DIR / "structure" / "variants_human1.json", OUT_DIR / "structure" / "variants_human2.json")
    current().output(out, so, OUT_DIR / "sterol" / f"sterol_summary{sfx}.json")
    print(json.dumps(summ, indent=1, default=str))
    return 0


if __name__ == "__main__":
    w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 4
    raise SystemExit(run_main(main, solver="glpk", workers=w))
