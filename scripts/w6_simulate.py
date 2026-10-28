"""W6 flux-state simulation: one model state per (arm, assay medium bounds), GLPK, <= 6 workers.

    .venv/bin/python scripts/w6_simulate.py --config provisional [--workers 6] [--pilot N] [--allow-dirty]

For every arm of the config and every W2 medium key used by an assay map (PRISM SOP medium, GDSC Screen Medium,
the line's culture medium) the bounds are built with the frozen W2 library (``media.uptake_bounds``), deduplicated
by bounds signature, and each distinct state gets:
  - WT pFBA flux (X_flux = asinh(v / mu));
  - the single-gene KO sweep (X_KO = e = 1 - ratio_r), same engine and status policy as W3 (``w3.ko_sweep``);
  - the KO pFBA flux of the curated drug-target genes (X_KOflux), with the exact zero-flux shortcut.
Every LP has the 300 s limit with one retry; a non-optimal status after the retry aborts the job (SolverInfraError).
A generic model has one state per medium: its features are line-invariant except through the medium.

Outputs: results/w6/<tag>/features/{states.tsv, flux_long.parquet, ko_ratio.parquet, ko_growth.tsv, checks.json};
results/w6/tmp/<tag>/koflux_long.parquet (bulky, regenerable).
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one BLAS thread per process (<= 6 worker processes in total; DESIGN §11.4)

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w6_common as C  # noqa: E402

from nextgem import paths, w6  # noqa: E402


def medium_keys_needed(lines_df: pd.DataFrame) -> pd.Series:
    """medium key -> number of cohort lines that use it in any assay map (descending)."""
    lines = list(lines_df.index)
    maps = {k: w6.assay_medium_map(k, lines, lines_df.medium_key, lines_df.gdsc_screen_medium) for k in w6.ASSAY_MAPS}
    return pd.concat(maps.values()).value_counts()


def build_states(cfg, keys, doc, pcs):
    from nextgem import media

    states = []  # (arm, signature, bounds, [medium keys])
    for arm in cfg.arms:
        by_sig: dict[str, dict] = {}
        for k in keys:
            b = media.uptake_bounds(None, arm.bounds_variant, doc[k], pcs[arm.parent], f=arm.f, gate_open=arm.gate_open)
            sig = media.bounds_signature(b)
            by_sig.setdefault(sig, {"bounds": b, "keys": []})["keys"].append(k)
        for sig, d in by_sig.items():
            states.append((arm, sig, d["bounds"], d["keys"]))
    return states


def main(a) -> int:
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from multiprocessing import get_context

    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import w2_engine

    cfg, cfg_path = C.load_config(a.config)
    run_tag = cfg.tag + ("_pilot" if a.pilot else "")
    out = C.run_dir(run_tag) / "features"
    tmp = C.tmp_dir(run_tag)
    tmp.mkdir(parents=True, exist_ok=True)
    lines_df = pd.read_parquet(C.tmp_dir("inputs") / "lines.parquet")
    table = pd.read_csv(C.OUT_DIR / "inputs" / "compounds.tsv", sep="\t", comment="#")
    targets = w6.target_genes(table)
    need = medium_keys_needed(lines_df)
    keys = list(need.index)
    if a.pilot:
        keys = keys[: a.pilot]
        targets = targets[: a.pilot_targets]
    _, _, doc, _, _, pcs, _, _ = w2_engine.sweep_context()
    states = build_states(cfg, keys, doc, pcs)
    parents = {arm.parent for arm in cfg.arms}
    if len(parents) != 1:
        raise ValueError("one parent per run")
    parent = parents.pop()
    genes_tsv = str(paths.data_root() / ("Human1/model/genes.tsv" if parent == "human1" else "Human-GEM-main/model/genes.tsv"))
    print(f"{len(states)} states ({len(keys)} medium keys x {len(cfg.arms)} arms), {len(targets)} target genes", flush=True)
    wd = C.CPUWatchdog(a.budget_cpu_h)
    results = {}
    jobs = [(f"{arm.name}|{sig}", b, arm.variant, targets, True, True) for arm, sig, b, _ in states]
    n_workers = min(a.workers, C.MAX_WORKERS, len(jobs))
    t0 = time.time()
    with ProcessPoolExecutor(n_workers, mp_context=get_context("spawn"), initializer=w6.init_worker,
                             initargs=(parent, genes_tsv)) as ex:
        meta = ex.submit(w6.worker_meta).result()
        futs = {ex.submit(w6.state_task, j): j[0] for j in jobs}
        for f in as_completed(futs):
            r = f.result()
            results[r["key"]] = r
            wd.add(r["cpu_s"])
            print(f"  {r['key']}: mu={r['wt_growth']:.4g} koflux LPs {r['n_lp_koflux']} (shortcut {r['n_shortcut_koflux']}) "
                  f"{r['wall_s']:.0f}s [{len(results)}/{len(jobs)}, {time.time() - t0:.0f}s]", flush=True)
    rxn_ids, genes = meta["rxn_ids"], meta["genes"]

    srows, flong, klong, kg, ratio = [], [], [], [], {}
    for arm, sig, b, mkeys in states:
        key = f"{arm.name}|{sig}"
        r = results[key]
        srows.append({"state": key, "arm": arm.name, "signature": sig, "medium_keys": ";".join(sorted(mkeys)),
                      "variant": (arm.variant or {}).get("name", "none"), "bounds_variant": arm.bounds_variant, "f": arm.f,
                      "gate_open": arm.gate_open, "n_open_exchanges": len(b), "wt_growth": r["wt_growth"],
                      "ko_wt_growth": r.get("ko_wt_growth"), "ko_n_solved": r.get("ko_n_solved"),
                      "shortcut_verify_max_dev": r.get("shortcut_verify_max_dev"), "n_lp_koflux": r["n_lp_koflux"],
                      "n_shortcut_koflux": r["n_shortcut_koflux"], "wall_s": round(r["wall_s"], 1), "cpu_s": round(r["cpu_s"], 1)})
        idx, val = r["flux"]
        flong.append(pd.DataFrame({"state": key, "rxn": np.asarray(rxn_ids)[idx], "v": val}))
        for s, (ki, kv) in r["koflux"].items():
            klong.append(pd.DataFrame({"state": key, "target": s, "rxn": np.asarray(rxn_ids)[ki], "v": kv}))
            kg.append({"state": key, "target": s, "ko_growth": r["ko_growth"][s], "wt_growth": r["wt_growth"]})
        ratio[key] = r["ratio_r"]
    states_df = pd.DataFrame(srows)
    hdr = [f"{cfg.label}; config {cfg_path.name}; parent {parent}; GLPK"]
    C.atomic_write_tsv(states_df, out / "states.tsv", header_lines=hdr)
    C.atomic_write_parquet(pd.concat(flong, ignore_index=True), out / "flux_long.parquet")
    C.atomic_write_parquet(pd.DataFrame(ratio, index=pd.Index(genes, name="gene")), out / "ko_ratio.parquet")  # genes x states
    C.atomic_write_tsv(pd.DataFrame(kg), out / "ko_growth.tsv", header_lines=hdr)
    if klong:
        C.atomic_write_parquet(pd.concat(klong, ignore_index=True), tmp / "koflux_long.parquet", register=False)
    medium_state = pd.DataFrame([{"arm": arm.name, "medium_key": k, "state": f"{arm.name}|{sig}"}
                                 for arm, sig, _, mk in states for k in mk])
    C.atomic_write_tsv(medium_state, out / "medium_state.tsv", header_lines=hdr)

    # known-answer check: the generic arm on M2 f=1 must equal the W2 KO sweeps (results/w2/ko) bitwise
    checks = {"label": cfg.label, "n_states": len(states), "n_medium_keys": len(keys), "pilot": bool(a.pilot),
              "targets": targets, "wall_s": round(time.time() - t0, 1)}
    w2 = pd.read_parquet(C.PROJECT / "results/w2/ko" / f"ratio_r_{parent}.parquet")
    lab = pd.read_csv(C.PROJECT / "results/w2/ko/sweep_labels.tsv", sep="\t")
    eq = {}
    for arm, sig, _, mk in states:
        if arm.variant is None and arm.bounds_variant == "M2" and arm.gate_open:
            sid = f"{parent}:{sig}"
            if sid in w2.index:
                a_ = w2.loc[sid, genes].to_numpy(float)  # w2: sweeps x genes
                b_ = np.asarray(ratio[f"{arm.name}|{sig}"], float)
                eq[sid] = bool(np.array_equal(a_, b_, equal_nan=True))
    checks["generic_ko_equals_w2"] = eq
    checks["generic_ko_equals_w2_all"] = bool(eq) and all(eq.values())
    checks["w2_labels_present"] = int(lab.parent.eq(parent).sum())
    C.atomic_write_json(checks, out / "checks.json")
    print(json.dumps({k: v for k, v in checks.items() if k != "targets"}, default=str), flush=True)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="provisional")
    ap.add_argument("--workers", type=int, default=C.MAX_WORKERS)
    ap.add_argument("--pilot", type=int, default=0, help="pilot: only the N most-used medium keys")
    ap.add_argument("--pilot-targets", type=int, default=8)
    ap.add_argument("--budget-cpu-h", type=float, default=4.0)
    ap.add_argument("--allow-dirty", action="store_true", help="label-free pilot only; recorded in the manifest")
    a = ap.parse_args()
    cfg, cfg_path = C.load_config(a.config)
    ins = [cfg_path, C.OUT_DIR / "inputs" / "compounds.tsv", C.PROJECT / "results/w2/media_library/compositions.tsv",
           C.PROJECT / "results/w2/media_library/exchange_map.tsv", C.PROJECT / "results/w2/floor/floor.tsv",
           C.PROJECT / "results/w2/floor/requirements.tsv", C.PROJECT / "results/w2/floor/floor_mult.tsv"]
    ins += [C.PROJECT / x for x in sorted({"results/w3/structure/variants_human1.json"})]
    sys.exit(C.run_main(lambda: main(a), inputs=ins, seeds=[11], solver="glpk", workers=min(a.workers, C.MAX_WORKERS),
                        allow_dirty=a.allow_dirty, budget_cpu_h=a.budget_cpu_h,
                        note={"step": "simulate", "config": cfg.tag, "label": cfg.label, "pilot": a.pilot}))
