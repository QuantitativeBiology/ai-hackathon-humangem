"""W2 label-free flip analysis (DESIGN §9 W2 "Metrics", domain O1): which generic (M0) binary calls a single nutrient
addition or removal can flip. Label-free: no label file and no labels.py.

Reference: the generic parent on M0 (Ham's GR inputs at 1000; binary call ratio_r < 0.01, §6.1).
Perturbations (one at a time, on M0):
  removal  — each M0 input exchange closed;
  addition — each exchange that any W2 library medium opens (documented, physiological and sensitivity media, the
             serum model and the micronutrient floor; ``w2_engine.build_sweeps``) and M0 lacks, opened at 1000.
WT viability uses the W2 precheck rule (WT > 1e-3 × M0 WT); non-viable perturbations are listed, not swept.
KO engine: ``w2_engine.ko_task`` (full sweep, exact shortcuts) for removals; for additions the exact pruning
``w2_score.addition_prune_safe`` restricts the KOs to the M0-lethal genes when WT_new/WT_M0 <= r_min/0.01 (else a
full sweep). Second, cache-only part: the M0 calls that the cached documented gate-free sweeps (M1, M2(f)) flip.

Outputs (results/w2/flips/): perturbations.tsv, flip_genes.tsv, cached_doc_flips.tsv, flips_summary.json.

GLPK LPs run under a time limit (--lp-timeout); a time-limited KO job is recorded as 'timeout', never as a ratio.
M4 run: additions only (--kinds addition, Human1). The removal half needs ~30 full sweeps per parent (~1 CPU-h each)
and is to be run from a pinned clean worktree (--kinds removal).

Usage: OMP_NUM_THREADS=1 .venv/bin/python scripts/w2_flips.py --parents human1 --workers 6 [--kinds addition removal]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
W2 = PROJECT / "results" / "w2"
OUT = W2 / "flips"
PRECHECK_REL = 1e-3


def main() -> int:
    from nextgem import paths, w2_engine as eng, w2_score as ws
    from nextgem.w2_runtime import Manifest, atomic_write_tsv, load_parent

    ap = argparse.ArgumentParser()
    ap.add_argument("--parents", nargs="+", default=["human1"])
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--verify", type=int, default=5)
    ap.add_argument("--kinds", nargs="+", default=["addition"], choices=("addition", "removal"),
                    help="perturbation kinds to sweep (removals need full sweeps: ~1 CPU-h per parent)")
    ap.add_argument("--lp-timeout", type=int, default=120, help="GLPK time limit per LP (s); timed-out jobs are recorded")
    a = ap.parse_args()
    if a.workers > 8:
        raise SystemExit("at most 8 workers")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with Manifest("w2_flips.py", solver="glpk", workers=a.workers,
                  inputs=[W2 / "ko" / f"ratio_r_{p}.parquet" for p in a.parents]
                  + [W2 / "controls" / "line_sweeps.tsv"]) as man:
        lines, sheets, doc, phys, sens, pcs, exchanges, hams = eng.sweep_context()
        sweeps = eng.build_sweeps(exchanges, hams, pcs, doc, phys, sens)
        pert_rows, flip_rows, cached_rows, summ = [], [], [], {}
        with eng.make_pool(a.workers, parents=tuple(a.parents), solver="glpk") as pool:
            for p in a.parents:
                genes = pool.submit(eng.genes_task, p).result()
                base = pd.read_parquet(paths.baseline_dir() / f"{p}_ko.parquet").set_index("gene")
                ref = base.ratio_r.reindex(genes).to_numpy(float)
                sym = base.symbol.reindex(genes).to_numpy()
                ref_call = ref < 0.01
                r_min = float(np.nanmin(ref[~ref_call]))
                m = load_parent(p, solver="glpk")
                rname = {r.id: (r.name or "") for r in m.exchanges}
                mname = {r.id: ";".join(x.name for x in r.metabolites) for r in m.exchanges}
                del m
                M0 = dict(hams[p])
                lib = set().union(*(set(s.bounds) for s in sweeps if s.parent == p and any(
                    lab["variant"] in ("M1", "M2", "M1g", "M2g") for lab in s.labels)))
                perts = []
                if "removal" in a.kinds:
                    perts += [("removal", e, {k: v for k, v in M0.items() if k != e}) for e in sorted(M0)]
                if "addition" in a.kinds:
                    perts += [("addition", e, {**M0, e: 1000.0}) for e in sorted(lib - set(M0))]
                wt0 = pool.submit(eng.wt_task, (p, "M0", M0)).result()["wt_growth"]
                wts = list(pool.map(eng.wt_task, [(p, f"{k}:{e}", b) for k, e, b in perts]))
                jobs = []
                for (kind, e, b), w in zip(perts, wts):
                    viable = w["wt_status"] == "optimal" and w["wt_growth"] > PRECHECK_REL * wt0
                    row = {"parent": p, "kind": kind, "exchange": e, "metabolite": mname.get(e, ""),
                           "wt_growth": w["wt_growth"], "wt_rel_M0": w["wt_growth"] / wt0, "viable": viable}
                    if viable:
                        if kind == "addition" and ws.addition_prune_safe(w["wt_growth"], wt0, r_min):
                            row["engine"] = "subset(M0-lethal genes; exact pruning)"
                            jobs.append((row, pool.submit(ws.ko_subset_task, (p, f"{kind}:{e}", b,
                                                                              [g for g, c in zip(genes, ref_call) if c],
                                                                              a.lp_timeout))))
                        else:
                            row["engine"] = "full sweep"
                            jobs.append((row, pool.submit(ws.ko_full_task, (p, f"{kind}:{e}", b, a.verify, a.lp_timeout))))
                    else:
                        row["engine"] = "none (WT non-viable)"
                        pert_rows.append(row)
                print(f"{p}: {len(perts)} perturbations, {len(jobs)} viable, WT {time.time() - t0:.0f} s", flush=True)
                for i, (row, fut) in enumerate(jobs, 1):
                    try:
                        r = fut.result()
                    except RuntimeError as err:
                        r = {"status": "timeout", "error": str(err)[:200]}
                    if r.get("status", "ok") != "ok":
                        row.update(status="timeout (LP time limit; not swept)", error=r.get("error", ""))
                        pert_rows.append(row)
                        continue
                    row["status"] = "ok"
                    new = r["ratio_r"]
                    known = ~np.isnan(new)
                    new_call = np.where(known, new < 0.01, ref_call)
                    to_l = np.flatnonzero(new_call & ~ref_call)
                    to_nl = np.flatnonzero(~new_call & ref_call)
                    row.update(n_flip_to_lethal=len(to_l), n_flip_to_nonlethal=len(to_nl),
                               genes_to_lethal=";".join(sorted(str(sym[j]) for j in to_l)),
                               genes_to_nonlethal=";".join(sorted(str(sym[j]) for j in to_nl)),
                               runtime_s=round(r.get("runtime_s", np.nan), 1) if "runtime_s" in r else np.nan)
                    pert_rows.append(row)
                    for j, d in [(j, "to_lethal") for j in to_l] + [(j, "to_nonlethal") for j in to_nl]:
                        flip_rows.append({"parent": p, "gene": genes[j], "symbol": sym[j],
                                          "entrez": base.entrez.reindex(genes).iloc[j], "direction": d,
                                          "kind": row["kind"], "exchange": row["exchange"], "metabolite": row["metabolite"]})
                    if i % 10 == 0:
                        print(f"  {p}: {i}/{len(jobs)} KO jobs, {time.time() - t0:.0f} s", flush=True)
                # cache-only: M0 calls flipped by the cached documented gate-free sweeps (M1, M2(f))
                rr = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
                ls = pd.read_csv(W2 / "controls" / "line_sweeps.tsv", sep="\t")
                ls = ls[(ls.parent == p) & ls.ladder_step.isin(ws.ELIGIBLE)]
                for step, g in ls.groupby("ladder_step"):
                    sw = g.sweep_id.unique()
                    calls = rr.loc[sw, genes].to_numpy() < 0.01
                    for j in np.flatnonzero((calls != ref_call).any(axis=0)):
                        cached_rows.append({"parent": p, "step": step, "symbol": sym[j], "gene": genes[j],
                                            "m0_call": bool(ref_call[j]), "n_media": len(sw),
                                            "n_media_flipped": int((calls[:, j] != ref_call[j]).sum())})
                summ[p] = {"kinds": a.kinds, "n_timeout": sum(1 for r in pert_rows if r.get("parent") == p
                                                                  and str(r.get("status", "")).startswith("timeout")),
                           "wt_M0": wt0, "r_min_nonlethal_M0": r_min, "n_removals_swept": sum(k == "removal" for k, _, _ in perts),
                           "n_additions": sum(k == "addition" for k, _, _ in perts), "n_viable": len(jobs)}
        P = pd.DataFrame(pert_rows)
        F = pd.DataFrame(flip_rows)
        C = pd.DataFrame(cached_rows)
        if len(F):
            fg = F.groupby(["parent", "symbol", "gene", "entrez", "direction"]).agg(
                n_perturbations=("exchange", "size"),
                perturbations=("exchange", lambda s: ";".join(f"{k}:{e}" for k, e in zip(F.loc[s.index, "kind"], s)))).reset_index()
        else:
            fg = pd.DataFrame(columns=["parent", "symbol", "gene", "entrez", "direction", "n_perturbations", "perturbations"])
        for p in a.parents:
            summ[p].update(n_genes_flippable=int(fg[fg.parent == p].gene.nunique()),
                           n_to_lethal=int(fg[(fg.parent == p) & (fg.direction == "to_lethal")].gene.nunique()),
                           n_to_nonlethal=int(fg[(fg.parent == p) & (fg.direction == "to_nonlethal")].gene.nunique()),
                           n_cached_doc_flipped=int(C[C.parent == p].gene.nunique()) if len(C) else 0)
        summ["wall_s"] = round(time.time() - t0, 1)
        outs = [atomic_write_tsv(P, OUT / "perturbations.tsv", float_format="%.6g"),
                atomic_write_tsv(fg, OUT / "flip_genes.tsv"),
                atomic_write_tsv(C, OUT / "cached_doc_flips.tsv")]
        (OUT / "flips_summary.json").write_text(json.dumps(summ, indent=1, default=float))
        outs.append(OUT / "flips_summary.json")
        man.add(*outs)
        man.extra.update(parents=a.parents, kinds=a.kinds, lp_timeout_s=a.lp_timeout, **{f"{p}_summary": summ[p] for p in a.parents})
        print(json.dumps(summ, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
