"""W2 step 4: generic single-gene KO sweeps per distinct medium variant, Human1 and Human2 (label-free outputs).

One sweep per distinct (parent x bounds signature) that passed the WT precheck (results/w2/prechecks/precheck_wt.tsv).
Engine: nextgem.w2_engine.ko_task (exact shortcuts; every sweep also re-solves a random sample of shortcut genes and
records the max deviation). Equality check: the M0 sweep must reproduce results/baseline/{parent}_ko.parquet ratio_r
for every gene (the frozen baseline engine nextgem.knockouts).

Outputs (results/w2/ko/):
  ratio_r_{parent}.parquet   sweeps x genes (Ensembl ids), float64 ratio_r
  sweeps.tsv                 per sweep: parent, wt growth, n LPs solved, infeasible KOs, lethal (<0.01) and partial
                             counts, runtime, shortcut verification
  sweep_labels.tsv           every ladder label -> sweep (with precheck and lethal counts)
  m0_equality.tsv            M0 sweep vs baseline parquet

Usage: .venv/bin/python scripts/w2_ko_sweeps.py [--workers 4] [--verify 20]
"""

from __future__ import annotations

import argparse
import pickle
import time
from concurrent.futures import as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from nextgem import w2_engine as eng
from nextgem.w2_runtime import PARENTS, Manifest, atomic_write_tsv, w2_dir


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--verify", type=int, default=20, help="shortcut genes re-solved per sweep")
    ap.add_argument("--m0-only", action="store_true", help="engine check: run only the M0 sweeps (to a tmp dir)")
    args = ap.parse_args()
    out = w2_dir() / ("ko" if not args.m0_only else "tmp/ko_m0_check")
    out.mkdir(parents=True, exist_ok=True)
    with Manifest("scripts/w2_ko_sweeps.py", solver="gurobi", workers=args.workers) as man:
        lines, sheets, doc, phys, sens, pcs, exchanges, hams = eng.sweep_context()
        sweeps = eng.build_sweeps(exchanges, hams, pcs, doc, phys, sens)
        pre = pd.read_csv(w2_dir() / "prechecks" / "precheck_wt.tsv", sep="\t").set_index("sweep_id")
        assert set(pre.index) == {s.sweep_id for s in sweeps}, "sweep set differs from the precheck run"
        run = [s for s in sweeps if pre.loc[s.sweep_id, "precheck_pass"]]
        # largest-first ordering is not known in advance; interleave parents
        run.sort(key=lambda s: (s.signature, s.parent))
        if args.m0_only:
            run = [s for s in run if any(lab["variant"] == "M0" for lab in s.labels)]
        print(f"{len(run)} sweeps pass the precheck ({len(sweeps) - len(run)} blocked)", flush=True)
        t0 = time.perf_counter()
        # per-sweep checkpoints (resume: sweeps already on disk are not re-run)
        parts = w2_dir() / "tmp" / ("ko_parts" if not args.m0_only else "ko_parts_m0")
        parts.mkdir(parents=True, exist_ok=True)

        def part(sid: str) -> Path:
            return parts / (sid.replace(":", "_") + ".pkl")

        todo = [s for s in run if not part(s.sweep_id).exists()]
        print(f"{len(run) - len(todo)} sweeps already checkpointed, {len(todo)} to run", flush=True)
        with eng.make_pool(args.workers) as pool:
            genes = {p: g for p, g in zip(PARENTS, pool.map(eng.genes_task, PARENTS))}
            futs = [pool.submit(eng.ko_task, (s.parent, s.sweep_id, s.bounds, args.verify)) for s in todo]
            for i, fut in enumerate(as_completed(futs)):
                r = fut.result()
                tmp = part(r["sweep_id"]).with_suffix(".tmp")
                tmp.write_bytes(pickle.dumps(r))
                tmp.replace(part(r["sweep_id"]))
                if (i + 1) % 10 == 0 or i + 1 == len(futs):
                    print(f"  {i + 1}/{len(futs)} sweeps, {time.perf_counter() - t0:.0f} s", flush=True)
        results = [pickle.loads(part(s.sweep_id).read_bytes()) for s in run]
        meta = []
        for p in PARENTS:
            rs = [r for r in results if r["parent"] == p]
            mat = pd.DataFrame(np.vstack([r["ratio_r"] for r in rs]).astype(np.float64),
                               index=pd.Index([r["sweep_id"] for r in rs], name="sweep_id"), columns=genes[p])
            f = out / f"ratio_r_{p}.parquet"
            mat.to_parquet(f)
            man.add(f)
            for r in rs:
                v = r["ratio_r"]
                meta.append(dict(sweep_id=r["sweep_id"], parent=p, wt_growth=r["wt_growth"], n_genes=len(v),
                                 n_lp_solved=r["n_solved"], n_genes_needing_lp=r.get("n_genes_lp"),
                                 n_infeasible=r["n_infeasible"],
                                 n_lethal=int((v < 0.01).sum()), n_partial=int(((v >= 0.01) & (v < 1)).sum()),
                                 n_ratio1=int((v == 1).sum()), n_nan=int(np.isnan(v).sum()),
                                 runtime_s=round(r["runtime_s"], 2), n_shortcut_verified=r["n_verified"],
                                 shortcut_verify_max_dev=r["shortcut_verify_max_dev"]))
        meta = pd.DataFrame(meta)
        lab = pd.read_csv(w2_dir() / "prechecks" / "sweep_labels.tsv", sep="\t")
        lab = lab.merge(meta[["sweep_id", "n_lethal", "n_partial"]], on="sweep_id", how="left")
        if args.m0_only:
            lab = lab[lab.sweep_id.isin(meta.sweep_id)]

        # equality with the frozen baseline engine on M0
        eq = []
        for p in PARENTS:
            sid = lab[(lab.parent == p) & (lab.variant == "M0")].sweep_id.iloc[0]
            base = pd.read_parquet(w2_dir().parent / "baseline" / f"{p}_ko.parquet").set_index("gene").ratio_r
            mine = pd.read_parquet(out / f"ratio_r_{p}.parquet").loc[sid].astype(float)
            d = (mine.reindex(base.index) - base).abs()
            eq.append(dict(parent=p, n_genes=len(base), n_equal=int((d == 0).sum()), max_abs_diff=float(d.max()),
                           n_diff_gt_1e6=int((d > 1e-6).sum()), genes_missing=int(mine.reindex(base.index).isna().sum())))
        eq = pd.DataFrame(eq)
        for n, df in {"sweeps.tsv": meta, "sweep_labels.tsv": lab, "m0_equality.tsv": eq}.items():
            man.add(atomic_write_tsv(df, out / n))
        man.extra = dict(n_sweeps=len(run), n_lp=int(meta.n_lp_solved.sum() + meta.n_shortcut_verified.sum()),
                         sum_task_s=round(meta.runtime_s.sum(), 1))
        print(eq.to_string())
        print(meta.describe().T[["mean", "min", "max"]].to_string())


if __name__ == "__main__":
    main()
