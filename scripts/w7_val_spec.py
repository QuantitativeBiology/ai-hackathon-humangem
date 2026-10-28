"""W7 validation-look preparation (NOT the look): label-free validation-line scores of promoted ideas + the shortlist
look spec for the reviewed generic harness (scripts/evaluate_val.py; DESIGN §W7 protocol 3, §17 A6(2), A11f-g).

Scores of an ML idea are model outputs: the idea is fitted on ``labels.train()`` (all training lines; DESIGN §4.2,
fitted components from the fitting tier only) with SEEDS[0:3], each fit asserted training-only
(``w7.assert_training_only``), and applied to the label-free features of the validation expression lines (same 198
lines as the M6 look file). The look score is the mean of the three seeds' predictions. No validation label is
read here. Methods: H1, H2 (generic, equal to the W1 scores), P1 (M6 cell (3), from the SHA-pinned M6 file),
W7:B6 (W7-01 = B6 as at M2, the current-best reference on training), and the candidate(s).

Usage: .venv/bin/python scripts/w7_val_spec.py --look W7_shortlist_1 --candidates W7-06 [--workers 4]
Outputs: results/w7/val/w7_val_scores_<look>.parquet, results/w7/val/look_spec_<look>.json
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import concurrent.futures as cf  # noqa: E402
import json  # noqa: E402
import multiprocessing  # noqa: E402
import resource  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w7_common as wc  # noqa: E402

_INP = {}
P1_NAME = "M6:H1+Bloo/M2(f=1)"


def _fit_predict(job):
    from nextgem import labels, w7

    inp = _INP["inp"]
    fit_long = labels.train(columns=w7.LABEL_COLS)
    w7.assert_training_only(fit_long.ModelID, inp.tier_of)
    idea = w7.IDEAS[job["idea"]]
    f = w7.fit(inp, fit_long, idea, job["seed"])
    lines = _INP["val_lines"]
    assert (inp.tier_of.reindex(lines) == "validation").all()
    p = w7.predict(inp, f, lines)
    print(f"  {job['idea']} seed {job['seed']}: fit on {f.n_rows} training pairs", flush=True)
    return job["idea"], job["seed"], p


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--look", required=True)
    ap.add_argument("--candidates", required=True, help="comma-separated idea ids, simplest first")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seeds", default="11,23,37")
    ap.add_argument("--notes", default="")
    a = ap.parse_args(argv)
    from nextgem import w1, w7
    from nextgem.run import provenance

    cands = a.candidates.split(",")
    for c in cands:
        wc.assert_logged(c)
    seeds = [int(s) for s in a.seeds.split(",")]
    with provenance(wc.MANIFEST, inputs=wc.input_files(), seeds=seeds, workers=a.workers, allow_dirty=True,
                    author=wc.AUTHOR) as prov:
        prov.note(label_files_frozen_sha256=wc.frozen_label_hashes(), look=a.look, candidates=cands,
                  step="validation look file (training fit, label-free validation features)")
        t0 = time.time()
        inp = wc.build_inputs(tiers=("train", "validation"), art=any("e_ART" in w7.IDEAS[c].extra for c in cands))
        m6v = pd.read_parquet(wc.PROJECT / "results" / "m6" / "val" / "m6_val_scores.parquet")
        val_lines = sorted(m6v.ModelID.unique())
        _INP.update(inp=inp, val_lines=val_lines)
        ideas = ["W7-01"] + [c for c in cands if c != "W7-01"]
        jobs = [{"idea": i, "seed": s} for i in ideas for s in seeds]
        preds = {}
        with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=multiprocessing.get_context("fork")) as pool:
            for i, s, p in pool.map(_fit_predict, jobs, chunksize=1):
                preds.setdefault(i, []).append(p)
        genes = inp.genes
        df = m6v[["ModelID", "entrez", "H1", "H2", P1_NAME]].copy()
        df = df.sort_values(["ModelID", "entrez"]).reset_index(drop=True)
        for i, ps in preds.items():
            M = sum(ps) / len(ps)
            long = M.stack().rename("s").reset_index()
            long.columns = ["ModelID", "entrez", "s"]
            long["entrez"] = long.entrez.astype(np.int64)
            col = "W7:B6" if i == "W7-01" else f"W7:{i}"
            df = df.merge(long.rename(columns={"s": col}), on=["ModelID", "entrez"], how="left", validate="1:1")
            if df[col].isna().any():
                raise RuntimeError(f"{col}: missing predictions")
            df[col] = df[col].astype(float)
        gen = w1.generic_table(pd.Series(genes))
        for s in ("H1", "H2"):
            g = gen[f"e_{s}"].reindex(df.entrez).to_numpy()
            if not np.array_equal(df[s].to_numpy(), g):
                raise RuntimeError(f"{s} differs from the W1 generic scores")
        sf = wc.atomic_write_parquet(df, wc.VAL / f"w7_val_scores_{a.look}.parquet")
        methods = ["H1", "H2", P1_NAME, "W7:B6"] + [f"W7:{c}" for c in cands if c != "W7-01"]
        cand_cols = [f"W7:{c}" if c != "W7-01" else "W7:B6" for c in cands]
        prov_files = ["results/m6/scores/SHA256SUMS.tsv",
                      "results/m6/val/m6_val_scores.parquet", "src/nextgem/w7.py", "scripts/w7_common.py",
                      "scripts/w7_screen.py", "scripts/w7_summarise.py", "scripts/w7_val_spec.py", "scripts/w1_b6.py",
                      "src/nextgem/w1.py"]
        prov_files += [str(p.relative_to(wc.PROJECT)) for p in sorted(wc.OUT.glob("screen_*.tsv"))]
        prov_files += [str(p.relative_to(wc.PROJECT)) for p in sorted(wc.OUT.glob("nested_*_beta.json"))]
        if (wc.OUT / "sweeps" / "ratio_r_human1.parquet").exists():
            prov_files += ["results/w7/sweeps/ratio_r_human1.parquet", "results/w7/sweeps/sweeps.json"]
        spec = {
            "look_id": a.look, "kind": "shortlist", "workstream": "W7", "milestone": "M7",
            "scores": {"path": str(sf.relative_to(wc.PROJECT)), "sha256": wc_sha(sf), "layout": "line_gene"},
            "methods": methods, "candidates": cand_cols, "reference": "H1", "labels": ["L_prob"],
            "genesets": ["all", "noCE"], "universes": ["U"], "strata": ["none", "expressed_pairs"],
            "B": 2000, "seed": 11,
            "decision_rule": ("§8 gate per candidate: eligible iff its paired validation Δ̄_lin vs G* (H1/M0, the W1 "
                              "look's G-star) has CI lower bound > 0 (L_prob, all genes, U, no stratum). Among eligible "
                              "candidates the paired-difference one-SE rule (§9) picks the first candidate in spec order "
                              "(simplest first) whose Δ̄_lin(best − candidate) <= its bootstrap SE. W7 protocol: one "
                              "validation look per promoted idea; scalar only. P2 (DESIGN §W7 step 7) = best validated "
                              "idea, the simpler of two within one paired SE. P1 and W7:B6 are descriptive methods "
                              "(C5 contrast P2 − B6 is a T2 test-tier estimate, not decided here)."),
            "provenance": {p: wc_sha(wc.PROJECT / p) for p in prov_files},
            "notes": a.notes + f" autoresearch/ideas.tsv SHA-256 at spec time: {wc_sha(wc.IDEAS_TSV)}.",
        }
        f2 = wc.atomic_write_text(json.dumps(spec, indent=1, ensure_ascii=False) + "\n", wc.VAL / f"look_spec_{a.look}.json")
        prov.output(sf, f2)
        ru = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_h = sum(r.ru_utime + r.ru_stime for r in ru) / 3600
        wc.log_budget(f"w7_val_spec {a.look}", cpu_h, (time.time() - t0) / 3600, f"{len(jobs)} full-training fits")
        print(json.dumps({"look": a.look, "scores": spec["scores"], "methods": methods, "cpu_h": round(cpu_h, 3)}))
    return 0


def wc_sha(p: Path) -> str:
    from nextgem.w7 import sha256

    return sha256(p)


if __name__ == "__main__":
    raise SystemExit(main())
