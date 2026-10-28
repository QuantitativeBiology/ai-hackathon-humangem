"""W7 screening run: nested 12-fold LOLO on the training tier for logged ideas (DESIGN §W7, §9 TG; A12a).

Labels: only through ``nextgem.cv.lolo`` (``labels.fold(k)`` to fit; held-out labels reach only the harness's score
callback, which returns per-line tAP tables). No validation or test label is read. Every fit asserts that its
lines are training-tier lines (``w7.assert_training_only``).

Each job is one (fold, idea, seed): fit on training minus the fold, predict the fold's expression lines, score.
Fixed comparators per fold: P1 (= H1 + W3 B-loo / M2(f=1), M6 cell (3)), G* (H1), B2 (fold-fitted gene prior), and
the W7-03 artefact-groups variant when its sweeps exist. Ideas of kind 'stack' choose beta by inner LOLO inside the
outer fitting set (scripts/w7_nested_pairs.py, SEEDS[0]); here B6 is refit on training minus the fold.

Outputs: results/w7/tmp/perline_<tag>.parquet (fold, method, seed, ModelID, per-line tAP table); one MANIFEST row;
a results/budget.tsv row.

Usage: .venv/bin/python scripts/w7_screen.py --ideas W7-01,W7-02 --tag r1 [--seeds 11,23,37] [--workers 8]
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import concurrent.futures as cf  # noqa: E402
import datetime as dt  # noqa: E402
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


def _job(job: dict) -> list[pd.DataFrame]:
    from nextgem import cv, w1, w7

    inp = _INP["inp"]
    genes = inp.genes
    t0 = time.process_time()
    kind, fold, seed = job["kind"], job["fold"], job["seed"]

    if kind == "fixed":
        def fp(fit_long, ids):
            w7.assert_training_only(fit_long.ModelID, inp.tier_of)
            ev = w7.eval_lines(ids, inp)
            fitted = w1.fit_components(fit_long, inp.line_class, labels_=(w7.LABEL,))
            out = {"P1": inp.gem_line["e_P1"].reindex(index=ev, columns=genes).astype(float),
                   "Gstar": inp.gem_line["e_Gstar"].reindex(index=ev, columns=genes).astype(float),
                   "B2": pd.DataFrame(np.tile(fitted.prior[w7.LABEL].reindex(genes).fillna(0.0).to_numpy(), (len(ev), 1)),
                                      index=ev, columns=genes)}
            if "e_ART" in inp.gem_line:
                out["ART"] = inp.gem_line["e_ART"].reindex(index=ev, columns=genes).astype(float)
            return out
    else:
        idea = w7.IDEAS[job["idea"]]

        def fp(fit_long, ids):
            ev = w7.eval_lines(ids, inp)
            if idea.kind == "stack":
                beta = _INP["beta"][idea.id][fold]  # nested inner LOLO choice (scripts/w7_nested_pairs.py)
                job["beta"] = beta
                f = w7.fit(inp, fit_long, w7.IDEAS["W7-01"], seed)
                p = w7.predict(inp, f, ev)
                return {idea.id: w7.stack_scores(p, inp.gem_line["e_P1"], beta)}
            f = w7.fit(inp, fit_long, idea, seed)
            return {idea.id: w7.predict(inp, f, ev)}

    res = cv.lolo(fp, lambda p, h: w7.score_lines(p, h, genes), columns=w7.LABEL_COLS, folds=[fold])
    out = []
    for method, lt in res[fold].items():
        t = lt.reset_index().rename(columns={"index": "ModelID"})
        t.insert(0, "seed", seed)
        t.insert(0, "method", method)
        t.insert(0, "fold", fold)
        if "beta" in job:
            t["beta"] = job["beta"]
        out.append(t)
    print(f"  {kind} {job.get('idea', '')} {fold} seed {seed}: {time.process_time() - t0:.0f} s CPU"
          + (f", beta {job['beta']}" if "beta" in job else ""), flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ideas", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--seeds", default="11,23,37")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--no-fixed", action="store_true")
    ap.add_argument("--with-art", action="store_true", help="include the W7-03 artefact variant (only after its "
                    "label-free W3-MEV check, results/w7/sweeps/sweeps.json)")
    ap.add_argument("--folds", default="")
    a = ap.parse_args(argv)
    if a.workers > 8:
        raise SystemExit("at most 8 workers")
    from nextgem import labels, w7
    from nextgem.run import provenance

    wc.register_ideas()
    ideas = [i for i in a.ideas.split(",") if i]
    for i in ideas:
        if i not in w7.IDEAS:
            raise SystemExit(f"{i} is not a registered idea")
        wc.assert_logged(i)
    seeds = [int(s) for s in a.seeds.split(",")]
    with provenance(wc.MANIFEST, inputs=wc.input_files(), seeds=seeds, workers=a.workers, allow_dirty=True,
                    author=wc.AUTHOR) as prov:
        prov.note(label_files_frozen_sha256=wc.frozen_label_hashes(), ideas=ideas, tag=a.tag)
        t0 = time.time()
        need_art = a.with_art
        if need_art and "mev_serum_sweeps_core_all_lethal" not in json.loads((wc.OUT / "sweeps" / "sweeps.json").read_text()):
            raise SystemExit("W7-03 label-free MEV check missing")
        if any("e_ART" in w7.IDEAS[i].extra for i in ideas) and not need_art:
            raise SystemExit("idea needs --with-art")
        _INP["inp"] = wc.build_inputs(tiers=("train",), art=need_art)
        _INP["beta"] = {}
        for i in ideas:
            if w7.IDEAS[i].kind == "stack":
                bf = wc.OUT / f"nested_{i}_beta.json"
                _INP["beta"][i] = json.loads(bf.read_text())["beta"]
                prov.input(bf)
        folds = [f for f in (a.folds.split(",") if a.folds else labels.folds())]
        jobs = []
        if not a.no_fixed:
            jobs += [{"kind": "fixed", "fold": f, "seed": 0} for f in folds]
        for i in ideas:
            for s in seeds:
                jobs += [{"kind": "idea", "idea": i, "fold": f, "seed": s} for f in folds]
        print(f"{len(jobs)} jobs on {a.workers} workers", flush=True)
        ctx = multiprocessing.get_context("fork")
        out = []
        with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx) as pool:
            for r in pool.map(_job, jobs, chunksize=1):
                out += r
        df = pd.concat(out, ignore_index=True)
        f = wc.atomic_write_parquet(df, wc.TMP / f"perline_{a.tag}.parquet")
        prov.output(f)
        ru = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_h = sum(r.ru_utime + r.ru_stime for r in ru) / 3600
        wall_h = (time.time() - t0) / 3600
        wc.log_budget(f"w7_screen {a.tag} ({','.join(ideas)})", cpu_h, wall_h,
                      f"{len(jobs)} jobs, seeds {seeds}, {a.workers} workers")
        prov.note(cpu_h=round(cpu_h, 3))
        print(json.dumps({"tag": a.tag, "rows": len(df), "cpu_h": round(cpu_h, 3), "wall_h": round(wall_h, 3)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
