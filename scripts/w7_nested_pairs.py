"""W7 nested inner LOLO for stacking ideas (DESIGN §11.1: nested inner LOLO for every choice; W7-04).

For every unordered pair {k, j} of training folds, one B6 model (idea W7-01, SEEDS[0]) is fitted on training minus
{k, j}. It serves two inner folds: inner j of outer k (``cv.lolo(..., outer=k, folds=[j])``) and inner k of outer j;
both calls hand the fit the same labels (training minus {k, j}, asserted). Held-out labels reach only the harness's
score callback. For each outer fold k, beta is the grid value with the highest lineage-equal mean inner tAP_line over
its 11 inner folds. Nothing about fold k's own labels enters its beta.

Outputs: results/w7/tmp/nested_<idea>.parquet (outer, inner, beta, per-line tAP), results/w7/nested_<idea>_beta.json.
Usage: .venv/bin/python scripts/w7_nested_pairs.py --idea W7-04 --workers 8
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import concurrent.futures as cf  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import multiprocessing  # noqa: E402
import resource  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w7_common as wc  # noqa: E402

_INP = {}


def _pair(job):
    from nextgem import cv, w7
    from nextgem.config import SEEDS

    inp, idea = _INP["inp"], w7.IDEAS[job["idea"]]
    k, j = job["pair"]
    cache = {}
    t0 = time.process_time()

    def fp(fit_long, ids):
        s = frozenset(fit_long.ModelID.unique())
        if "f" not in cache:
            cache["f"], cache["set"] = w7.fit(inp, fit_long, w7.IDEAS["W7-01"], SEEDS[0]), s
        elif cache["set"] != s:
            raise AssertionError("pair fit sets differ")
        p = w7.predict(inp, cache["f"], w7.eval_lines(ids, inp))
        return {str(b): w7.stack_scores(p, inp.gem_line["e_P1"], b) for b in idea.beta_grid}

    out = []
    for outer, inner in ((k, j), (j, k)):
        res = cv.lolo(fp, lambda p, h: w7.score_lines(p, h, inp.genes), columns=w7.LABEL_COLS, outer=outer, folds=[inner])
        for b, lt in res[inner].items():
            t = lt.reset_index().rename(columns={"index": "ModelID"})
            t.insert(0, "beta", float(b))
            t.insert(0, "inner", inner)
            t.insert(0, "outer", outer)
            out.append(t)
    print(f"  pair {k} | {j}: {time.process_time() - t0:.0f} s CPU", flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--idea", required=True)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args(argv)
    if a.workers > 8:
        raise SystemExit("at most 8 workers")
    from nextgem import labels, w7
    from nextgem.config import SEEDS
    from nextgem.run import provenance

    wc.assert_logged(a.idea)
    if w7.IDEAS[a.idea].kind != "stack":
        raise SystemExit("not a stacking idea")
    with provenance(wc.MANIFEST, inputs=wc.input_files(), seeds=[SEEDS[0]], workers=a.workers, allow_dirty=True,
                    author=wc.AUTHOR) as prov:
        prov.note(label_files_frozen_sha256=wc.frozen_label_hashes(), idea=a.idea, step="nested inner LOLO")
        t0 = time.time()
        _INP["inp"] = wc.build_inputs(tiers=("train",), art=False)
        folds = labels.folds()
        jobs = [{"idea": a.idea, "pair": p} for p in itertools.combinations(folds, 2)]
        print(f"{len(jobs)} pair fits on {a.workers} workers", flush=True)
        out = []
        with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=multiprocessing.get_context("fork")) as pool:
            for r in pool.map(_pair, jobs, chunksize=1):
                out += r
        df = pd.concat(out, ignore_index=True)
        f1 = wc.atomic_write_parquet(df, wc.TMP / f"nested_{a.idea}.parquet")
        d = df.dropna(subset=["tap"])
        sc = d.groupby(["outer", "beta", "inner"]).tap.mean().groupby(["outer", "beta"]).mean()
        choice = {o: float(sc.loc[o].idxmax()) for o in folds}
        rec = {"idea": a.idea, "seed": SEEDS[0], "beta_grid": list(w7.IDEAS[a.idea].beta_grid), "beta": choice,
               "n_pair_fits": len(jobs)}
        f2 = wc.atomic_write_text(json.dumps(rec, indent=1) + "\n", wc.OUT / f"nested_{a.idea}_beta.json")
        prov.output(f1, f2)
        ru = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_h = sum(r.ru_utime + r.ru_stime for r in ru) / 3600
        wc.log_budget(f"w7_nested_pairs {a.idea}", cpu_h, (time.time() - t0) / 3600, f"{len(jobs)} pair fits")
        prov.note(cpu_h=round(cpu_h, 3))
        print(json.dumps(rec), f"cpu_h {cpu_h:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
