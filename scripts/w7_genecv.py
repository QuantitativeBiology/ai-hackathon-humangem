"""W7 gene-grouped CV for promoted ideas (DESIGN §W7 protocol 4; §12 W7 row), training tier only.

Lineage x gene-group cross-fitting: genes of U are split into 5 groups (SeedSequence(SEEDS[0]) permutation). For
held-out lineage k and gene group g, the idea is fitted on training minus k restricted to genes outside g (the
labels of g's genes never enter the fit, so their B2/B2m features are missing, NaN) and predicts lineage k's genes
in g. A held-out line's full gene vector is assembled from its 5 group predictions and scored by tAP_line. This
measures how much of an idea's LOLO score rests on gene identity (gene priors) rather than transferable features.

Usage: .venv/bin/python scripts/w7_genecv.py --ideas W7-01,W7-06 --tag g1 [--workers 8]
Outputs: results/w7/tmp/perline_<tag>.parquet (methods suffixed '|genecv'); budget row.
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

N_GROUPS = 5
_INP = {}


def gene_groups(genes, seed: int, n: int = N_GROUPS) -> list[list[int]]:
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    perm = rng.permutation(np.asarray(genes))
    return [sorted(int(g) for g in perm[i::n]) for i in range(n)]


def _job(job):
    from nextgem import cv, w7

    inp, idea = _INP["inp"], w7.IDEAS[job["idea"]]
    groups = _INP["groups"]
    t0 = time.process_time()

    def fp(fit_long, ids):
        ev = w7.eval_lines(ids, inp)
        parts = []
        for g in groups:
            sub = fit_long[~fit_long.entrez.isin(g)]
            f = w7.fit(inp, sub, idea, job["seed"])
            if set(f.fitted.prior[w7.LABEL].dropna().index) & set(g):
                raise AssertionError("held-out genes entered the fit")
            parts.append(w7.predict(inp, f, ev)[g])
        return {f"{idea.id}|genecv": pd.concat(parts, axis=1)[inp.genes]}

    res = cv.lolo(fp, lambda p, h: w7.score_lines(p, h, inp.genes), columns=w7.LABEL_COLS, folds=[job["fold"]])
    out = []
    for m, lt in res[job["fold"]].items():
        t = lt.reset_index().rename(columns={"index": "ModelID"})
        t.insert(0, "seed", job["seed"])
        t.insert(0, "method", m)
        t.insert(0, "fold", job["fold"])
        out.append(t)
    print(f"  {idea.id} {job['fold']}: {time.process_time() - t0:.0f} s CPU", flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ideas", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args(argv)
    if a.workers > 8:
        raise SystemExit("at most 8 workers")
    from nextgem import labels
    from nextgem.config import SEEDS
    from nextgem.run import provenance

    ideas = a.ideas.split(",")
    for i in ideas:
        wc.assert_logged(i)
    with provenance(wc.MANIFEST, inputs=wc.input_files(), seeds=[SEEDS[0]], workers=a.workers, allow_dirty=True,
                    author=wc.AUTHOR) as prov:
        prov.note(label_files_frozen_sha256=wc.frozen_label_hashes(), ideas=ideas, step="gene-grouped CV")
        t0 = time.time()
        inp = wc.build_inputs(tiers=("train",), art=True)
        _INP.update(inp=inp, groups=gene_groups(inp.genes, SEEDS[0]))
        jobs = [{"idea": i, "fold": f, "seed": SEEDS[0]} for i in ideas for f in labels.folds()]
        out = []
        with cf.ProcessPoolExecutor(max_workers=a.workers, mp_context=multiprocessing.get_context("fork")) as pool:
            for r in pool.map(_job, jobs, chunksize=1):
                out += r
        f = wc.atomic_write_parquet(pd.concat(out, ignore_index=True), wc.TMP / f"perline_{a.tag}.parquet")
        prov.output(f)
        ru = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_h = sum(r.ru_utime + r.ru_stime for r in ru) / 3600
        wc.log_budget(f"w7_genecv {a.tag} ({','.join(ideas)})", cpu_h, (time.time() - t0) / 3600,
                      f"{len(jobs)} jobs x {N_GROUPS} gene groups, seed {SEEDS[0]}")
        print(json.dumps({"tag": a.tag, "cpu_h": round(cpu_h, 3)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
