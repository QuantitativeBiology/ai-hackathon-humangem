"""W7 screening summary (training tier only): seed-averaged per-line tAP -> tAP_lin per method and the §9 TG of each
idea against the current best; descriptive deltas vs P1 and B2. Reads only the per-line tables written by
scripts/w7_screen.py (no label file).

Usage: .venv/bin/python scripts/w7_summarise.py --tags r1_W7-01,r2_W7-02 --best W7-01 --ideas W7-02
Outputs: results/w7/screen_<name>.tsv; prints a compact table.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w7_common as wc  # noqa: E402


def load(tags) -> pd.DataFrame:
    df = pd.concat([pd.read_parquet(wc.TMP / f"perline_{t}.parquet") for t in tags], ignore_index=True)
    return df.drop_duplicates(["method", "seed", "ModelID"], keep="last")


def table(df: pd.DataFrame, method: str, seed=None) -> pd.DataFrame:
    d = df[df.method == method]
    if seed is not None:
        d = d[d.seed == seed]
    g = d.groupby("ModelID")
    t = pd.DataFrame({"n": g.n.first(), "n_pos": g.n_pos.first(), "estimable": g.estimable.first(),
                      "tap": g.tap.mean(), "n_seeds": g.seed.nunique()})
    return t


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", required=True)
    ap.add_argument("--best", default="auto", help="current best method; 'auto' = stronger of P1 and W7-01")
    ap.add_argument("--ideas", default="")
    ap.add_argument("--name", default="")
    a = ap.parse_args(argv)
    from nextgem import labels, metrics, w7
    from nextgem.config import SEEDS

    df = load(a.tags.split(","))
    split = labels.load_split().set_index("ModelID")
    lineage = split.lineage
    folds = labels.folds()
    methods = sorted(df.method.unique())
    T = {m: table(df, m) for m in methods}
    lines = sorted(set.intersection(*[set(t.index) for t in T.values()]))
    T = {m: t.loc[lines] for m, t in T.items()}
    rows = []
    for m in methods:
        t = T[m]
        est = t[t.estimable]
        bt = metrics.bootstrap_mean(est.tap, lineage.reindex(est.index), B=metrics.B_OTHER, seed=SEEDS[0])
        r = {"method": m, "n_lines": len(t), "n_estimable": len(est), "tap_lin": bt.estimate, "tap_lo": bt.ci_low,
             "tap_hi": bt.ci_high, "n_seeds": int(t.n_seeds.max())}
        seeds = sorted(df[df.method == m].seed.unique())
        if len(seeds) > 1:
            per = [w7.lin_mean(table(df, m, s).loc[lines].tap[t.estimable], lineage) for s in seeds]
            r["tap_lin_by_seed"] = json.dumps([round(x, 5) for x in per])
        rows.append(r)
    S = pd.DataFrame(rows).set_index("method")
    best = a.best
    if best == "auto":
        best = "W7-01" if S.loc["W7-01", "tap_lin"] >= S.loc["P1", "tap_lin"] else "P1"
    out = []
    for m in methods:
        r = S.loc[m].to_dict()
        r["method"] = m
        for ref in (best, "P1", "B2", "Gstar"):
            if ref == m or ref not in T:
                continue
            d = metrics.paired_line_delta(T[m], T[ref])
            res = w7.tg(d, lineage, folds, B=metrics.B_OTHER, seed=SEEDS[0])
            tag = "best" if ref == best else ref
            r.update({f"d_{tag}": res["delta_lin"], f"d_{tag}_lo": res["delta_lo"], f"d_{tag}_hi": res["delta_hi"],
                      f"d_{tag}_se": res["delta_se"], f"folds_{tag}": res["folds_positive"]})
            if ref == best:
                r["tg_vs_best"] = res["tg_pass"]
                r["fold_means_vs_best"] = res["fold_means"]
        r["current_best"] = best
        out.append(r)
    O = pd.DataFrame(out)
    name = a.name or a.tags.replace(",", "+")
    f = wc.atomic_write_tsv(O, wc.OUT / f"screen_{name}.tsv")
    cols = [c for c in ("method", "n_estimable", "tap_lin", "tap_lin_by_seed", "d_best", "d_best_lo", "d_best_hi",
                        "folds_best", "tg_vs_best", "d_P1", "d_B2") if c in O]
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.float_format", "{:.4f}".format):
        print(f"current best = {best}")
        print(O[cols].to_string(index=False))
    print(f"wrote {f.relative_to(wc.PROJECT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
