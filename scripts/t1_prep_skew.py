"""T1 prep (7): the §6.5 skew check that fixes the P1 spec's ``bootstrap_method`` for C1 / C1b.

Rule (DESIGN §6.5; ``nextgem.metrics.skew_triggers_bca``): BCa only if the training-tier distribution of line ΔtAP
(P1 − G*) has |skew| > 1; percentile otherwise. The spec has one bootstrap_method for C1 and C1b: it is decided by
the C1 distribution (primary label, all U genes); the C1b distribution (CE_fit removed, fold-scoped by LOLO as in
the M6 training scoring) is reported alongside, and a disagreement is flagged for the PI.

Labels only through ``nextgem.labels.train()`` and ``labels.fold(k)`` (via ``nextgem.cv.lolo``); no validation or
test label. Lines: the W1 primary training set (tier train, non-rare lineage, expression profile; 686 lines,
12 lineages). Scores: results/m6/scores/scores_train.large.parquet (SHA-256 checked), P1 = M6 cell (3)
``M6:H1+Bloo/M2(f=1)``, G* = ``H1``. Cross-check: per-line tAP equals results/m6/train/per_line.tsv (6-digit TSV).

Output: results/m6/t1/skew.json (scalars only: skew, n_lines, lineage-equal mean ΔtAP, decision).
``--synthetic-labels``: Bernoulli(0.12) labels, no label file opened; output under results/m6/tmp/t1/.

Usage: .venv/bin/python scripts/t1_prep_skew.py [--synthetic-labels]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as mc  # noqa: E402
import t1_prep_common as tc  # noqa: E402

SC = mc.SCDIR
P1_COL, G_COL = tc.SOURCE_COLUMN[tc.P1_H1], "H1"
SYNTH_RATE = 0.12


def primary_label() -> str:
    from nextgem import frozen

    for r in frozen.read_rows():
        if r["item"] == "label.primary" and r["kind"] == "value":
            return r["sha256_or_value"]
    raise RuntimeError("FROZEN.tsv has no label.primary row")


def training_lines(split: pd.DataFrame, lf: pd.DataFrame) -> list[str]:
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    return sorted(l for l in elig.index if l in lf.index and bool(lf.loc[l, "has_expression_matrix"]))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--synthetic-labels", action="store_true")
    a = ap.parse_args(argv)
    from nextgem import cv, labels, metrics, w1
    from nextgem.config import SEEDS

    t0 = time.time()
    src = SC / "scores_train.large.parquet"
    pins = pd.read_csv(SC / "SHA256SUMS.tsv", sep="\t").set_index("file").sha256
    ins = [src, SC / "SHA256SUMS.tsv", mc.TRDIR / "per_line.tsv"]
    out_f = (tc.T1_TMP / "skew_synthetic.json") if a.synthetic_labels else (tc.T1 / "skew.json")
    with tc.provenance(inputs=ins, seeds=[SEEDS[0]]) as prov:
        if tc.sha256(src) != pins[tc.rel(src)]:
            raise RuntimeError("scores_train.large.parquet: SHA-256 differs from SHA256SUMS.tsv")
        split = labels.load_split().set_index("ModelID")
        lf = w1.line_features()
        U = list(labels.gene_universe().entrez.astype(int))
        E = training_lines(split, lf)
        lab = primary_label()
        sc = pd.read_parquet(src, columns=["ModelID", "entrez", P1_COL, G_COL])
        sc = sc[sc.ModelID.isin(set(E))]
        S = {m: sc.pivot(index="ModelID", columns="entrez", values=c).reindex(index=E, columns=U)
             for m, c in (("P1", P1_COL), ("G*", G_COL))}
        if a.synthetic_labels:
            rng = np.random.default_rng(SEEDS[0])
            Y = pd.DataFrame((rng.random((len(E), len(U))) < SYNTH_RATE).astype(float), index=E, columns=U)
            line_ce = {l: set() for l in E}
        else:
            def fit_predict(fit_long, ids):
                return w1.fit_components(fit_long, lf.base_medium_class, labels_=(lab,)).ce[lab], ids

            res = cv.lolo(fit_predict, lambda pr, held: pr, columns=["gene_effect", lab])
            line_ce = {l: ce for ce, ids in res.values() for l in ids}
            tl = labels.train(columns=[lab])
            Y = w1.wide(tl, lab, E, U)
            del tl
        Ys = {"C1": Y, "C1b": Y.where(w1.ce_mask(E, U, line_ce))}
        lineage = split.lineage
        out = {"rule": "DESIGN §6.5: BCa for C1/C1b iff |skew| of the training-tier line ΔtAP (P1 − G*) > 1 "
                       "(metrics.skew_triggers_bca); decided on C1, C1b reported",
               "label": lab, "lines": "W1 primary training set", "p1_column": P1_COL, "gstar_column": G_COL,
               "synthetic_labels": bool(a.synthetic_labels)}
        ref = None if a.synthetic_labels else pd.read_csv(mc.TRDIR / "per_line.tsv", sep="\t")
        for k, Yk in Ys.items():
            lt = {m: metrics.per_line(S[m], Yk) for m in S}
            ok = lt["P1"].estimable & lt["G*"].estimable
            d = (lt["P1"].tap - lt["G*"].tap)[ok]
            dec = tc.skew_decision(d)
            dec["n_lineages"] = int(lineage.reindex(d.index).nunique())
            dec["lineage_equal_mean_delta"] = float(metrics.lineage_equal_mean(d, lineage.reindex(d.index)))
            if ref is not None:
                lc = "all" if k == "C1" else "noCE"
                mx = 0.0
                for m, cell in (("P1", "H1+Bloo/M2"), ("G*", "H1")):
                    r = ref[(ref.cell == cell) & (ref.label_cell == lc)].set_index("ModelID").tap
                    mx = max(mx, float((lt[m].tap[ok] - r.reindex(d.index)).abs().max()))
                dec["max_abs_diff_vs_m6_per_line_tsv"] = mx
                if not mx <= 1e-6:
                    raise RuntimeError(f"{k}: per-line tAP differs from results/m6/train/per_line.tsv by {mx}")
            out[k] = dec
        out["bootstrap_method"] = out["C1"]["bootstrap_method"]
        out["c1b_disagrees"] = out["C1b"]["bootstrap_method"] != out["C1"]["bootstrap_method"]
        out["wall_s"] = round(time.time() - t0, 1)
        tc.atomic_write_json(out, out_f)
        prov.output(out_f)
        prov.note(label_api="nextgem.labels.train()/fold(k) only" if not a.synthetic_labels else "synthetic",
                  bootstrap_method=out["bootstrap_method"])
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
