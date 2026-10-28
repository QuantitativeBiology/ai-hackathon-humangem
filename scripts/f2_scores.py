"""F2 label-free scores and checks (A15 post hoc; RESEARCH_TREE D18). No label is read.

(1) Gene-map check: e = 1 − ratio_r of M6's P1 ratio sweeps (results/m6/sweeps) through this script's own gene
    map must equal M6's sweep-keyed e table (results/m6/scores/sweep_e.parquet) bitwise.
(2) W3-MEV check for F2 (W7-03 sweeps): the mevalonate core (HMGCR, MVK, PMVK, MVD) must be lethal
    (ratio_r < 0.01) on every M2(f=1) serum signature. Re-derived from the ratio parquet and compared with
    results/w7/sweeps/mev.tsv; P1's core status on the same signatures is shown alongside. Line counts are over
    all 1,361 assigned lines (no tier split).
(3) Validation scores for the look M10_F2: results/m10/val/f2_val_scores.parquet, the 198 validation expression
    lines x U, columns H1, H2, P1 (M6:H1+Bloo/M2(f=1)) and F2 (F2:H1+artBloo/M2(f=1)); H1/H2/P1 asserted bitwise
    equal to results/m6/val/m6_val_scores.parquet.

Outputs: results/f2/mev_check.tsv, results/f2/scores_checks.json, results/m10/val/f2_val_scores.parquet.
Usage: .venv/bin/python scripts/f2_scores.py [--allow-dirty]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f2_common as fc  # noqa: E402

VAL_SCORES = fc.M10VAL / "f2_val_scores.parquet"
VAL_COLUMNS = ["H1", "H2", fc.P1_METHOD, fc.F2_METHOD]


def check_gene_map(U, gmap) -> dict:
    R = pd.read_parquet(fc.M6_SWEEPS).set_index("sweep")
    mine = fc.e_rows(R, gmap, U)
    mine.index = [k.split("|", 1)[1] for k in mine.index]
    ref = fc.p1_sweep_e(U).reindex(mine.index)
    eq = bool(np.array_equal(mine.to_numpy(), ref.to_numpy()))
    if not eq:
        raise AssertionError("own e-from-ratio differs from M6 sweep_e on P1 sweeps")
    return {"gene_map_check_sweeps": int(len(mine)), "gene_map_check_bitwise_equal": eq}


def mev_check(gmap) -> tuple[pd.DataFrame, dict]:
    R = pd.read_parquet(fc.W7_SWEEPS).set_index("sweep")
    sym = fc.symbol_map()
    core = {g: s for g, s in sym.items() if s in fc.MEV_CORE and g in R.columns}
    if sorted(core.values()) != sorted(fc.MEV_CORE):
        raise AssertionError(f"mevalonate core genes not all in the sweeps: {sorted(core.values())}")
    meta = pd.read_csv(fc.W7_SWEEPS_TSV, sep="\t").set_index("key")
    p1e = fc.p1_sweep_e(fc.universe())
    ent = gmap.reindex(list(core)).to_dict()
    sweep_of = fc.line_m2_sweep()
    nlines = sweep_of.value_counts()
    rows = []
    for key in R.index:
        sid = key.split("|", 1)[1]
        m = meta.loc[key]
        if m.step != fc.M2:
            continue
        r = R.loc[key, list(core)].astype(float)
        p1_ratio = 1.0 - p1e.loc[sid, [ent[g] for g in core]].to_numpy()
        rows.append({"w2_sweep_id": sid, "medium_keys": m.medium_keys, "serum": bool(m.serum),
                     "n_lines": int(nlines.get(sid, 0)), "F2_core_max_ratio": float(r.max()),
                     "F2_core_all_lethal": bool((r < fc.LETHAL).all()),
                     "P1_core_min_ratio": float(p1_ratio.min()),
                     "P1_core_all_lethal": bool((p1_ratio < fc.LETHAL).all())})
    t = pd.DataFrame(rows).sort_values(["serum", "medium_keys"], ascending=[False, True], ignore_index=True)
    # agreement with W7's own MEV table
    w7 = pd.read_csv(fc.W7_MEV, sep="\t")
    w7 = w7[w7.step == fc.M2]
    mine = {(k, sym[g]): float(R.loc[k, g]) for k in R.index for g in core if meta.loc[k, "step"] == fc.M2}
    agree = all(abs(mine[(k, g)] - v) < 1e-9 for k, g, v in zip(w7.key, w7.gene, w7.ratio_r)) and len(w7) == len(mine)
    ser = t[t.serum]
    summ = {"mev_serum_signatures": int(len(ser)),
            "mev_F2_serum_signatures_core_lethal": int(ser.F2_core_all_lethal.sum()),
            "mev_F2_serum_core_max_ratio": float(ser.F2_core_max_ratio.max()),
            "mev_F2_serum_free_signatures_core_lethal": f"{int(t[~t.serum].F2_core_all_lethal.sum())}/{int((~t.serum).sum())}",
            "mev_lines_on_serum_signatures": int(ser.n_lines.sum()),
            "mev_lines_core_lethal_F2": int(t.loc[t.F2_core_all_lethal, "n_lines"].sum()),
            "mev_lines_core_lethal_P1": int(t.loc[t.P1_core_all_lethal, "n_lines"].sum()),
            "mev_n_lines_total": int(t.n_lines.sum()),
            "mev_P1_serum_signatures_core_lethal": int(ser.P1_core_all_lethal.sum()),
            "mev_agrees_with_w7_mev_tsv": bool(agree)}
    if not agree:
        raise AssertionError("re-derived core ratios disagree with results/w7/sweeps/mev.tsv")
    return t, summ


def val_scores(U, gmap) -> pd.DataFrame:
    from nextgem import labels, w1

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    V = sorted(l for l in split.index[split.tier == "validation"] if bool(lf.loc[l, "has_expression_matrix"]))
    sweep_of = fc.line_m2_sweep()
    df = pd.DataFrame({"ModelID": np.repeat(V, len(U)), "entrez": np.tile(np.asarray(U, dtype=np.int64), len(V))})
    for s, p in (("H1", "human1"), ("H2", "human2")):
        df[s] = np.tile(fc.generic_e(U, p).to_numpy(float), len(V))
    df[fc.P1_METHOD] = fc.line_matrix(fc.p1_sweep_e(U), sweep_of, V).to_numpy(float).ravel()
    df[fc.F2_METHOD] = fc.line_matrix(fc.f2_sweep_e(U, gmap), sweep_of, V).to_numpy(float).ravel()
    m6 = pd.read_parquet(fc.M6_VAL_SCORES, columns=["ModelID", "entrez", "H1", "H2", fc.P1_METHOD])
    if sorted(m6.ModelID.unique()) != V:
        raise AssertionError("validation expression lines differ from the M6 look's 198")
    m6 = m6.set_index(["ModelID", "entrez"]).reindex(pd.MultiIndex.from_frame(df[["ModelID", "entrez"]]))
    for c in ("H1", "H2", fc.P1_METHOD):
        if not np.array_equal(m6[c].to_numpy(), df[c].to_numpy()):
            raise AssertionError(f"{c}: differs from results/m6/val/m6_val_scores.parquet")
    if not np.isfinite(df[VAL_COLUMNS].to_numpy()).all():
        raise AssertionError("non-finite validation scores")
    return df[["ModelID", "entrez"] + VAL_COLUMNS]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    inputs = [fc.W7_SWEEPS, fc.W7_SWEEPS_TSV, fc.W7_MEV, fc.M6_SWEEPS, fc.M6_SWEEP_E, fc.M6_ASSIGN, fc.M6_VAL_SCORES,
              fc.W2_LINE_SWEEPS]
    with fc.provenance(inputs=inputs, allow_dirty=a.allow_dirty) as prov:
        U, gmap = fc.universe(), fc.gene_map()
        checks = check_gene_map(U, gmap)
        mev, ms = mev_check(gmap)
        checks.update(ms)
        df = val_scores(U, gmap)
        F2v = df[fc.F2_METHOD].to_numpy()
        P1v = df[fc.P1_METHOD].to_numpy()
        checks.update({"val_lines": int(df.ModelID.nunique()), "val_genes": len(U), "val_rows": int(len(df)),
                       "val_H1_H2_P1_equal_M6_bitwise": True,
                       "val_pairs_F2_ne_P1": int((F2v != P1v).sum()),
                       "val_pairs_called_lethal_F2": int((F2v >= 1 - fc.LETHAL).sum()),
                       "val_pairs_called_lethal_P1": int((P1v >= 1 - fc.LETHAL).sum())})
        outs = [fc.write_tsv(mev, fc.F2DIR / "mev_check.tsv"), fc.write_parquet(df, VAL_SCORES)]
        checks["val_scores_sha256"] = fc.sha256(VAL_SCORES)
        outs.append(fc.write_json(checks, fc.F2DIR / "scores_checks.json"))
        prov.output(*outs)
        prov.note(**{k: v for k, v in checks.items() if not isinstance(v, (list, dict))})
        print(json.dumps(checks, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
