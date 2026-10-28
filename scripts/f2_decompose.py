"""Gene-group decomposition of P1's and F2's training gain over G* (A15 post hoc / exploratory; coordinator request
after D18). Question: does part of P1's gain come from its own W3-MEV artefact? Reads TRAINING labels only, through
nextgem.labels.train() (via f2_score_train.real_labels); aggregates only are written.

Groups (label-free, fixed here before the first label read; disjoint, assigned in this priority order):
  (a) mev_core      HMGCR, MVK, PMVK, MVD (DESIGN mevalonate positive control; W3-MEV)
  (b) sterol        the W3 sterol panel, Human1 (results/w3/sterol/sterol_summary.tsv, panel "sterol": CYP51A1,
                    DHCR24, DHCR7, EBP, FDFT1, HSD17B7, LSS, MSMO1, NSDHL, SC5D, SQLE, TM7SF2)
  (c) cofactor_pool the 35 Human1 calls removed by B-loo:cofactor_pool_biomass on M0
                    (results/w3/review/rv_changed_calls.tsv), minus (a)-(b)
  (d) medium_flip   the W2 flip lists, Human1: results/w2/flips/flip_genes.tsv (single-addition flips on M0) and
                    results/w2/flips/cached_doc_flips.tsv (documented-medium flips, all steps), minus (a)-(c)
  (e) other         every other U gene

Method: for each line and method M in {P1, F2}, the hybrid scorer takes G*'s scores everywhere except the group's
genes, which take M's scores; the group's contribution is Δ̄_lin(hybrid − G*) (lineage-equal, L_prob, all U genes,
B = 2,000 lineage-stratified bootstrap from SEEDS[0], folds positive as in TG). AP is not additive, so the sum of
the five contributions is reported next to the total Δ̄_lin(M − G*), with the residual (total − sum).

Output: results/f2/decomposition.tsv (+ decomposition.json with the group gene lists and counts).
Usage: .venv/bin/python scripts/f2_decompose.py [--synthetic-labels] [--allow-dirty]
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
import f2_score_train as st  # noqa: E402

STEROL_SUMMARY = fc.RES / "w3" / "sterol" / "sterol_summary.tsv"
W3_CALLS = fc.RES / "w3" / "review" / "rv_changed_calls.tsv"
W2_FLIPS = fc.RES / "w2" / "flips" / "flip_genes.tsv"
W2_DOC_FLIPS = fc.RES / "w2" / "flips" / "cached_doc_flips.tsv"
GROUPS = ("mev_core", "sterol", "cofactor_pool", "medium_flip", "other")
METHODS = ("P1", "F2")


def gene_groups(U: list[int]) -> tuple[dict, dict]:
    """Label-free disjoint groups of U Entrez IDs, and their symbols."""
    gmap = fc.gene_map()
    sym = fc.symbol_map()
    sym2ent = {s: int(gmap[g]) for g, s in sym.items() if g in gmap.index and isinstance(s, str)}
    Uset = set(U)

    def ents(symbols):
        miss = [s for s in symbols if s not in sym2ent]
        if miss:
            raise AssertionError(f"symbols not in Human1: {miss}")
        return {sym2ent[s] for s in symbols} & Uset

    st_ = pd.read_csv(STEROL_SUMMARY, sep="\t")
    sterol_syms = sorted(st_.loc[(st_.parent == "human1") & (st_.panel == "sterol"), "gene"].unique())
    calls = pd.read_csv(W3_CALLS, sep="\t")
    cof = calls[(calls.parent == "human1") & (calls.variant == fc.P1_VARIANT) & (calls.change == "removed")]
    if len(cof) != 35:
        raise AssertionError(f"expected the 35 W3 cofactor-pool calls, found {len(cof)}")
    fl = pd.read_csv(W2_FLIPS, sep="\t")
    fl = set(fl.loc[fl.parent == "human1", "entrez"].astype(int))
    dfl = pd.read_csv(W2_DOC_FLIPS, sep="\t")
    dfl = set(gmap.reindex(dfl.loc[dfl.parent == "human1", "gene"]).dropna().astype(int))
    raw = {"mev_core": ents(fc.MEV_CORE), "sterol": ents(sterol_syms),
           "cofactor_pool": set(cof.entrez.astype(int)) & Uset, "medium_flip": (fl | dfl) & Uset}
    out, taken = {}, set()
    for g in GROUPS[:-1]:
        out[g] = sorted(raw[g] - taken)
        taken |= set(out[g])
    out["other"] = sorted(Uset - taken)
    ent2sym = {e: s for s, e in sym2ent.items()}
    return out, {g: sorted(ent2sym.get(e, str(e)) for e in v) for g, v in out.items() if g != "other"}


def run(a) -> tuple[pd.DataFrame, dict]:
    from nextgem.config import SEEDS

    E, lineage, folds, U, S, ver = st.score_inputs()
    groups, syms = gene_groups(U)
    Y, _ = (st.synthetic_labels if a.synthetic_labels else st.real_labels)(E, U, folds, lineage)
    G = S["Gstar"]
    lt_g = fc.per_line_tap(G, Y)
    rows = []
    for m in METHODS:
        total = fc.contrast(fc.per_line_tap(S[m], Y), lt_g, lineage, folds, B=2000, seed=SEEDS[0])
        parts = {}
        for g in GROUPS:
            H = G.copy()
            cols = groups[g]
            H[cols] = S[m][cols].to_numpy()
            r = fc.contrast(fc.per_line_tap(H, Y), lt_g, lineage, folds, B=2000, seed=SEEDS[0])
            parts[g] = r
            changed = int((S[m][cols].to_numpy() != G[cols].to_numpy()).sum())
            rows.append({"method": m, "group": g, "n_genes": len(cols), "n_pairs_changed_vs_Gstar": changed,
                         "contribution": r["delta_lin"], "ci_lo": r["ci"][0], "ci_hi": r["ci"][1],
                         "folds_positive": r["folds_positive"], "share_of_total": r["delta_lin"] / total["delta_lin"]})
        s = sum(p["delta_lin"] for p in parts.values())
        rows.append({"method": m, "group": "sum_of_groups", "n_genes": len(U), "contribution": s,
                     "share_of_total": s / total["delta_lin"]})
        rows.append({"method": m, "group": "total", "n_genes": len(U),
                     "n_pairs_changed_vs_Gstar": int((S[m].to_numpy() != G.to_numpy()).sum()),
                     "contribution": total["delta_lin"], "ci_lo": total["ci"][0], "ci_hi": total["ci"][1],
                     "folds_positive": total["folds_positive"], "share_of_total": 1.0})
        rows.append({"method": m, "group": "residual_total_minus_sum", "n_genes": len(U),
                     "contribution": total["delta_lin"] - s, "share_of_total": (total["delta_lin"] - s) / total["delta_lin"]})
    t = pd.DataFrame(rows)
    # F2 − P1 per group (difference of the two contributions; same lines)
    p = t[t.method == "P1"].set_index("group").contribution
    f = t[t.method == "F2"].set_index("group").contribution
    for g in list(GROUPS) + ["sum_of_groups", "total", "residual_total_minus_sum"]:
        t = pd.concat([t, pd.DataFrame([{"method": "F2-P1", "group": g, "contribution": f[g] - p[g]}])],
                      ignore_index=True)
    meta = {"synthetic_labels": bool(a.synthetic_labels), "label": st.LABEL, "n_lines": len(E),
            "n_lineages": int(lineage.reindex(E).nunique()), "B": 2000, "seed": SEEDS[0],
            "group_sizes": {g: len(v) for g, v in groups.items()}, "group_symbols": syms,
            "verification": {k: v for k, v in ver.items() if "equal" in k}}
    return t, meta


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--synthetic-labels", action="store_true", help="Bernoulli(0.12) labels; no label is read")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    from nextgem.config import SEEDS

    outdir = fc.F2DIR if not a.synthetic_labels else fc.F2DIR / "tmp" / "decompose_synthetic"
    inputs = [fc.W7_SWEEPS, fc.M6_SWEEP_E, fc.M6_ASSIGN, fc.M6_TRAIN_SCORES, fc.F1_DELTA, STEROL_SUMMARY, W3_CALLS,
              W2_FLIPS, W2_DOC_FLIPS]
    with fc.provenance(inputs=inputs, seeds=[SEEDS[0]], allow_dirty=a.allow_dirty) as prov:
        t, meta = run(a)
        outs = [fc.write_tsv(t, outdir / "decomposition.tsv"), fc.write_json(meta, outdir / "decomposition.json")]
        prov.output(*outs)
        prov.note(synthetic_labels=meta["synthetic_labels"], **{f"n_{g}": n for g, n in meta["group_sizes"].items()})
        with pd.option_context("display.width", 200):
            print(t.to_string(index=False))
        print(json.dumps(meta["group_sizes"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
