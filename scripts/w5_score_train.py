"""W5 training-tier scoring (DESIGN §9 W5, §6, TG): generic vs W5a vs W5b on both parents, LOLO, reusing
``nextgem.w1`` / ``nextgem.metrics``. Labels only through ``nextgem.cv.lolo`` (``labels.fold(k)``) and
``labels.train()``; no validation or test label is read.

Methods (score e = 1 - ratio_r on M0; U genes, Entrez dedup as W1):
  H1, H2                 generic parents (results/baseline; equal to W1's)
  H1_TYMS                H1 with MAR04644 irreversible (the H1-only label-free curation item)
  H1_W5a, H2_W5a         W5a GPRs (H1 on the TYMS base);  H1_W5a_gpr_only: W5a GPRs on the released H1
  H?_W5a_pid30           sensitivity: W5a with paralogs at >= 30% identity (declared before scoring; not shortlist-
                         eligible without a PI decision)
  H1_W5b, H2_W5b         W5a + training-fold co-essentiality: an OR node of a target group whose cross-branch pairs are all
                         non-paralog and co-complex becomes AND if every pair has Chronos r >= 0.4 over >= 30 fitting
                         lines (gene_effect of labels.fold(k)); recomputed inside every LOLO fold
  H?_W5b_loose           sensitivity: co-essentiality + non-paralog only (no complex requirement)
  H?_neg_s<seed>         negative controls (random matched edits, SEEDS[0:5])
Cells: labels L_prob (primary), L_ge05, L_ge10 x gene sets {all, minus CE_fit (LOLO)} x universes {U, U∩}, expression
lines of the 12 eligible training lineages (686). Paired Δ̄_lin with lineage-stratified bootstrap CI (B = 2,000,
SEEDS[0]); TG (§9): Δ̄_lin >= 0.005, CI lower bound > 0, Δ > 0 in >= 9 of 12 folds (< 5 estimable lines = not positive).
Fixed-cutoff binary calls (ratio_r < 0.01) x L_prob MCC on U ∩ parent (secondary).

Usage: OMP_NUM_THREADS=1 .venv/bin/python scripts/w5_score_train.py [--workers 6]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys

import cobra

cobra.Configuration().solver = "glpk"

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from nextgem import metrics, paths, w5  # noqa: E402
from nextgem import w5_evidence as ev  # noqa: E402
from nextgem.config import SEEDS  # noqa: E402

LABEL_COLS = ["gene_effect", "dep_prob", "L_prob", "L_ge05", "L_ge06", "L_ge10"]
LABELS = ("L_prob", "L_ge05", "L_ge10")
GENESETS = ("all", "noCE")
TG_MIN_DELTA = 0.005
TG_MIN_FOLDS = 9
G: dict = {}


def coess_ok(ge: pd.DataFrame, pairs, ent: dict) -> tuple[bool, list]:
    rs = []
    for a, b in pairs:
        ea, eb = ent.get(a), ent.get(b)
        if ea is None or eb is None or ea not in ge.columns or eb not in ge.columns:
            return False, rs
        x = ge[[ea, eb]].dropna()
        if len(x) < ev.COESS_MIN_N:
            return False, rs
        r = float(np.corrcoef(x[ea], x[eb])[0, 1])
        rs.append(round(r, 4))
        if not r >= ev.COESS_R:
            return False, rs
    return bool(pairs), rs


def build(workers: int):
    from nextgem import cv, labels, w1

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    ucap = [g for g, a, b in zip(U.entrez, U.in_human1, U.in_human2) if a and b]
    gen = w1.generic_table(U.entrez)
    line_class = lf.base_medium_class
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    E_all = sorted(elig.index)
    E = [l for l in E_all if lf.loc[l, "has_expression_matrix"]]
    lineage = split.lineage
    # ---- label-free arms
    const, ko, eng, gps, bases, w5a_trees, ent = {}, {}, {}, {}, {}, {}, {}
    for p in w5.PARENTS:
        s = w5.SHORT[p]
        k = pd.read_parquet(w5.w5_dir() / f"ko_{p}.parquet")
        ko[p] = k
        for col in [c for c in k.columns if c.startswith("ratio_")]:
            arm = col[6:]
            e = w5.e_by_entrez(k.gene, k.entrez, k[col], genes)
            name = s if arm == "generic" else (f"{s}_TYMS" if arm == "base" and p == "human1" else f"{s}_{arm}")
            if arm == "base" and p == "human2":
                continue  # H2 base = generic
            const[name] = e
        assert (const[s] == gen[f"e_{s}"].reindex(genes)).all(), f"generic {s} differs from W1"
        m = w5.load_parent(p, tyms_irreversible=(p == "human1"))
        eng[p] = w5.KOEngine(m)
        gps[p] = w5.model_gprs(m)
        bases[p] = k.set_index("gene").ratio_base
        cur = pd.read_csv(w5.w5_dir() / f"curated_W5a_{p}.tsv", sep="\t")
        w5a_trees[p] = {rx: w5.parse_rule(r.curated_rule) for r in cur[cur.edited].itertuples() for rx in r.reactions.split(";")}
        ent[p] = w5.gene_table(p).entrez.dropna().astype(int).to_dict()
    cands = json.loads((w5.w5_dir() / "w5b_candidates.json").read_text())
    need_entrez = sorted({ent[p][g] for p in w5.PARENTS for v in cands[p].values() for c in v["candidates"]
                          for pr in c["pairs"] for g in pr if g in ent[p]})
    fold_log = []

    def fit_predict(fit_long, ids):
        fitted = w1.fit_components(fit_long, line_class)
        ge = fit_long[fit_long.entrez.isin(need_entrez)].pivot(index="ModelID", columns="entrez", values="gene_effect")
        out = {}
        for p in w5.PARENTS:
            s = w5.SHORT[p]
            for variant in ("W5b", "W5b_loose"):
                new = dict(w5a_trees[p])
                for grp, v in cands[p].items():
                    paths_ = []
                    for c in v["candidates"]:
                        if variant == "W5b" and not c["anchored"]:
                            continue
                        ok, rs = coess_ok(ge, [tuple(x) for x in c["pairs"]], ent[p])
                        if ok:
                            paths_.append(tuple(c["path"]))
                            fold_log.append({"fold": lineage.get(ids[0]), "parent": p, "variant": variant, "group": int(grp), "path": c["path"],
                                             "pairs": ";".join("|".join(x) for x in c["pairs"]), "r": rs,
                                             "fit_lines": int(ge.shape[0])})
                    if paths_:
                        t = ev.convert_or_nodes(w5.parse_rule(v["curated_rule"]), paths_)
                        for rx in v["reactions"]:
                            new[rx] = t
                rr = w5.local_sweep(eng[p], gps[p], bases[p], new, workers)
                kk = ko[p]
                out[f"{s}_{variant}"] = w5.e_by_entrez(kk.gene, kk.entrez, rr.reindex(kk.gene).to_numpy(), genes)
        return fitted, ids, out

    folds = cv.lolo(fit_predict, lambda pr, held: pr, columns=LABEL_COLS)
    scores: dict[str, list] = {}
    line_ce = {lab: {} for lab in LABELS}
    for k, (fitted, ids, w5b) in folds.items():
        for l in ids:
            for lab in LABELS:
                line_ce[lab][l] = fitted.ce[lab]
        ids = [i for i in ids if i in set(E_all)]
        for name, e in {**const, **w5b}.items():
            scores.setdefault(name, []).append(pd.DataFrame(np.tile(e.reindex(genes).to_numpy(), (len(ids), 1)),
                                                            index=ids, columns=genes))
    scores = {n: pd.concat(v).loc[E] for n, v in scores.items()}
    train_long = labels.train(columns=LABEL_COLS)
    Y = {lab: w1.wide(train_long, lab, E, genes) for lab in LABELS}
    G.update(dict(scores=scores, Y=Y, genes=genes, ucap=ucap, lineage=lineage, line_ce=line_ce, E=E, gen=gen,
                  fold_log=fold_log, train_long=train_long, U=U))


def ref_of(name: str) -> str:
    return name.split("_")[0]


def run_cell(job):
    from nextgem import w1

    label, geneset, universe = job
    Y = G["Y"][label]
    Y = Y.loc[:, G["genes"] if universe == "U" else G["ucap"]]
    if geneset == "noCE":
        Y = Y.where(w1.ce_mask(Y.index, Y.columns, G["line_ce"][label]))
    rows, lts = [], {}
    for n, S in G["scores"].items():
        row, lt, _ = w1.method_summary(S.loc[Y.index, Y.columns], Y, G["lineage"], gene_unit=False)
        row.update({"label": label, "geneset": geneset, "universe": universe, "method": n})
        rows.append(row)
        lts[n] = lt
    deltas = []
    for n in lts:
        refs = {ref_of(n), "H1"} - {n}
        if n.endswith("_W5b") or n.endswith("_W5b_loose"):
            refs.add(f"{ref_of(n)}_W5a")
        for ref in sorted(refs):
            d = metrics.paired_line_delta(lts[n], lts[ref])
            ds = w1.delta_summary(lts[n], lts[ref], G["lineage"])
            lin = G["lineage"].reindex(d.index)
            per_fold = d.groupby(lin).mean()
            n_est = d.notna().groupby(lin).sum()
            pos = int(((per_fold > 0) & (n_est >= 5)).sum())
            tg = bool(ds["delta_lin"] >= TG_MIN_DELTA and ds["delta_lo"] > 0 and pos >= TG_MIN_FOLDS)
            deltas.append({"label": label, "geneset": geneset, "universe": universe, "method": n, "ref": ref, **ds,
                           "folds_positive": pos, "n_folds": int(per_fold.size), "TG_pass": tg,
                           "n_lines_changed": int((d.abs() > 0).sum()), "max_abs_line_delta": float(d.abs().max())})
    per_line = None
    if label == "L_prob" and universe == "U":
        per_line = pd.concat([lt.assign(method=n) for n, lt in lts.items()]).rename_axis("ModelID").reset_index()
        per_line["geneset"] = geneset
    return rows, deltas, per_line


def binary_rows():
    rows, Y = [], G["Y"]["L_prob"]
    for p, s in (("human1", "H1"), ("human2", "H2")):
        gs = [g for g in G["gen"].index[G["gen"][f"in_{s}"]] if g in Y.columns]
        T = Y.loc[:, gs]
        pls = {}
        for n, S in G["scores"].items():
            if ref_of(n) != s:
                continue
            C = (1.0 - S.loc[T.index, gs]).round(6) < metrics.FBA_CUTOFF
            pl = metrics.binary_per_line(C, T)
            pls[n] = pl
            lin = G["lineage"].reindex(pl.index)
            for m in ("mcc", "f1", "precision", "recall"):
                v = pl[m].dropna()
                bt = metrics.bootstrap_mean(v, lin.reindex(v.index), B=metrics.B_OTHER, seed=SEEDS[0])
                rows.append({"parent": s, "method": n, "metric": m, "lineage_equal_mean": bt.estimate, "ci_lo": bt.ci_low,
                             "ci_hi": bt.ci_high, "n_lines_used": len(v), "n_nan": int(pl[m].isna().sum())})
        for n, pl in pls.items():
            if n == s:
                continue
            d = (pl.mcc - pls[s].mcc).dropna()
            bt = metrics.bootstrap_mean(d, G["lineage"].reindex(d.index), B=metrics.B_OTHER, seed=SEEDS[0])
            rows.append({"parent": s, "method": n, "metric": f"mcc_minus_{s}", "lineage_equal_mean": bt.estimate,
                         "ci_lo": bt.ci_low, "ci_hi": bt.ci_high, "n_lines_used": len(d), "n_nan": 0})
    return pd.DataFrame(rows)


def gene_table_out(gc: pd.DataFrame) -> pd.DataFrame:
    """Training per-gene context for genes whose score changes under a curated arm (W5a, TYMS, W5b in any fold)."""
    Y = G["Y"]["L_prob"]
    frac = Y.mean()
    n = Y.notna().sum()
    S = G["scores"]
    rows = []
    for name in [k for k in S if any(k.endswith(x) for x in ("_W5a", "_W5a_pid30", "_W5a_gpr_only", "_TYMS", "_W5b", "_W5b_loose"))]:
        ref = ref_of(name)
        diff = (S[name] != S[ref]).any()
        for g in diff.index[diff]:
            e_ref = float(S[ref][g].iloc[0])
            e_new = S[name][g]
            rows.append({"method": name, "entrez": int(g), "e_generic": e_ref, "e_new_min": float(e_new.min()),
                         "e_new_max": float(e_new.max()), "n_lines_changed": int((e_new != e_ref).sum()),
                         "frac_pos_train_L_prob": float(frac.get(g, np.nan)), "n_lines": int(n.get(g, 0))})
    df = pd.DataFrame(rows)
    sym = gc.drop_duplicates("entrez").set_index("entrez").symbol if len(gc) else pd.Series(dtype=str)
    if len(df):
        df["symbol"] = df.entrez.map(sym)
        df["direction"] = np.where(df.e_new_max > df.e_generic, "more essential", "less essential")
        df["call_agrees_majority_after"] = np.where(df.direction == "more essential", df.frac_pos_train_L_prob >= 0.5,
                                                    df.frac_pos_train_L_prob < 0.5)
    return df


def main(prov, a) -> int:
    t = w5.Timer()
    build(a.workers)
    print(f"built {len(G['scores'])} methods on {len(G['E'])} lines, {t.lap()}", flush=True)
    jobs = [(lab, gs, "U") for lab in LABELS for gs in GENESETS] + [("L_prob", "all", "Ucap")]
    rows, deltas, pls = [], [], []
    with mp.get_context("fork").Pool(min(a.workers, len(jobs))) as pool:
        for r, d, pl in pool.imap_unordered(run_cell, jobs):
            rows += r
            deltas += d
            if pl is not None:
                pls.append(pl)
    out = w5.w5_dir()
    grid = pd.DataFrame(rows).sort_values(["label", "geneset", "universe", "method"])
    dl = pd.DataFrame(deltas).sort_values(["label", "geneset", "universe", "method", "ref"])
    pl = pd.concat(pls, ignore_index=True)
    pl["lineage"] = pl.ModelID.map(G["lineage"])
    # per-lineage Δ for the primary contrasts
    lin_rows = []
    for gs, sub in pl.groupby("geneset"):
        piv = sub.pivot(index="ModelID", columns="method", values="tap")
        for n in piv.columns:
            for ref in sorted({ref_of(n), "H1"} - {n}):
                d = (piv[n] - piv[ref]).dropna()
                lm = d.groupby(G["lineage"].reindex(d.index)).agg(["mean", "size"])
                for lin, r in lm.iterrows():
                    lin_rows.append({"geneset": gs, "method": n, "ref": ref, "lineage": lin, "delta_mean": r["mean"], "n_lines": int(r["size"])})
    lin_df = pd.DataFrame(lin_rows)
    gc = pd.read_csv(out / "gene_changes.tsv", sep="\t")
    gt = gene_table_out(gc)
    bn = binary_rows()
    fl = pd.DataFrame(G["fold_log"])
    files = [w5.atomic_tsv(grid, out / "train_grid.tsv"), w5.atomic_tsv(dl, out / "train_deltas.tsv"),
             w5.atomic_tsv(lin_df, out / "train_per_lineage_delta.tsv"), w5.atomic_tsv(gt, out / "train_gene_changes.tsv"),
             w5.atomic_tsv(bn, out / "train_binary.tsv"), w5.atomic_tsv(fl, out / "w5b_fold_edits.tsv")]
    prim = dl[(dl.label == "L_prob") & (dl.universe == "U")]
    neg = prim[prim.method.str.contains("_neg_s") & (prim.ref == prim.method.str.split("_").str[0])]
    summ = {
        "n_lines": len(G["E"]), "methods": sorted(G["scores"]),
        "primary": prim[prim.geneset == "all"].set_index(["method", "ref"])[["delta_lin", "delta_lo", "delta_hi", "folds_positive", "TG_pass", "n_lines_changed"]].reset_index().to_dict(orient="records"),
        "noCE": prim[prim.geneset == "noCE"].set_index(["method", "ref"])[["delta_lin", "delta_lo", "delta_hi", "folds_positive", "TG_pass"]].reset_index().to_dict(orient="records"),
        "neg_controls_mean_delta": neg.groupby([neg.method.str.split("_").str[0], "geneset"]).delta_lin.agg(["mean", "min", "max"]).reset_index().to_dict(orient="records"),
        "w5b_fold_edits": int(len(fl)), **t.lap(),
    }
    (out / "train_summary.json").write_text(json.dumps(summ, indent=1, default=float))
    files.append(out / "train_summary.json")
    prov.input(paths.splits_dir() / "lineage_split.tsv", out / "ko_human1.parquet", out / "ko_human2.parquet",
               w5.w5_dir() / "w5b_candidates.json", out / "curated_W5a_human1.tsv", out / "curated_W5a_human2.tsv",
               paths.results_dir() / "w1" / "line_features.tsv")
    prov.output(*files)
    print(json.dumps(summ, indent=1, default=float))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    sys.exit(w5.run_main(lambda prov: main(prov, args), solver="glpk", workers=args.workers, seeds=[SEEDS[0]],
                         reads_labels=True))
