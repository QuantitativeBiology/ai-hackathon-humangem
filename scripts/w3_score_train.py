"""W3 on the training tier (DESIGN §6, §9 W3, §9 TG; M4 training work). Labels only through ``nextgem.labels``
(``train()`` and ``fold(k)`` via ``nextgem.cv.lolo``). No validation or test label is read.

Reuses ``nextgem.w1`` (fit_components, method_summary, delta_summary, b8, binary_endpoints, ce_mask) and
``nextgem.metrics`` (tAP, bootstrap, sign-flip); nothing is re-implemented.

Evaluation lines: training lines of the 12 eligible lineages with an expression profile (686; W1 primary set).
Primary cell: L_prob, all U genes, no stratum; secondary: minus CE_fit (LOLO).
Reference: T0 = generic parent, biomass only, M0 (e = 1 − ratio_r; equals W1's generic score, asserted).

Modes (DESIGN W3 table; final a-priori list §5 sets as ``L-*``; README §"modes"):
  criterion: T1 (GR ∪ ER), T2 (GR ∪ ER ∪ IC ∪ BS), T3 (all 57 = P57; equals W1's task score, asserted),
             T4_k (GR ∪ functional group k, 11 groups), T5 (GR ∪ CellFie tasks active in the line; thresholds from
             fitting-line expression, per LOLO fold), L-T3 (8 BS nucleotide), L-T4_j (12 IC singly),
             L-T5 / L-T5+ (membrane), L-T6 / L-T6+ (cofactors), L-T8 (Full amino-acid tasks), L-T10 (T10 biomass +
             L-T3 constraints)
  sweeps:    D-BS, D-ER, O-NGAM_k, B-loo_j, B-add_c (@1e-4; @1e-3 sensitivity), O-cap, O-MOMA (if kept)
  medium:    every mode also on the RPMI M2(f = 1) proxy (secondary; Δ vs T0 on that medium)
Controls: pseudo-task twins (SEEDS[0:5]) for every Essential-task criterion mode; matched call rate (the same
number of added / removed calls given to the uncalled genes with the highest / called genes with the lowest
fold-scoped B2 [L_prob]); ETF capability-only positive control (label-free).
TG (§9): Δ̄_lin vs T0 >= 0.005, paired bootstrap CI lower bound > 0, Δ > 0 in >= 9 of 12 folds (fold with < 5
estimable lines = not positive). Shortlist (W3 decision): <= 2 candidates that pass TG and beat both negative
controls, chosen with the paired-difference one-SE rule (simplest first).
Context claims: line-level contexts (base medium, serum-free, background groups) per mode: one-sided sign-flip
p (10,000 flips, SEEDS[0]) on the lineage-equal Δ̄_lin, BH q < 0.05 across all tests (training only; a claim
also needs the same sign on validation, which is NOT looked at here).
A8 (exploratory): gene-pFBA KO-cost tie-break vs ratio_r (T0) and vs B8± on the same parent.

Outputs (results/w3/train/): modes.tsv, mode_scores.tsv (summary + Δ + TG + controls), per_line.tsv (primary cell),
per_fold.tsv, context_claims.tsv, per_task.tsv, binary.tsv, gpfba.tsv, shortlist.json, train_summary.json.

Usage: OMP_NUM_THREADS=1 .venv/bin/python scripts/w3_score_train.py --workers 7
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w3_common import OUT_DIR, atomic_write_tsv, current, run_main  # noqa: E402

TRDIR = OUT_DIR / "train"
LABEL_COLS = ["gene_effect", "dep_prob", "L_prob", "L_ge05", "L_ge06", "L_ge10"]
PARENTS = ("human1", "human2")
SHORT = {"human1": "H1", "human2": "H2"}
TG_MIN_DELTA = 0.005
TG_MIN_FOLDS = 9
CALL = 0.99  # e >= 0.99 <=> ratio_r < 0.01 or a task failure (the §6.1 binary call)
CONTEXT_MIN_LINES = 10
G: dict = {}


# ============================================================================= label-free mode construction
def entrez_map(parent: str) -> pd.Series:
    """Model gene -> Entrez as in w1.generic_table (baseline parquet order, NaN dropped, first per Entrez)."""
    from nextgem import paths

    ko = pd.read_parquet(paths.baseline_dir() / f"{parent}_ko.parquet", columns=["gene", "entrez"]).dropna(subset=["entrez"])
    ko = ko.assign(entrez=ko.entrez.astype(int)).drop_duplicates("entrez")
    return ko.set_index("gene").entrez


def to_U(s: pd.Series, emap: pd.Series, U: list, fill=np.nan) -> pd.Series:
    """Gene-indexed series -> Entrez-indexed over U (genes without an Entrez / duplicates dropped as in W1)."""
    v = s.reindex(emap.index)
    out = pd.Series(v.to_numpy(), index=emap.to_numpy()).reindex(U)
    return out.fillna(fill) if fill is not None else out


def fail_sets(parent: str, which: str, seed: int = -1) -> tuple[dict, pd.DataFrame]:
    """task key -> set of model genes whose KO makes the task infeasible; plus the task meta table."""
    d = OUT_DIR / "tasks"
    long = pd.read_parquet(d / f"{which}_{parent}.parquet", columns=["task", "seed", "gene", "feasible", "capacity"])
    long = long[long.seed == seed]
    meta = pd.read_csv(d / f"{which}_{parent}_tasks.tsv", sep="\t")
    meta = meta[meta.seed == seed] if "seed" in meta else meta
    fs = {t: set(g.gene[~g.feasible]) for t, g in long.groupby("task")}
    return fs, meta, long


def build_modes(U: list, ucap: list) -> tuple[dict, pd.DataFrame, dict]:
    """All line-constant modes per parent: {(parent, mode): e Series over U}; mode table; extras for T5/pseudo."""
    from nextgem import w3

    E, rows, extra = {}, [], {}
    for p in PARENTS:
        s = SHORT[p]
        em = entrez_map(p)
        sw = pd.read_parquet(OUT_DIR / "sweeps" / f"ratio_r_{p}.parquet").set_index("sweep")
        fs_e, meta_e, long_e = fail_sets(p, "essential")
        fs_f, meta_f, _ = fail_sets(p, "full_subset")
        meta_e = meta_e.set_index("task")
        gr = meta_e.index[meta_e.category == "GR"][0]
        cap = long_e.pivot(index="gene", columns="task", values="capacity")
        in_p = to_U(pd.Series(True, index=em.index), em, U, fill=False).astype(bool)
        extra[p] = {"emap": em, "in_p": in_p, "gr": gr, "fs_e": fs_e, "meta_e": meta_e, "fs_f": fs_f}

        def add(name, e_gene, family, desc, n_items, medium="M0", fails=None):
            e = w3.mode_e(e_gene, fails)
            eu = to_U(e, em, U, fill=0.0)
            E[(p, name if medium == "M0" else f"{name}@RPMI")] = eu.round(6)
            rows.append({"parent": p, "mode": name if medium == "M0" else f"{name}@RPMI", "family": family,
                         "medium": medium, "description": desc, "n_items": n_items})

        def fails_of(keys, fs):
            out = set()
            for k in keys:
                out |= fs[k]
            return pd.Series(True, index=sorted(out), dtype=bool)

        cats = meta_e.category
        groups = meta_e.group
        ess_keys = list(meta_e.index)
        nongr = [k for k in ess_keys if k != gr]
        crit = {
            "T1": ([gr] + [k for k in nongr if cats[k] == "ER"], "biomass + ER"),
            "T2": ([gr] + [k for k in nongr if cats[k] in ("ER", "IC", "BS")], "biomass + ER + IC + BS"),
            "T3": ([gr] + nongr, "biomass + all 56 non-GR tasks (= P57)"),
        }
        for gi, gname in enumerate([g for g in dict.fromkeys(groups[nongr])], 1):
            crit[f"T4_{gi}"] = ([gr] + [k for k in nongr if groups[k] == gname], f"biomass + group {gi}: {gname}")
        crit["L-T3"] = ([gr] + [k for k in nongr if cats[k] == "BS" and "de novo synthesis" in meta_e.description[k]
                            and meta_e.description[k].split()[0] in ("ATP", "CTP", "GTP", "UTP", "dATP", "dCTP", "dGTP", "dTTP")],
                        "list T3: biomass + 8 BS nucleotide tasks")
        ic = [k for k in nongr if cats[k] == "IC"]
        for j, k in enumerate(ic, 1):
            crit[f"L-T4_{j}"] = ([gr, k], f"list T4_{j}: biomass + IC task {meta_e.description[k]}")
        memb = [k for k in nongr if (cats[k] == "BS" and any(x in meta_e.description[k] for x in
                ("Phosphatidylcholine", "Phosphatidylethanolamine", "Phosphatidylserine", "Phosphatidylinositol", "Cholesterol")))
                or (cats[k] == "SU" and meta_e.description[k] in ("Choline uptake", "Inositol uptake"))]
        crit["L-T5"] = ([gr] + memb, "list T5: biomass + PC/PE/PS/PI/cholesterol + choline/inositol uptake")
        cof = [k for k in nongr if cats[k] == "BS" and any(x in meta_e.description[k] for x in ("Thiamin", "Coenzyme A", "FAD", "Heme"))]
        crit["L-T6"] = ([gr] + cof, "list T6: biomass + TPP/CoA/FAD/heme")
        for name, (keys, desc) in crit.items():
            add(name, sw.loc["M0|T0"], "criterion", desc, len(keys) - 1, fails=fails_of(keys, fs_e))
            add(name, sw.loc["RPMI_M2f1|T0"], "criterion", desc, len(keys) - 1, medium="RPMI", fails=fails_of(keys, fs_e))
        # Full-list sets
        meta_f = pd.read_csv(OUT_DIR / "tasks" / f"full_subset_{p}_tasks.tsv", sep="\t").set_index("task")
        fid = {k: int(meta_f.category[k]) for k in meta_f.index if meta_f.status[k] == "ok"}
        t5p = [k for k, i in fid.items() if 104 <= i <= 111]
        t6p = [k for k, i in fid.items() if 112 <= i <= 118]
        t8 = [k for k, i in fid.items() if 42 <= i <= 99]
        extra[p]["full_dropped"] = [k for k in meta_f.index if meta_f.status[k] != "ok"]
        for name, ek, fk, desc in (("L-T5+", crit["L-T5"][0], t5p, "list T5 + Full 104-111"),
                                   ("L-T6+", crit["L-T6"][0], t6p, "list T6 + Full 112-118"),
                                   ("L-T8", [gr], t8, "list T8: biomass + Full 42-50, 52-68, 77-99 (non-should-fail)")):
            f = pd.concat([fails_of(ek, fs_e), fails_of(fk, fs_f)])
            f = f[~f.index.duplicated()]
            add(name, sw.loc["M0|T0"], "criterion", desc, len(ek) - 1 + len(fk), fails=f)
            add(name, sw.loc["RPMI_M2f1|T0"], "criterion", desc, len(ek) - 1 + len(fk), medium="RPMI", fails=f)
        extra[p]["criterion"] = {k: v[0] for k, v in crit.items()}
        # T0 and sweep modes
        for med, tag in (("M0", "M0"), ("RPMI", "RPMI_M2f1")):
            add("T0", sw.loc[f"{tag}|T0"], "reference", "biomass only", 0, medium=med)
            for sid in sw.index:
                mm, v = sid.split("|", 1)
                if mm != tag or v == "T0":
                    continue
                if v.startswith("L-T10"):
                    add("L-T10", sw.loc[sid], "criterion", "list T10: de-artefacted biomass + list-T3 constraints", 9,
                        medium=med, fails=fails_of(crit["L-T3"][0], fs_e))
                    continue
                fam = v.split(":")[0].split("_")[0]
                if v.endswith("@1e-3"):
                    fam = "B-add-sens"  # DESIGN: 1e-3 is a sensitivity, never a shortlist candidate
                add(v, sw.loc[sid], fam, v, 1, medium=med)
            # O-cap: e = 1 - min(ratio_biomass, min over T3 tasks of the capacity ratio)
            r0 = sw.loc[f"{tag}|T0"]
            capmin = cap[nongr].min(axis=1).reindex(r0.index).fillna(1.0)
            add("O-cap", np.minimum(r0, capmin), "O-cap", "1 - min(ratio, min_t capacity_t)", 56,
                medium=med)
        # O-MOMA (M0 only, if run and kept)
        mf = OUT_DIR / "moma" / f"moma_{p}.parquet"
        if mf.exists():
            mo = pd.read_parquet(mf).set_index("gene").ratio_r
            add("O-MOMA", mo, "O-MOMA", "QP MOMA biomass ratio", 0)
        # A8 exploratory
        gf = OUT_DIR / "gpfba" / f"gpfba_{p}.parquet"
        if gf.exists():
            sc = pd.read_parquet(gf).set_index("gene").score
            add("A8-gpFBA", sc, "exploratory", "gene-pFBA KO-cost tie-break (A8; exploratory)", 0)
        # pseudo twins (Essential criterion modes, M0)
        pf = OUT_DIR / "tasks" / f"pseudo_{p}.parquet"
        if pf.exists():
            pl = pd.read_parquet(pf, columns=["task", "seed", "gene", "feasible"])
            pmeta = pd.read_csv(OUT_DIR / "tasks" / f"pseudo_{p}_tasks.tsv", sep="\t")
            extra[p]["pseudo_meta"] = pmeta
            for sd, g in pl.groupby("seed"):
                pfs = {t: set(x.gene[~x.feasible]) for t, x in g.groupby("task")}
                for name, keys in extra[p]["criterion"].items():
                    twin = [k for k in keys if k != gr]
                    f = fails_of([gr], fs_e)
                    have = [k for k in twin if k in pfs]
                    if have:
                        f = pd.concat([f, fails_of(have, pfs)])
                        f = f[~f.index.duplicated()]
                    e = w3.mode_e(sw.loc["M0|T0"], f)
                    E[(p, f"PSEUDO[{sd}]:{name}")] = to_U(e, em, U, fill=0.0).round(6)
                    rows.append({"parent": p, "mode": f"PSEUDO[{sd}]:{name}", "family": "control_pseudo", "medium": "M0",
                                 "description": f"pseudo twin of {name}, seed {sd}", "n_items": len(have)})
    return E, pd.DataFrame(rows), extra


def cellfie_T5_rows(p: str, lines, fit_lines, U, extra, expr) -> pd.DataFrame:
    """T5 score rows for ``lines``: GR ∪ CellFie tasks active in the line (thresholds from ``fit_lines``)."""
    from nextgem import context, w3

    ex = extra[p]
    if "cellfie" not in ex:
        fs_c, meta_c, _ = fail_sets(p, "cellfie")
        ok = meta_c[meta_c.status == "ok"].set_index("task")
        ess = {t: json.loads(v) for t, v in ok.ess_rxns.items() if isinstance(v, str)}
        import cobra

        cobra.Configuration().solver = "glpk"
        from nextgem.model import load_model

        m = load_model(p)
        gprs = context.compile_gprs(m)
        ex["cellfie"] = (fs_c, ess, gprs)
    fs_c, ess, gprs = ex["cellfie"]
    em = ex["emap"]
    thr = w3.cellfie_thresholds(expr, [l for l in fit_lines if l in expr.index])
    gs_all = w3.cellfie_gene_scores(expr, lines, thr)
    r0 = pd.read_parquet(OUT_DIR / "sweeps" / f"ratio_r_{p}.parquet").set_index("sweep").loc["M0|T0"]
    gr_fail = ex["fs_e"][ex["gr"]]
    rows = []
    for l in lines:
        gsc = {}
        v = gs_all.loc[l]
        for g, e in em.items():
            x = v.get(e, np.nan)
            if pd.notna(x):
                gsc[g] = float(x)
        rs = {}
        active = []
        for t, rx in ess.items():
            vals = []
            for r in rx:
                if r in gprs:
                    if r not in rs:
                        rs[r] = w3.eval_gpr_score(gprs[r], gsc)
                    if not np.isnan(rs[r]):
                        vals.append(rs[r])
            if vals and np.mean(vals) >= w3.CELLFIE_ACTIVE:
                active.append(t)
        f = set(gr_fail)
        for t in active:
            f |= fs_c[t]
        e = w3.mode_e(r0, pd.Series(True, index=sorted(f), dtype=bool))
        rows.append(to_U(e, em, U, fill=0.0).round(6).rename(l))
        G.setdefault("t5_active_counts", {}).setdefault(p, {})[l] = len(active)
    return pd.DataFrame(rows)


# ============================================================================= label-side build
def build():
    from nextgem import cv, labels, w1

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    Udf = labels.gene_universe()
    U = list(Udf.entrez.astype(int))
    ucap = [g for g, a, b in zip(Udf.entrez, Udf.in_human1, Udf.in_human2) if a and b]
    omics = w1.load_omics()
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    E_all = sorted(elig.index)
    E = [l for l in E_all if lf.loc[l, "has_expression_matrix"]]
    lineage = split.lineage
    modes, mode_table, extra = build_modes(U, ucap)
    gen = w1.generic_table(Udf.entrez)
    # equality: T0 == W1 generic e; T3 == W1 task-augmented e
    checks = {}
    for p in PARENTS:
        s = SHORT[p]
        checks[f"T0_equals_w1_e_{s}"] = bool(np.allclose(modes[(p, "T0")].to_numpy(), gen[f"e_{s}"].reindex(U).to_numpy(), atol=0, rtol=0))
        checks[f"T3_equals_w1_task_e_{s}"] = bool(np.allclose(modes[(p, "T3")].to_numpy(), gen[f"et_{s}"].reindex(U).to_numpy(), atol=0, rtol=0))
    if not all(checks.values()):
        raise RuntimeError(f"equality with W1 failed: {checks}")
    # LOLO: fold-scoped components (B2 prior for the matched-call-rate control, CE_fit, T5 thresholds)
    line_class = lf.base_medium_class

    def fit_predict(fit_long, ids):
        return w1.fit_components(fit_long, line_class), sorted(fit_long.ModelID.unique()), ids

    folds = cv.lolo(fit_predict, lambda pr, held: pr, columns=LABEL_COLS)
    fold_of, prior, ce, t5 = {}, {}, {}, {p: [] for p in PARENTS}
    for k, (fitted, fit_lines, ids) in folds.items():
        ids = [i for i in ids if i in set(E)]
        for l in ids:
            fold_of[l] = k
        prior[k] = fitted.prior["L_prob"]
        ce[k] = fitted.ce["L_prob"]
        for p in PARENTS:
            t5[p].append(cellfie_T5_rows(p, ids, fit_lines, U, extra, omics["expr"]))
    for p in PARENTS:
        modes_line = pd.concat(t5[p]).loc[E]
        G.setdefault("line_modes", {})[(p, "T5")] = modes_line
        mode_table.loc[len(mode_table)] = {"parent": p, "mode": "T5", "family": "criterion", "medium": "M0",
                                           "description": "biomass + CellFie-inspired active tasks (line-specific, LOLO thresholds)",
                                           "n_items": np.nan}
    train_long = labels.train(columns=LABEL_COLS)
    Y = {lab: w1.wide(train_long, lab, E, U) for lab in ("L_prob", "L_ge06")}
    line_ce = {l: ce[fold_of[l]] for l in E}
    Y["L_prob_noCE"] = Y["L_prob"].where(w1.ce_mask(E, U, line_ce))
    G.update(dict(split=split, lf=lf, U=U, ucap=ucap, E=E, lineage=lineage, modes=modes, mode_table=mode_table,
                  extra=extra, Y=Y, fold_of=fold_of, prior=prior, folds=list(folds), checks=checks, gen=gen,
                  expr=omics["expr"]))


# ============================================================================= scoring
def score_matrix(key):
    p, mode = key
    if key in G.get("line_modes", {}):
        return G["line_modes"][key]
    if mode.startswith("MCR:"):
        return mcr_matrix(p, mode[4:], symmetric=False)
    if mode.startswith("MCRS:"):
        return mcr_matrix(p, mode[5:], symmetric=True)
    if mode.startswith("B8"):
        from nextgem import w1

        X = w1._expr(G["expr"], G["E"], G["U"])
        return w1.b8(G["modes"][(p, "T0")], X, +1 if mode.startswith("B8+") else -1)
    e = G["modes"][key]
    return pd.DataFrame(np.tile(e.to_numpy(), (len(G["E"]), 1)), index=G["E"], columns=G["U"])


def mcr_matrix(p: str, mode: str, symmetric: bool = False) -> pd.DataFrame:
    """Matched call rate (DESIGN W3, literal = primary): T0 plus the same number of *extra* calls as the mode, given
    to the uncalled parent genes with the highest fold-scoped B2 (L_prob). ``symmetric`` (ours, sensitivity): also
    remove as many calls as the mode removes, from the called genes with the lowest B2."""
    U, E = G["U"], G["E"]
    t0 = G["modes"][(p, "T0")]
    in_p = G["extra"][p]["in_p"]
    rows = {}
    if (p, mode) in G.get("line_modes", {}):
        M = G["line_modes"][(p, mode)]
        per_line = {l: M.loc[l] for l in E}
    else:
        per_line = None
        x = G["modes"][(p, mode)]
    cache = {}
    for l in E:
        x_l = per_line[l] if per_line is not None else x
        c0, c1 = t0 >= CALL, x_l >= CALL
        n_add, n_drop = int((c1 & ~c0).sum()), (int((c0 & ~c1).sum()) if symmetric else 0)
        k = G["fold_of"][l]
        key = (k, n_add, n_drop)
        if key not in cache:
            b2 = G["prior"][k].reindex(U).fillna(0.0)
            e = t0.copy()
            if n_add:
                cand = b2[(~c0) & in_p].sort_values(ascending=False, kind="mergesort")
                e[cand.index[:n_add]] = 1.0
            if n_drop:
                cand = b2[c0 & in_p].sort_values(ascending=True, kind="mergesort")
                e[cand.index[:n_drop]] = 0.0
            cache[key] = e.to_numpy()
        rows[l] = cache[key]
    return pd.DataFrame.from_dict(rows, orient="index", columns=U).loc[E]


def run_job(job):
    """(parent, mode, ycell) -> per-line table (and gene-unit summary for main modes)."""
    from nextgem import metrics, w1

    p, mode, ycell, full = job
    S = score_matrix((p, mode))
    Y = G["Y"][ycell]
    S = S.loc[Y.index, Y.columns]
    if full:
        row, lt, gt = w1.method_summary(S, Y, G["lineage"], gene_unit=True)
    else:
        lt = metrics.per_line(S, Y)
        row = {}
    return job, row, lt


def fold_stats(lt_a, lt_b) -> dict:
    from nextgem import metrics

    d = metrics.paired_line_delta(lt_a, lt_b)
    lin = G["lineage"].reindex(d.index)
    per = d.groupby(lin).agg(["mean", "size"])
    folds = G["folds"]
    n_pos = sum(1 for f in folds if f in per.index and per.loc[f, "size"] >= 5 and per.loc[f, "mean"] > 0)
    return {"folds_positive": n_pos, "n_folds": len(folds),
            "fold_means": json.dumps({f: (round(float(per.loc[f, "mean"]), 6) if f in per.index else None) for f in folds})}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    a = ap.parse_args(argv)
    t0 = time.time()
    from nextgem import metrics, w1
    from nextgem.config import SEEDS

    build()
    print(f"built: {len(G['E'])} lines, {len(G['modes'])} line-constant modes, {time.time() - t0:.0f} s", flush=True)
    mt = G["mode_table"]
    main_modes = [(r.parent, r.mode) for r in mt.itertuples() if r.family != "control_pseudo"]
    jobs = [(p, m, "L_prob", True) for p, m in main_modes]
    jobs += [(p, m, "L_prob_noCE", False) for p, m in main_modes]
    jobs += [(r.parent, r.mode, "L_prob", False) for r in mt.itertuples() if r.family == "control_pseudo"]
    m0_modes = [(p, m) for p, m in main_modes if not m.endswith("@RPMI") and m not in ("T0", "A8-gpFBA")]
    jobs += [(p, f"MCR:{m}", "L_prob", False) for p, m in m0_modes]
    jobs += [(p, f"MCRS:{m}", "L_prob", False) for p, m in m0_modes]
    jobs += [(p, b, "L_prob", False) for p in PARENTS for b in ("B8+", "B8-")]
    res = {}
    with mp.get_context("fork").Pool(a.workers) as pool:
        for i, (job, row, lt) in enumerate(pool.imap_unordered(run_job, jobs, chunksize=2), 1):
            res[job[:3]] = (row, lt)
            if i % 50 == 0:
                print(f"{i}/{len(jobs)} jobs, {time.time() - t0:.0f} s", flush=True)
    lineage = G["lineage"]
    out = []
    for p, m in main_modes:
        row, lt = res[(p, m, "L_prob")]
        ref = "T0@RPMI" if m.endswith("@RPMI") else "T0"
        r = {"parent": p, "mode": m, **{k: v for k, v in row.items()}}
        r.update(w1.delta_summary(lt, res[(p, ref, "L_prob")][1], lineage))
        r.update(fold_stats(lt, res[(p, ref, "L_prob")][1]))
        r["delta_noCE"] = w1.delta_summary(res[(p, m, "L_prob_noCE")][1], res[(p, ref, "L_prob_noCE")][1], lineage)["delta_lin"]
        r["tg_pass"] = bool(r["delta_lin"] >= TG_MIN_DELTA and r["delta_lo"] > 0 and r["folds_positive"] >= TG_MIN_FOLDS) if m != ref else False
        if (p, f"MCR:{m}", "L_prob") in res:
            d = w1.delta_summary(lt, res[(p, f"MCR:{m}", "L_prob")][1], lineage)
            r.update({"vs_mcr_delta": d["delta_lin"], "vs_mcr_lo": d["delta_lo"], "vs_mcr_hi": d["delta_hi"],
                      "mcr_delta_vs_T0": w1.delta_summary(res[(p, f"MCR:{m}", "L_prob")][1], res[(p, "T0", "L_prob")][1], lineage)["delta_lin"]})
        if (p, f"MCRS:{m}", "L_prob") in res:
            d = w1.delta_summary(lt, res[(p, f"MCRS:{m}", "L_prob")][1], lineage)
            r.update({"vs_mcr_sym_delta": d["delta_lin"], "vs_mcr_sym_lo": d["delta_lo"], "vs_mcr_sym_hi": d["delta_hi"]})
        ps = [res[k][1] for k in res if k[0] == p and k[1].startswith("PSEUDO[") and k[1].split(":", 1)[1] == m and k[2] == "L_prob"]
        if ps:
            pm = ps[0].copy()
            pm["tap"] = np.mean([x["tap"].to_numpy() for x in ps], axis=0)
            d = w1.delta_summary(lt, pm, lineage)
            pd_vs_t0 = [w1.delta_summary(x, res[(p, "T0", "L_prob")][1], lineage)["delta_lin"] for x in ps]
            r.update({"vs_pseudo_delta": d["delta_lin"], "vs_pseudo_lo": d["delta_lo"], "vs_pseudo_hi": d["delta_hi"],
                      "pseudo_delta_vs_T0_mean": float(np.mean(pd_vs_t0)), "pseudo_delta_vs_T0_max": float(np.max(pd_vs_t0)),
                      "n_pseudo_seeds": len(ps)})
        if m == "A8-gpFBA":
            for b in ("B8+", "B8-"):
                d = w1.delta_summary(lt, res[(p, b, "L_prob")][1], lineage)
                r.update({f"vs_{b}_delta": d["delta_lin"], f"vs_{b}_lo": d["delta_lo"], f"vs_{b}_hi": d["delta_hi"]})
        out.append(r)
    scores = pd.DataFrame(out).merge(mt, on=["parent", "mode"], how="left")
    scores["beats_mcr"] = scores.vs_mcr_lo > 0
    scores["beats_pseudo"] = np.where(scores.vs_pseudo_delta.notna(), scores.vs_pseudo_lo > 0, np.nan)
    scores["beats_controls"] = scores.beats_mcr & scores.beats_pseudo.fillna(True).astype(bool)
    TRDIR.mkdir(parents=True, exist_ok=True)
    atomic_write_tsv(mt, TRDIR / "modes.tsv")
    atomic_write_tsv(scores, TRDIR / "mode_scores.tsv")
    pl = pd.concat([res[(p, m, "L_prob")][1].assign(parent=p, mode=m) for p, m in main_modes]).rename_axis("ModelID").reset_index()
    pl["lineage"] = pl.ModelID.map(lineage)
    atomic_write_tsv(pl[["parent", "mode", "ModelID", "lineage", "n", "n_pos", "estimable", "tap", "e_random"]], TRDIR / "per_line.tsv")
    # ---------------------------------------------------------------- context claims (BH over all tests)
    lm = pd.read_csv(OUT_DIR.parent / "w2" / "media_library" / "line_media.tsv", sep="\t").set_index("ModelID")
    bg = pd.read_csv(OUT_DIR.parent / "w1" / "background_members.tsv", sep="\t").set_index("ModelID")
    E = G["E"]
    ctxs = {}
    for c, s in lm.reindex(E).base_class.dropna().groupby(lm.reindex(E).base_class.dropna()):
        ctxs[f"base_medium={c}"] = list(s.index)
    sf = lm.reindex(E).serum_free
    ctxs["serum_free=True"] = list(sf.index[sf == True])  # noqa: E712
    ctxs["serum_free=False"] = list(sf.index[sf == False])  # noqa: E712
    for gcol in bg.columns:
        mem = bg.reindex(E)[gcol]
        ctxs[f"background={gcol}"] = list(mem.index[mem == 1])
    claims = []
    for (p, m) in main_modes:
        if m in ("T0", "T0@RPMI"):
            continue
        ref = "T0@RPMI" if m.endswith("@RPMI") else "T0"
        d = metrics.paired_line_delta(res[(p, m, "L_prob")][1], res[(p, ref, "L_prob")][1])
        for cname, members in ctxs.items():
            dd = d.reindex([l for l in members if l in d.index]).dropna()
            lin = lineage.reindex(dd.index)
            if len(dd) < CONTEXT_MIN_LINES or lin.nunique() < 2:
                continue
            if np.allclose(dd.to_numpy(), 0):
                obs, pv = 0.0, 1.0
            else:
                obs, pv = metrics.sign_flip_test(dd, lin, seed=SEEDS[0])
            claims.append({"parent": p, "mode": m, "context": cname, "n_lines": len(dd), "n_lineages": lin.nunique(),
                           "delta_lin": obs, "p_one_sided": pv})
    cl = pd.DataFrame(claims)
    if len(cl):
        from scipy.stats import false_discovery_control

        cl["q_bh"] = false_discovery_control(cl.p_one_sided.to_numpy(), method="bh")
        cl["claim_training"] = cl.q_bh < 0.05
    atomic_write_tsv(cl, TRDIR / "context_claims.tsv")
    # ---------------------------------------------------------------- per-task table (57 Essential tasks)
    per_task = per_task_table(res)
    atomic_write_tsv(per_task, TRDIR / "per_task.tsv")
    # ---------------------------------------------------------------- shortlist
    sl = shortlist(scores, res)
    # ---------------------------------------------------------------- binary endpoints (T0, T3, shortlist) on L_ge06 / L_prob
    bin_rows = binary_rows(scores, sl)
    atomic_write_tsv(bin_rows, TRDIR / "binary.tsv")
    # ---------------------------------------------------------------- positive controls (label-free)
    pos = positive_controls()
    summary = {"n_lines": len(E), "n_modes_scored": len(main_modes), "n_jobs": len(jobs), "checks": G["checks"],
               "n_context_tests": int(len(cl)), "n_context_claims": int(cl.claim_training.sum()) if len(cl) else 0,
               "shortlist": sl, "positive_controls": pos, "wall_s": round(time.time() - t0, 1),
               "t5_active_tasks_per_line": {p: pd.Series(v).describe().to_dict() for p, v in G.get("t5_active_counts", {}).items()},
               "n_candidates_screened": int((scores.family != "reference").sum())}
    (TRDIR / "shortlist.json").write_text(json.dumps(sl, indent=1, default=float))
    (TRDIR / "train_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    outs = [TRDIR / f for f in ("modes.tsv", "mode_scores.tsv", "per_line.tsv", "context_claims.tsv", "per_task.tsv",
                                "binary.tsv", "shortlist.json", "train_summary.json")]
    current().output(*outs)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("shortlist",)}, indent=1, default=float))
    print(json.dumps(sl, indent=1, default=float))
    return 0


def per_task_table(res) -> pd.DataFrame:
    """Per Essential task and parent: KO-failing genes (label-free), training precision of those calls (per-line,
    lineage-equal) and tAP / Δ̄_lin of 'T0 + this task' vs T0 (exploratory)."""
    from nextgem import metrics, w1

    rows = []
    Y = G["Y"]["L_prob"]
    lineage = G["lineage"]
    for p in PARENTS:
        ex = G["extra"][p]
        em = ex["emap"]
        t0 = G["modes"][(p, "T0")]
        lt0 = res[(p, "T0", "L_prob")][1]
        sym = pd.read_parquet(OUT_DIR.parent / "baseline" / f"{p}_ko.parquet", columns=["gene", "symbol"]).set_index("gene").symbol
        for t, row in ex["meta_e"].iterrows():
            genes = sorted(ex["fs_e"][t])
            fails = to_U(pd.Series(True, index=genes, dtype=bool), em, G["U"], fill=False).astype(bool) if genes else pd.Series(False, index=G["U"])
            e = t0.where(~fails, 1.0)
            S = pd.DataFrame(np.tile(e.to_numpy(), (len(Y), 1)), index=Y.index, columns=Y.columns)
            lt = metrics.per_line(S, Y)
            d = w1.delta_summary(lt, lt0, lineage)
            called = fails[fails].index
            prec = Y.loc[:, called].mean(axis=1) if len(called) else pd.Series(np.nan, index=Y.index)
            new = [g for g in called if t0[g] < CALL]
            rows.append({"parent": p, "task": t, "category": row.category, "group": row.group,
                         "description": row.description, "n_genes_failing_model": len(genes),
                         "n_failing_in_U": int(len(called)), "n_failing_not_T0_lethal": len(new),
                         "failing_symbols": ";".join(sorted({str(sym.get(g, g)) for g in genes}))[:400],
                         "precision_lin": metrics.lineage_equal_mean(prec.dropna(), lineage) if prec.notna().any() else np.nan,
                         "precision_new_calls_lin": (metrics.lineage_equal_mean(Y.loc[:, new].mean(axis=1).dropna(), lineage)
                                                     if new else np.nan),
                         "tap_lin_T0_plus_task": metrics.lineage_equal_mean(lt.tap[lt.estimable], lineage),
                         "delta_lin_vs_T0": d["delta_lin"], "delta_lo": d["delta_lo"], "delta_hi": d["delta_hi"]})
    return pd.DataFrame(rows)


# DESIGN §17 A11 (PI rulings after notes/reviews/M4_W3_review.md), applied to the shortlist:
# (a) the matched-call-rate control is read literally (extra calls only), so biomass leave-outs have no applicable
#     pre-registered negative control and qualify on TG alone; (b) flags carried verbatim; (c) withdrawals;
#     the look's gate candidates are the Human1 (G*) variants only, Human2 variants are descriptive.
A11_WITHDRAWN = {
    "B-loo:metabolite_pool_biomass": ("A11c: withdrawn before the validation look (TG-floor effect within the winner's-curse "
                                      "band, gain only in common essentials, removes genuine demands); training result only"),
}
A11_FLAGS = {
    "B-loo:cofactor_pool_biomass": {
        "W3-D1": ("passes TG; no pre-registered negative control applies to biomass leave-outs (the matched-call-rate control "
                  "is defined for added calls only; pseudo-tasks n/a). The control reading was fixed after the training "
                  "aggregates were seen. Loses to the B2-informed symmetric removal comparator; beats a post-hoc label-free "
                  "random-removal control. No claim may rest on 'beats controls'."),
        "W3-MEV": ("under serum media (RPMI M2 f=1) the variant makes HMGCR/MVK/PMVK/MVD non-lethal in both parents, a "
                   "mevalonate positive-control failure."),
    },
}

COMPLEXITY_ORDER = {"reference": 0, "B-loo": 1, "B-add": 1, "O-NGAM": 1, "D-BS": 2, "D-ER": 2, "O-cap": 3,
                    "criterion": 3, "O-MOMA": 3}


def shortlist(scores: pd.DataFrame, res) -> dict:
    """W3 decision: <= 2 candidates that pass TG and beat both negative controls; primary parent = G* (Human1).
    Among eligible candidates within one paired-difference SE of the best, the simplest (fewest items, then
    family order) is taken first; then the next best not within its SE band."""
    from nextgem import metrics

    out = {"rule": ("TG pass & beats the literal matched-call-rate control & pseudo-tasks where they apply (A11a: neither "
                    "applies to biomass leave-outs, which qualify on TG alone); A11c withdrawals; paired one-SE rule; <= 2 "
                    "(DESIGN W3)"),
           "primary_parent": "human1 (G*)", "gate_candidates": "Human1 (G*) variants only; Human2 variants are descriptive (A11c)",
           "candidates": [], "eligible": {}, "withdrawn": [], "amendment": "DESIGN §17 A11"}
    for p in PARENTS:
        el = scores[(scores.parent == p) & scores.tg_pass & scores.beats_controls & ~scores["mode"].str.endswith("@RPMI")
                    & ~scores.family.isin(["exploratory", "B-add-sens", "reference"])].copy()
        out["eligible"][p] = el["mode"].tolist()
        for m in el["mode"]:
            if m in A11_WITHDRAWN:
                r = el[el["mode"] == m].iloc[0]
                out["withdrawn"].append({"parent": p, "mode": m, "delta_lin": r.delta_lin, "ci": [r.delta_lo, r.delta_hi],
                                         "folds_positive": int(r.folds_positive), "reason": A11_WITHDRAWN[m]})
        el = el[~el["mode"].isin(list(A11_WITHDRAWN))]
        if p != "human1" or el.empty:
            continue
        el = el.sort_values("delta_lin", ascending=False)
        best = el.iloc[0]
        lt_best = res[(p, best["mode"], "L_prob")][1]
        se = {}
        for m in el["mode"]:
            d = metrics.paired_line_delta(lt_best, res[(p, m, "L_prob")][1])
            bt = metrics.bootstrap_mean(d, G["lineage"].reindex(d.index), B=metrics.B_OTHER, seed=11)
            se[m] = (bt.estimate, float(np.std(bt.replicates, ddof=1)))
        el["diff_vs_best"] = el["mode"].map(lambda m: se[m][0])
        el["se_diff"] = el["mode"].map(lambda m: se[m][1])
        el["within_1se"] = el.diff_vs_best <= el.se_diff
        el["complexity"] = el.family.map(lambda f: COMPLEXITY_ORDER.get(f, 3)) * 1000 + el.n_items.fillna(0)
        pick = el[el.within_1se].sort_values(["complexity", "delta_lin"], ascending=[True, False]).head(1)
        rest = el[~el["mode"].isin(pick["mode"])].sort_values("delta_lin", ascending=False).head(1)
        for _, r in pd.concat([pick, rest]).iterrows():
            out["candidates"].append({"parent": p, "mode": r["mode"], "delta_lin": r.delta_lin, "ci": [r.delta_lo, r.delta_hi],
                                      "folds_positive": int(r.folds_positive), "vs_mcr": r.vs_mcr_delta,
                                      "vs_pseudo": r.get("vs_pseudo_delta"), "vs_mcr_symmetric": r.get("vs_mcr_sym_delta"),
                                      "within_1se_of_best": bool(r.within_1se), "flags": A11_FLAGS.get(r["mode"], {})})
        n = len(el)
        sd = float(np.std([v for v in scores[(scores.parent == p) & (scores.family != "reference")].delta_lin.dropna()], ddof=1))
        wc = sd * float(np.sqrt(2 * np.log(max(2, int((scores.parent == p).sum())))))
        out["winners_curse_sigma_sqrt_2lnm"] = wc
        tg = scores[(scores.parent == p) & scores.tg_pass & ~scores["mode"].str.endswith("@RPMI")]
        out["winners_curse_note"] = (f"sigma*sqrt(2 ln m) = {wc:.4f} exceeds the training effect of every TG-passing Human1 "
                                     f"mode (max {tg.delta_lin.max():.4f})" if len(tg) and wc > tg.delta_lin.max() else
                                     f"sigma*sqrt(2 ln m) = {wc:.4f}")
        out["n_eligible_h1"] = n
    return out


def binary_rows(scores: pd.DataFrame, sl: dict) -> pd.DataFrame:
    from nextgem import w1

    Y06 = G["Y"]["L_ge06"]
    rows = []
    modes = {("human1", "T0"), ("human2", "T0"), ("human1", "T3"), ("human2", "T3")}
    modes |= {(c["parent"], c["mode"]) for c in sl.get("candidates", []) + sl.get("withdrawn", [])}
    for p, m in sorted(modes):
        if (p, m) not in G["modes"]:
            continue
        calls = G["modes"][(p, m)] >= CALL
        genes = [g for g in G["U"] if G["extra"][p]["in_p"][g]]
        summ, _ = w1.binary_endpoints(calls, Y06, G["lineage"], genes)
        rows.append(summ.assign(parent=p, mode=m, label="L_ge06"))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def positive_controls() -> dict:
    """Label-free W3 positive controls: ETF capability-only FPs under T3; mevalonate core under serum + B-add."""
    out = {}
    for p in PARENTS:
        em = G["extra"][p]["emap"]
        sym = pd.read_parquet(OUT_DIR.parent / "baseline" / f"{p}_ko.parquet", columns=["gene", "symbol", "entrez"])
        e_of = lambda m, s: float(G["modes"][(p, m)].get(int(sym.entrez[sym.symbol == s].iloc[0]), np.nan)) if (sym.symbol == s).any() and pd.notna(sym.entrez[sym.symbol == s].iloc[0]) else None  # noqa: E731
        out[f"{p}:ETF_T3"] = {s: {"T0": e_of("T0", s), "T3": e_of("T3", s)} for s in ("ETFA", "ETFB", "ETFDH")}
        out[f"{p}:ETF_capability_only_FP"] = all(v["T0"] is not None and v["T0"] < 0.01 and v["T3"] == 1.0
                                                for v in out[f"{p}:ETF_T3"].values())
        sw = pd.read_parquet(OUT_DIR / "sweeps" / f"ratio_r_{p}.parquet").set_index("sweep")
        mev = {}
        for s in ("HMGCR", "MVK", "PMVK", "MVD", "HMGCS1", "IDI1", "FDPS", "GGPS1"):
            g = sym.gene[sym.symbol == s]
            if g.empty:
                continue
            g = g.iloc[0]
            mev[s] = {k: float(sw.loc[k, g]) for k in sw.index if k.startswith("RPMI_M2f1|B-add") or k in ("M0|T0", "RPMI_M2f1|T0")}
        out[f"{p}:mevalonate_serum_badd"] = mev
        out[f"{p}:mevalonate_core_no_loss"] = all(
            all(v < 0.01 for v in mev[s].values()) for s in ("HMGCR", "MVK", "PMVK", "MVD") if s in mev)
        # A11b / W3-MEV: the same core under the B-loo variants that pass TG, on M0 and on the serum medium (RPMI M2 f=1)
        for v in ("B-loo:cofactor_pool_biomass", "B-loo:metabolite_pool_biomass"):
            core = {}
            for s_ in ("HMGCR", "MVK", "PMVK", "MVD"):
                g = sym.gene[sym.symbol == s_].iloc[0]
                core[s_] = {med: float(sw.loc[f"{med}|{v}", g]) for med in ("M0", "RPMI_M2f1")}
            out[f"{p}:mevalonate_core_under_{v}"] = core
            out[f"{p}:mevalonate_core_lethal_under_{v}|M0"] = all(x["M0"] < 0.01 for x in core.values())
            out[f"{p}:mevalonate_core_lethal_under_{v}|RPMI_M2f1"] = all(x["RPMI_M2f1"] < 0.01 for x in core.values())
    return out


if __name__ == "__main__":
    w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 7
    raise SystemExit(run_main(main, seeds=[11, 23, 37, 41, 53], solver="", workers=w, reads_labels=True))
