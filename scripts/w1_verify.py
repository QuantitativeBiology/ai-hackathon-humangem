"""W1 verification (DESIGN.md §12): controls and a second code path. Labels only via ``nextgem.labels``.

Checks (results/w1/verification.tsv, one row each):
A. Second code path (training tier, both parents, primary cell L_prob / all genes / U / expression lines):
   per-line tAP recomputed from scratch — KO parquet read with pyarrow, Entrez join written here, labels
   pivoted with numpy, tAP from an independent implementation (expected precision of each positive over a
   uniformly random within-block order, summed position by position) — compared with
   results/w1/train_per_line.tsv and the grid summaries (lineage-equal mean, lift, pooled tAP); a Monte-Carlo
   check of tAP (random tie-breaking + sklearn AP, 2,000 orders) on 5 lines; and per-line MCC for
   H2 P57 x L-ge06 via sklearn.matthews_corrcoef vs results/w1/train_binary.tsv.
B. Positive controls: oracle scorer tAP = 1 on every estimable line; oracle caller MCC = 1 (from w1_binary);
   CAD/UMPS lethal in H2, CAD lethal in H1 (ratio_r = 0).
C. Negative controls: constant scorer B0 tAP = E[AP_random] exactly; B1 random lift ~ 1 (CI covers 1 for
   >= 4 of 5 seeds is not required: reported); labels shuffled across genes within each line (SEEDS[0]) give
   H2 lift ~ 1; constant caller MCC NaN (primary) / 0 (sensitivity); a fold-scoped memorising scorer gets
   exactly E[AP_random] on every held-out fold (nextgem.cv.lolo); LAP3 (H2) / TGDS (H1) ratio_r = 1.
D. §12 §6.8 row: P57 calls equal Human-GEM ``taskEssentialGenes.find_task_essential_categories`` on both
   generic parents (``--task-equality``; Gurobi, 7 workers + parent = 8 sessions), and
   ``metabolicTasks_Essential.txt`` vs Robinson 2020 data file S3 (sheet HumanGEM).

Usage: .venv/bin/python scripts/w1_verify.py [--task-equality]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w1_common import OUT_DIR, PROJECT, atomic_write_tsv, current, run_main  # noqa: E402

CHECKS: list[dict] = []
TOL = 1e-6  # result files are written with 6 significant digits (%.6g), so file-vs-recompute agreement is ~5e-7


def check(name, passed, detail):
    CHECKS.append({"check": name, "passed": bool(passed), "detail": str(detail)})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def tap_independent(score: np.ndarray, y: np.ndarray) -> float:
    """Expected AP over uniform random orders within tied score blocks, written from the definition:
    a positive in a block of size n with p positives, preceded by N items (P positive), is at in-block
    position r (1..n) with probability 1/n; given r, the expected number of other positives before it in
    the block is (r-1)(p-1)/(n-1); E[precision] = (P + 1 + that) / (N + r) (exact, since the denominator
    is fixed given r)."""
    order = np.argsort(-score, kind="stable")
    s, t = score[order], y[order]
    R = t.sum()
    if R == 0:
        return float("nan")
    total, N, P = 0.0, 0, 0
    i = 0
    while i < len(s):
        j = i
        while j < len(s) and s[j] == s[i]:
            j += 1
        n, p = j - i, int(t[i:j].sum())
        if p:
            r = np.arange(1, n + 1, dtype=float)
            extra = (r - 1) * (p - 1) / (n - 1) if n > 1 else np.zeros(1)
            total += p * np.mean((P + 1 + extra) / (N + r))
        N += n
        P += p
        i = j
    return total / R


def section_a(train_long, E, lineage, genes):
    from sklearn.metrics import average_precision_score, matthews_corrcoef

    pl = pd.read_csv(OUT_DIR / "train_per_line.tsv", sep="\t")
    grid = pd.read_csv(OUT_DIR / "train_grid.tsv", sep="\t")
    prim = grid[(grid.label == "L_prob") & (grid.geneset == "all") & (grid.universe == "U") & (grid.stratum == "none")
                & (grid.lineset == "expr")].set_index("method")
    # labels: own numpy pivot
    tl = train_long[train_long.ModelID.isin(E)]
    li = {l: i for i, l in enumerate(E)}
    gi = {g: i for i, g in enumerate(genes)}
    Y = np.full((len(E), len(genes)), np.nan)
    Y[tl.ModelID.map(li).to_numpy(), tl.entrez.map(gi).to_numpy()] = tl.L_prob.astype(float).to_numpy()
    for parent, short in (("human1", "H1"), ("human2", "H2")):
        t = pq.read_table(PROJECT / "results" / "baseline" / f"{parent}_ko.parquet", columns=["entrez", "ratio_r"]).to_pandas()
        t = t[t.entrez.notna()]
        e = np.zeros(len(genes))
        for ez, r in zip(t.entrez.astype(int), t.ratio_r.astype(float)):
            if ez in gi:
                e[gi[ez]] = round(1.0 - r, 6)
        mine = []
        for l in E:
            y = Y[li[l]]
            m = ~np.isnan(y)
            npos = int(y[m].sum())
            mine.append(tap_independent(e[m], y[m]) if npos >= 5 else np.nan)
        mine = pd.Series(mine, index=E)
        ref = pl[(pl.method == short) & (pl.geneset == "all")].set_index("ModelID").tap.reindex(E)
        diff = (mine - ref).abs().max()
        same_na = bool((mine.isna() == ref.isna()).all())
        check(f"A1_second_path_tap_line_{short}", diff < TOL and same_na,
              f"{int(mine.notna().sum())} estimable lines; max |diff| = {diff:.2e}; NA pattern identical = {same_na}")
        lin = lineage.reindex(E)
        lem = mine.groupby(lin).mean().mean()
        check(f"A2_second_path_tap_lin_{short}", abs(lem - prim.loc[short, "tap_lin"]) < TOL,
              f"lineage-equal mean tAP {lem:.6f} vs grid {prim.loc[short, 'tap_lin']:.6f}")
        # pooled
        m = ~np.isnan(Y)
        pooled = tap_independent(np.tile(e, (len(E), 1))[m], Y[m])
        check(f"A3_second_path_pooled_{short}", abs(pooled - prim.loc[short, "pooled_tap"]) < TOL,
              f"pooled tAP {pooled:.6f} vs grid {prim.loc[short, 'pooled_tap']:.6f}")
        # Monte-Carlo: random tie-breaking + sklearn AP
        rng = np.random.default_rng(11)
        for l in list(mine.dropna().index[:5]):
            y = Y[li[l]]
            mk = ~np.isnan(y)
            vals = [average_precision_score(y[mk], e[mk] + 1e-9 * rng.random(mk.sum())) for _ in range(2000)]
            se = np.std(vals) / np.sqrt(len(vals))
            check(f"A4_monte_carlo_{short}_{l}", abs(np.mean(vals) - mine[l]) < 4 * se + 1e-12,
                  f"MC mean {np.mean(vals):.5f} (SE {se:.1e}) vs tAP {mine[l]:.5f}")
    # binary second path: H2 P57 x L-ge06, parent genes, expression lines
    p57 = pq.read_table(OUT_DIR / "p57_human2.parquet", columns=["entrez", "p57"]).to_pandas().dropna(subset=["entrez"])
    call = {int(a): bool(b) for a, b in zip(p57.entrez, p57.p57)}
    pg = [g for g in genes if g in call]
    Y6 = np.full((len(E), len(genes)), np.nan)
    Y6[tl.ModelID.map(li).to_numpy(), tl.entrez.map(gi).to_numpy()] = tl.L_ge06.astype(float).to_numpy()
    c = np.array([call[g] for g in pg])
    cols = [gi[g] for g in pg]
    mcc = []
    for l in E:
        y = Y6[li[l], cols]
        m = ~np.isnan(y)
        yt, cp = y[m].astype(int), c[m].astype(int)
        den_ok = len(set(yt)) > 1 and len(set(cp)) > 1
        mcc.append(matthews_corrcoef(yt, cp) if den_ok else np.nan)
    mcc = pd.Series(mcc, index=E)
    lem = mcc.groupby(lineage.reindex(E)).mean().mean()
    tb = pd.read_csv(OUT_DIR / "train_binary.tsv", sep="\t")
    r = tb[(tb.protocol == "Human2-paper") & (tb.method == "H2") & (tb.call_rule == "P57") & (tb.gene_set == "parent")
           & (tb.lineset == "expr") & (tb.metric == "mcc") & (tb.nan_policy == "exclude")].iloc[0]
    check("A5_second_path_mcc_H2_P57", abs(lem - r.lineage_equal_mean) < TOL,
          f"sklearn lineage-equal MCC {lem:.6f} vs train_binary {r.lineage_equal_mean:.6f}")


def section_bc(train_long, E, lineage, genes):
    from nextgem import cv, metrics, w1
    from nextgem.config import SEEDS

    grid = pd.read_csv(OUT_DIR / "train_grid.tsv", sep="\t")
    prim = grid[(grid.label == "L_prob") & (grid.geneset == "all") & (grid.universe == "U") & (grid.stratum == "none")
                & (grid.lineset == "expr")].set_index("method")
    Y = w1.wide(train_long, "L_prob", E, genes)
    # B: oracle
    lt = metrics.per_line(Y.fillna(0.0), Y)
    est = lt[lt.estimable]
    check("B1_oracle_tap_1", bool(np.allclose(est.tap, 1.0)), f"{len(est)} estimable lines, min tAP {est.tap.min():.6f}")
    tb = pd.read_csv(OUT_DIR / "train_binary.tsv", sep="\t")
    o = tb[(tb.method == "oracle") & (tb.metric == "mcc") & (tb.nan_policy == "exclude")].lineage_equal_mean.iloc[0]
    check("B2_oracle_mcc_1", abs(o - 1) < 1e-12, f"oracle MCC {o}")
    ko2 = pd.read_parquet(PROJECT / "results" / "baseline" / "human2_ko.parquet").set_index("symbol").ratio_r
    ko1 = pd.read_parquet(PROJECT / "results" / "baseline" / "human1_ko.parquet").set_index("symbol").ratio_r
    check("B3_CAD_UMPS_lethal_H2", ko2["CAD"] == 0 and ko2["UMPS"] == 0, f"H2 CAD {ko2['CAD']}, UMPS {ko2['UMPS']}")
    check("B4_CAD_lethal_H1_UMPS_known_issue", ko1["CAD"] == 0, f"H1 CAD {ko1['CAD']}, UMPS {ko1['UMPS']} (known: reversible TYMS)")
    # C: negative controls
    b0 = prim.loc["B0"]
    check("C1_constant_equals_null", abs(b0.tap_lin - b0.e_random_lin) < 1e-12,
          f"B0 tAP_lin {b0.tap_lin:.6f} vs E[AP_random] {b0.e_random_lin:.6f}")
    for sd in SEEDS[0:5]:
        r = prim.loc[f"B1_s{sd}"]
        check(f"C2_random_lift_{sd}", abs(r.lift_lin - 1) < 0.05,
              f"lift {r.lift_lin:.4f}; tAP {r.tap_lin:.5f} [{r.tap_lin_lo:.5f}, {r.tap_lin_hi:.5f}] vs null {r.e_random_lin:.5f}")
    gen = w1.generic_table(pd.Series(genes))
    S = pd.DataFrame(np.tile(gen.e_H2.reindex(genes).to_numpy(), (len(E), 1)), index=E, columns=genes)
    rng = np.random.default_rng(SEEDS[0])
    Ys = Y.copy()
    arr = Ys.to_numpy()
    for i in range(arr.shape[0]):
        arr[i] = arr[i][rng.permutation(arr.shape[1])]
    Ys = pd.DataFrame(arr, index=Y.index, columns=Y.columns)
    lt = metrics.per_line(S, Ys)
    e = lt[lt.estimable]
    lift = metrics.lineage_equal_mean(e.tap, lineage) / metrics.lineage_equal_mean(e.e_random, lineage)
    check("C3_gene_shuffled_labels_H2_lift_1", abs(lift - 1) < 0.05, f"H2 lift with labels shuffled across genes: {lift:.4f} "
          f"(real {prim.loc['H2', 'lift_lin']:.4f})")
    c = tb[(tb.method == "constant_nonessential") & (tb.metric == "mcc")]
    ex = c[c.nan_policy == "exclude"].iloc[0]
    zr = c[c.nan_policy == "zero"].iloc[0]
    check("C4_constant_caller_mcc", np.isnan(ex.lineage_equal_mean) and zr.lineage_equal_mean == 0,
          f"primary {ex.lineage_equal_mean} (n_nan {ex.n_nan}); sensitivity {zr.lineage_equal_mean}")

    # fold-scoped memoriser: remembers every fitting (line, gene) label; unseen lines get 0 for every gene
    def fit_predict(fit_long, ids):
        memo = {(m, g): v for m, g, v in zip(fit_long.ModelID, fit_long.entrez, fit_long.L_prob) if pd.notna(v)}
        return pd.DataFrame([[memo.get((l, g), 0.0) for g in genes] for l in ids], index=ids, columns=genes).astype(float)

    def score(pred, held):
        Yh = w1.wide(held, "L_prob", list(pred.index), genes)
        lt = metrics.per_line(pred, Yh)
        e = lt[lt.estimable]
        return float((e.tap - e.e_random).abs().max()) if len(e) else 0.0

    res = cv.lolo(fit_predict, score, columns=["L_prob"])
    worst = max(res.values())
    check("C5_fold_scoped_memoriser_null", worst < 1e-12, f"max |tAP - E[AP_random]| over 12 held-out folds = {worst:.1e}")
    check("C6_nonlethal_controls", ko2["LAP3"] == 1 and ko1["TGDS"] == 1, f"H2 LAP3 {ko2['LAP3']}, H1 TGDS {ko1['TGDS']}")
    # orientation audit (label-free direction check on training): lethal genes are more often positive
    frac = Y.mean()
    leth = gen.fba_call_H2.reindex(genes)
    a, b = frac[leth.to_numpy()].mean(), frac[~leth.to_numpy()].mean()
    check("C7_orientation_lethal_more_positive", a > b, f"mean L_prob prevalence: H2-lethal {a:.3f} vs non-lethal {b:.3f}")


def section_d(task_equality: bool):
    import cobra

    from nextgem.model import GEMS
    from raven_toolbox.tasks.tasklist import parse_task_list

    # task file vs Robinson data file S3
    if task_equality:
        return section_d2()
    s3 = pd.read_excel(PROJECT / "papers" / "aaz1482_data_file_s3.xlsx", sheet_name="HumanGEM")
    s3_tasks, cur = [], None
    for _, r in s3.iterrows():
        if str(r["Unnamed: 0"]) == "#":
            cur = None
            continue
        if pd.notna(r["DESCRIPTION"]):
            cur = {"id": str(r["ID"]), "desc": str(r["DESCRIPTION"]).strip(), "in": set(), "out": set()}
            s3_tasks.append(cur)
        if cur is not None:
            for k, col in (("in", "IN"), ("out", "OUT")):
                if pd.notna(r[col]):
                    cur[k] |= {x.strip().replace("[s]", "[e]") for x in str(r[col]).split(";")}  # S3 extracellular = [s]
    for name in ("human1", "human2"):
        tasks = parse_task_list(GEMS[name].essential_tasks)
        n_same_id = sum(t.id == s["id"] and t.description.strip() == s["desc"] for t, s in zip(tasks, s3_tasks))
        n_same_io = sum({i[0] for i in t.inputs} == s["in"] and {o[0] for o in t.outputs} == s["out"]
                        for t, s in zip(tasks, s3_tasks))
        check(f"D1_tasks_vs_Robinson_S3_{name}", len(tasks) == len(s3_tasks) == 57 and n_same_id == 57,
              f"{len(tasks)} tasks vs S3 {len(s3_tasks)}; identical (id, description) {n_same_id}/57; "
              f"identical IN and OUT metabolite sets {n_same_io}/57 (S3 compartment [s] read as [e])")



def section_d2():
    import cobra

    from nextgem.model import GEMS

    sys.path.insert(0, str(PROJECT / "data" / "Human-GEM-main" / "code" / "test"))
    from taskEssentialGenes import find_task_essential_categories  # Human-GEM code (reuse, DESIGN §2)

    from nextgem.model import load_model

    for name in ("human1", "human2"):
        t0 = time.time()
        cobra.Configuration().solver = "gurobi"
        m = load_model(name)
        ref = find_task_essential_categories(m, GEMS[name].essential_tasks, processes=7)
        del m
        mine = pd.read_parquet(OUT_DIR / f"p57_{name}.parquet").set_index("gene")
        a = set(mine.index[mine.p57])
        b = set(ref)
        cats_mine = {g: set(c.split(",")) for g, c in mine.broken_categories.items() if c}
        cat_eq = sum(cats_mine.get(g, set()) == ref[g] for g in b & a)
        check(f"D2_p57_equals_taskEssentialGenes_{name}", a == b,
              f"ours {len(a)} vs Human-GEM {len(b)}; only ours {sorted(a - b)[:5]}; only theirs {sorted(b - a)[:5]}; "
              f"category sets equal for {cat_eq}/{len(a & b)} (theirs stop once all 5 categories are hit); {time.time() - t0:.0f} s")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-equality", action="store_true")
    a = ap.parse_args(argv)
    t0, c0 = time.time(), time.process_time()
    from nextgem import labels

    split = labels.load_split().set_index("ModelID")
    lf = pd.read_csv(OUT_DIR / "line_features.tsv", sep="\t").set_index("ModelID")
    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    E = [l for l in sorted(elig.index) if lf.loc[l, "has_expression_matrix"]]
    if not a.task_equality:
        tl = labels.train(columns=["L_prob", "L_ge06"])
        section_a(tl, E, split.lineage, genes)
        section_bc(tl, E, split.lineage, genes)
    section_d(a.task_equality)
    out = OUT_DIR / ("verification_tasks.tsv" if a.task_equality else "verification.tsv")
    atomic_write_tsv(pd.DataFrame(CHECKS), out, float_format=None)
    current().output(out)
    if not a.task_equality:
        current().input(*[OUT_DIR / f for f in ("train_per_line.tsv", "train_grid.tsv", "train_binary.tsv",
                                                 "p57_human2.parquet")], PROJECT / "results" / "baseline" / "human1_ko.parquet",
                        PROJECT / "results" / "baseline" / "human2_ko.parquet", PROJECT / "papers" / "aaz1482_data_file_s3.xlsx")
    else:
        current().input(OUT_DIR / "p57_human1.parquet", OUT_DIR / "p57_human2.parquet")
    n_fail = sum(not c["passed"] for c in CHECKS)
    print(f"{len(CHECKS)} checks, {n_fail} failed, {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(run_main(main, seeds=[11], solver="gurobi" if "--task-equality" in sys.argv else "",
                              workers=7 if "--task-equality" in sys.argv else 1, reads_labels="--task-equality" not in sys.argv))
