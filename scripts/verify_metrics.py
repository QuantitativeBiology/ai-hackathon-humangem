"""Synthetic verification of src/nextgem/metrics.py beyond the unit tests (DESIGN.md v3 §6.1, §12 'tAP and metrics').

Synthetic data only; no label file is opened. Runs on <= 8 cores. Writes results/m1_metrics/verify_metrics.json.

1. Monte-Carlo check (§6.1): mean AP of score + 1e-9 noise over 10^5 jitters vs tAP, on a generic-FBA-shaped line.
2. Bootstrap calibration (§12 negative control; stats r1 SF5): 500 synthetic Δ = 0 datasets per regime with a
   lineage random effect on Δ, near-binary score blocks and training-like per-line prevalence. Reports the coverage
   of the lineage-stratified bootstrap CI for the conditional estimand (mean of the 7 lineages' true Δ) and of the
   lineage t-interval for the population estimand (0), plus sign-flip rejection rates.
3. Positive control (§12): a known gap δ = 0.01 added to every line is recovered.
4. Runtime of the production-size calls.

Usage: .venv/bin/python scripts/verify_metrics.py [--sims 500] [--workers 8]
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from nextgem import metrics as M  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "m1_metrics" / "verify_metrics.json"
SEEDS0 = 11  # SEEDS[0], DESIGN §7
# Illustrative 7-lineage test tier of ~350 lines, sizes 15-126 (notes/design_review_r2/statistics.md, R2-10).
LINEAGE_SIZES = [126, 75, 50, 35, 28, 20, 15]
# Training-like per-line prevalence (notes/design_review_r1/statistics.md: median 320 positives per line at
# L-ge05 on all genes, 37 without common essentials). Gene counts: |U| = 2,866; U minus CE taken as 2,450.
REGIMES = {
    "all_genes": dict(G=2866, median_pos=320),
    "no_common_essentials": dict(G=2450, median_pos=37),
}
TAUS = [0.0, 0.01]  # lineage random-effect SD on Δ (0.01 = SESOI)
GAP = 0.01  # positive control: known gap


def lineage_index():
    idx = [f"S{i:04d}" for i in range(sum(LINEAGE_SIZES))]
    return pd.Series(np.repeat([f"K{k}" for k in range(len(LINEAGE_SIZES))], LINEAGE_SIZES), index=idx)


def simulate_line_delta(rng, G, median_pos, n_top=200, n_swap=60, beta=3.0):
    """Per-line ΔtAP of two exchangeable near-binary scorers (E[Δ | lineage] = 0 exactly)."""
    L = sum(LINEAGE_SIZES)
    logw = rng.normal(0.0, 1.5, G)  # gene-level dependency propensity (common vs rare)
    npos = np.clip(np.rint(median_pos * np.exp(rng.normal(0, 0.3, L))), 5, G // 2).astype(int)
    keys = logw + rng.gumbel(size=(L, G))
    rank = np.argsort(-keys, axis=1)
    Y = np.zeros((L, G))
    rows = np.arange(L)[:, None]
    Y[rows, rank] = (np.arange(G)[None, :] < npos[:, None]).astype(float)
    # generic-model-like base scorer: 150 lethal (e = 1), 50 partial, rest tied at 0
    e0 = np.argsort(-(0.7 * logw + rng.gumbel(size=G)))[:n_top]
    rest = np.setdiff1d(np.arange(G), e0)
    base = np.zeros(G)
    base[e0[:150]] = 1.0
    base[e0[150:]] = rng.choice([0.25, 0.5, 0.75], n_top - 150)

    def method():
        S = np.tile(base, (L, 1))
        drop = np.argpartition(rng.random((L, n_top)), n_swap, axis=1)[:, :n_swap]
        S[rows, e0[drop]] = 0.0
        k_in = beta * Y[:, rest] + rng.gumbel(size=(L, rest.size))
        add = np.argpartition(-k_in, n_swap, axis=1)[:, :n_swap]
        S[rows, rest[add]] = 1.0
        return S

    SA, SB = method(), method()
    ta, _, pa = M.tap_columns(SA.T, Y.T)
    tb, _, _ = M.tap_columns(SB.T, Y.T)
    tie_share = float(np.mean((SA == 0).sum(1) / G))
    return ta - tb, ta, npos, tie_share


def one_sim(args):
    regime, sim_ss = args
    p = REGIMES[regime]
    ss_data, ss_u, ss_rest = sim_ss.spawn(3)
    delta_raw, ta, npos, tie_share = simulate_line_delta(np.random.default_rng(ss_data), p["G"], p["median_pos"])
    lin = lineage_index()
    K = len(LINEAGE_SIZES)
    rng_u = np.random.default_rng(ss_u)
    z = rng_u.normal(0.0, 1.0, K)
    s_within = float(pd.Series(delta_raw, index=lin.index).groupby(lin).std().mean())
    out = []
    scen = [("null", tau, 0.0) for tau in TAUS] + [("gap", TAUS[-1], GAP)]
    for (name, tau, gap), ss in zip(scen, ss_rest.spawn(len(scen))):
        u = tau * z
        d = pd.Series(delta_raw + u[np.repeat(np.arange(K), LINEAGE_SIZES)] + gap, index=lin.index)
        theta_cond = float(u.mean() + gap)
        theta_pop = gap
        ss_b, ss_f = ss.spawn(2)
        boot = M.bootstrap_mean(d, lin, B=M.B_CONFIRMATORY, seed=ss_b)
        _, p_flip = M.sign_flip_test(d, lin, n_flip=M.N_SIGN_FLIPS, seed=ss_f)
        tr = M.transfer_summary(d, lin)
        out.append(
            dict(
                regime=regime, scenario=name, tau=tau, gap=gap, estimate=boot.estimate,
                theta_cond=theta_cond, theta_pop=theta_pop,
                cover_cond=boot.ci_low <= theta_cond <= boot.ci_high,
                cover_pop_boot=boot.ci_low <= theta_pop <= boot.ci_high,
                cover_t_pop=tr["t_low"] <= theta_pop <= tr["t_high"],
                cover_t_cond=tr["t_low"] <= theta_cond <= tr["t_high"],
                boot_width=boot.ci_high - boot.ci_low, t_width=tr["t_high"] - tr["t_low"],
                reject_flip=p_flip <= M.ALPHA, t_low_pos=tr["t_low"] > 0, generalises=tr["generalises"],
                s_within=s_within, mean_tap=float(np.mean(ta)), median_pos=float(np.median(npos)),
                tie_share=tie_share,
            )
        )
    return out


def mc_check(n_jitter=100_000, seed=SEEDS0):
    """§6.1: mean AP over 10^5 jitters (score + 1e-9 noise) vs tAP on a generic-FBA-shaped line."""
    rng = np.random.default_rng(seed)
    G = 2866
    s = np.zeros(G)
    idx = rng.permutation(G)
    s[idx[:150]] = 1.0
    s[idx[150:190]] = rng.choice([0.25, 0.5, 0.75], 40)
    y = np.zeros(G)
    y[idx[:75]] = 1
    y[rng.choice(idx[150:], 245, replace=False)] = 1
    target = M.tap(s, y)
    aps = np.empty(n_jitter)
    bs = 500
    for b0 in range(0, n_jitter, bs):
        noisy = s[None, :] + 1e-9 * rng.random((bs, G))
        order = np.argsort(-noisy, axis=1)
        ys = y[order]
        prec = np.cumsum(ys, 1) / np.arange(1, G + 1)
        aps[b0 : b0 + bs] = (prec * ys).sum(1) / y.sum()
    se = aps.std(ddof=1) / np.sqrt(n_jitter)
    return dict(
        n_jitter=n_jitter, tap=target, mc_mean=float(aps.mean()), mc_se=float(se),
        z=float((aps.mean() - target) / se), sklearn_step_ap=M.ap_step_sklearn(s, y),
        e_random=float(M.expected_ap_random(G, int(y.sum()))),
    )


def runtime_benchmarks():
    rng = np.random.default_rng(SEEDS0)
    lin = lineage_index()
    L, G = len(lin), 2866
    Y = pd.DataFrame((rng.random((L, G)) < 0.11).astype(float), index=lin.index, columns=[f"g{j}" for j in range(G)])
    SA = pd.DataFrame(np.round(rng.random((L, G)), 2), index=Y.index, columns=Y.columns)
    SB = pd.DataFrame(rng.integers(0, 3, (L, G)) / 2.0, index=Y.index, columns=Y.columns)
    res = {"shape_lines_x_genes": [L, G]}
    t0 = time.perf_counter(); M.per_line(SA, Y); res["per_line_s"] = time.perf_counter() - t0
    t0 = time.perf_counter(); g = M.per_gene(SA, Y); res["per_gene_s"] = time.perf_counter() - t0
    res["n_genes_eligible"] = int(g["estimable"].sum())
    d = M.paired_line_delta(M.per_line(SA, Y), M.per_line(SB, Y))
    t0 = time.perf_counter(); M.confirmatory_contrast(d, lin, contrast="C1", sesoi=M.SESOI, seed=SEEDS0); res["confirmatory_contrast_B1e4_flip1e4_s"] = time.perf_counter() - t0
    B = 200
    t0 = time.perf_counter(); M.gene_bootstrap(SA, SB, Y, lin, B=B, seed=SEEDS0)
    dt = time.perf_counter() - t0
    res["gene_bootstrap_B200_s"] = dt
    res["gene_bootstrap_B1e4_extrapolated_s"] = dt * 10_000 / B
    n_perm = 50
    t0 = time.perf_counter(); M.c3_permutation_null(SA, Y, lin, n_perm=n_perm, seed=SEEDS0)
    dt = time.perf_counter() - t0
    res["c3_null_50perm_s"] = dt
    res["c3_null_2000perm_extrapolated_s"] = dt * 2000 / n_perm
    return res


def summarise(rows):
    df = pd.DataFrame(rows)
    out = []
    for (regime, scen, tau), g in df.groupby(["regime", "scenario", "tau"], sort=False):
        n = len(g)

        def cov(col):
            c = float(g[col].mean())
            return {"value": c, "mc_se": float(np.sqrt(c * (1 - c) / n))}

        out.append(
            dict(
                regime=regime, scenario=scen, tau=tau, gap=float(g["gap"].iloc[0]), n_sims=n,
                coverage_boot_conditional=cov("cover_cond"),
                coverage_boot_population=cov("cover_pop_boot"),
                coverage_t_population=cov("cover_t_pop"),
                coverage_t_conditional=cov("cover_t_cond"),
                reject_rate_sign_flip=float(g["reject_flip"].mean()),
                rate_t_low_positive=float(g["t_low_pos"].mean()),
                rate_generalises=float(g["generalises"].mean()),
                mean_estimate=float(g["estimate"].mean()),
                sd_estimate=float(g["estimate"].std()),
                mean_boot_ci_width=float(g["boot_width"].mean()),
                mean_t_ci_width=float(g["t_width"].mean()),
                mean_within_lineage_sd_delta=float(g["s_within"].mean()),
                mean_line_tap=float(g["mean_tap"].mean()),
                median_positives_per_line=float(g["median_pos"].median()),
                mean_share_genes_tied_at_zero=float(g["tie_share"].mean()),
            )
        )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", type=int, default=500)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    t_start = time.perf_counter()
    report = {"seed": SEEDS0, "lineage_sizes": LINEAGE_SIZES, "regimes": REGIMES, "taus": TAUS, "gap": GAP,
              "B": M.B_CONFIRMATORY, "n_flip": M.N_SIGN_FLIPS}

    t0 = time.perf_counter()
    report["mc_check"] = mc_check()
    report["mc_check"]["seconds"] = time.perf_counter() - t0
    print("MC check:", json.dumps(report["mc_check"]), flush=True)

    t0 = time.perf_counter()
    master = np.random.SeedSequence(SEEDS0)
    tasks = []
    for regime, rss in zip(REGIMES, master.spawn(len(REGIMES))):
        tasks += [(regime, ss) for ss in rss.spawn(args.sims)]
    with Pool(min(args.workers, 8)) as pool:
        rows = [r for chunk in pool.imap(one_sim, tasks, chunksize=4) for r in chunk]
    report["calibration"] = summarise(rows)
    report["calibration_seconds"] = time.perf_counter() - t0
    for r in report["calibration"]:
        print(json.dumps(r), flush=True)

    t0 = time.perf_counter()
    report["runtime"] = runtime_benchmarks()
    report["runtime"]["seconds"] = time.perf_counter() - t0
    print("Runtime:", json.dumps(report["runtime"]), flush=True)
    report["total_seconds"] = time.perf_counter() - t_start
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix(".tmp")
    tmp.write_text(json.dumps(report, indent=2, default=float) + "\n")
    os.replace(tmp, OUT)
    print(f"wrote {OUT} in {report['total_seconds']:.0f} s")


if __name__ == "__main__":
    main()
