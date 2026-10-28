"""W2 M4 scoring helpers (DESIGN §9 W2 "Metrics", "Controls", "Decision"; §6; §9 TG). Pure functions; no label
file is opened here. Label frames come from ``labels.train()`` / ``labels.fold(k)`` in the calling script
(``scripts/w2_score_train.py``).

Score convention (as W1/W3): e = 1 − ratio_r of the line's sweep for a ladder step, mapped to the gene universe U
(Entrez); a gene absent from a parent scores 0; rounded to 6 decimals. Binary call (§6.1): ratio_r < 0.01, i.e.
e >= 0.99 (``CALL``).

Ladder steps (``results/w2/controls/line_sweeps.tsv``): every line has one sweep per step and parent. Steps whose
name contains ``g(`` or equals ``M1g`` carry the transporter gate (+g), a *context* lever (PI decision r2-8): they
never enter "medium alone" or the W2 shortlist. Pre-registered eligibility for the W2 shortlist (fixed before any
W2 score was computed): the documented gate-free variants M1 and M2(f) for f in {0.001, 0.01, 0.1, 1}. M3 is
uncalibrated and excluded (A9i; no M3 sweep exists). M-open (negative control), M0 (reference) and the
physiological counterfactuals M-HPLM / M-Plasmax (applied to all lines; not the documented medium) are scored and
reported but are not eligible.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

CALL = 0.99  # e >= 0.99 <=> ratio_r < 0.01 (§6.1 fixed cutoff)
TG_MIN_DELTA = 0.005  # §9 TG (1)
TG_MIN_FOLDS = 9  # §9 TG (3): >= 9 of the 12 training LOLO folds
FOLD_MIN_LINES = 5  # §9 TG: a training fold with < 5 estimable lines counts as "not positive"
F_LEVELS = (0.001, 0.01, 0.1, 1.0)
SHUFFLE_MIN_LINES = 8  # DESIGN W2 negative control: power only with >= 8 attribution lines in each of two media
ELIGIBLE = ("M1",) + tuple(f"M2(f={f:g})" for f in F_LEVELS)  # documented, gate-free (M3 excluded, A9i)
DOC_LADDER = ("M0",) + ELIGIBLE  # oracle ladder (gate-free documented, incl. the reference)
STRATA = {
    "attribution": ("attribution",),
    "onboarded": ("onboarded_sensitivity",),
    "other": ("imputed", "fallback_sheet_missing", "unmappable", "mixture_on_disk"),
}


@dataclass(frozen=True)
class Step:
    name: str
    family: str  # negative_control | reference | documented | physiological
    gated: bool
    f: float | None
    eligible: bool

    @property
    def gate_free_twin(self) -> str:
        """The same step without the gate (M1g -> M1, M2g(f=x) -> M2(f=x), M-HPLM:M2g(...) -> M-HPLM:M2(...))."""
        return self.name.replace("M1g", "M1").replace("M2g(", "M2(")


def parse_step(name: str) -> Step:
    gated = name == "M1g" or "M2g(" in name
    m = re.search(r"f=([0-9.]+)", name)
    f = float(m.group(1)) if m else None
    if name == "M-open":
        fam = "negative_control"
    elif name == "M0":
        fam = "reference"
    elif name.startswith("M-HPLM") or name.startswith("M-Plasmax"):
        fam = "physiological"
    else:
        fam = "documented"
    return Step(name, fam, gated, f, (name in ELIGIBLE))


def step_order(names) -> list[str]:
    """Display order: M-open, M0, M1, M2(f) ascending, gated twins, then physiological."""
    def key(n):
        s = parse_step(n)
        fam = {"negative_control": 0, "reference": 1, "documented": 2, "physiological": 3}[s.family]
        phys = 0 if not n.startswith("M-") or n == "M-open" else (1 if "HPLM" in n else 2)
        return (fam, phys, s.gated, 0 if "M1" in n else 1, s.f or 0.0)

    return sorted(set(names), key=key)


# ----------------------------------------------------------------------------------------------- scores
def gene_to_entrez(baseline_ko: pd.DataFrame) -> pd.Series:
    """Model gene -> Entrez exactly as ``w1.generic_table`` (baseline parquet order, NaN dropped, first per Entrez)."""
    ko = baseline_ko[["gene", "entrez"]].dropna(subset=["entrez"])
    ko = ko.assign(entrez=ko.entrez.astype(int)).drop_duplicates("entrez")
    return ko.set_index("gene").entrez


def sweep_e_matrix(ratio: pd.DataFrame, emap: pd.Series, U) -> pd.DataFrame:
    """sweeps x U matrix of e = 1 − ratio_r (absent gene -> 0), rounded to 6 decimals."""
    r = ratio.reindex(columns=emap.index)
    e = (1.0 - r).fillna(0.0)
    e.columns = emap.to_numpy()
    return e.reindex(columns=list(U), fill_value=0.0).round(6)


def line_matrix(sweep_of_line: pd.Series, E_sweeps: pd.DataFrame, lines) -> pd.DataFrame:
    """lines x U score matrix: each line takes the row of its sweep (``sweep_of_line``: ModelID -> sweep_id)."""
    sid = sweep_of_line.reindex(list(lines))
    if sid.isna().any():
        raise KeyError(f"{int(sid.isna().sum())} lines have no sweep")
    M = E_sweeps.loc[sid.to_numpy()]
    M.index = list(lines)
    return M


def shuffled_sweeps(shuffled_key: pd.Series, step: str, parent: str, labels_tab: pd.DataFrame) -> pd.Series:
    """ModelID -> sweep_id of ``step`` under a shuffled medium key (``sweep_labels.tsv`` lookup, gate open)."""
    s = parse_step(step)
    lt = labels_tab[(labels_tab.parent == parent) & (labels_tab.variant == ("M1" if s.name == "M1" else "M2"))
                    & (labels_tab.gate_open.astype(str) == "True")]
    if s.f is not None:
        lt = lt[np.isclose(lt.f.astype(float), s.f)]
    lookup = lt.drop_duplicates("medium_key").set_index("medium_key").sweep_id
    out = shuffled_key.map(lookup)
    if out.isna().any():
        raise KeyError(f"no {step} sweep for medium keys {sorted(set(shuffled_key[out.isna()]))}")
    return out


# ----------------------------------------------------------------------------------------------- gates
def folds_positive(d: pd.Series, lineage: pd.Series, folds, min_lines: int = FOLD_MIN_LINES) -> dict:
    """§9 TG (3) on a per-line Δ series: a fold is positive iff it has >= ``min_lines`` estimable lines and a mean
    Δ > 0; a fold with fewer lines is 'not positive' (never dropped; the denominator stays len(folds))."""
    lin = lineage.reindex(d.index)
    g = d.groupby(lin).agg(["mean", "size"])
    pos, means, sizes = 0, {}, {}
    for f in folds:
        if f in g.index:
            means[f], sizes[f] = float(g.loc[f, "mean"]), int(g.loc[f, "size"])
            if sizes[f] >= min_lines and means[f] > 0:
                pos += 1
        else:
            means[f], sizes[f] = None, 0
    return {"folds_positive": pos, "n_folds": len(folds),
            "fold_means": json.dumps({k: (round(v, 6) if v is not None else None) for k, v in means.items()}),
            "fold_sizes": json.dumps(sizes)}


def tg_pass(delta_lin: float, delta_lo: float, n_pos_folds: int, min_folds: int = TG_MIN_FOLDS) -> bool:
    return bool(np.isfinite(delta_lin) and delta_lin >= TG_MIN_DELTA and delta_lo > 0 and n_pos_folds >= min_folds)


def pick_shortlist(rows: pd.DataFrame, se_of: dict) -> dict:
    """W2 decision (≤ 1): among TG-passing eligible steps, the paired one-SE rule (§9): candidates within one
    paired-difference SE of the best; simplest first (M1 < M2(f); among M2(f) the larger Δ̄_lin). ``rows``: one row
    per eligible step with step, delta_lin, tg_pass. ``se_of``: step -> (Δ̄_lin(best − step), SE)."""
    el = rows[rows.tg_pass].copy()
    if el.empty:
        return {"shortlist": [], "decision": "W2 negative: no eligible gate-free ladder variant passes TG; M0 is kept"}
    best = el.sort_values("delta_lin", ascending=False).iloc[0].step
    el["diff_vs_best"] = el.step.map(lambda s: se_of[s][0])
    el["se_diff"] = el.step.map(lambda s: se_of[s][1])
    el["within_1se"] = el.diff_vs_best <= el.se_diff
    el["complexity"] = el.step.map(lambda s: 0 if s == "M1" else 1)
    pick = el[el.within_1se].sort_values(["complexity", "delta_lin"], ascending=[True, False]).iloc[0]
    return {"shortlist": [pick.step], "best_by_delta": best,
            "within_1se": el.loc[el.within_1se, "step"].tolist(),
            "decision": f"W2 shortlist: {pick.step} (one-SE rule among {el.step.tolist()})"}


# ----------------------------------------------------------------------------------------------- attribution
def flip_counts(calls_ref: pd.DataFrame, calls_new: pd.DataFrame, Y: pd.DataFrame) -> dict:
    """Attribution (ii): generic (reference) binary errors flipped to correct, net of correct calls flipped to wrong,
    over the non-NaN label pairs of ``Y`` (lines x genes, identical index/columns to the call frames)."""
    m = Y.notna().to_numpy()
    y = np.where(m, Y.to_numpy(dtype=float), 0.0).astype(bool)
    a = calls_ref.loc[Y.index, Y.columns].to_numpy(bool)
    b = calls_new.loc[Y.index, Y.columns].to_numpy(bool)
    wrong_a, wrong_b = (a != y) & m, (b != y) & m
    fixed = wrong_a & ~wrong_b
    broken = ~wrong_a & wrong_b & m
    out = {"n_pairs": int(m.sum()), "errors_ref": int(wrong_a.sum()), "fp_ref": int((a & ~y & m).sum()),
           "fn_ref": int((~a & y & m).sum()), "fixed": int(fixed.sum()), "broken": int(broken.sum()),
           "fixed_fp": int((fixed & a).sum()), "fixed_fn": int((fixed & ~a).sum()),
           "n_flipped_pairs": int(((a != b) & m).sum())}
    out["net_flip_fraction"] = (out["fixed"] - out["broken"]) / out["errors_ref"] if out["errors_ref"] else np.nan
    out["per_line_fixed"] = pd.Series(fixed.sum(1), index=Y.index)
    out["per_line_broken"] = pd.Series(broken.sum(1), index=Y.index)
    out["per_line_errors"] = pd.Series(wrong_a.sum(1), index=Y.index)
    return out


def oracle_ratio(d_doc: pd.Series, d_oracle: pd.Series, lineage: pd.Series, *, B: int, seed) -> dict:
    """Attribution (iii): (doc − M0)/(oracle − M0) on lineage-equal means, with a joint lineage-stratified line
    bootstrap (the same resampling of lines for numerator and denominator; percentile CI)."""
    from nextgem import metrics

    idx = d_doc.dropna().index.intersection(d_oracle.dropna().index)
    lin = lineage.reindex(idx)
    w = metrics.lineage_weights(lin).to_numpy()
    num_v, den_v = d_doc.loc[idx].to_numpy(), d_oracle.loc[idx].to_numpy()
    num, den = float(w @ num_v), float(w @ den_v)
    C = metrics.stratified_bootstrap_counts(lin, B, seed)
    rn, rd = C @ (w * num_v), C @ (w * den_v)
    with np.errstate(divide="ignore", invalid="ignore"):
        rr = np.where(rd > 0, rn / rd, np.nan)
    rr = rr[np.isfinite(rr)]
    lo, hi = (np.quantile(rr, [0.025, 0.975]) if rr.size else (np.nan, np.nan))
    return {"doc_minus_m0": num, "oracle_minus_m0": den, "ratio": num / den if den > 0 else np.nan,
            "ratio_lo": float(lo), "ratio_hi": float(hi), "n_lines": int(len(idx)), "n_boot_finite": int(rr.size)}


def shuffle_eligible_lineages(media_key: pd.Series, lineage: pd.Series, min_lines: int = SHUFFLE_MIN_LINES) -> list:
    """Lineages with >= ``min_lines`` lines in each of at least two media (DESIGN W2 shuffled-medium power rule).
    'Medium' = base medium class of the key (text before '+' / '|')."""
    base = media_key.map(lambda k: re.split(r"[+|]", str(k))[0])
    tab = pd.crosstab(lineage.reindex(base.index), base)
    return sorted(tab.index[(tab >= min_lines).sum(axis=1) >= 2])


# ----------------------------------------------------------------------------------------------- flip analysis
def addition_prune_safe(wt_new: float, wt_ref: float, r_min_nonlethal: float, cutoff: float = 0.01) -> bool:
    """Exact pruning for a nutrient *addition* (bounds only relax, so KO_new >= KO_ref and WT_new >= WT_ref):
    a reference non-lethal gene (ratio_ref >= r_min_nonlethal >= cutoff) keeps ratio_new >= ratio_ref·WT_ref/WT_new
    >= cutoff whenever WT_new/WT_ref <= r_min_nonlethal/cutoff. Then only the reference-lethal genes can flip."""
    return wt_ref > 0 and wt_new / wt_ref <= r_min_nonlethal / cutoff


def ko_subset_task(args) -> dict:
    """Single-gene KO ratio_r for a subset of a parent's genes under ``bounds`` (worker of ``w2_engine.make_pool``;
    same semantics and rounding as ``w2_engine.ko_task``, incl. the identical-closure grouping and the exact
    zero-flux shortcut). Genes not listed are returned as NaN."""
    import math

    from optlang.interface import INFEASIBLE, OPTIMAL

    from nextgem import w2_engine as eng
    from nextgem.knockouts import RATIO_DECIMALS, WT_GROWTH_TOL

    parent, tag, bounds, subset = args[:4]
    m = eng._MODELS[parent]
    if len(args) > 4 and args[4]:
        m.solver.configuration.timeout = args[4]  # LP time limit (s); a time-limited LP raises (never silently 0)
    genes = eng._GENES[parent]
    clos = eng._CLOSURES[parent]
    pos = {g: i for i, g in enumerate(genes)}
    idx = [pos[g] for g in subset]
    ko = np.full(len(genes), np.nan)
    with m:
        m.medium = bounds
        wt, wst = eng._solve(m)
        if wst != OPTIMAL:
            raise RuntimeError(f"{tag}: WT status {wst}")
        pv = m.solver.primal_values
        carrying = {rid for rid, fwd, rev in eng._VARNAMES[parent] if abs(pv[fwd]) + abs(pv[rev]) > eng.ZERO_FLUX_TOL}
        groups: dict[frozenset, list[int]] = {}
        for i in idx:
            c = clos[genes[i]]
            if not c or not any(r in carrying for r in c):
                ko[i] = 1.0
            else:
                groups.setdefault(frozenset(c), []).append(i)
        for members in groups.values():
            with m:
                m.genes.get_by_id(genes[members[0]]).knock_out()
                g, st = eng._solve(m)
            if st not in (OPTIMAL, INFEASIBLE):
                raise RuntimeError(f"{tag}: KO {genes[members[0]]} status {st}")
            for i in members:
                ko[i] = g / wt if wt > WT_GROWTH_TOL else math.nan
    ratio_r = np.clip(np.round(ko, RATIO_DECIMALS), 0.0, 1.0) + 0.0
    return dict(tag=tag, parent=parent, wt_growth=wt, ratio_r=ratio_r, n_solved=len(groups))


def ko_full_task(args) -> dict:
    """``w2_engine.ko_task`` with an LP time limit on the worker's model (GLPK stalls on rare degenerate LPs). A
    time-limited LP makes ko_task raise; the error is returned as ``status = 'timeout'`` instead of a ratio vector."""
    from nextgem import w2_engine as eng

    parent, tag, bounds, verify, timeout = args
    eng._MODELS[parent].solver.configuration.timeout = timeout
    try:
        r = eng.ko_task((parent, tag, bounds, verify))
        r["status"] = "ok"
        return r
    except RuntimeError as e:
        return {"sweep_id": tag, "parent": parent, "status": "timeout", "error": str(e)[:200]}
