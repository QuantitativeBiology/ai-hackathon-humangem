"""W1 baseline benchmark library (DESIGN.md §5–§6, §8, W1). Pure functions over label frames passed in.

Shared by ``scripts/w1_score_train.py`` (training tier, LOLO) and ``scripts/evaluate_val.py`` (the single W1
validation look). This module never opens a label file: label frames come from ``labels.train()`` /
``labels.fold(k)`` / ``labels.validation()`` in the calling entry point.

Conventions
-----------
- Score frames are lines x genes (columns = Entrez int, the gene universe U), higher = more essential.
- Generic FBA score e = 1 - ratio_r (§6.1); a gene absent from a parent scores 0 (bottom tied block, §3.3).
  Task-augmented (§6.1): e = 1 if the KO fails a P57 task, else 1 - ratio_r.
- Fitted baselines (B2, B2m, B2c, CE_fit, selective genes, B6) are fitted on the *fitting lines* passed in
  (training minus fold for LOLO; training for validation looks).
- B8+/B8-: rank by e of the parent, ties broken by higher / lower expression (dense lexicographic rank, exact).
  B8': e with genes at log2(TPM+1) < 1 set to 0.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from nextgem import metrics, paths
from nextgem.config import CE_FIT_FRACTION, SEEDS

GRID_LABELS = ("L_prob", "L_ge05", "L_ge10")
PARENTS = ("human1", "human2")
SHORT = {"human1": "H1", "human2": "H2"}
B2M_MIN_CLASS_LINES = 10  # §8: B2m falls back to B2 if the class has < 10 fitting lines
EXPRESSED = 1.0  # log2(TPM+1) >= 1 (TPM >= 1), §5.4
SELECTIVE = (0.05, 0.95)  # §5.4
B_OTHER = metrics.B_OTHER


def w1_dir() -> Path:
    return paths.results_dir() / "w1"


# ----------------------------------------------------------------------------- label-free inputs
def gene_universe() -> pd.DataFrame:
    from nextgem import labels

    return labels.gene_universe()


def generic_table(U: pd.Series) -> pd.DataFrame:
    """Per U gene (Entrez): ratio_r, e, task-augmented e and calls for both parents (label-free)."""
    out = pd.DataFrame(index=pd.Index(sorted(U.astype(int)), name="entrez"))
    for p in PARENTS:
        ko = pd.read_parquet(paths.baseline_dir() / f"{p}_ko.parquet", columns=["entrez", "ratio_r"]).dropna(subset=["entrez"])
        ko = ko.assign(entrez=ko.entrez.astype(int)).drop_duplicates("entrez").set_index("entrez").ratio_r
        s = SHORT[p]
        out[f"in_{s}"] = out.index.isin(ko.index)
        out[f"ratio_{s}"] = ko.reindex(out.index)
        out[f"e_{s}"] = (1.0 - out[f"ratio_{s}"]).fillna(0.0).round(6)
        out[f"fba_call_{s}"] = out[f"ratio_{s}"] < metrics.FBA_CUTOFF  # NaN (absent) -> False
        pf = w1_dir() / f"p57_{p}.parquet"
        if pf.exists():
            t = pd.read_parquet(pf, columns=["entrez", "p57", "pgr"]).dropna(subset=["entrez"])
            t = t.assign(entrez=t.entrez.astype(int)).drop_duplicates("entrez").set_index("entrez")
            out[f"p57_{s}"] = t.p57.astype(bool).reindex(out.index, fill_value=False)
            out[f"pgr_{s}"] = t.pgr.astype(bool).reindex(out.index, fill_value=False)
            out[f"et_{s}"] = np.where(out[f"p57_{s}"], 1.0, out[f"e_{s}"])
    return out


def load_omics() -> dict[str, pd.DataFrame]:
    """Label-free matrices written by scripts/w1_inputs.py (lines x Entrez)."""
    d = w1_dir() / "tmp"
    out = {}
    for k in ("expr", "cn", "mut_dam", "mut_hot"):
        df = pd.read_parquet(d / f"{k}.parquet").set_index("ModelID")
        df.columns = df.columns.astype(int)
        out[k] = df.astype("float32")
    return out


def line_features() -> pd.DataFrame:
    return pd.read_csv(w1_dir() / "line_features.tsv", sep="\t").set_index("ModelID")


def gene_features() -> pd.DataFrame:
    return pd.read_csv(w1_dir() / "gene_features.tsv", sep="\t").set_index("entrez")


def gene_panels() -> pd.DataFrame:
    return pd.read_csv(paths.data_root() / "controls" / "gene_panels.tsv", sep="\t")


# ----------------------------------------------------------------------------- label frames
def wide(long: pd.DataFrame, col: str, lines, genes) -> pd.DataFrame:
    """Long label frame -> lines x genes float (NaN = excluded pair or missing)."""
    w = long.pivot(index="ModelID", columns="entrez", values=col)
    w = w.reindex(index=list(lines), columns=list(genes))
    return w.astype(float)


# ----------------------------------------------------------------------------- fitted components
@dataclass
class Fitted:
    """Components fitted on one set of fitting lines (one LOLO fold, or all training for validation)."""

    fit_lines: list
    prior: dict = field(default_factory=dict)  # label -> Series(entrez): fraction positive (B2)
    prior_class: dict = field(default_factory=dict)  # label -> DataFrame(class x entrez)
    class_n: pd.Series | None = None  # fitting lines per medium class
    b2c: pd.Series | None = None  # -mean gene effect (B2c)
    ce: dict = field(default_factory=dict)  # label -> set(entrez): CE_fit (>= 90% of fitting lines)
    rate: dict = field(default_factory=dict)  # label -> fitting-line positive rate (pairs)


def fit_components(fit_long: pd.DataFrame, line_class: pd.Series, labels_=GRID_LABELS + ("L_ge06",)) -> Fitted:
    f = Fitted(fit_lines=sorted(fit_long.ModelID.unique()))
    cls = fit_long.ModelID.map(line_class).fillna("unknown")
    f.class_n = line_class.reindex(f.fit_lines).fillna("unknown").value_counts()
    for lab in labels_:
        d = fit_long[["entrez", lab]].assign(cls=cls.to_numpy())
        d = d[d[lab].notna()]
        y = d[lab].astype(float)
        f.prior[lab] = y.groupby(d.entrez).mean()
        f.prior_class[lab] = y.groupby([d.cls, d.entrez]).mean().unstack("entrez")
        f.ce[lab] = set(f.prior[lab][f.prior[lab] >= CE_FIT_FRACTION].index.astype(int))
        f.rate[lab] = float(y.mean())
    ge = fit_long[["entrez", "gene_effect"]].dropna()
    f.b2c = -ge.groupby("entrez").gene_effect.mean()
    return f


def _const(series: pd.Series, lines, genes, fill=0.0) -> pd.DataFrame:
    v = series.reindex(list(genes)).fillna(fill).to_numpy(dtype=float)
    return pd.DataFrame(np.tile(v, (len(lines), 1)), index=list(lines), columns=list(genes))


def _expr(expr: pd.DataFrame, lines, genes) -> pd.DataFrame:
    return expr.reindex(index=list(lines), columns=list(genes)).astype(float).fillna(0.0)


def b8(e: pd.Series, expr_m: pd.DataFrame, direction: int) -> pd.DataFrame:
    """Exact lexicographic dense rank of (e, ±expression) per line (§8 B8+/B8-)."""
    E = np.rint(e.reindex(expr_m.columns).fillna(0.0).to_numpy() * 1e6).astype(np.int64)
    out = np.empty(expr_m.shape)
    X = expr_m.to_numpy()
    for i in range(X.shape[0]):
        r = pd.Series(direction * X[i]).rank(method="dense").to_numpy().astype(np.int64)
        out[i] = E * 100_000 + r
    return pd.DataFrame(out, index=expr_m.index, columns=expr_m.columns)


def score_matrices(lines, genes, fitted: Fitted, gen: pd.DataFrame, line_class: pd.Series, expr: pd.DataFrame,
                   b6: pd.DataFrame | None = None, random_seeds=SEEDS[0:5]) -> dict[str, pd.DataFrame]:
    """All W1 methods for ``lines`` x ``genes``. Label-specific ones are keyed ``B2[L]`` / ``B2m[L]``."""
    lines, genes = list(lines), list(genes)
    S: dict[str, pd.DataFrame] = {}
    for s in ("H1", "H2"):
        S[s] = _const(gen[f"e_{s}"], lines, genes)
        if f"et_{s}" in gen:
            S[f"{s}_task"] = _const(gen[f"et_{s}"], lines, genes)
    S["B0"] = pd.DataFrame(0.0, index=lines, columns=genes)
    for sd in random_seeds:  # seeded per (seed, line), so a line's draw does not depend on its fold
        S[f"B1_s{sd}"] = pd.DataFrame(np.vstack([np.random.default_rng([sd, _line_num(l)]).random(len(genes)) for l in lines]),
                                      index=lines, columns=genes)
    X = _expr(expr, lines, genes)
    S["B3"] = X
    for lab, prior in fitted.prior.items():
        if lab == "L_ge06":
            continue
        S[f"B2[{lab}]"] = _const(prior, lines, genes)
        pc = fitted.prior_class[lab]
        rows = []
        for l in lines:
            c = line_class.get(l, "unknown")
            if fitted.class_n.get(c, 0) >= B2M_MIN_CLASS_LINES and c in pc.index:
                rows.append(pc.loc[c].reindex(genes).fillna(prior.reindex(genes)).fillna(0.0).to_numpy())
            else:
                rows.append(prior.reindex(genes).fillna(0.0).to_numpy())
        S[f"B2m[{lab}]"] = pd.DataFrame(np.vstack(rows), index=lines, columns=genes)
    S["B2c"] = _const(fitted.b2c, lines, genes, fill=float(fitted.b2c.min()))
    for s in ("H1", "H2"):
        e = gen[f"e_{s}"]
        S[f"B8+_{s}"] = b8(e, X, +1)
        S[f"B8-_{s}"] = b8(e, X, -1)
        S[f"B8'_{s}"] = S[s].where(X >= EXPRESSED, 0.0)
    if b6 is not None:
        S["B6[L_prob]"] = b6.reindex(index=lines, columns=genes).astype(float)
        if S["B6[L_prob]"].isna().any().any():
            raise ValueError("B6 predictions missing for some pairs")
    return S


def _line_num(model_id: str) -> int:
    return int(str(model_id).split("-")[-1])


# ----------------------------------------------------------------------------- B6 features (§8)
B6_FEATURES = ["b2", "b2m", "expr", "cn", "mut_dam", "mut_hot", "neigh_mean_expr", "max_paralog_expr",
               "n_isozyme_genes_h2", "medium_class", "serum"]


def expression_aggregates(expr: pd.DataFrame, genes, neighbours: dict, paralogs: dict):
    """Per line x gene: mean expression of H2 enzyme neighbours; max expression of Ensembl paralogs."""
    cols = {c: i for i, c in enumerate(expr.columns)}
    X = expr.to_numpy(dtype="float32")
    nm = np.full((X.shape[0], len(genes)), np.nan, dtype="float32")
    pm = np.full((X.shape[0], len(genes)), np.nan, dtype="float32")
    for j, g in enumerate(genes):
        idx = [cols[e] for e in neighbours.get(str(g), []) if e in cols]
        if idx:
            nm[:, j] = np.nanmean(X[:, idx], axis=1)
        idx = [cols[e] for e in paralogs.get(str(g), []) if e in cols]
        if idx:
            pm[:, j] = np.nanmax(X[:, idx], axis=1)
    return (pd.DataFrame(nm, index=expr.index, columns=list(genes)), pd.DataFrame(pm, index=expr.index, columns=list(genes)))


def b6_matrix(lines, genes, b2_rows: pd.DataFrame, b2m_rows: pd.DataFrame, omics: dict, aggs, gfeat: pd.DataFrame,
              lf: pd.DataFrame, class_codes: dict) -> np.ndarray:
    """(n_lines * n_genes, len(B6_FEATURES)) float32, row-major over lines then genes. NaN = missing."""
    lines, genes = list(lines), list(genes)
    n, m = len(lines), len(genes)

    def mat(df):
        return df.reindex(index=lines, columns=genes).to_numpy(dtype="float32")

    iso = gfeat["n_isozyme_genes_h2"].reindex(genes).to_numpy(dtype="float32")
    med = np.array([class_codes.get(lf.loc[l, "base_medium_class"], -1) if l in lf.index else -1 for l in lines], dtype="float32")
    ser = np.array([lf.loc[l, "serum"] if l in lf.index else np.nan for l in lines], dtype="float32")
    cols = [mat(b2_rows), mat(b2m_rows), mat(omics["expr"]), mat(omics["cn"]), mat(omics["mut_dam"]), mat(omics["mut_hot"]),
            mat(aggs[0]), mat(aggs[1]), np.tile(iso, (n, 1)), np.repeat(med[:, None], m, 1), np.repeat(ser[:, None], m, 1)]
    return np.stack([c.reshape(-1) for c in cols], axis=1)


def heldout_priors(fitted, lines, line_class, genes, label="L_prob"):
    """B6 prior features of held-out lines: B2 / B2m from all fitting lines (fallback B2 if class < 10 lines)."""
    p = fitted.prior[label].reindex(genes)
    pc = fitted.prior_class[label]
    b2 = pd.DataFrame(np.tile(p.to_numpy(), (len(lines), 1)), index=lines, columns=genes)
    rows = []
    for l in lines:
        c = line_class.get(l, "unknown")
        if fitted.class_n.get(c, 0) >= B2M_MIN_CLASS_LINES and c in pc.index:
            rows.append(pc.loc[c].reindex(genes).fillna(p).to_numpy())
        else:
            rows.append(p.to_numpy())
    return b2, pd.DataFrame(np.vstack(rows), index=lines, columns=genes)


def methods_for(label: str, names) -> list[str]:
    out = []
    for n in names:
        if "[" in n and not n.endswith(f"[{label}]"):
            continue
        out.append(n)
    return out


# ----------------------------------------------------------------------------- masks / strata
def stratum_specs(panels: pd.DataFrame, gfeat: pd.DataFrame) -> dict:
    """§5.4 reporting strata: name -> (kind, payload). kind: 'rows' (line filter on line features),
    'genes_in' / 'genes_ex' (gene sets), 'pairs' (callable on (expr, cn) frames), 'selective'."""
    def pset(section_sub_prefix: str):
        sec, sub = section_sub_prefix.split(":", 1)
        m = (panels.section.astype(str) == sec) & panels.subpanel.astype(str).str.startswith(sub)
        return set(panels.loc[m, "entrez"].astype(int))

    multicopy = pset("1.6:b")
    housekeeping = pset("1.6:g")
    moonlighting = pset("1.6:f")
    pseudo = pset("1.4:d")
    addr = {s: set(gfeat.index[gfeat[f"addressable_{p}"]].astype(int)) for p, s in (("human1", "H1"), ("human2", "H2"))}
    return {
        "none": ("none", None),
        "expressed_pairs": ("pairs", "expressed"),
        "selective_genes": ("selective", None),
        "addressable_H2": ("genes_in", addr["H2"]),
        "addressable_H1": ("genes_in", addr["H1"]),
        "lib_Avana": ("rows", ("library_class", "Avana-containing")),
        "lib_Humagne": ("rows", ("library_class", "Humagne-only")),
        "lib_KY": ("rows", ("library_class", "KY-only")),
        "organoid": ("rows", ("model_type", "Organoid")),
        "cell_line": ("rows", ("model_type", "Cell Line")),
        "cn_amp4_in": ("pairs", "cn_amp4_in"),
        "cn_amp4_ex": ("pairs", "cn_amp4_ex"),
        "multicopy_in": ("genes_in", multicopy),
        "multicopy_ex": ("genes_ex", multicopy),
        "housekeeping_in": ("genes_in", housekeeping),
        "housekeeping_ex": ("genes_ex", housekeeping),
        "moonlighting_in": ("genes_in", moonlighting),
        "moonlighting_ex": ("genes_ex", moonlighting),
        "pseudoreaction_ex": ("genes_ex", pseudo),
    }


def apply_stratum(Y: pd.DataFrame, spec, lf: pd.DataFrame, expr, cn, sel_mask: pd.DataFrame | None) -> pd.DataFrame:
    kind, payload = spec
    if kind == "none":
        return Y
    if kind == "rows":
        col, val = payload
        keep = [l for l in Y.index if lf.loc[l, col] == val]
        return Y.loc[keep]
    if kind == "genes_in":
        return Y.where(pd.DataFrame(np.broadcast_to(Y.columns.isin(list(payload)), Y.shape), index=Y.index, columns=Y.columns))
    if kind == "genes_ex":
        return Y.where(pd.DataFrame(np.broadcast_to(~Y.columns.isin(list(payload)), Y.shape), index=Y.index, columns=Y.columns))
    if kind == "selective":
        return Y.where(sel_mask.reindex(index=Y.index, columns=Y.columns).fillna(False))
    if kind == "pairs":
        if payload == "expressed":
            X = expr.reindex(index=Y.index, columns=Y.columns)
            return Y.where(X >= EXPRESSED)
        C = cn.reindex(index=Y.index, columns=Y.columns)
        amp = C >= 4.0
        return Y.where(amp) if payload == "cn_amp4_in" else Y.where(~amp)
    raise ValueError(kind)


# ----------------------------------------------------------------------------- cell summaries
def method_summary(S: pd.DataFrame, Y: pd.DataFrame, lineage: pd.Series, *, seed=SEEDS[0], B=B_OTHER,
                   gene_unit=True, auroc=False):
    lt = metrics.per_line(S, Y)
    gt = metrics.per_gene(S, Y) if gene_unit else None
    pl = metrics.pooled(S, Y)
    lin = lineage.reindex(lt.index)
    est = lt[lt.estimable]
    row = {"n_lines": len(lt), "n_lines_estimable": int(lt.estimable.sum()),
           "na": metrics.cell_is_na(lt, gt), "n_lineages": int(lin[lt.estimable].nunique())}
    if len(est):
        bt = metrics.bootstrap_mean(est.tap, lin, B=B, seed=seed)
        e_lin = metrics.lineage_equal_mean(est.e_random, lin)
        row.update({"tap_lin": bt.estimate, "tap_lin_lo": bt.ci_low, "tap_lin_hi": bt.ci_high,
                    "e_random_lin": e_lin, "lift_lin": bt.estimate / e_lin,
                    "delta_null_lin": metrics.lineage_equal_mean(est.delta_null, lin),
                    "tap_line": float(est.tap.mean()), "e_random_line": float(est.e_random.mean()),
                    "lift_line": float(est.tap.mean() / est.e_random.mean()),
                    "prevalence_lin": metrics.lineage_equal_mean(est.prevalence, lin)})
    if gt is not None:
        row.update({"n_genes_estimable": int(gt.estimable.sum()), "m_delta_ap_gene": metrics.m_delta_ap_gene(gt)})
    row.update({"pooled_tap": pl.tap, "pooled_e_random": pl.e_random, "pooled_lift": pl.lift, "pooled_n": pl.n,
                "pooled_n_pos": pl.n_pos})
    if auroc:
        S2, Y2 = S.loc[est.index], Y.loc[est.index]
        a = [metrics.auroc(S2.loc[l][Y2.loc[l].notna()].to_numpy(), Y2.loc[l].dropna().to_numpy()) for l in est.index]
        row["auroc_lin"] = metrics.lineage_equal_mean(pd.Series(a, index=est.index), lin)
        ap = [metrics.ap_step_sklearn(S2.loc[l][Y2.loc[l].notna()].to_numpy(), Y2.loc[l].dropna().to_numpy()) for l in est.index]
        row["ap_sklearn_lin"] = metrics.lineage_equal_mean(pd.Series(ap, index=est.index), lin)  # sensitivity row only (§6.1)
    return row, lt, gt


def delta_summary(lt_a: pd.DataFrame, lt_b: pd.DataFrame, lineage: pd.Series, *, seed=SEEDS[0], B=B_OTHER) -> dict:
    d = metrics.paired_line_delta(lt_a, lt_b)
    if d.dropna().empty:
        return {"delta_lin": np.nan, "delta_lo": np.nan, "delta_hi": np.nan, "delta_line": np.nan}
    lin = lineage.reindex(d.index)
    bt = metrics.bootstrap_mean(d, lin, B=B, seed=seed)
    return {"delta_lin": bt.estimate, "delta_lo": bt.ci_low, "delta_hi": bt.ci_high,
            "delta_line": metrics.line_weighted_mean(d)}


def ce_mask(lines, genes, line_ce: dict) -> pd.DataFrame:
    """True where the pair is kept (gene not in the line's CE_fit)."""
    g = np.asarray(list(genes))
    return pd.DataFrame(np.vstack([~np.isin(g, list(line_ce[l])) for l in lines]), index=list(lines), columns=list(genes))


def selective_mask(lines, genes, line_prior: dict) -> pd.DataFrame:
    rows = []
    for l in lines:
        p = line_prior[l].reindex(list(genes))
        rows.append(((p >= SELECTIVE[0]) & (p <= SELECTIVE[1])).fillna(False).to_numpy())
    return pd.DataFrame(np.vstack(rows), index=list(lines), columns=list(genes))


# ----------------------------------------------------------------------------- §6.8 binary endpoints
def parent_gene_set(gen: pd.DataFrame, parent_short: str) -> list[int]:
    """§6.8: the parent model's genes that have a CRISPR column (= U ∩ parent)."""
    return list(gen.index[gen[f"in_{parent_short}"]])


def binary_endpoints(calls: pd.Series, truth: pd.DataFrame, lineage: pd.Series, genes, *, seed=SEEDS[0], B=B_OTHER):
    """Per-line confusion + metrics on ``genes`` (calls: Series over Entrez, line-constant for generic models)."""
    T = truth.loc[:, list(genes)]
    C = pd.DataFrame(np.tile(calls.reindex(list(genes)).fillna(False).to_numpy(bool), (len(T), 1)),
                     index=T.index, columns=T.columns)
    pl = metrics.binary_per_line(C, T)
    lin = lineage.reindex(pl.index)
    summ = metrics.summarise_binary(pl, lin)
    cis = {}
    for m in metrics.BINARY_METRICS:
        v = pl[m].dropna()
        if len(v):
            bt = metrics.bootstrap_mean(v, lin, B=B, seed=seed)
            cis[m] = (bt.ci_low, bt.ci_high)
    summ["ci_lo"] = [cis.get(m, (np.nan, np.nan))[0] if p == "exclude" else np.nan for m, p in zip(summ.metric, summ.nan_policy)]
    summ["ci_hi"] = [cis.get(m, (np.nan, np.nan))[1] if p == "exclude" else np.nan for m, p in zip(summ.metric, summ.nan_policy)]
    return summ, pl


def config_sha(obj) -> str:
    import hashlib

    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()
