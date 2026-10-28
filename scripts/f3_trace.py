"""F3 (A15 follow-up; post hoc / descriptive): pFBA traces of the 20 worst T1 lineage x gene cells of P1
(DESIGN.md §10.3 "The 20 worst lineage x gene cells are traced by pFBA"; §17 A15). Fits no rule, changes no model
or score, opens no label file: inputs are the unsealed T1 worst-cell list (results/m6/test/unsealed/T1/), the
label-free P1 build inputs (line assignment, sweep_e, the W3 B-loo spec, the W2 media library) and label-free
gene / line features.

P1 = Human1 + W3 ``B-loo:cofactor_pool_biomass`` on documented medium M2(f=1), gate-free (RESEARCH_TREE D11). For
each cell the lines are the lineage's test-tier lines with an expression matrix (P1's evaluation lines); each line
has one M2(f=1) bounds signature (results/m6/scores/line_assignment.tsv), and a cell is traced once per distinct
signature (bounds rebuilt by ``scripts/m6_sweeps.bounds_by_id``; KO engine settings as ``nextgem.w3.ko_sweep``:
GLPK, biomass MAR13082 objective, ``Gene.knock_out``).

Per (cell, medium):
  gene KO           growth ratio KO / WT (recomputed; checked against 1 - e of P1's sweep_e for that sweep);
  reaction KO       counterfactual: every reaction the gene appears in is closed (all isozymes / alternative
                    complexes off) -> tells whether the GPR OR alone hides the lethality;
  pFBA              WT and, where the gene's reactions carry WT pFBA flux and the KO stays viable, KO; top
                    rerouted reactions and exchanges by |v_KO - v_WT|;
  FVA               the gene's reactions at >= 99.9 % of WT growth and with growth free (blocked or not);
  bypass test       each exchange whose uptake rises in the KO pFBA (top 5) is closed in the KO and in the WT; an
                    exchange is a medium bypass if KO growth without it is < 1 % of WT growth without it.
Mechanism class of a false negative (first match, see ``classify_fn``):
  outside_model            gene not in Human1 (outside the addressable set)
  isozyme_in_gpr           KO closes no reaction / not all, and the reaction KO is lethal (ratio < 0.01)
  medium_uptake_bypass     reaction KO viable; WT pFBA uses the reactions; the KO relies on an extra uptake
  alternative_pathway      reaction KO viable; WT pFBA uses the reactions; no single-exchange bypass
  function_not_represented reaction KO viable and the reactions carry no WT pFBA flux (blocked or unused):
                           nothing the biomass objective needs depends on the gene's modelled function
False positives (none among the 20 cells) are traced by ``blocked_biomass_components`` (which biomass precursor
the KO cannot make).

Outputs (results/f3/): cells.tsv (one row per cell), cells_by_medium.tsv, reroutes.tsv, crosstab.tsv,
crosstab_long.tsv, f3_summary.json, MANIFEST.tsv. README.md is written by hand from these tables.

Usage: OMP_NUM_THREADS=1 .venv/bin/python scripts/f3_trace.py --workers 6 [--allow-dirty]
"""

from __future__ import annotations

import argparse
import ast
import concurrent.futures as cf
import json
import multiprocessing
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "results" / "f3"
UNSEALED = PROJECT / "results" / "m6" / "test" / "unsealed" / "T1"
WORST = UNSEALED / "worst_cells.tsv"
FAILURE_ACCOUNT = UNSEALED / "failure_account.tsv"
LINE_ASSIGNMENT = PROJECT / "results" / "m6" / "scores" / "line_assignment.tsv"
SWEEP_E = PROJECT / "results" / "m6" / "scores" / "sweep_e.parquet"
ANNOTATIONS = PROJECT / "results" / "m6" / "t1" / "p1_gene_annotations.parquet"
GENE_MAP = PROJECT / "results" / "baseline" / "human1_gene_mapping.tsv"
AUTHOR = "F3 engineer agent (Claude Opus)"
MAX_WORKERS = 6

PARENT = "human1"
STEP = "M2(f=1)"
BIOMASS = "MAR13082"
LETHAL = 0.01  # P1's fixed FBA call: ratio_r < 0.01
FLUX_TOL = 1e-6
FVA_FRACTION = 0.999
TOP_REACTIONS = 10
TOP_EXCHANGES = 8
BYPASS_CANDIDATES = 5
EXPRESSED = 1.0  # failure-account "unexpressed" = expr < 1.0 (log2(TPM + 1))
CN_AMP = 4.0
FA_GENE_CLASSES = ("paralog_buffered", "medium_sensitive", "complex_subunit", "housekeeping_adjacent",
                   "biomass_artefact_lethal", "expected_fba_failure")
FN_CLASSES = ("outside_model", "isozyme_in_gpr", "medium_uptake_bypass", "alternative_pathway",
              "function_not_represented")


# ============================================================================= pure helpers (tested on toys)
def gpr_branches(rule: str) -> list[frozenset]:
    """Top-level OR branches of a GPR rule, each as the set of genes it names ([] for an empty rule)."""
    rule = (rule or "").strip()
    if not rule:
        return []
    body = ast.parse(rule, mode="eval").body

    def names(node) -> frozenset:
        return frozenset(n.id for n in ast.walk(node) if isinstance(n, ast.Name))

    def flatten_or(node) -> list:
        if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.Or):
            out = []
            for v in node.values:
                out.extend(flatten_or(v))
            return out
        return [node]

    return [names(b) for b in flatten_or(body)]


def gpr_structure(rule: str, gene: str) -> dict:
    """AND/OR structure of ``rule`` as seen from ``gene``: kind, branches, whether the gene is a complex subunit,
    and the genes of the OR branches that do not need it (the alternatives that survive its KO)."""
    br = gpr_branches(rule)
    if not br:
        kind = "none"
    elif len(br) == 1:
        kind = "single" if len(br[0]) == 1 else "AND"
    elif all(len(b) == 1 for b in br):
        kind = "OR"
    else:
        kind = "OR-of-AND"
    with_g = [b for b in br if gene in b]
    without = [b for b in br if gene not in b]
    alt = sorted(set().union(*without)) if without else []
    return {"kind": kind, "n_branches": len(br), "n_genes": len(set().union(*br)) if br else 0,
            "n_branches_with_gene": len(with_g), "subunit": any(len(b) > 1 for b in with_g),
            "closed_by_ko": bool(with_g) and not without, "alternatives": alt}


def flux_diff(wt: pd.Series, ko: pd.Series, *, exclude=(), tol: float = FLUX_TOL) -> pd.DataFrame:
    """Reactions ranked by |v_KO - v_WT| (descending), |delta| > tol, ``exclude`` dropped."""
    ix = wt.index.union(ko.index)
    w, k = wt.reindex(ix, fill_value=0.0), ko.reindex(ix, fill_value=0.0)
    d = pd.DataFrame({"v_wt": w, "v_ko": k, "delta": k - w})
    d = d[(d.delta.abs() > tol) & ~d.index.isin(list(exclude))]
    return d.reindex(d.delta.abs().sort_values(ascending=False).index)


def uptake_increase(wt: pd.Series, ko: pd.Series, exchanges, tol: float = FLUX_TOL) -> list[str]:
    """Exchanges (uptake = negative flux) whose uptake is larger in the KO than in the WT, largest first."""
    ex = [e for e in exchanges if e in wt.index or e in ko.index]
    w = -wt.reindex(ex, fill_value=0.0).clip(upper=0.0)
    k = -ko.reindex(ex, fill_value=0.0).clip(upper=0.0)
    inc = (k - w)[(k - w) > tol]
    return list(inc.sort_values(ascending=False).index)


def classify_fn(in_model: bool, rxn_ko_ratio: float, wt_flux_on_reactions: float, bypass: list[str]) -> str:
    """Mechanism class of a false negative (module docstring)."""
    if not in_model:
        return "outside_model"
    if np.isfinite(rxn_ko_ratio) and rxn_ko_ratio < LETHAL:
        return "isozyme_in_gpr"
    if wt_flux_on_reactions <= FLUX_TOL:
        return "function_not_represented"
    if bypass:
        return "medium_uptake_bypass"
    return "alternative_pathway"


def modal_class(classes: pd.Series, weights: pd.Series) -> str:
    """Line-weighted modal class over media (ties -> FN_CLASSES order)."""
    w = weights.groupby(classes.to_numpy()).sum()
    best = w.max()
    return next(c for c in (*FN_CLASSES, *sorted(w.index)) if c in w.index and w[c] == best)


def crosstab(cells: pd.DataFrame, fa_cols) -> pd.DataFrame:
    """mechanism class x failure-account class: number of cells carrying each (a cell counts in every class it has)."""
    rows = []
    for mech, g in cells.groupby("mechanism_class"):
        r = {"mechanism_class": mech, "n_cells": len(g)}
        for c in fa_cols:
            r[c] = int(g[c].astype(bool).sum())
        rows.append(r)
    t = pd.DataFrame(rows)
    tot = {"mechanism_class": "all", "n_cells": len(cells), **{c: int(cells[c].astype(bool).sum()) for c in fa_cols}}
    return pd.concat([t, pd.DataFrame([tot])], ignore_index=True)


def blocked_biomass_components(model, biomass_id: str = BIOMASS, skip=("MAM01371c", "MAM02040c")) -> dict:
    """Max production of each biomass reactant (demand objective, biomass free); {met id: max}. For FP traces:
    a component whose max is ~0 in the KO but not in the WT is what the KO cannot make."""
    import cobra

    b = model.reactions.get_by_id(biomass_id)
    out = {}
    for m, c in b.metabolites.items():
        if c >= 0 or m.id in skip:
            continue
        with model:
            d = cobra.Reaction(f"F3_DM_{m.id}", lower_bound=0.0, upper_bound=1000.0)
            d.add_metabolites({m: -1.0})
            model.add_reactions([d])
            b.bounds = (0.0, 1000.0)
            model.objective = d.id
            v = model.slim_optimize(error_value=0.0)
        out[m.id] = float(v) if np.isfinite(v) else 0.0
    return out


# ============================================================================= model work (one worker per medium)
_W: dict = {}


def _init():
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem.model import load_model

    _W["model"] = load_model(PARENT)


def _growth(model) -> float:
    from nextgem import w3

    g, st = w3.solve_growth(model)
    return float(g) if st == "optimal" else 0.0


def _pfba(model) -> pd.Series | None:
    from cobra.flux_analysis import pfba

    try:
        return pfba(model).fluxes
    except Exception:  # noqa: BLE001 - infeasible KO
        return None


def _close(model, rids) -> None:
    for r in rids:
        model.reactions.get_by_id(r).bounds = (0.0, 0.0)


def trace_medium(job: dict) -> dict:
    """All per-gene traces for one medium signature (bounds + B-loo)."""
    from cobra.flux_analysis import flux_variability_analysis as fva
    from nextgem import w3

    m = _W["model"]
    t0 = time.time()
    v = w3.Variant.from_json(job["variant_json"])
    genes = job["genes"]  # [(entrez, ensembl or None, fp: bool)]
    exch = [r.id for r in m.exchanges]
    rows, rer, ex_rows = [], [], []
    with m:
        m.medium = job["bounds"]
        v.apply(m)
        m.objective = BIOMASS
        wt = _growth(m)
        wt_flux = _pfba(m)
        all_r = sorted({r.id for _, g, _ in genes if g for r in m.genes.get_by_id(g).reactions})
        f1 = fva(m, all_r, fraction_of_optimum=FVA_FRACTION, processes=1)
        f0 = fva(m, all_r, fraction_of_optimum=0.0, processes=1)
        wt_bio = None
        for entrez, g, is_fp in genes:
            row = {"medium_sweep": job["key"], "entrez": entrez, "ensembl": g, "wt_growth": wt}
            if not g:
                row.update(in_model=False, mechanism_class=classify_fn(False, np.nan, 0.0, []))
                rows.append(row)
                continue
            gene = m.genes.get_by_id(g)
            R = sorted(r.id for r in gene.reactions)
            with m:
                gene.knock_out()
                closed = [r for r in R if m.reactions.get_by_id(r).bounds == (0.0, 0.0)]
                gko = _growth(m)
                gko_flux = _pfba(m) if (closed and gko > FLUX_TOL and
                                        wt_flux[closed].abs().sum() > FLUX_TOL) else None
            with m:
                _close(m, R)
                rko = _growth(m)
                rko_flux = _pfba(m) if (rko > FLUX_TOL and wt_flux[R].abs().sum() > FLUX_TOL) else None
                bypass, tested = [], []
                if rko_flux is not None:
                    for e in uptake_increase(wt_flux, rko_flux, exch)[:BYPASS_CANDIDATES]:
                        with m:
                            m.reactions.get_by_id(e).lower_bound = 0.0
                            k_no = _growth(m)
                        tested.append((e, k_no))
            # WT growth without each tested exchange (the KO context is closed above)
            for e, k_no in tested:
                with m:
                    m.reactions.get_by_id(e).lower_bound = 0.0
                    w_no = _growth(m)
                rel = k_no / w_no if w_no > FLUX_TOL else np.nan
                ex_rows.append({"medium_sweep": job["key"], "entrez": entrez, "exchange": e,
                                "exchange_name": m.reactions.get_by_id(e).name,
                                "wt_growth_without": w_no, "rxnko_growth_without": k_no, "rel": rel,
                                "bypass": bool(np.isfinite(rel) and rel < LETHAL)})
                if np.isfinite(rel) and rel < LETHAL:
                    bypass.append(e)
            wt_on_R = float(wt_flux[R].abs().sum())
            row.update(in_model=True, reactions=";".join(R), closed_by_gene_ko=";".join(closed),
                       gene_ko_ratio=round(gko / wt, 6) if wt > 0 else np.nan,
                       rxn_ko_ratio=round(rko / wt, 6) if wt > 0 else np.nan,
                       wt_pfba_flux_on_reactions=wt_on_R,
                       wt_pfba_flux_on_closed=float(wt_flux[closed].abs().sum()) if closed else 0.0,
                       fva_opt=";".join(f"{r}:[{f1.loc[r, 'minimum']:.3g},{f1.loc[r, 'maximum']:.3g}]" for r in R),
                       blocked=bool(all(max(abs(f0.loc[r, "minimum"]), abs(f0.loc[r, "maximum"])) <= FLUX_TOL
                                        for r in R)),
                       required_at_opt=bool(any(f1.loc[r, "minimum"] > FLUX_TOL or f1.loc[r, "maximum"] < -FLUX_TOL
                                                for r in R)),
                       bypass_exchanges=";".join(bypass),
                       mechanism_class=classify_fn(True, rko / wt if wt > 0 else np.nan, wt_on_R, bypass))
            if is_fp:
                with m:
                    if wt_bio is None:
                        wt_bio = blocked_biomass_components(m)
                with m:
                    gene.knock_out()
                    kb = blocked_biomass_components(m)
                row["fp_blocked_components"] = ";".join(k for k in kb if kb[k] <= FLUX_TOL < wt_bio.get(k, 0.0))
            rows.append(row)
            for kind, fl in (("gene_ko", gko_flux), ("reaction_ko", rko_flux)):
                if fl is None:
                    continue
                d = flux_diff(wt_flux, fl, exclude=R)
                is_ex = d.index.isin(exch)
                for part, top in ((d[~is_ex], TOP_REACTIONS), (d[is_ex], TOP_EXCHANGES)):
                    for rank, (rid, r) in enumerate(part.head(top).iterrows(), 1):
                        rx = m.reactions.get_by_id(rid)
                        rer.append({"medium_sweep": job["key"], "entrez": entrez, "ko": kind,
                                    "type": "exchange" if rid in exch else "reaction", "rank": rank,
                                    "reaction": rid, "name": rx.name, "subsystem": rx.subsystem,
                                    "v_wt": r.v_wt, "v_ko": r.v_ko, "delta": r.delta})
    return {"key": job["key"], "rows": rows, "reroutes": rer, "bypass_tests": ex_rows, "wall_s": time.time() - t0}


# ============================================================================= main
def load_inputs() -> dict:
    from nextgem import w1

    worst = pd.read_csv(WORST, sep="\t")
    lf = w1.line_features()
    la = pd.read_csv(LINE_ASSIGNMENT, sep="\t")
    la = la[(la.parent == PARENT) & (la.step == STEP)].set_index("ModelID")
    gm = pd.read_csv(GENE_MAP, sep="\t", dtype={"entrez": str})
    gm = gm.dropna(subset=["entrez"])
    e2g = gm.groupby(gm.entrez.astype(int)).gene.apply(list).to_dict()
    g2sym = dict(zip(gm.gene, gm.symbol))
    g2ent = dict(zip(gm.gene, gm.entrez.astype(int)))
    return {"worst": worst, "lf": lf, "la": la, "e2g": e2g, "g2sym": g2sym, "g2ent": g2ent}


def cell_lines(lf: pd.DataFrame, lineage: str) -> list[str]:
    s = lf[(lf.tier == "test") & (lf.lineage == lineage) & lf.has_expression_matrix.astype(bool)]
    return sorted(s.index)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    if a.workers > MAX_WORKERS:
        raise SystemExit(f"at most {MAX_WORKERS} worker processes")
    from nextgem.run import provenance

    inputs = [WORST, FAILURE_ACCOUNT, LINE_ASSIGNMENT, SWEEP_E, ANNOTATIONS, GENE_MAP,
              PROJECT / "results" / "w3" / "structure" / f"variants_{PARENT}.json",
              PROJECT / "results" / "w1" / "line_features.tsv", PROJECT / "results" / "w1" / "gene_features.tsv"]
    OUT.mkdir(parents=True, exist_ok=True)
    with provenance(OUT / "MANIFEST.tsv", inputs=inputs, solver="glpk", workers=a.workers,
                    allow_dirty=a.allow_dirty, author=AUTHOR) as prov:
        outs, summ = run(a)
        prov.output(*outs)
        prov.note(**summ)
    return 0


def run(a):
    import m6_common as mc
    import m6_sweeps as ms
    from nextgem import w1

    t_start = time.time()
    I = load_inputs()
    worst, la = I["worst"], I["la"]
    # lines and medium signatures per lineage
    lin_lines = {lin: cell_lines(I["lf"], lin) for lin in worst.lineage.unique()}
    sig_lines = {}
    for lin, lines in lin_lines.items():
        sid = la.w2_sweep_id.reindex(lines)
        if sid.isna().any():
            raise RuntimeError(f"{lin}: {int(sid.isna().sum())} lines without an {STEP} sweep")
        sig_lines[lin] = sid.value_counts()
    sids = sorted({s for v in sig_lines.values() for s in v.index})
    bounds, info = ms.bounds_by_id({PARENT: {STEP: sids}})
    vj = ms.variant_json(PARENT, mc.BLOO)
    genes = []
    for r in worst.itertuples():
        gl = I["e2g"].get(int(r.entrez), [])
        if len(gl) > 1:
            raise RuntimeError(f"{r.entrez}: maps to {len(gl)} Human1 genes")
        genes.append((int(r.entrez), gl[0] if gl else None, bool(r.n_fp > 0)))
    genes = list(dict.fromkeys(genes))
    jobs = [{"key": mc.sweep_key(mc.BLOO, s), "w2_sweep_id": s, "bounds": bounds[s], "variant_json": vj,
             "genes": genes} for s in sids]
    print(f"{len(worst)} cells, {len(genes)} genes, {len(sids)} medium signatures; {a.workers} workers", flush=True)
    res = []
    ctx = multiprocessing.get_context("spawn")
    with cf.ProcessPoolExecutor(max_workers=min(a.workers, len(jobs)), mp_context=ctx, initializer=_init) as pool:
        for r in pool.map(trace_medium, jobs, chunksize=1):
            print(f"  {r['key']}: {r['wall_s']:.0f} s", flush=True)
            res.append(r)
    long = pd.DataFrame([x for r in res for x in r["rows"]])
    rer = pd.DataFrame([x for r in res for x in r["reroutes"]])
    byp = pd.DataFrame([x for r in res for x in r["bypass_tests"]])
    long["w2_sweep_id"] = long.medium_sweep.str.split("|").str[1]
    long["medium_keys"] = long.w2_sweep_id.map(lambda s: _mk(info, s))
    # cross-check gene-KO ratio against P1's sweep_e (e = 1 - ratio_r)
    se = pd.read_parquet(SWEEP_E, columns=["key", "entrez", "e"])
    se = se[se.key.isin(set(long.medium_sweep))].set_index(["key", "entrez"]).e
    long["p1_ratio_r"] = [1.0 - se.get((k, e), np.nan) for k, e in zip(long.medium_sweep, long.entrez)]
    inm = long.in_model.astype(bool)
    dev = (long.gene_ko_ratio[inm] - long.p1_ratio_r[inm]).abs()
    check_ok = bool((dev <= 1e-6).all())
    # cell table
    ann = pd.read_parquet(ANNOTATIONS).set_index("entrez")
    gf = w1.gene_features()
    omics = w1.load_omics()
    expr, cn = omics["expr"], omics["cn"]
    cells = []
    for r in worst.itertuples():
        lines = lin_lines[r.lineage]
        w = sig_lines[r.lineage]
        sub = long[long.entrez == r.entrez].set_index("w2_sweep_id").reindex(w.index)
        g = I["e2g"].get(int(r.entrez), [None])[0]
        ex = expr.reindex(index=lines)
        cnl = cn.reindex(index=lines)
        row = {"lineage": r.lineage, "entrez": int(r.entrez),
               "symbol": I["g2sym"].get(g) if g else _symbol_from_features(int(r.entrez)),
               "n_lines": int(r.n_lines), "n_lines_traced": len(lines), "n_fn": int(r.n_fn), "n_fp": int(r.n_fp),
               "error_type": "FN" if r.n_fn >= r.n_fp else "FP", "ensembl": g, "in_human1": g is not None,
               "addressable_human1": bool(gf.addressable_human1.get(int(r.entrez), False)),
               "n_media": len(w), "media": " | ".join(f"{_mk(info, s)}={n}" for s, n in w.items())}
        if g:
            row.update(_gpr_cols(g, I))
            row.update(gene_ko_ratio_min=float(sub.gene_ko_ratio.min()), gene_ko_ratio_max=float(sub.gene_ko_ratio.max()),
                       p1_ratio_r_min=float(sub.p1_ratio_r.min()),
                       rxn_ko_ratio_min=float(sub.rxn_ko_ratio.min()), rxn_ko_ratio_max=float(sub.rxn_ko_ratio.max()),
                       wt_pfba_flux_on_reactions_max=float(sub.wt_pfba_flux_on_reactions.max()),
                       blocked_all_media=bool(sub.blocked.astype(bool).all()),
                       required_at_opt_any=bool(sub.required_at_opt.astype(bool).any()),
                       bypass_exchanges=";".join(sorted({x for s in sub.bypass_exchanges.fillna("") for x in s.split(";") if x})))
        row["mechanism_class"] = modal_class(sub.mechanism_class, pd.Series(w.to_numpy(), index=w.index))
        row["mechanism_class_by_medium"] = " | ".join(f"{_mk(info, s)}:{c}" for s, c in sub.mechanism_class.items())
        row["class_constant_across_media"] = bool(sub.mechanism_class.nunique() == 1)
        # failure-account classes (label-free ones; common_essential needs labels and is not reproduced)
        row["fa_missing_from_model"] = g is None
        for c in FA_GENE_CLASSES:
            row[f"fa_{c}"] = bool(ann[c].get(int(r.entrez), False))
        col = int(r.entrez)
        row["frac_lines_unexpressed"] = float((ex[col] < EXPRESSED).mean()) if col in ex.columns else np.nan
        row["frac_lines_cn_amplified"] = float((cnl[col] >= CN_AMP).mean()) if col in cnl.columns else np.nan
        row["fa_unexpressed"] = bool(row["frac_lines_unexpressed"] >= 0.5) if np.isfinite(row["frac_lines_unexpressed"]) else False
        row["fa_cn_amplified"] = bool(row["frac_lines_cn_amplified"] >= 0.5) if np.isfinite(row["frac_lines_cn_amplified"]) else False
        row["median_expr"] = float(ex[col].median()) if col in ex.columns else np.nan
        if g:
            alts = [x for x in row["alternatives"].split(";") if x]
            ae = {s: float(ex[I["sym2ent"][s]].median()) for s in alts if I["sym2ent"].get(s) in ex.columns}
            row["alternatives_median_expr"] = ";".join(f"{s}={v:.2f}" for s, v in sorted(ae.items(), key=lambda t: -t[1])[:8])
            row["n_alternatives_expressed"] = int(sum(v >= EXPRESSED for v in ae.values()))
        cells.append(row)
    cells = pd.DataFrame(cells)
    cells["mechanism_note"] = [_note(r) for r in cells.itertuples()]
    fa_cols = ["fa_missing_from_model", "fa_unexpressed", "fa_cn_amplified"] + [f"fa_{c}" for c in FA_GENE_CLASSES]
    ct = crosstab(cells, fa_cols)
    fa = pd.read_csv(FAILURE_ACCOUNT, sep="\t").set_index("class")
    ctl = []
    for c in ["missing_from_model", "unexpressed", "cn_amplified", *FA_GENE_CLASSES, "common_essential"]:
        col = f"fa_{c}"
        ctl.append({"failure_account_class": c,
                    "n_worst_cells": int(cells[col].sum()) if col in cells else np.nan,
                    "worst_cells_mechanisms": ";".join(f"{k}={v}" for k, v in cells.loc[cells[col], "mechanism_class"]
                                                       .value_counts().items()) if col in cells else "not reproduced (label-derived)",
                    "frac_of_all_fn_T1": float(fa.loc[c, "frac_of_fn"]) if c in fa.index else np.nan})
    ctl = pd.DataFrame(ctl)
    long_out = long.drop(columns=[c for c in ("fp_blocked_components",) if c in long and long[c].isna().all()])
    outs = [mc.atomic_write_tsv(cells, OUT / "cells.tsv"), mc.atomic_write_tsv(long_out, OUT / "cells_by_medium.tsv"),
            mc.atomic_write_tsv(rer, OUT / "reroutes.tsv"), mc.atomic_write_tsv(byp, OUT / "bypass_tests.tsv"),
            mc.atomic_write_tsv(ct, OUT / "crosstab.tsv"), mc.atomic_write_tsv(ctl, OUT / "crosstab_long.tsv")]
    summ = {"n_cells": int(len(cells)), "n_fn_cells": int((cells.error_type == "FN").sum()),
            "n_fp_cells": int((cells.error_type == "FP").sum()), "n_media": len(sids),
            "classes": cells.mechanism_class.value_counts().to_dict(),
            "gene_ko_ratio_matches_p1_sweep_e": check_ok, "max_abs_dev_vs_p1": float(dev.max()) if len(dev) else 0.0,
            "lines_traced_equal_n_lines": bool((cells.n_lines_traced == cells.n_lines).all()),
            "wall_s": round(time.time() - t_start, 1), "solver": "glpk", "workers": a.workers}
    outs.append(mc.atomic_write_text(json.dumps(summ, indent=1) + "\n", OUT / "f3_summary.json"))
    print(json.dumps(summ, indent=1), flush=True)
    if not check_ok:
        raise RuntimeError("recomputed gene-KO ratios differ from P1's sweep_e")
    return outs, summ


def _mk(info: dict, sid: str) -> str:
    """Medium key(s) of a bounds signature, ',' between keys (';' and '|' are field separators)."""
    return info[sid]["medium_keys"].replace(";", ",")


def _symbol_from_features(entrez: int) -> str:
    from nextgem import paths

    f = paths.baseline_dir() / "human2_gene_mapping.tsv"
    t = pd.read_csv(f, sep="\t", dtype={"entrez": str}).dropna(subset=["entrez"])
    hit = t[t.entrez.astype(int) == entrez]
    return str(hit.symbol.iloc[0]) if len(hit) else ""


_MODEL_GPR: dict = {}


def _gpr_cols(g: str, I: dict) -> dict:
    """GPR description of gene ``g`` from the SBML (cached reaction rules, no solver)."""
    if not _MODEL_GPR:
        import cobra

        cobra.Configuration().solver = "glpk"
        from nextgem.model import load_model

        m = load_model(PARENT)
        for r in m.reactions:
            _MODEL_GPR[r.id] = (r.gene_reaction_rule, r.name, r.subsystem, [x.id for x in r.genes])
        I["sym2ent"] = {s: e for gg, s in I["g2sym"].items() for e in [I["g2ent"].get(gg)] if e is not None}
    rx = sorted(r for r, v in _MODEL_GPR.items() if g in v[3])
    parts, kinds, alts, closed = [], [], set(), []
    for r in rx:
        s = gpr_structure(_MODEL_GPR[r][0], g)
        kinds.append(s["kind"])
        alts.update(s["alternatives"])
        if s["closed_by_ko"]:
            closed.append(r)
        parts.append(f"{r}({_MODEL_GPR[r][2]}; {_MODEL_GPR[r][1] or '-'}): {s['kind']}, {s['n_genes']} genes, "
                     f"{s['n_branches']} branches{', subunit' if s['subunit'] else ''}")
    return {"n_reactions": len(rx), "reactions": ";".join(rx), "gpr_kinds": ";".join(sorted(set(kinds))),
            "gpr_detail": " | ".join(parts), "closed_by_gene_ko": ";".join(closed),
            "alternatives": ";".join(sorted(I["g2sym"].get(x, x) for x in alts))}


def _note(r) -> str:
    if not r.in_human1:
        return "Not in Human1: P1 cannot score the gene (e = 0 by the W1 convention)."
    alt = f"{len(r.alternatives.split(';'))} alternative genes" if isinstance(r.alternatives, str) and r.alternatives else "no alternative gene"
    s = (f"Gene KO ratio {r.gene_ko_ratio_min:.3g}-{r.gene_ko_ratio_max:.3g}; KO closes "
         f"{len([x for x in str(r.closed_by_gene_ko).split(';') if x and x != 'nan'])}/{r.n_reactions} reactions "
         f"({alt} in the GPR); all-reactions KO ratio {r.rxn_ko_ratio_min:.3g}-{r.rxn_ko_ratio_max:.3g}; WT pFBA flux "
         f"on its reactions {r.wt_pfba_flux_on_reactions_max:.3g}{'; reactions blocked on every medium' if r.blocked_all_media else ''}.")
    if r.bypass_exchanges:
        s += f" Bypass via {r.bypass_exchanges}."
    return s


if __name__ == "__main__":
    raise SystemExit(main())
