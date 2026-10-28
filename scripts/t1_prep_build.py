"""T1 prep (2), (3), (5), (6): P1's label-free input files for scripts/evaluate_test.py (DESIGN §6.8, §10.1-10.3).
No label file is opened (split, line features, scores, sweeps, task matrices, panels, paralogs and CORUM only).

(2) Scores, layout line_gene, exactly the tier's expression lines x U, columns H1, H2, P1[H1], P1[H2], M2[H1]:
    taken from the M6 score files (results/m6/scores/scores_<tier>.large.parquet, SHA-256 checked against
    results/m6/scores/SHA256SUMS.tsv), renamed from the M6 cell names, and re-derived independently from the tracked
    sweep_e.parquet + line_assignment.tsv (must agree bitwise). H1 / H2 must equal the W1 generic e (as the harness
    checks). validation -> results/m6/t1/p1_scores_validation.parquet (tracked; also checked equal to the M6 look's
    results/m6/val/m6_val_scores.parquet); test -> results/m6/tmp/t1/p1_scores_test.parquet (git-ignored).
(3) §6.8 calls of the arms P1[H1], P1[H2], long (ModelID, entrez, method, p57, pgr, fba), layout line_gene:
    p57 / pgr from the B-loo Essential-task matrix (scripts/t1_prep_tasks.py; paper protocol, GR = growth on Ham's,
    so line-constant per parent); fba = ratio_r < 0.01 of the arm's line score (line-specific: M2(f=1) per line).
    Genes of U absent from the parent are False. Cross-checks: the M6 construction (generic non-GR task failures
    from results/w3 ∪ GR from the B-loo|M0 sweep, KO growth < 1) gives the same p57 / pgr; the fba rule reproduces
    W1's generic fba_call_H1 / H2 from e_H1 / e_H2. validation -> results/m6/t1/, test -> results/m6/tmp/t1/.
(5) Per-task calls of the primary arm's model (H1 + B-loo): results/m6/t1/p1_task_calls.parquet (entrez, task, fails).
(6) Failure-account gene annotations over U (label-free; results/m6/t1/p1_gene_annotations.parquet):
    paralog_buffered      : Human1 gene whose single KO closes no reaction (an OR alternative remains on every
                            reaction it catalyses) and which shares a reaction (hence a compartment) with an
                            Ensembl 116 within-species paralog (data/paralogs/ensembl_human_paralogs.tsv)
    medium_sensitive      : a-priori panel §1.2 (medium-dependent essentials) ∪ genes whose P1(H1) FBA<0.01 call
                            differs between any two H1 + B-loo sweeps (M0 and the 20 M2(f=1) signatures)
    complex_subunit       : a subunit (by gene name) of a CORUM 5.1 human complex with >= 2 subunits
    housekeeping_adjacent : a-priori panel §1.6g
    biomass_artefact_lethal: a-priori panel §1.4a-c (cofactor/vitamin and small-metabolite pool dead-ends, glycogen)
    expected_fba_failure  : the pre-registered list of DESIGN §10.3: MAT2A, PRMT5 (MTAP collateral), V-ATPase
                            (ATP6V* by family), GPX4, PGD
    Not annotated (column omitted, so the harness reports them as such): gene_less_bypass (no label-free rule fixed
    in advance) and model_contradicted_loss (a context-family repair log; P1 has no context family).

Also writes results/m6/t1/build_summary.json (checks, counts, input SHA-256 including the data/ inputs that the
P1 spec's provenance cannot pin). Usage: .venv/bin/python scripts/t1_prep_build.py
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as mc  # noqa: E402
import t1_prep_common as tc  # noqa: E402

SC = mc.SCDIR
W3 = tc.PROJECT / "results" / "w3"
EXPECTED_FBA_FAILURE_SYMBOLS = ("MAT2A", "PRMT5", "GPX4", "PGD")
EXPECTED_FBA_FAILURE_FAMILY = r"^ATP6V\w+$"
OUT = {
    "scores_validation": tc.T1 / "p1_scores_validation.parquet",
    "scores_test": tc.T1_TMP / "p1_scores_test.parquet",
    "calls_validation": tc.T1 / "p1_calls_validation.parquet",
    "calls_test": tc.T1_TMP / "p1_calls_test.parquet",
    "task_calls": tc.T1 / "p1_task_calls.parquet",
    "gene_annotations": tc.T1 / "p1_gene_annotations.parquet",
    "summary": tc.T1 / "build_summary.json",
}
G: dict = {}


# ============================================================================= context
def load_context() -> None:
    from nextgem import labels, paths, w1

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    U = list(labels.gene_universe().entrez.astype(int))
    gen = w1.generic_table(pd.Series(U))
    emap, sym = {}, {}
    for p in mc.PARENTS:
        b = pd.read_parquet(paths.baseline_dir() / f"{p}_ko.parquet", columns=["gene", "entrez", "symbol"])
        b = b.dropna(subset=["entrez"])
        emap[p] = pd.Series(b.entrez.astype(int).to_numpy(), index=b.gene.to_numpy())
        sym.update(dict(zip(b.symbol, b.entrez.astype(int))))
    G.update(split=split, lf=lf, U=U, gen=gen, emap=emap, sym=sym)


def eval_lines(tier: str) -> list[str]:
    s, lf = G["split"], G["lf"]
    tl = sorted(s.index[s.tier == tier])
    return [l for l in tl if l in lf.index and bool(lf.loc[l, "has_expression_matrix"])]


def m6_sha_pins() -> dict:
    t = pd.read_csv(SC / "SHA256SUMS.tsv", sep="\t")
    return dict(zip(t.file, t.sha256))


# ============================================================================= (2) scores
def rederive_from_sweeps(lines, method: str) -> pd.DataFrame:
    """lines x U e-matrix of an M6 cell from the tracked sweep_e + line_assignment (independent of the .large files)."""
    cell = next(c for c in mc.CELLS if c.method == tc.SOURCE_COLUMN[method])
    la = pd.read_csv(SC / "line_assignment.tsv", sep="\t")
    la = la[(la.parent == cell.parent) & (la.step == cell.step)].set_index("ModelID")
    keys = la.w2_sweep_id.reindex(list(lines)).map(lambda s: mc.sweep_key(cell.variant, s))
    if keys.isna().any():
        raise RuntimeError(f"{method}: lines without a sweep")
    se = G.setdefault("sweep_e", pd.read_parquet(SC / "sweep_e.parquet", columns=["key", "entrez", "e"]))
    E = se[se.key.isin(set(keys))].pivot(index="key", columns="entrez", values="e").reindex(columns=G["U"])
    M = E.loc[keys.to_numpy()]
    M.index = list(lines)
    return M


def build_scores(tier: str, checks: list) -> pd.DataFrame:
    src = SC / f"scores_{tier}.large.parquet"
    pin = m6_sha_pins()[tc.rel(src)]
    if tc.sha256(src) != pin:
        raise RuntimeError(f"{src.name}: SHA-256 differs from results/m6/scores/SHA256SUMS.tsv")
    lines = eval_lines(tier)
    df = tc.rename_scores(pd.read_parquet(src), lines)
    if set(df.ModelID) != set(lines) or set(df.entrez) != set(G["U"]) or len(df) != len(lines) * len(G["U"]):
        raise RuntimeError(f"{tier}: scores do not cover exactly the {len(lines)} expression lines x U")
    for m in tc.METHODS:
        if not np.isfinite(df[m].to_numpy()).all():
            raise RuntimeError(f"{tier}: non-finite {m} (a failed build would need the failed_builds file)")
    for s in ("H1", "H2"):
        v = G["gen"][f"e_{s}"].reindex(G["U"]).to_numpy()
        W = df.pivot(index="ModelID", columns="entrez", values=s).reindex(index=lines, columns=G["U"]).to_numpy()
        checks.append({"check": f"{tier}: {s} equals the W1 generic e_{s} (atol 1e-9)",
                       "pass": bool(np.allclose(W, np.tile(v, (len(lines), 1)), rtol=0, atol=1e-9))})
    for m in (tc.P1_H1, tc.P1_H2, tc.MED_H1):
        W = df.pivot(index="ModelID", columns="entrez", values=m).reindex(index=lines, columns=G["U"])
        R = rederive_from_sweeps(lines, m)
        checks.append({"check": f"{tier}: {m} equals its re-derivation from sweep_e + line_assignment (bitwise)",
                       "pass": bool(np.array_equal(W.to_numpy(), R.to_numpy()))})
    if tier == "validation":
        v = pd.read_parquet(mc.VALDIR / "m6_val_scores.parquet")
        v = tc.rename_scores(v, lines)
        checks.append({"check": "validation: equals the M6 look's m6_val_scores.parquet (renamed)",
                       "pass": bool(v.equals(df))})
    checks.append({"check": f"{tier}: n_lines = {len(lines)}, n_genes = {len(G['U'])}", "pass": True})
    return df


# ============================================================================= (3) calls
def task_matrix(parent: str) -> tuple[pd.DataFrame, str]:
    long = pd.read_parquet(tc.T1 / f"tasks_bloo_{parent}.parquet")
    meta = pd.read_csv(tc.T1 / f"tasks_bloo_{parent}_tasks.tsv", sep="\t")
    gr = meta.task[meta.category == "GR"]
    if len(gr) != 1 or meta.task.nunique() != 57 or (meta.status != "ok").any():
        raise RuntimeError(f"tasks_bloo_{parent}: expected 57 ok tasks with one GR task")
    return long, gr.iloc[0]


def gene_calls(arm: str, checks: list) -> pd.DataFrame:
    """Line-constant p57 / pgr of an arm over U, and the M6-construction cross-check."""
    p = tc.ARM_PARENT[arm]
    long, gr = task_matrix(p)
    gc = tc.p57_pgr(long, gr, G["U"], G["emap"][p])
    # M6 construction: generic non-GR task failures (W3) ∪ GR from the B-loo|M0 biomass sweep (KO growth < 1)
    gl = pd.read_parquet(W3 / "tasks" / f"essential_{p}.parquet", columns=["task", "seed", "gene", "feasible"])
    gl = gl[gl.seed == -1]
    gm = pd.read_csv(W3 / "tasks" / f"essential_{p}_tasks.tsv", sep="\t")
    ggr = gm.task[gm.category == "GR"].iloc[0]
    nongr = tc.p57_pgr(gl[gl.task != ggr], "__none__", G["U"], G["emap"][p]).p57
    cat = pd.read_csv(mc.SWDIR / "catalog.tsv", sep="\t")
    m0 = cat[(cat.parent == p) & (cat.variant == mc.BLOO) & (cat.step == mc.M0)]
    if len(m0) != 1:
        raise RuntimeError(f"{p}: expected one B-loo|M0 sweep in the M6 catalog")
    se = G.setdefault("sweep_e", pd.read_parquet(SC / "sweep_e.parquet", columns=["key", "entrez", "e"]))
    e0 = se[se.key == m0.key.iloc[0]].set_index("entrez").e.reindex(G["U"])
    in_p = G["gen"][f"in_{mc.SHORT[p]}"].reindex(G["U"]).to_numpy(bool)
    gr_m6 = pd.Series(mc.gr_fail(1.0 - e0, float(m0.wt_growth.iloc[0])).to_numpy(bool) & in_p, index=gc.index)
    checks.append({"check": f"{arm}: PGR from the B-loo task run equals the M6 construction (B-loo|M0 sweep, KO "
                            f"growth < 1)", "pass": bool((gc.pgr == gr_m6).all())})
    checks.append({"check": f"{arm}: P57 from the B-loo task run equals the M6 construction (generic non-GR ∪ GR)",
                   "pass": bool((gc.p57 == (nongr | gr_m6)).all())})
    s = mc.SHORT[p]
    checks.append({"check": f"{arm}: absent-from-parent genes are False", "pass": bool(not (gc.p57 & ~in_p).any())})
    G.setdefault("counts", {})[arm] = {"n_p57": int(gc.p57.sum()), "n_pgr": int(gc.pgr.sum()),
                                       "n_p57_generic": int(G["gen"][f"p57_{s}"].sum()),
                                       "n_pgr_generic": int(G["gen"][f"pgr_{s}"].sum())}
    return gc


def build_calls(tier: str, scores: pd.DataFrame, gcs: dict, checks: list) -> pd.DataFrame:
    lines = sorted(scores.ModelID.unique())
    parts = []
    for arm in tc.ARMS:
        E = scores.pivot(index="ModelID", columns="entrez", values=arm).reindex(index=lines, columns=G["U"])
        parts.append(tc.calls_long(arm, gcs[arm], E))
        G["counts"][arm][f"n_fba_pairs_{tier}"] = int(parts[-1].fba.sum())
        G["counts"][arm][f"n_ratio_exactly_0.01_{tier}"] = int((np.round(1 - E.to_numpy(), 6) == 0.01).sum())
    df = pd.concat(parts, ignore_index=True)
    checks.append({"check": f"{tier}: calls cover exactly lines x U per arm",
                   "pass": bool(all(len(df[df.method == a]) == len(lines) * len(G["U"]) for a in tc.ARMS))})
    return df


def fba_rule_check(checks: list) -> None:
    for s in ("H1", "H2"):
        e = G["gen"][f"e_{s}"].reindex(G["U"]).to_numpy()
        checks.append({"check": f"fba rule on e_{s} reproduces W1 fba_call_{s}",
                       "pass": bool((tc.fba_call_from_e(e) == G["gen"][f"fba_call_{s}"].reindex(G["U"]).to_numpy(bool)).all())})


# ============================================================================= (6) annotations
def _panel(section: str, sub_prefixes=None) -> set[int]:
    from nextgem import w1

    p = w1.gene_panels()
    p = p[p.section.astype(str) == section]
    if sub_prefixes:
        p = p[p.subpanel.astype(str).str[:2].isin([f"{x}:" for x in sub_prefixes])]
    return set(p.entrez.dropna().astype(int))


def paralog_buffered_h1() -> set[int]:
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import paths, w3
    from nextgem.model import load_model

    m = load_model("human1")
    gd = w3.gene_disabled_reactions(m)
    model_genes = {g.id for g in m.genes}
    par = pd.read_csv(paths.data_root() / "paralogs" / "ensembl_human_paralogs.tsv", sep="\t",
                      usecols=["Gene stable ID", "Human paralogue gene stable ID"]).dropna()
    par = par[par["Gene stable ID"].isin(model_genes) & par["Human paralogue gene stable ID"].isin(model_genes)]
    pmap = par.groupby("Gene stable ID")["Human paralogue gene stable ID"].apply(set).to_dict()
    out = set()
    for g in m.genes:
        if gd.get(g.id) or not g.reactions:  # the KO closes a reaction: not buffered
            continue
        partners = set().union(*[{x.id for x in r.genes} for r in g.reactions]) - {g.id}
        if partners & pmap.get(g.id, set()):
            out.add(g.id)
    em = G["emap"]["human1"]
    return {int(em[g]) for g in out if g in em.index}


def medium_sensitive_sweeps() -> set[int]:
    se = G.setdefault("sweep_e", pd.read_parquet(SC / "sweep_e.parquet", columns=["key", "entrez", "e"]))
    cat = pd.read_csv(mc.SWDIR / "catalog.tsv", sep="\t")
    keys = set(cat.key[(cat.parent == "human1") & (cat.variant == mc.BLOO)])
    W = se[se.key.isin(keys)].pivot(index="key", columns="entrez", values="e")
    C = pd.DataFrame(tc.fba_call_from_e(W.to_numpy()), index=W.index, columns=W.columns)
    G["n_bloo_h1_sweeps"] = len(W)
    return {int(g) for g in C.columns[C.any() & ~C.all()]}


def complex_subunits() -> set[int]:
    from nextgem import paths

    c = pd.read_csv(paths.data_root() / "complexes" / "corum_humanComplexes.txt", sep="\t",
                    usecols=["organism", "subunits_gene_name"], dtype=str)
    c = c[c.organism == "Human"].dropna(subset=["subunits_gene_name"])
    names = set()
    for s in c.subunits_gene_name:
        sub = {x.strip() for x in s.split(";") if x.strip()}
        if len(sub) >= 2:
            names |= sub
    return {G["sym"][n] for n in names if n in G["sym"]}


def expected_fba_failure() -> set[int]:
    fam = re.compile(EXPECTED_FBA_FAILURE_FAMILY)
    return {e for s, e in G["sym"].items() if isinstance(s, str) and (s in EXPECTED_FBA_FAILURE_SYMBOLS or fam.match(s))}


def build_annotations() -> tuple[pd.DataFrame, dict]:
    med_panel, med_sweep = _panel("1.2"), medium_sensitive_sweeps()
    sets = {
        "paralog_buffered": paralog_buffered_h1(),
        "medium_sensitive": med_panel | med_sweep,
        "complex_subunit": complex_subunits(),
        "housekeeping_adjacent": _panel("1.6", ["g"]),
        "biomass_artefact_lethal": _panel("1.4", ["a", "b", "c"]),
        "expected_fba_failure": expected_fba_failure(),
    }
    U = pd.Index(G["U"], name="entrez")
    df = pd.DataFrame({"entrez": U.to_numpy(dtype="int64")})
    for k, s in sets.items():
        df[k] = U.isin(list(s))
    info = {k: int(df[k].sum()) for k in sets}
    info.update(medium_sensitive_panel_in_U=int(U.isin(list(med_panel)).sum()),
                medium_sensitive_sweeps_in_U=int(U.isin(list(med_sweep)).sum()),
                n_bloo_h1_sweeps=G["n_bloo_h1_sweeps"],
                not_annotated={"gene_less_bypass": "no label-free rule fixed in advance",
                               "model_contradicted_loss": "context-family repair log; P1 has no context family"})
    return df, info


# ============================================================================= main
def inputs() -> list[Path]:
    from nextgem import paths

    d = paths.data_root()
    return ([SC / "SHA256SUMS.tsv", SC / "sweep_e.parquet", SC / "line_assignment.tsv", mc.SWDIR / "catalog.tsv",
             mc.VALDIR / "m6_val_scores.parquet"] + [SC / f"scores_{t}.large.parquet" for t in ("validation", "test")]
            + [tc.T1 / f"tasks_bloo_{p}{s}" for p in mc.PARENTS for s in (".parquet", "_tasks.tsv")]
            + [W3 / "tasks" / f"essential_{p}{s}" for p in mc.PARENTS for s in (".parquet", "_tasks.tsv")]
            + [paths.baseline_dir() / f"{p}_ko.parquet" for p in mc.PARENTS]
            + [d / "controls" / "gene_panels.tsv", d / "paralogs" / "ensembl_human_paralogs.tsv",
               d / "complexes" / "corum_humanComplexes.txt"])


def main(argv=None) -> int:
    t0 = time.time()
    ins = inputs()
    with tc.provenance(inputs=ins) as prov:
        load_context()
        checks: list = []
        fba_rule_check(checks)
        gcs = {arm: gene_calls(arm, checks) for arm in tc.ARMS}
        written = []
        for tier in ("validation", "test"):
            sc = build_scores(tier, checks)
            written.append(tc.atomic_write_parquet(sc, OUT[f"scores_{tier}"]))
            ca = build_calls(tier, sc, gcs, checks)
            written.append(tc.atomic_write_parquet(ca, OUT[f"calls_{tier}"]))
            del sc, ca
        long, _ = task_matrix("human1")
        written.append(tc.atomic_write_parquet(tc.task_calls_table(long, G["emap"]["human1"]), OUT["task_calls"]))
        ann, ann_info = build_annotations()
        written.append(tc.atomic_write_parquet(ann, OUT["gene_annotations"]))
        failed = [c for c in checks if not c["pass"]]
        summary = {"checks": checks, "all_checks_pass": not failed, "counts": G["counts"], "annotations": ann_info,
                   "n_eval_lines": {t: len(eval_lines(t)) for t in ("validation", "test")}, "n_genes_U": len(G["U"]),
                   "outputs": {tc.rel(p): tc.sha256(p) for p in written},
                   "inputs_sha256": {str(p.relative_to(tc.PROJECT)) if str(p).startswith(str(tc.PROJECT)) else str(p):
                                     tc.sha256(p) for p in ins},
                   "wall_s": round(time.time() - t0, 1)}
        tc.atomic_write_json(summary, OUT["summary"])
        prov.output(*written, OUT["summary"])
        prov.note(all_checks_pass=not failed)
    print(json.dumps({k: summary[k] for k in ("all_checks_pass", "counts", "annotations", "n_eval_lines", "outputs")},
                     indent=1))
    if failed:
        print(f"FAILED CHECKS: {failed}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
