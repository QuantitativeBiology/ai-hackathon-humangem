"""W4 toy tests with known answers (DESIGN.md §12 row W4): one toy network per family (F0, GIMME,
FASTCORE, E-Flux, CORDA), the repair rule, the KO gene-state convention, the loss-rule drop under
shuffled input order, and the KO-sweep zero-flux shortcut against brute force. Every test runs on
GLPK and on Gurobi (GLPK = Gurobi per family on toys). Label-free: no DepMap file is read.
"""

from __future__ import annotations

import math

import cobra
import numpy as np
import pandas as pd
import pytest

from nextgem import context as cx
from nextgem.model import BIOMASS_REACTION as BIO

SOLVERS = ["glpk", "gurobi"]


def make_toy(spec, solver, objective=BIO):
    """spec: (id, equation, gpr, lb, ub)."""
    m = cobra.Model("toy")
    m.solver = solver
    for rid, eq, gpr, lb, ub in spec:
        r = cobra.Reaction(rid)
        m.add_reactions([r])
        r.build_reaction_from_string(eq, verbose=False)
        r.gene_reaction_rule = gpr
        r.bounds = (lb, ub)
    m.objective = objective
    return cx.parent_from_model(m)


def scores(**kw):
    return pd.Series(kw, dtype="float64")


def build(parent, family, s, lost=None):
    lm = cx.extract(parent, family, s.reindex(parent.genes))
    lm.line = "toy"
    return cx.finish(parent, lm, lost)


def brute_force(parent, lm):
    """KO ratio of every gene by solving each KO (no shortcut)."""
    close, rebound = cx.closure_map(parent, lm)
    m = parent.model
    out = {}
    with m:
        cx._apply_bounds(m, lm.bounds())
        wt = m.slim_optimize()
        for g in parent.genes:
            with m:
                cx._apply_bounds(m, {**rebound.get(g, {}), **{r: (0.0, 0.0) for r in close.get(g, set())}})
                v = m.slim_optimize(error_value=0.0)
            out[g] = round(min(max(v / wt, 0.0), 1.0), 6)
    return pd.Series(out)


# --------------------------------------------------------------------------- toy networks
TOY1 = [  # two routes glc -> b; F0 closes both
    ("EX_glc", "--> glc", "", 0, 10),
    ("R1", "glc --> a", "gA", 0, 1000),
    ("R2", "a --> b", "gB1 or gB2", 0, 1000),
    ("R3", "glc --> b", "gC and gD", 0, 1000),
    ("R4", "b --> p", "gE", 0, 1000),
    (BIO, "p --> ", "", 0, 1000),
]


@pytest.mark.parametrize("solver", SOLVERS)
def test_f0_prunes_and_repair_reopens_cheapest(solver):
    p = make_toy(TOY1, solver)
    s = scores(gA=5, gB1=-2, gB2=-3, gC=2, gD=-0.5, gE=1)
    lm = cx.extract(p, "F0", s.reindex(p.genes))
    assert lm.closed == {"R2", "R3"}  # R2: both isozymes off; R3: gD off
    lm = cx.finish(p, lm)
    # weights: R2 = 1 + 2 = 3, R3 = 1 + 0.5 = 1.5 -> the L1 LP re-opens R3 only
    assert lm.reopened == {"R3"}
    assert lm.info["wt_growth"] == pytest.approx(10.0)
    assert not lm.info["failed"]
    ko = cx.ko_sweep(p, lm).set_index("gene")
    assert ko.loc["gC", "ratio_r"] == 0 and ko.loc["gD", "ratio_r"] == 0  # re-opened R3 keeps its full GPR
    assert ko.loc["gA", "ratio_r"] == 1 and ko.loc["gA", "status"] == "zero_flux"
    assert ko.loc["gB1", "status"] == "no_effect"
    assert ko.loc["gE", "ratio_r"] == 0
    pd.testing.assert_series_equal(ko.ratio_r.sort_index(), brute_force(p, lm).sort_index(), check_names=False)


@pytest.mark.parametrize("solver", SOLVERS)
def test_off_isozyme_does_not_buffer(solver):
    """Gene-state convention: an open reaction's pruned GPR is used, so KO of the only expressed
    isozyme is lethal; the generic model (all genes functional) is buffered."""
    p = make_toy(TOY1, solver)
    s = scores(gA=5, gB1=1, gB2=-1, gC=-1, gD=-1, gE=1)
    lm = build(p, "F0", s)
    assert lm.closed == {"R3"} and lm.reopened == set()
    ko = cx.ko_sweep(p, lm).set_index("gene")
    assert ko.loc["gB1", "ratio_r"] == 0
    assert ko.loc["gB2", "status"] == "no_effect"
    generic = cx.LineModel(line="g", family="F0", closed=set())
    kg = cx.ko_sweep(p, generic).set_index("gene")
    assert kg.loc["gB1", "ratio_r"] == 1  # buffered by gB2 (and R3) in the generic model


TOY_GIMME = [
    ("EX_glc", "--> glc", "", 0, 10),
    ("Rh", "glc --> b", "gH", 0, 1000),
    ("Rl", "glc --> b", "gL", 0, 1000),
    ("Rq", "glc --> q", "gQ", 0, 1000),
    ("Rx", "b --> x", "gX", 0, 1000),
    ("EX_x", "x --> ", "", 0, 1000),
    (BIO, "b + q --> ", "", 0, 1000),
]


@pytest.mark.parametrize("solver", SOLVERS)
def test_gimme_known_answer(solver):
    p = make_toy(TOY_GIMME, solver)
    assert p.wt_generic == pytest.approx(5.0)
    s = scores(gH=2, gL=-1, gQ=-1, gX=-2)
    lm = cx.extract(p, "F2", s.reindex(p.genes))
    # biomass >= 4.5 forces Rq (low, used, kept); b via the expressed Rh; Rl and Rx unused -> closed
    assert lm.closed == {"Rl", "Rx"}
    assert lm.info["gimme_objective"] == pytest.approx(4.5)
    lm = cx.finish(p, lm)
    assert lm.reopened == set() and lm.info["wt_growth"] == pytest.approx(5.0)


TOY_FASTCORE = [
    ("EX_glc", "--> glc", "", 0, 10),
    ("Ra", "glc --> a", "gRa", 0, 1000),
    ("Rb1", "glc --> m1", "gRb1", 0, 1000),
    ("Rb2", "m1 --> a", "gRb2", 0, 1000),
    ("Rc", "a --> c", "gRc", 0, 1000),
    ("Rrev", "d <=> c", "gRr", -1000, 1000),
    ("Rd2", "d --> d2", "", 0, 1000),
    ("EX_d2", "d2 --> ", "", 0, 1000),
    ("Rdead", "a --> z", "gZ", 0, 1000),
    (BIO, "c --> ", "", 0, 1000),
]


@pytest.mark.parametrize("solver", SOLVERS)
def test_fastcore_known_answer(solver):
    p = make_toy(TOY_FASTCORE, solver)
    cons = cx.consistent_reactions(p)
    assert "Rdead" not in cons and len(cons) == len(p.rxn_ids) - 1
    s = scores(gRa=-1, gRb1=-1, gRb2=-1, gRc=1, gRr=1, gZ=-1)
    lm = cx.extract(p, "F3", s.reindex(p.genes))
    support = set(p.rxn_ids) - lm.closed
    # core {Rc, Rrev}; Rrev only in reverse (flip), shortest non-core paths: Ra (not Rb1+Rb2), biomass, Rd2/EX_d2
    assert support == {"EX_glc", "Ra", "Rc", BIO, "Rrev", "Rd2", "EX_d2"}
    assert lm.info["unreachable_core"] == 0


TOY_EFLUX = [
    ("EX_glc", "--> glc", "", 0, 1000),
    ("R1", "glc --> a", "gA", 0, 1000),
    ("R2", "a --> b", "gB1 or gB2", 0, 1000),
    (BIO, "b --> ", "", 0, 1000),
]


@pytest.mark.parametrize("solver", SOLVERS)
def test_eflux_known_answer(solver):
    p = make_toy(TOY_EFLUX, solver)
    # linear values 2**s: gA 16, gB1 8, gB2 4 -> cap R1 16, R2 12 (OR = sum); bounds 1000 and 750
    lm = build(p, "F4", scores(gA=4, gB1=3, gB2=2))
    assert lm.eflux["R2"] == (0, 750.0) and lm.eflux["R1"] == (0, 1000.0)
    assert lm.info["wt_growth"] == pytest.approx(750.0) and lm.reopened == set()
    ko = cx.ko_sweep(p, lm).set_index("gene")
    assert ko.loc["gB1", "ratio_r"] == pytest.approx(round(1 / 3, 6))  # cap 4 -> bound 250
    assert ko.loc["gB2", "ratio_r"] == pytest.approx(round(2 / 3, 6))  # cap 8 -> bound 500
    assert ko.loc["gA", "ratio_r"] == 0
    pd.testing.assert_series_equal(ko.ratio_r.sort_index(), brute_force(p, lm).sort_index(), check_names=False)
    # floor of 10: gA 256, gB1 off (0), gB2 at threshold (1) -> 1000/256 < 10 -> bound 10
    lm2 = build(p, "F4", scores(gA=8, gB1=-1, gB2=0))
    assert lm2.eflux["R2"] == (0, 10.0) and lm2.info["wt_growth"] == pytest.approx(10.0)


TOY_CORDA = [
    ("EX_glc", "--> glc", "", 0, 10),
    ("Rlow", "glc --> a", "gL", 0, 1000),
    ("Rneg", "glc --> a", "gN", 0, 1000),
    ("Rhigh", "a --> b", "gH", 0, 1000),
    ("Runk", "b --> u", "", 0, 1000),
    ("EX_u", "u --> ", "", 0, 1000),
    (BIO, "b --> ", "", 0, 1000),
]


@pytest.mark.parametrize("solver", SOLVERS)
def test_corda_known_answer(solver):
    p = make_toy(TOY_CORDA, solver)
    s = scores(gL=-0.5, gN=-3, gH=2)  # confidences: low 1, negative -1, high 3
    assert cx.corda_confidence_genes(s) == {"gL": 1, "gN": -1, "gH": 3}
    lm = cx.extract(p, "F3b", s.reindex(p.genes))
    included = set(p.rxn_ids) - lm.closed
    assert {"EX_glc", "Rlow", "Rhigh"} <= included  # the low path (penalty 1) beats the negative one (100)
    assert "Rneg" not in included
    assert len(included & {BIO, "Runk"}) >= 1  # exactly which zero-confidence outlet is kept is degenerate


TOY_LOSS = [
    ("EX_glc", "--> glc", "", 0, 10),
    ("Rx", "glc --> a", "gX", 0, 1000),
    ("Ryz", "a --> b", "gY or gZ", 0, 1000),
    ("Rw", "a --> w", "gW", 0, 1000),
    ("EX_w", "w --> ", "", 0, 1000),
    (BIO, "b --> ", "", 0, 1000),
]


@pytest.mark.parametrize("solver", SOLVERS)
def test_loss_rule_order_independent(solver):
    p = make_toy(TOY_LOSS, solver)
    p.entrez = {"gX": "300", "gY": "200", "gZ": "100", "gW": "400"}
    expr = pd.Series({"gX": 3.0, "gY": 5.0, "gZ": 5.0, "gW": 1.0})
    expected = {"lost": ["gW", "gX", "gY", "gZ"], "contradicted": ["gX"], "dropped_joint": ["gZ"],
                "accepted": ["gW", "gY"]}
    rng = np.random.default_rng(11)
    for _ in range(20):
        order = list(rng.permutation(["gX", "gY", "gZ", "gW"]))
        res = cx.accept_losses(p, order, expr.loc[order])
        assert {k: res[k] for k in expected} == expected  # gY/gZ tie on expression -> lower Entrez (gZ) dropped
        assert res["joint_growth"] == pytest.approx(10.0)
    # accepted loss closures are never re-opened by the repair
    lm = build(p, "F0", scores(gX=1, gY=1, gZ=1, gW=1), lost={"gW", "gY", "gZ"})
    assert lm.info["failed"] and "Ryz" in lm.final_closed and "Ryz" not in lm.reopened


def test_presence_score_matches_boolean_pruning():
    rng = np.random.default_rng(23)
    nodes = [("or", (("and", ("a", "b")), "c")), ("and", ("a", ("or", ("b", "c")))), ("or", ("a", "b", "c"))]
    for _ in range(300):
        vals = {g: (np.nan if rng.random() < 0.3 else float(rng.normal())) for g in "abc"}
        off = {g for g, v in vals.items() if not np.isnan(v) and v < 0}
        for n in nodes:
            sc = cx.presence_score(n, {g: v for g, v in vals.items() if not np.isnan(v)})
            assert (not math.isnan(sc) and sc < 0) == (not cx.gpr_true(n, off))


# --------------------------------------------------------------------------- inputs (no solver)
def test_loss_rule_definition_known_answer():
    from nextgem import context_inputs as ci

    genes = ["hom", "het_loh", "het_silent", "het_only", "deep_del", "shallow", "none"]
    dmg = pd.DataFrame([[2, 1, 1, 1, 0, 0, 0]], index=["L"], columns=genes, dtype=float)
    cn = pd.DataFrame([[1.0, 0.4, 1.0, 1.0, 0.2, 0.4, 1.0]], index=["L"], columns=genes)
    rna = pd.DataFrame([[5, 5, 0.5, 5, 5, 5, 5]], index=["L"], columns=genes, dtype=float)
    lost, flags = ci.lost_genes("L", dmg, cn, rna)
    assert lost == sorted(["hom", "het_loh", "het_silent", "deep_del"]) and not flags["cn_missing"]
    lost_v1, _ = ci.lost_genes("L", dmg, cn, rna, rule="v1")
    assert lost_v1 == sorted(["hom", "het_loh", "het_silent", "het_only", "deep_del"])
    lost_05, _ = ci.lost_genes("L", dmg, cn, rna, cn_del=0.5)  # CN sensitivity 0.5 adds the shallow loss
    assert "shallow" in lost_05
    lost_nocn, flags = ci.lost_genes("L", dmg, None, rna)  # mutations only, flagged
    assert flags["cn_missing"] and lost_nocn == sorted(["hom", "het_silent"])


def test_panel_training_only_and_deterministic():
    from nextgem import context_inputs as ci

    rng = np.random.default_rng(5)
    rows = []
    for k, L in enumerate([f"L{i:02d}" for i in range(14)]):
        tier = "train" if k < 12 else "test"
        for j in range(15 + k):
            rows.append((f"ACH-{k:02d}{j:03d}", tier, L, bool(rng.random() > 0.1), True, False, "x"))
    split = pd.DataFrame(rows, columns=["ModelID", "tier", "lineage", "has_expression", "has_cn", "has_procan", "base_medium_class"])
    p1, p2 = ci.stage_a_panel(split), ci.stage_a_panel(split)
    assert p1.equals(p2) and len(p1) == 30 and p1.lineage.nunique() == 10
    assert set(p1.ModelID) <= set(split[(split.tier == "train") & split.has_expression].ModelID)


def test_permutation_preserves_each_gene_distribution():
    from nextgem import context_inputs as ci

    x = pd.DataFrame(np.arange(40, dtype=float).reshape(8, 5), index=[f"l{i}" for i in range(8)], columns=list("abcde"))
    p = ci.permute_across_lines(x, seed=11)
    assert p.equals(ci.permute_across_lines(x, seed=11))
    for c in x.columns:
        assert sorted(p[c]) == sorted(x[c])
    assert not p.equals(x)


def test_lcp_scores_known_answer():
    """A5 (amendment A4 d): s = w_r z_RNA + w_p z_prot + w_c loss - tau; z against training statistics;
    z_prot only for soluble ProCan genes of lines with proteomics; loss = -1 for accepted losses."""
    from nextgem import context_inputs as ci

    genes = pd.DataFrame({"soluble": [True, True, False, True]}, index=pd.Index(["a", "b", "c", "d"], name="gene"))
    rna = pd.DataFrame([[3.0, 1.0, 2.0, 5.0]], index=["L"], columns=list("abcd"))
    prot = pd.DataFrame([[4.0, np.nan, 9.0]], index=["L"], columns=list("abc"))  # b not detected; d not in ProCan
    zs = ci.ZStats(rna_mean=pd.Series({"a": 1.0, "b": 1.0, "c": 2.0, "d": 5.0}),
                   rna_sd=pd.Series({"a": 2.0, "b": 0.5, "c": 1.0, "d": 0.0}),
                   prot_mean=pd.Series({"a": 2.0, "b": 0.0, "c": 3.0}), prot_sd=pd.Series({"a": 1.0, "b": 2.0, "c": 1.0}),
                   prot_floor=-4.0)
    s = ci.lcp_scores("L", 1.0, 0.5, 1.0, -0.5, genes, rna, prot, zs, ["c"])
    # a: z_rna 1, z_prot 2 -> 1 + 1 + 0 + 0.5;  b: z_rna 0, z_prot (-4-0)/2 = -2 -> -1 + 0.5
    # c: non-soluble (z_prot 0), lost -> 0 - 1 + 0.5;  d: zero training SD -> z 0, not in ProCan -> 0.5
    assert s.to_dict() == pytest.approx({"a": 2.5, "b": -0.5, "c": -0.5, "d": 0.5})
    s_noprot = ci.lcp_scores("L", 1.0, 1.0, 0.0, 0.0, genes, rna, None, zs, [])
    assert s_noprot.to_dict() == pytest.approx({"a": 1.0, "b": 0.0, "c": 0.0, "d": 0.0})
    assert len(ci.LCP_GRID) == 24 and len(set(ci.LCP_GRID)) == 24


def _stage_a_module():
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location("w4_stage_a", Path(__file__).parents[1] / "scripts" / "w4_stage_a.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_stage_a_refuses_undeclared_solver():
    from types import SimpleNamespace

    from nextgem import context_run as cr

    mod = _stage_a_module()
    assert cx.DECLARED_SOLVER == {"F0": "glpk", "F2": "glpk", "F4": "glpk", "F5": "glpk", "F3": "gurobi", "F3b": "gurobi"}
    m = cobra.Model("t")
    m.solver = "glpk"
    old = dict(cr._W)
    try:
        cr._W["parent"] = SimpleNamespace(model=m)
        for fam, meta in (("F3", {}), ("F0", {"solver": "gurobi"})):  # declared, and declared-after-override
            job = cr.Job(build_id="x", line="L", family=fam, rule="R1", modality="rna", variant="real",
                         scores=pd.Series(dtype=float), meta=meta)
            with pytest.raises(RuntimeError, match="declared"):
                mod.worker_run(job, None)
    finally:
        cr._W.clear()
        cr._W.update(old)
    a = mod.parse_args(["--pilot", "--out", "x"])  # --pilot alone keeps the declared solvers
    assert a.solver_override == []
    with pytest.raises(SystemExit):
        mod.main(["--solver-override", "F0=gurobi", "--dry-run"])  # override needs --pilot


def test_glpk_worker_never_touches_gurobi():
    """Review B1: a GLPK worker (and the Stage-A parent's model read) must not create any gurobipy Model/Env."""
    import subprocess
    import sys

    code = (
        "import gurobipy\n"
        "def boom(*a, **k): raise SystemExit('gurobipy used')\n"
        "gurobipy.Model = boom; gurobipy.Env = boom\n"
        "from nextgem import context_run as cr\n"
        "cr.init_worker('human2', 'glpk', None)\n"
        "p = cr._W['parent']; assert cr.cx.solver_name(p.model) == 'glpk'\n"
        "import cobra; from nextgem import paths\n"
        "m = cobra.io.read_sbml_model(str(paths.data_root() / 'Human-GEM-main/model/Human-GEM.xml'))\n"
        "assert 'glpk' in m.solver.interface.__name__\n"
        "print('OK')\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0 and r.stdout.strip().endswith("OK"), r.stdout[-2000:] + r.stderr[-2000:]
    assert "Academic license" not in r.stdout + r.stderr  # no WLS banner


def test_build_independent_of_job_order():
    """Review B2: the same build gives a bitwise-identical open mask and KO table whether it runs first or after
    another build in the same lane (one fresh process per build), for a GLPK family (F0, repair) and for F3 on
    Gurobi (FASTCORE alternate optima)."""
    import json as _json

    from nextgem import context_inputs as ci
    from nextgem import context_run as cr
    from nextgem import paths

    mod = _stage_a_module()
    gm = pd.read_parquet(paths.results_dir() / "w4" / "threshold_gene_means.parquet").set_index("gene")
    st = _json.loads((paths.results_dir() / "w4" / "threshold_stats.json").read_text())
    import cobra as _cobra

    _cobra.Configuration().solver = "glpk"
    genes = ci.gene_table(_cobra.io.read_sbml_model(str(paths.data_root() / "Human-GEM-main/model/Human-GEM.xml")))
    rs = ci.ThresholdStats("rna", st["rna"]["p25"], st["rna"]["p75"], gm.rna_mean.dropna(), st["rna"]["floor_value"], 714)
    lines = ["ACH-002024", "ACH-000957"]
    rna = ci.read_expression(genes, lines)
    sub = tuple(sorted(genes.index)[::10])  # 287-gene KO subset keeps the test short

    def job(line, fam):
        lane = cx.DECLARED_SOLVER[fam]
        return cr.Job(build_id=f"t|{line}|{fam}", line=line, family=fam, rule="R1", modality="rna", variant="real",
                      scores=ci.signed_scores(line, "R1", "rna", genes, rna, rs), tasks=False, genes=sub,
                      meta={"solver": lane}), lane

    results = {}

    def on_done(key, status, res, payload, wall):
        assert status == "ok", res
        results[key] = res

    items = []
    for fam in ("F0", "F3"):
        x, lane = job("ACH-002024", fam)
        y, _ = job("ACH-000957", "F2" if fam == "F0" else "F3")
        items += [(f"{fam}:first", lane, x, None, None), (f"{fam}:other", lane, y, None, None),
                  (f"{fam}:after", lane, x, None, None)]
    mod.Runner({"glpk": 1, "gurobi": 1}, job_timeout=1800, fastcc_cache=str(mod.FASTCC)).run_all(items, on_done)
    for fam in ("F0", "F3"):
        a, b = results[f"{fam}:first"], results[f"{fam}:after"]
        assert a["rows"][0]["open_mask"] == b["rows"][0]["open_mask"]
        pd.testing.assert_frame_equal(a["ko"]["none"].reset_index(drop=True), b["ko"]["none"].reset_index(drop=True))


def test_a6prime_line_choice():
    """Amendment A7: per panel lineage, the panel line closest (Euclidean) to the lineage median profile over all
    its training lines; ties by ModelID; never a non-training line."""
    from nextgem import context_inputs as ci

    split = pd.DataFrame({"ModelID": ["a1", "a2", "a3", "a4", "b1", "b2", "b3", "t1"],
                          "tier": ["train"] * 7 + ["test"], "lineage": ["A"] * 4 + ["B"] * 3 + ["A"],
                          "has_expression": [True] * 8})
    rna = pd.DataFrame({"g1": [0, 2, 4, 10, 1, 1, 5, 100], "g2": [0, 2, 4, 10, 3, 3, 5, 100]},
                       index=split.ModelID, dtype=float)
    panel = pd.DataFrame({"ModelID": ["a1", "a3", "a4", "b2", "b1", "b3"], "lineage": ["A"] * 3 + ["B"] * 3})
    out = ci.a6prime_lines(panel, rna, split).set_index("lineage")
    # A: median over training a1..a4 = (3, 3); distances a1 18**.5, a3 2**.5, a4 98**.5 -> a3 (test line t1 ignored)
    assert out.loc["A", "ModelID"] == "a3" and out.loc["A", "distance"] == pytest.approx(2 ** 0.5)
    assert out.loc["A", "n_lineage_lines"] == 4
    # B: median (1, 3); b1 and b2 tie at 0 -> lower ModelID b1
    assert out.loc["B", "ModelID"] == "b1" and out.loc["B", "distance"] == 0
