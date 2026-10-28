"""Tests for the M5 / A12(b) context pre-check (scripts/m5_build.py, scripts/m5_score_train.py). No label file is read."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
import pytest

from nextgem import context as cx
from nextgem.model import BIOMASS_REACTION as BIO

ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def build_mod():
    return _script("m5_build")


@pytest.fixture(scope="module")
def score_mod():
    return _script("m5_score_train")


def _toy(solver="glpk"):
    """Two isozyme routes to p; MAM10012c (the cofactor pool) is a biomass reactant made only via gF."""
    m = cobra.Model("toy")
    m.solver = solver
    spec = [("EX_glc", "--> glc", "", 0, 10), ("R1", "glc --> a", "gA", 0, 1000), ("R2", "a --> p", "gB1 or gB2", 0, 1000),
            ("R3", "glc --> MAM10012c", "gF", 0, 1000), (BIO, "p + 0.1 MAM10012c --> ", "", 0, 1000)]
    for rid, eq, gpr, lb, ub in spec:
        r = cobra.Reaction(rid)
        m.add_reactions([r])
        r.build_reaction_from_string(eq, verbose=False)
        r.gene_reaction_rule = gpr
        r.bounds = (lb, ub)
    m.objective = BIO
    return m


def test_m5_variant_spec_is_cofactor_pool(build_mod):
    v = build_mod.load_variant()
    assert v.name == "B-loo:cofactor_pool_biomass" and v.kind == "remove" and v.params == {"mets": ["MAM10012c"]}


def test_m5_variant_applies_permanently_and_rescues_toy(build_mod):
    m = _toy()
    assert m.slim_optimize() > 0
    with m:  # gF is lethal with the pool in the biomass
        m.genes.gF.knock_out()
        assert m.slim_optimize(error_value=0.0) == pytest.approx(0.0)
    build_mod.load_variant().apply(m)  # outside any context: permanent
    assert "MAM10012c" not in {x.id for x in m.reactions.get_by_id(BIO).metabolites}
    with m:
        m.genes.gF.knock_out()
        assert m.slim_optimize() == pytest.approx(10.0)


def test_m5_generic_build_equals_brute_force_toy(build_mod):
    parent = cx.parent_from_model(_toy())
    res = build_mod.generic_build(parent, tasks=False)
    ko = res["ko"]["none"].set_index("gene").ratio_r
    expect = {}
    for g in parent.genes:
        with parent.model as m:
            m.genes.get_by_id(g).knock_out()
            expect[g] = round(m.slim_optimize(error_value=0.0) / parent.wt_generic, 6)
    assert ko.to_dict() == pytest.approx(expect)
    assert ko["gA"] == 0.0 and ko["gB1"] == 1.0 and ko["gF"] == 0.0
    assert res["rows"][0]["family"] == "generic" and res["rows"][0]["n_reopened"] == 0


def test_m5_enumerate_jobs(build_mod, score_mod):
    lines = [f"ACH-{i:06d}" for i in range(30)]
    specs = build_mod.enumerate_jobs(lines)
    assert len(specs) == 61 and len({s["build_id"] for s in specs}) == 61
    assert specs[0]["family"] == "generic"
    rest = specs[1:]
    assert sorted({s["family"] for s in rest}) == ["F2", "F4"]
    assert all(s["rule"] == "R1" and s["modality"] == "rna" and s["variant"] == "real" and s["genomic"] == ("none",)
               for s in rest)
    ids = {s["build_id"] + "|none" for s in rest}
    for fam in ("F2", "F4"):
        assert set(score_mod.arm_build_ids(lines, fam).values()) <= ids


def test_m5_task_lp_statuses(build_mod):
    assert build_mod.task_lp_statuses("a=3;b=0;c=infeasible;d=error: x") == {"solved": 2, "infeasible": 1, "error: x": 1}
    assert build_mod.task_lp_statuses("") == {}


def test_m5_cli_guards(build_mod, score_mod):
    for argv in (["--lines", "2"], ["--workers", "9", "--dry-run"], ["--job-timeout", "100", "--dry-run"],
                 ["--pilot", "--lines", "2"], ["--pilot", "--out", "results/m5/x", "--lines", "1"]):
        with pytest.raises(SystemExit):
            build_mod.main(argv)
    with pytest.raises(SystemExit, match="tmp"):
        score_mod.main(["--synthetic-labels", "--out", "results/m5/train_dry"])
    with pytest.raises(SystemExit, match="pilot"):  # refused before any label is read
        score_mod.main(["--builds", "results/m5/tmp/pilot", "--out", "results/m5/tmp/x"])


def test_m5_stage_b_verdict_threshold(score_mod):
    assert score_mod.A12B_MIN_DELTA == 0.005
    assert score_mod.stage_b_verdict(0.005) and score_mod.stage_b_verdict(0.0123)
    assert not score_mod.stage_b_verdict(0.00499) and not score_mod.stage_b_verdict(-0.01)
    assert not score_mod.stage_b_verdict(float("nan"))


def test_m5_single_config_contrast_known_answer(score_mod):
    lines = [f"L{i}" for i in range(9)]
    lineage = pd.Series(["A"] * 3 + ["B"] * 3 + ["C"] * 3, index=lines)
    delta = pd.Series([0.01, 0.02, 0.03, -0.01, -0.01, 0.04, 0.0, np.nan, 0.06], index=lines)
    r = score_mod.single_config_contrast(delta, lineage, "x", B=200)
    expect = np.mean([0.02, (-0.01 - 0.01 + 0.04) / 3, 0.03])
    assert r["delta_lin"] == pytest.approx(expect)
    assert r["folds_positive"] == 3 and r["n_folds"] == 3 and r["n_lines_estimable"] == 8
    assert r["fold_means"]["C"] == pytest.approx(0.03)
    assert r["lo"] <= r["delta_lin"] <= r["hi"]


def test_m5_h1_entrez_map_matches_generic_table(score_mod):
    from nextgem import labels, paths, w1
    from nextgem import w4_analysis as wa

    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    gen = w1.generic_table(U.entrez)
    h1 = pd.read_parquet(paths.baseline_dir() / "human1_ko.parquet", columns=["gene", "ratio_r"]).set_index("gene").ratio_r
    S = wa.score_frame(pd.DataFrame([h1.to_numpy()], index=["x"], columns=h1.index), genes, score_mod.h1_entrez_map())
    assert np.array_equal(S.iloc[0].to_numpy(), gen["e_H1"].reindex(genes).fillna(0).to_numpy())


def test_m5_worker_parent_is_h1_bloo_on_glpk():
    """The worker init loads Human1 on GLPK (no gurobipy Model/Env, no WLS banner) with the B-loo biomass applied."""
    code = (
        "import sys, gurobipy\n"
        "def boom(*a, **k): raise SystemExit('gurobipy used')\n"
        "gurobipy.Model = boom; gurobipy.Env = boom\n"
        f"sys.path.insert(0, {str(ROOT / 'scripts')!r})\n"
        "import m5_build as mb\n"
        "from nextgem import context_run as cr\n"
        "mb.init_parent('glpk')\n"
        "p = cr._W['parent']\n"
        "assert cr.cx.solver_name(p.model) == 'glpk' and p.name == mb.PARENT_TAG\n"
        "b = p.model.reactions.get_by_id('MAR13082')\n"
        "assert 'MAM10012c' not in {m.id for m in b.metabolites}, 'B-loo not applied'\n"
        "assert len(p.tasks) == 10 and p.model.notes.get('nextgem_release') == 'human1'\n"
        "print('WT', round(p.wt_generic, 4)); print('OK')\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=900, cwd=str(ROOT))
    assert r.returncode == 0 and r.stdout.strip().endswith("OK"), r.stdout[-2000:] + r.stderr[-2000:]
    assert "Academic license" not in r.stdout + r.stderr
