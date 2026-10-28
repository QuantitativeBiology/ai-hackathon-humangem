"""Label-free toy tests for the T1 preparation scripts (scripts/t1_prep_*.py): spec shape vs the test harness's
static checks, call construction (P57/PGR line-constant, FBA<0.01 line-specific), score renaming, the per-task table,
the §6.5 skew decision, and the task-matrix helpers. No label file and no real model is touched."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import t1_prep_common as tc  # noqa: E402

_S = importlib.util.spec_from_file_location("evaluate_test", SCRIPTS / "evaluate_test.py")
et = importlib.util.module_from_spec(_S)
_S.loader.exec_module(et)
ev = et.ev
HEX = "a" * 64


def _toy_spec(**over):
    spec = {
        "spec_id": "P1_M6", "kind": "p1", "milestone": "M6", "primary": tc.PRIMARY, "methods": list(tc.METHODS),
        "arms": dict(tc.ARMS),
        "scores": {"validation": {"path": "results/m6/t1/p1_scores_validation.parquet", "sha256": HEX, "layout": "line_gene"},
                   "test": {"path": "results/m6/tmp/t1/p1_scores_test.parquet", "sha256": HEX, "layout": "line_gene"}},
        "calls": {"validation": {"path": "results/m6/t1/p1_calls_validation.parquet", "sha256": HEX, "layout": "line_gene"},
                  "test": {"path": "results/m6/tmp/t1/p1_calls_test.parquet", "sha256": HEX, "layout": "line_gene"}},
        "failed_builds": {"validation": None, "test": None},
        "attribution": {k: (list(v) if isinstance(v, list) else v) for k, v in tc.ATTRIBUTION.items()},
        "task_calls": {"path": "results/m6/t1/p1_task_calls.parquet", "sha256": HEX},
        "gene_annotations": {"path": "results/m6/t1/p1_gene_annotations.parquet", "sha256": HEX},
        "bootstrap_method": "percentile", "fitted_on_labels": False,
        "provenance": {"scripts/t1_prep_build.py": HEX, "results/m6/t1/skew.json": HEX}, "notes": "",
    }
    spec.update(over)
    return spec


# ============================================================================= spec vs harness
def test_method_names_pass_harness_checks():
    et._check_method_names(list(tc.METHODS), "toy")
    for m in tc.METHODS:
        assert not ev.LABEL_LIKE_RE.search(m) and not ev.ID_RE.match(m) and m != "G*"
    assert set(tc.SOURCE_COLUMN) == set(tc.METHODS)
    assert tc.SOURCE_COLUMN[tc.P1_H1] == "M6:H1+Bloo/M2(f=1)" and tc.SOURCE_COLUMN[tc.MED_H1] == "M6:H1/M2(f=1)"


def test_attribution_follows_template_rules():
    a = tc.ATTRIBUTION
    assert a["i"] == [tc.PRIMARY, "G*"] and a["iv"] == [tc.PRIMARY, "G*"]  # no curation in P1: P1 = P1 w/ generic GPRs
    assert a["iii"] == a["v"] == et.LEVER_NOT_SELECTED == tc.LEVER_NOT_SELECTED
    assert a["ii"][0] == tc.MED_H1 and a["ii"][1] in ("H1", "G*")
    assert set(a) == set(et.ATTRIBUTION)


def test_toy_spec_passes_check_p1_spec():
    et.check_p1_spec(_toy_spec())


def test_toy_spec_refusals():
    with pytest.raises(et.PreflightError):
        et.check_p1_spec(_toy_spec(bootstrap_method="bootstrap"))
    with pytest.raises(et.PreflightError):
        et.check_p1_spec(_toy_spec(attribution={**tc.ATTRIBUTION, "i": ["H1", "G*"]}))


def test_test_files_are_under_tmp_and_gitignored_pattern():
    import t1_prep_build as tb

    for k in ("scores_test", "calls_test"):
        r = tc.rel(tb.OUT[k])
        assert r.startswith("results/m6/tmp/") and r.endswith(".parquet")
    for k in ("scores_validation", "calls_validation", "task_calls", "gene_annotations"):
        assert "/tmp/" not in tc.rel(tb.OUT[k]) and ".large." not in tc.rel(tb.OUT[k])


def test_annotation_classes_known_to_harness():
    import t1_prep_build as tb

    for c in ("paralog_buffered", "medium_sensitive", "complex_subunit", "housekeeping_adjacent",
              "biomass_artefact_lethal", "expected_fba_failure"):
        assert c in et.FAILURE_GENE_CLASSES
    assert set(tb.EXPECTED_FBA_FAILURE_SYMBOLS) == {"MAT2A", "PRMT5", "GPX4", "PGD"}


# ============================================================================= calls
def test_fba_call_from_e_cutoff_edges():
    e = np.array([0.0, 0.99, 0.990001, 1.0, 0.5, 0.989999])
    assert tc.fba_call_from_e(e).tolist() == [False, False, True, True, False, False]


def test_p57_pgr_absent_genes_false_and_gr_subset():
    long = pd.DataFrame({"task": ["T1", "T1", "GR", "GR", "T2"], "gene": ["g1", "g2", "g2", "g3", "gX"],
                         "feasible": [False, True, False, False, False]})
    emap = pd.Series({"g1": 10, "g2": 20, "g3": 30})  # gX has no Entrez
    out = tc.p57_pgr(long, "GR", [10, 20, 30, 40], emap)
    assert out.p57.tolist() == [True, True, True, False]
    assert out.pgr.tolist() == [False, True, True, False]
    assert (out.pgr <= out.p57).all()


def test_calls_long_shape_and_values():
    gc = pd.DataFrame({"p57": [True, False, False], "pgr": [True, False, False]}, index=pd.Index([1, 2, 3], name="entrez"))
    E = pd.DataFrame([[1.0, 0.0, 0.995], [0.2, 1.0, 0.0]], index=["ACH-000001", "ACH-000002"], columns=[1, 2, 3])
    df = tc.calls_long(tc.P1_H1, gc, E)
    assert len(df) == 6 and set(df.method) == {tc.P1_H1}
    w = df.set_index(["ModelID", "entrez"])
    assert w.loc[("ACH-000001", 3), "fba"] and not w.loc[("ACH-000002", 3), "fba"]
    assert w.loc[("ACH-000002", 2), "fba"] and w.loc[("ACH-000002", 1), "p57"]  # p57 line-constant
    assert sorted(df.columns) == sorted(["ModelID", "entrez", "method", "p57", "pgr", "fba"])
    assert all(pd.api.types.is_bool_dtype(df[c]) for c in ("p57", "pgr", "fba"))
    assert pd.api.types.is_integer_dtype(df.entrez)


def test_calls_long_requires_gene_coverage():
    gc = pd.DataFrame({"p57": [True], "pgr": [False]}, index=pd.Index([1], name="entrez"))
    E = pd.DataFrame([[1.0, 0.0]], index=["ACH-000001"], columns=[1, 2])
    with pytest.raises(ValueError):
        tc.calls_long(tc.P1_H1, gc, E)


def test_task_calls_table_and_duplicates():
    long = pd.DataFrame({"task": ["A", "A", "B"], "gene": ["g1", "g2", "g1"], "feasible": [False, True, True]})
    out = tc.task_calls_table(long, pd.Series({"g1": 5, "g2": 6}))
    assert out.to_dict("records") == [{"entrez": 5, "task": "A", "fails": True}, {"entrez": 6, "task": "A", "fails": False},
                                      {"entrez": 5, "task": "B", "fails": False}]
    with pytest.raises(ValueError):
        tc.task_calls_table(long, pd.Series({"g1": 5, "g2": 5}))


# ============================================================================= scores
def test_rename_scores_subsets_and_renames():
    rows = []
    for l in ("ACH-000001", "ACH-000002", "ACH-000003"):
        for g in (1, 2):
            r = {"ModelID": l, "entrez": g}
            for c in set(tc.SOURCE_COLUMN.values()) | {"M6:H1+Bloo/M0", "M6:H2/M2(f=1)"}:
                r[c] = 0.25
            rows.append(r)
    df = pd.DataFrame(rows)
    out = tc.rename_scores(df, ["ACH-000003", "ACH-000001"])
    assert list(out.columns) == ["ModelID", "entrez", *tc.METHODS]
    assert sorted(out.ModelID.unique()) == ["ACH-000001", "ACH-000003"] and len(out) == 4
    assert all(pd.api.types.is_float_dtype(out[m]) for m in tc.METHODS)


# ============================================================================= skew
def test_skew_decision_threshold():
    rng = np.random.default_rng(0)
    sym = rng.normal(size=2000)
    d = tc.skew_decision(sym)
    assert d["bootstrap_method"] == "percentile" and d["n_lines"] == 2000 and abs(d["skew"]) < 1
    sk = rng.exponential(size=2000) ** 2  # strongly right-skewed
    d = tc.skew_decision(np.r_[sk, np.nan])
    assert d["bootstrap_method"] == "bca" and d["n_lines"] == 2000 and d["skew"] > 1


# ============================================================================= task-matrix helpers
def test_gene_long_maps_ko_groups_and_compare_generic():
    import t1_prep_tasks as tt

    tasks = pd.DataFrame({"key": ["ESS001", "ESS057"]})
    S = {"gd": {"g1": frozenset({"r1"}), "g2": frozenset({"r1"}), "g3": frozenset()}, "kos": [frozenset({"r1"})],
         "tasks": tasks}
    F = {0: np.array([False]), 1: np.array([True])}
    gm = pd.DataFrame({"gene": ["g1", "g2", "g3"], "entrez": ["1", "2", "3"]})
    long = tt.gene_long(S, F, gm)
    t0 = long[long.task == "ESS001"].set_index("gene").feasible
    assert t0.to_dict() == {"g1": False, "g2": False, "g3": True}  # g1/g2 share the KO set; g3 closes nothing
    gen = long.assign(seed=-1).copy()
    gen.loc[(gen.task == "ESS057") & (gen.gene == "g3"), "feasible"] = False
    ch = tt.compare_generic(long, gen)
    assert ch.to_dict() == {"ESS001": 0, "ESS057": 1}
