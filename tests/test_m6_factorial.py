"""Label-free toy tests for the M6 integration factorial (scripts/m6_*.py)."""

from __future__ import annotations

import json
import signal
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import m6_common as mc  # noqa: E402


def test_cells_fixed_by_d10():
    """D10 cells; candidate order per DESIGN §9 (medium < tasks < combination), D10 as corrected (review F-B1)."""
    assert [c.key for c in mc.CANDIDATES] == ["H1/M2", "H1+Bloo/M0", "H1+Bloo/M2"]
    assert mc.REFERENCE.method == "H1" and mc.CELL["H2"].method == "H2"
    assert len(mc.CELLS) == 8 and len({c.method for c in mc.CELLS}) == 8
    for c in mc.CANDIDATES:
        assert c.parent == "human1" and c.role == "candidate"
    for c in mc.CELLS:
        if c.parent == "human2" and c.key != "H2":
            assert c.role == "descriptive"
    assert (mc.CELL["H1+Bloo/M2"].variant, mc.CELL["H1+Bloo/M2"].step) == (mc.BLOO, mc.M2)
    assert mc.CELL["H1/M2"].variant == mc.T0 and mc.CELL["H1+Bloo/M0"].step == mc.M0


def test_method_names_acceptable_to_harness():
    import evaluate_val as ev

    for c in mc.CELLS:
        m = c.method
        assert not ev.LABEL_LIKE_RE.search(m) and not ev.ID_RE.match(m) and m.strip() == m


def _toy_line_sweeps():
    rows = []
    for p in mc.PARENTS:
        for l, med in (("ACH-000001", "a"), ("ACH-000002", "b"), ("ACH-000003", "a")):
            rows.append({"ModelID": l, "parent": p, "ladder_step": "M0", "sweep_id": f"{p}:m0"})
            rows.append({"ModelID": l, "parent": p, "ladder_step": "M2(f=1)", "sweep_id": f"{p}:{med}"})
            rows.append({"ModelID": l, "parent": p, "ladder_step": "M1", "sweep_id": f"{p}:x"})
    return pd.DataFrame(rows)


def test_cell_matrix_maps_lines_to_their_sweep_rows():
    ls = _toy_line_sweeps()
    U = [10, 20]
    Esw = pd.DataFrame([[0.0, 1.0], [0.5, 0.0], [1.0, 1.0], [0.2, 0.3]], columns=U,
                       index=[mc.sweep_key(mc.T0, "human1:m0"), mc.sweep_key(mc.BLOO, "human1:a"),
                              mc.sweep_key(mc.BLOO, "human1:b"), mc.sweep_key(mc.BLOO, "human1:m0")])
    M = mc.cell_matrix(mc.CELL["H1+Bloo/M2"], ls, Esw, ["ACH-000002", "ACH-000001", "ACH-000003"])
    assert M.loc["ACH-000002"].tolist() == [1.0, 1.0]
    assert M.loc["ACH-000001"].tolist() == [0.5, 0.0] == M.loc["ACH-000003"].tolist()
    M0 = mc.cell_matrix(mc.CELL["H1+Bloo/M0"], ls, Esw, ["ACH-000001"])
    assert M0.iloc[0].tolist() == [0.2, 0.3]
    with pytest.raises(KeyError):
        mc.cell_matrix(mc.CELL["H1/M2"], ls, Esw, ["ACH-000001"])  # T0|human1:a missing
    with pytest.raises(KeyError):
        mc.cell_matrix(mc.CELL["H1"], ls, Esw, ["ACH-999999"])


def test_line_sweep_ids_refuses_duplicates():
    ls = _toy_line_sweeps()
    ls = pd.concat([ls, ls.iloc[[0]]])
    with pytest.raises(ValueError):
        mc.line_sweep_ids(ls, "human1", "M0")


def test_e_from_ratio_absent_genes_and_failed_rows():
    emap = pd.Series([10, 20], index=["g1", "g2"])
    R = pd.DataFrame([[0.0, 0.4999996], [np.nan, np.nan]], columns=["g1", "g2"], index=["ok", "failed"])
    E = mc.e_from_ratio(R, emap, [10, 20, 30])
    assert E.loc["ok"].tolist() == [1.0, 0.5, 0.0]  # absent gene 30 -> 0; 6 decimals
    assert E.loc["failed"].isna().all()


def test_gr_fail_uses_the_task_lb():
    r = pd.Series([0.0, 0.01, 0.0152, 0.0153, 1.0], index=list("abcde"))
    out = mc.gr_fail(r, 65.5734)  # 1/65.5734 = 0.01525
    assert out.tolist() == [True, True, True, False, False]


def test_interaction_toy():
    idx = ["l1", "l2", "l3"]
    t0 = pd.Series([0.30, 0.30, 0.30], index=idx)
    t1 = pd.Series([0.34, 0.34, np.nan], index=idx)
    t2 = pd.Series([0.32, 0.32, 0.32], index=idx)
    t3 = pd.Series([0.36, 0.35, 0.36], index=idx)
    I = mc.interaction(t0, t1, t2, t3)
    assert list(I.index) == ["l1", "l2"]  # l3 not finite in every cell
    assert I.round(10).tolist() == [0.0, -0.01]  # additive -> 0; sub-additive -> negative


def test_tg_rule():
    assert mc.tg_pass(0.005, 0.001, 9)
    assert not mc.tg_pass(0.0049, 0.001, 12)
    assert not mc.tg_pass(0.02, 0.0, 12)
    assert not mc.tg_pass(0.02, 0.01, 8)
    assert not mc.tg_pass(np.nan, 0.01, 12)


def test_cpu_limit_signal_raises_timeout():
    import m6_sweeps as ms

    old = signal.signal(signal.SIGPROF, ms._on_prof)
    try:
        signal.setitimer(signal.ITIMER_PROF, 0.05)
        with pytest.raises(ms.JobCPUTimeout):
            x = 0
            while True:
                x += 1
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0.0)
        signal.signal(signal.SIGPROF, old)


@pytest.mark.parametrize("behaviour,expected", [
    ("ok", "ok"), ("timeout", "failed-timeout"), ("lp_time_limit", "failed-timeout"),
    ("nonoptimal", "failed-nonoptimal"), ("infeasible", "failed-infeasible"), ("zero_wt", "failed-wt")])
def test_job_status_mapping(monkeypatch, behaviour, expected):
    import m6_sweeps as ms
    from nextgem import w3

    genes = ["g1", "g2"]
    monkeypatch.setattr(ms, "_W", {"parent": "human1", "model": None, "genes": genes, "closures": {}})

    def fake(model, g, clos, bounds, v):
        if behaviour == "timeout":
            raise ms.JobCPUTimeout()
        if behaviour == "lp_time_limit":
            raise w3.SolverInfraError("non-optimal LP status 'time_limit' then 'time_limit'")
        if behaviour == "nonoptimal":
            raise w3.SolverInfraError("non-optimal LP status 'undefined' then 'undefined'")
        wt, st = {"infeasible": (0.0, "infeasible"), "zero_wt": (0.0, "optimal")}.get(behaviour, (50.0, "optimal"))
        return {"wt_growth": wt, "wt_status": st, "ratio_r": np.array([0.0, 1.0]), "n_solved": 1, "n_genes_lp": 1,
                "n_infeasible": 0, "shortcut_verify_max_dev": 0.0}

    monkeypatch.setattr(w3, "ko_sweep", fake)
    job = {"parent": "human1", "key": "k", "variant": mc.BLOO, "w2_sweep_id": "human1:a", "purpose": "new",
           "bounds": {}, "variant_json": w3.Variant(mc.T0).to_json(), "cpu_limit": 100}
    r = ms._job(job)
    assert r["status"] == expected
    if expected == "ok":
        assert r["ratio_r"].tolist() == [0.0, 1.0]
    else:
        assert np.isnan(r["ratio_r"]).all()


def test_look_spec_structure_and_workstream_rule(tmp_path):
    """The M6 look spec passes every static harness check, as a milestone-level look (A13)."""
    import evaluate_val as ev
    import m6_val_spec as vs

    spec = vs.build_spec(scores_sha="0" * 64, provenance={"scripts/m6_common.py": "1" * 64}, train_summary={})
    assert spec["candidates"] == [c.method for c in mc.CANDIDATES]
    assert spec["methods"][:2] == ["H1", "H2"] and set(spec["methods"]) == {c.method for c in mc.CELLS}
    assert spec["reference"] == "H1" and spec["scores"]["layout"] == "line_gene"
    ev.check_spec(spec, "M6_factorial", "m6")
    alt = json.loads(json.dumps(spec))
    alt["milestone"] = "M4"  # an M<n> look must be that milestone's
    with pytest.raises(ev.PreflightError, match="workstream"):
        ev.check_spec(alt, "M6_factorial", "m6")
