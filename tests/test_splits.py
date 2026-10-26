"""§12 "Splits": determinism; reproduces the simulated draws in design_numbers; no lineage in two tiers;
all lines assigned; refuses to overwrite; refuses any S whose attempt range overlaps the retired range;
refuses code, DESIGN.md or FROZEN.tsv differing from the freeze commit, and dirty code.

Only simulated seeds are used (1..1000, as in design_numbers.py, and >= 500000 for sandbox runs).
The production seed is never run here."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess

import pandas as pd
import pytest

import toy
from nextgem import splits

FREEZE = "61b0beb"
SANDBOX_SEED = 500001


def _ref():
    out = subprocess.run(["git", "-C", str(toy.PROJECT), "show", f"{FREEZE}:results/design/design_numbers.json"],
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out)


def _make_splits_module():
    spec = importlib.util.spec_from_file_location("make_splits", toy.SCRIPTS / "make_splits.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def cohort():
    from nextgem import cohort as c

    return c.load()


@pytest.fixture(scope="module")
def ref():
    return _ref()


# ----------------------------------------------------------------------------- pure draw function
def test_seed_range_refusals():
    for s in (20260927, 20261000, 20261026, 20261125, 0, -1, True):
        with pytest.raises(splits.SeedRefused):
            splits.check_seed(s)
    for s in (1, 38 + 1000, 20260926, 20261126):
        splits.check_seed(s)
    with pytest.raises(splits.SeedRefused):  # per-attempt check inside the loop as well
        splits.run_attempts(20261026, {"S1": [], "S2": [], "S3": []}, {}, {})


def test_retired_seed_refused_by_cli(tmp_path):
    r = subprocess.run([toy.PY, str(toy.SCRIPTS / "make_splits.py"), "--seed", "20261026", "--freeze-commit", FREEZE,
                        "--sandbox", str(tmp_path / "sb")], capture_output=True, text=True)
    assert r.returncode == 2 and "retired" in r.stderr
    assert not (tmp_path / "sb" / "results" / "splits" / "SPLIT_SHA256").exists()


def test_determinism(ref):
    strata = ref["split_simulation"]["strata"]
    assert splits.draw(12345, strata) == splits.draw(12345, strata)
    assert splits.draw(12345, strata) != splits.draw(12346, strata)


def test_reproduces_reference_draws(ref):
    sim = ref["split_simulation"]
    for sd, d in sim["reference_draws_seeds_1_to_20"].items():
        t, v = splits.draw(int(sd), sim["strata"])
        assert (list(t), list(v)) == (d["test"], d["val"]), sd


def test_cohort_matches_design_numbers(cohort, ref):
    assert cohort.lineage_counts() == ref["cohort"]["lineage_counts"]
    assert len(cohort.screened) == ref["cohort"]["n_screened_lines"]
    assert {k: int(v) for k, v in cohort.model_type.value_counts().items()} == ref["cohort"]["model_type_screened"]
    assert {k: int(v) for k, v in cohort.library_class.value_counts().items()} == ref["screens"]["library_class_counts"]
    assert {k: int(v) for k, v in cohort.medium_source.value_counts().items()} == ref["medium"]["tier_counts"]
    assert {k: int(v) for k, v in cohort.base_medium_class.value_counts().items()} == ref["medium"]["base_class_counts"]
    S = set(cohort.screened)
    cov = ref["coverage_screened"]
    assert len(S & cohort.expr) == cov["expression_default"] and len(S - cohort.cn) == cov["no_cn"]
    assert len(S & cohort.nusinow) == cov["nusinow_screened"] and len(S & cohort.procan) == cov["procan_screened"]
    assert len(S & cohort.proteomics_any) == cov["any_proteomics_screened"]
    assert len(S & cohort.metabolomics) == cov["metabolomics_screened"]
    assert len(S & cohort.prism) == cov["prism_24q2_or_20q2_screened"] and len(S & cohort.ctrp) == cov["ctrp_screened"]
    _, _, strata = splits.eligible_and_strata(cohort.lineage_counts())
    assert strata == ref["split_simulation"]["strata"]


EXCLUDED = range(38, 138)  # simulated seeds overlapping the production attempt range are never evaluated here


def test_frozen_function_equality_and_pinned_inputs(cohort, ref):
    """M3: per-lineage expression/proteomics counts equal coverage_by_lineage of the frozen JSON; the draw
    and checks equal the verbatim 61b0beb functions (seeds 1-1000 except 38-137, never evaluated here);
    the split.* pins in FROZEN.tsv match. (The 1-1000 acceptance statistics include 38-137, so that
    cross-check runs only in the production draw.)"""
    ms = _make_splits_module()
    sim = ref["split_simulation"]
    assert ms.check_coverage(cohort, ref) == []
    e, p = cohort.by_lineage(cohort.expr), cohort.by_lineage(cohort.proteomics_any)
    fdraw, fchecks = splits.frozen_reference_functions(toy.PROJECT, FREEZE, sim["strata"], e, p)
    for sd in (x for x in range(1, 1001) if x not in EXCLUDED):
        t, v = splits.draw(sd, sim["strata"])
        assert (t, v) == tuple(map(list, fdraw(sd)))
        assert splits.checks(t, v, e, p) == fchecks(t, v)
    pins, inputs = ms.compute_pins(cohort, toy.PROJECT / "data")
    recorded = ms.frozen_rows(toy.PROJECT, None, "value")
    assert ms.check_pins(pins, inputs, ms.resolve(toy.PROJECT, FREEZE), recorded, ref) == []
    tampered = {**recorded, "split.expr_by_lineage_sha256": "0" * 64}
    assert ms.check_pins(pins, inputs, ms.resolve(toy.PROJECT, FREEZE), tampered, ref)
    bad = dict(inputs, **{"depmap/Model.csv": "0" * 64})
    assert any("Model.csv" in f for f in ms.check_pins(pins, bad, ms.resolve(toy.PROJECT, FREEZE), recorded, ref))


def test_attempt_loop_matches_frozen_loop(cohort, ref):
    """M3: run_attempts accepts the same attempt as a loop on the frozen functions (simulated seeds)."""
    ms = _make_splits_module()
    sim = ref["split_simulation"]
    e, p = cohort.by_lineage(cohort.expr), cohort.by_lineage(cohort.proteomics_any)
    fdraw, fchecks = splits.frozen_reference_functions(toy.PROJECT, FREEZE, sim["strata"], e, p)
    redraws = 0
    for s in range(600000, 600200):
        t, v, log = splits.run_attempts(s, sim["strata"], e, p)
        assert ms.frozen_attempt_loop(s, fdraw, fchecks) == (len(log) - 1, list(t), list(v))
        redraws += len(log) > 1
    assert redraws > 0  # the redraw path is exercised


def _redraw_seed(cohort, ref):
    """First simulated seed >= 500000 whose attempt 0 fails a check (found in memory)."""
    sim = ref["split_simulation"]
    e, p = cohort.by_lineage(cohort.expr), cohort.by_lineage(cohort.proteomics_any)
    for s in range(500100, 501000):
        t, v = splits.draw(s, sim["strata"])
        if not all(splits.checks(t, v, e, p).values()):
            return s
    raise AssertionError("no redraw seed found")


# ----------------------------------------------------------------------------- sandbox end-to-end
def _run_sandbox(sb, seed, *extra):
    return subprocess.run([toy.PY, str(toy.SCRIPTS / "make_splits.py"), "--seed", str(seed), "--freeze-commit", FREEZE,
                           "--sandbox", str(sb), *extra], capture_output=True, text=True)


def test_sandbox_draw_end_to_end(tmp_path, cohort, ref):
    sb = tmp_path / "sb"
    r = _run_sandbox(sb, SANDBOX_SEED)
    assert r.returncode == 0, r.stderr[-3000:]
    out = sb / "results" / "splits"
    lt = pd.read_csv(out / "lineage_split.tsv", sep="\t")
    assert len(lt) == 1361 and lt.ModelID.is_unique and set(lt.ModelID) == set(cohort.screened)
    assert lt.groupby("lineage").tier.nunique().max() == 1  # no lineage in two tiers
    elig = lt[lt.stratum != "rare"]
    assert elig.groupby("tier").lineage.nunique().to_dict() == {"test": 7, "train": 12, "validation": 4}
    assert set(lt[lt.stratum == "rare"].tier) == {"train"} and (lt.stratum == "rare").sum() == 30
    for s, (nt, nv) in {"S1": (2, 1), "S2": (2, 1), "S3": (3, 2)}.items():
        g = elig[elig.stratum == s].groupby("tier").lineage.nunique()
        assert (g.get("test", 0), g.get("validation", 0)) == (nt, nv)
    t_expected, v_expected = splits.draw(SANDBOX_SEED, ref["split_simulation"]["strata"])
    meta = json.loads((out / "split_meta.json").read_text())
    if meta["n_redraws"] == 0:
        assert set(lt[lt.tier == "test"].lineage) == set(t_expected)
    assert (out / "SPLIT_SHA256").read_text().split()[0] == hashlib.sha256((out / "lineage_split.tsv").read_bytes()).hexdigest()
    assert meta["seed"] == SANDBOX_SEED and meta["freeze_commit"].startswith(FREEZE) and meta["sandbox"]
    assert set(meta["script_sha256"]) == {"scripts/make_splits.py", "src/nextgem/splits.py", "src/nextgem/cohort.py"}
    att = pd.read_csv(out / "attempt_log.tsv", sep="\t")
    assert len(att) == meta["n_attempts"] and att.accepted.iloc[-1]
    for c in ("has_expression", "has_cn", "has_nusinow", "has_procan", "has_metabolomics", "has_prism", "has_ctrp",
              "library_class", "n_screens", "screen_type", "model_type", "medium_source", "base_medium_class"):
        assert c in lt.columns
    # ledger rows in the sandbox (uncommitted), chain intact
    from nextgem import ledger

    rows = ledger.parse((sb / "results" / "test_ledger.tsv").read_text())
    assert [x["event"] for x in rows] == ["split_started", "split_completed"]
    assert rows[1]["output_sha"] == meta["split_sha256"] and ledger.verify_chain("test", rows) == []
    # refuses to overwrite
    r2 = _run_sandbox(sb, SANDBOX_SEED + 1)
    assert r2.returncode == 2 and "never redrawn" in r2.stderr
    assert (out / "SPLIT_SHA256").read_text().split()[0] == meta["split_sha256"]


def test_sandbox_records_redraws(tmp_path, cohort, ref):
    seed = _redraw_seed(cohort, ref)
    sb = tmp_path / "sb"
    r = _run_sandbox(sb, seed)
    assert r.returncode == 0, r.stderr[-3000:]
    out = sb / "results" / "splits"
    att = pd.read_csv(out / "attempt_log.tsv", sep="\t")
    meta = json.loads((out / "split_meta.json").read_text())
    assert meta["n_redraws"] >= 1 and len(att) == meta["n_redraws"] + 1
    assert (~att.accepted.iloc[:-1]).all() and att.accepted.iloc[-1]
    assert att.failed_checks.iloc[0] == "d_not_both_heme_heldout"
    assert list(att.rng_seed) == list(range(seed, seed + len(att)))
    assert meta["rejected_attempts"] == list(range(len(att) - 1))


def test_sandbox_refuses_inside_project(tmp_path):
    r = _run_sandbox(toy.PROJECT / "results" / "sandbox_should_not_exist", SANDBOX_SEED)
    assert r.returncode == 2 and not (toy.PROJECT / "results" / "sandbox_should_not_exist").exists()


# ----------------------------------------------------------------------------- git checks (toy repos)
TOY_DESIGN = """# design v3
line 2

### Amendment log (after the M1 freeze)
| Date | Section | Change | Reason | Before/after the affected analysis |
|---|---|---|---|---|
| — | — | — | — | — |

### Revision history
x
"""


@pytest.fixture()
def repo(tmp_path):
    """Toy repo: freeze commit (DESIGN.md) -> code commit (src, scripts, tests, FROZEN.tsv) -> seed source."""
    r = tmp_path / "repo"
    toy.init_repo(r)
    (r / "DESIGN.md").write_text(TOY_DESIGN)
    freeze = toy.commit_all(r, "freeze")
    for rel, text in {"src/nextgem/splits.py": "x = 1\n", "scripts/make_splits.py": "y = 2\n", "src/nextgem/cohort.py": "c = 1\n",
                      "tests/t.py": "z = 3\n"}.items():
        (r / rel).parent.mkdir(parents=True, exist_ok=True)
        (r / rel).write_text(text)
    sha = lambda rel: hashlib.sha256((r / rel).read_bytes()).hexdigest()  # noqa: E731
    (r / "FROZEN.tsv").write_text("item\tkind\tpath\tsha256_or_value\tfrozen_on\tamendment\tnote\n" + "".join(
        f"{p}\tfile\t{p}\t{sha(p)}\t2026-10-27\t\t\n" for p in ("scripts/make_splits.py", "src/nextgem/splits.py", "src/nextgem/cohort.py")))
    code = toy.commit_all(r, "code")
    (r / "results" / "splits").mkdir(parents=True)
    (r / "results" / "splits" / "SEED_SOURCE.md").write_text(f"Split seed supplied by the user: 777 (freeze commit {freeze})\n")
    toy.commit_all(r, "seed source")
    return r, freeze, code


def test_git_checks_pass_when_clean(repo):
    ms = _make_splits_module()
    r, freeze, code = repo
    assert ms.git_checks(r, freeze, code, 777) == []
    (r / "results" / "test_ledger.tsv").write_text("ledger rows\n")  # ledger rows may be committed in between
    toy.commit_all(r, "ledger row")
    assert ms.git_checks(r, freeze, code, 777) == []


def test_git_checks_literal_freeze_rule_refuses_code_added_after_freeze(repo):
    ms = _make_splits_module()
    r, freeze, code = repo
    fails = ms.git_checks(r, freeze, freeze, 777)
    assert any("differ between" in f for f in fails)


@pytest.mark.parametrize("change", ["src_commit", "design_commit", "frozen_commit", "dirty_src", "untracked_tests",
                                    "dirty_design", "design_lines_removed", "design_sentence_added_before_code", "seed_missing", "seed_source_absent",
                                    "freeze_not_ancestor"])
def test_git_checks_refuse(repo, change):
    ms = _make_splits_module()
    r, freeze, code = repo
    seed = 777
    if change == "src_commit":
        (r / "src/nextgem/splits.py").write_text("x = 2\n")
        toy.commit_all(r, "edit")
    elif change == "design_commit":
        (r / "DESIGN.md").write_text(TOY_DESIGN + "new\n")
        toy.commit_all(r, "edit")
    elif change == "frozen_commit":
        with open(r / "FROZEN.tsv", "a") as f:
            f.write("x\tvalue\t\tv\t\t\t\n")
        toy.commit_all(r, "edit")
    elif change == "dirty_src":
        (r / "src/nextgem/splits.py").write_text("x = 3\n")
    elif change == "untracked_tests":
        (r / "tests" / "new.py").write_text("")
    elif change == "dirty_design":
        (r / "DESIGN.md").write_text("changed\n")
    elif change == "design_lines_removed":
        # code commit rebuilt on top of a DESIGN.md that lost a line since the freeze commit
        (r / "DESIGN.md").write_text(TOY_DESIGN.replace("line 2\n", ""))
        code = toy.commit_all(r, "code2")
    elif change == "design_sentence_added_before_code":
        (r / "DESIGN.md").write_text(TOY_DESIGN.replace("line 2\n", "line 2\nan added sentence\n"))
        code = toy.commit_all(r, "code2")
    elif change == "seed_missing":
        seed = 778
    elif change == "seed_source_absent":
        toy.git(r, "rm", "-q", "results/splits/SEED_SOURCE.md")
        toy.commit_all(r, "rm seed")
    elif change == "freeze_not_ancestor":
        toy.git(r, "checkout", "-q", "--orphan", "other")
        toy.commit_all(r, "orphan")
        code = toy.git(r, "rev-parse", "HEAD").stdout.strip()
    fails = ms.git_checks(r, freeze, code, seed)
    assert fails, change


def test_cli_enforces_git_checks_before_any_draw(repo, tmp_path):
    r, freeze, code = repo
    (r / "src/nextgem/splits.py").write_text("dirty\n")
    sb = tmp_path / "sb"
    res = subprocess.run([toy.PY, str(toy.SCRIPTS / "make_splits.py"), "--seed", str(SANDBOX_SEED), "--freeze-commit", freeze,
                          "--code-commit", code, "--reference-commit", FREEZE, "--sandbox", str(sb), "--enforce-git",
                          "--repo", str(r)], capture_output=True, text=True)
    assert res.returncode == 2 and "uncommitted" in res.stderr, res.stderr[-2000:]
    assert not (sb / "results" / "splits" / "SPLIT_SHA256").exists()


def test_amendment_row_in_section17_allowed(repo):
    """A1: an amendment row added to the §17 table (replacing the placeholder) before C passes."""
    ms = _make_splits_module()
    r, freeze, _ = repo
    (r / "DESIGN.md").write_text(TOY_DESIGN.replace("| — | — | — | — | — |", "| 2026-10-27 | 4.1 | A1: code commit | why | before |"))
    code = toy.commit_all(r, "amendment A1 + code")
    assert ms.git_checks(r, freeze, code, 777) == []


def test_production_overrides_refused(tmp_path, monkeypatch):
    """M1: --repo, --reference-commit and NEXTGEM_* overrides are refused without --sandbox (before any draw)."""
    base = [toy.PY, str(toy.SCRIPTS / "make_splits.py"), "--seed", str(SANDBOX_SEED), "--freeze-commit", FREEZE]
    cases = [(base + ["--repo", str(tmp_path)], {}), (base + ["--reference-commit", "HEAD"], {}),
             (base, {"NEXTGEM_DATA_ROOT": str(tmp_path)}), (base, {"NEXTGEM_STATE_ROOT": str(tmp_path)})]
    import os

    marker = toy.PROJECT / "results" / "splits" / "SPLIT_SHA256"
    before = marker.read_bytes() if marker.exists() else None
    for args, env in cases:
        r = subprocess.run(args, capture_output=True, text=True, env={**os.environ, **env})
        assert r.returncode == 2 and "REFUSED" in r.stderr, (args, env, r.stderr[-500:])
    assert (marker.read_bytes() if marker.exists() else None) == before  # the real split is untouched


def test_sandbox_requires_large_seeds(tmp_path):
    """M2: sandbox seeds below 500000 are refused before anything is drawn."""
    for seed in (1, 499999):
        r = _run_sandbox(tmp_path / f"sb{seed}", seed)
        assert r.returncode == 2 and "500000" in r.stderr
        assert not (tmp_path / f"sb{seed}" / "results" / "splits").exists()
