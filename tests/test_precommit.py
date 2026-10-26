"""§11.1 pre-commit hook (installed by scripts/install_git_hooks.py) in a toy repository: rejects files
under data/ or papers/, files > 50 MB, scanner hits, test-keyed tables until the T2 unseal row exists,
FROZEN.tsv violations, and a failing test edited together with the code fix."""

from __future__ import annotations

import os
import subprocess

import numpy as np
import pandas as pd
import pytest

import toy


@pytest.fixture(scope="module")
def hooked(tmp_path_factory):
    t = toy.make_toy(tmp_path_factory.mktemp("pc"))
    r = toy.build(t)
    assert r.returncode == 0, r.stderr[-1500:]
    toy.commit_all(t.state, "labels built: ledger/FROZEN/reports")
    r = subprocess.run([toy.PY, str(toy.SCRIPTS / "install_git_hooks.py"), "--repo", str(t.state)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return t


def _commit(t, *files, msg="x"):
    env = {**os.environ, "NEXTGEM_DATA_ROOT": str(t.data)}
    toy.git(t.state, "add", "-f", *files)
    r = toy.git(t.state, "commit", "-q", "-m", msg, check=False, env=env)
    if r.returncode != 0:
        toy.git(t.state, "reset", "-q", "HEAD", "--", *files, check=False)
    return r


def _write(t, rel, text="x\n"):
    p = t.state / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return rel


def test_hooks_installed(hooked):
    for h in ("pre-commit", "pre-push"):
        assert "nextgem guard hook" in (hooked.state / ".git" / "hooks" / h).read_text()


def test_data_and_papers_rejected(hooked):
    r = _commit(hooked, _write(hooked, "data/some.csv"))
    assert r.returncode != 0 and "data/ or papers/" in r.stderr
    r = _commit(hooked, _write(hooked, "papers/x.txt"))
    assert r.returncode != 0
    r = _commit(hooked, _write(hooked, "notes/ok.md"))
    assert r.returncode == 0, r.stderr


def test_large_file_rejected(hooked):
    p = hooked.state / "results" / "big.bin"
    with open(p, "wb") as f:
        f.truncate(51 * 1024 * 1024)
    r = _commit(hooked, "results/big.bin")
    p.unlink()
    assert r.returncode != 0 and "50 MB" in r.stderr


def test_label_copy_rejected_expression_accepted(hooked):
    t = hooked
    ge = t.ge[toy.U_ENTREZ].iloc[:15]
    ge.columns = [f"G{e - 1000} ({e})" for e in ge.columns]
    ge.index.name = "ModelID"
    ge.to_csv(t.state / "results" / "innocent.csv")
    r = _commit(t, "results/innocent.csv")
    assert r.returncode != 0 and "fingerprint" in r.stderr
    (t.state / "results" / "innocent.csv").unlink()
    rng = np.random.default_rng(3)
    expr = pd.DataFrame(rng.gamma(2, 1.5, (len(t.ids_in("train")), 60)), index=t.ids_in("train"),
                        columns=[f"G{i} ({1000 + i})" for i in range(60)])
    expr.index.name = "ModelID"
    expr.to_csv(t.state / "results" / "expr_train.csv")
    r = _commit(t, "results/expr_train.csv")
    assert r.returncode == 0, r.stderr


def test_frozen_violation_rejected(hooked):
    t = hooked
    _write(t, "results/frozen_thing.txt", "v1\n")
    r = t.run("check_frozen.py", "--add", "results/frozen_thing.txt")
    assert r.returncode == 0, r.stdout + r.stderr
    r = _commit(t, "results/frozen_thing.txt", "FROZEN.tsv")
    assert r.returncode == 0, r.stderr
    _write(t, "results/frozen_thing.txt", "v2\n")
    r = _commit(t, "results/frozen_thing.txt")
    assert r.returncode != 0 and "frozen SHA-256" in r.stderr
    toy.git(t.state, "checkout", "--", "results/frozen_thing.txt")


def test_failing_test_not_edited_with_code_fix(hooked):
    t = hooked
    of = t.state / "results" / "test_runs" / "fake.tsv"
    of.parent.mkdir(parents=True, exist_ok=True)
    of.write_text("nodeid\toutcome\ntests/test_x.py::test_a\tfailed\n")
    (t.state / "results" / "test_runs.tsv").write_text(
        "timestamp\tcommit\tfull_run\tn_passed\tn_failed\tn_xfailed\tn_skipped\tn_error\ttestids_sha256\toutcomes_file\n"
        "2026-10-27T00:00:00\tc\tFalse\t0\t1\t0\t0\t0\tx\tresults/test_runs/fake.tsv\n")
    r = _commit(t, _write(t, "tests/test_x.py", "def test_a(): pass\n"), _write(t, "src/y.py", "y = 1\n"))
    assert r.returncode != 0 and "same commit" in r.stderr
    r = _commit(t, "src/y.py")
    assert r.returncode == 0, r.stderr
    (t.state / "results" / "test_runs.tsv").unlink()


def test_test_keyed_table_until_unseal(hooked, monkeypatch):
    t = hooked
    rng = np.random.default_rng(5)
    pred = pd.DataFrame(rng.uniform(0, 1, (5, 60)), index=t.ids_in("test"), columns=[f"G{i} ({1000 + i})" for i in range(60)])
    pred.index.name = "ModelID"
    pred.to_csv(t.state / "results" / "test_predictions.csv")
    r = _commit(t, "results/test_predictions.csv")
    assert r.returncode != 0 and "test-keyed" in r.stderr
    # T2 completed + unseal: ledger appends commit through the same hook
    monkeypatch.setenv("NEXTGEM_STATE_ROOT", str(t.state))
    monkeypatch.setenv("NEXTGEM_DATA_ROOT", str(t.data))
    from nextgem import ledger

    with ledger.look("test", "T2", workstream="M7", config_sha="c"):
        pass
    ledger.append("test", "unseal")
    r = _commit(t, "results/test_predictions.csv")
    assert r.returncode == 0, r.stderr


def test_record_commit_not_blocked_by_failing_full_run(hooked):
    """Review S1: a test failure blocks only commits that touch code, tests or FROZEN.tsv."""
    t = hooked
    d = t.state / "results" / "test_runs"
    d.mkdir(parents=True, exist_ok=True)
    (d / "r1.tsv").write_text("nodeid\toutcome\ntests/test_x.py::test_a\tpassed\n")
    (d / "r2.tsv").write_text("nodeid\toutcome\ntests/test_x.py::test_a\tfailed\n")
    hdr = "timestamp\tcommit\tfull_run\tn_passed\tn_failed\tn_xfailed\tn_skipped\tn_error\ttestids_sha256\toutcomes_file\n"
    (t.state / "results" / "test_runs.tsv").write_text(
        hdr + "t1\tc\tTrue\t1\t0\t0\t0\t0\tx\tresults/test_runs/r1.tsv\n" + "t2\tc\tTrue\t0\t1\t0\t0\t0\tx\tresults/test_runs/r2.tsv\n")
    try:
        r = _commit(t, _write(t, "notes/record.md", "a record\n"))
        assert r.returncode == 0, r.stderr
        r = _commit(t, _write(t, "src/z.py", "z = 1\n"))
        assert r.returncode != 0 and "earlier full run" in r.stderr
    finally:
        (t.state / "results" / "test_runs.tsv").unlink()


def test_pre_push_scans_pushed_history(hooked, tmp_path):
    """Review S4: a label copy committed and deleted again is still caught in the pushed commits."""
    t = hooked
    remote = tmp_path / "remote.git"
    subprocess.run(["git", "init", "-q", "--bare", str(remote)], check=True)
    toy.git(t.state, "remote", "add", "origin", str(remote), check=False)
    ge = t.ge[toy.U_ENTREZ].loc[t.ids_in("test")]
    ge.columns = [f"G{e - 1000} ({e})" for e in ge.columns]
    ge.index.name = "ModelID"
    ge.to_csv(t.state / "results" / "hist_copy.csv")
    toy.git(t.state, "add", "results/hist_copy.csv")
    toy.git(t.state, "commit", "-q", "--no-verify", "-m", "copy")  # simulates a bypassed pre-commit
    toy.git(t.state, "rm", "-q", "results/hist_copy.csv")
    toy.git(t.state, "commit", "-q", "--no-verify", "-m", "delete it again")
    env = {**os.environ, "NEXTGEM_DATA_ROOT": str(t.data), "NEXTGEM_PREPUSH_SKIP_TESTS": "1"}
    r = toy.git(t.state, "push", "-q", "origin", "HEAD:refs/heads/main", check=False, env=env)
    assert r.returncode != 0 and "hist_copy.csv" in r.stderr and "fingerprint" in r.stderr, r.stderr[-1500:]
