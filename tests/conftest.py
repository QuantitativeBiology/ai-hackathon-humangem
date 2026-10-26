"""pytest configuration: toy fixtures, and the test-run ledger (DESIGN.md §11.3: results/test_runs.tsv).

Run the suite (from the project root):
    .venv/bin/python -m pytest tests -q --basetemp=<scratch dir>
A run over the whole ``tests`` directory without -k/-m selection is recorded as a full run.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import toy as toy_mod  # noqa: E402

_OUTCOMES: dict[str, str] = {}


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    """Tests start with the real roots; a test that redirects them uses monkeypatch."""
    monkeypatch.delenv("NEXTGEM_DATA_ROOT", raising=False)
    monkeypatch.delenv("NEXTGEM_STATE_ROOT", raising=False)


@pytest.fixture(scope="module")
def built_toy(tmp_path_factory):
    """A toy project with labels built through scripts/build_labels.py (kappa 'high')."""
    old = {k: os.environ.pop(k, None) for k in ("NEXTGEM_DATA_ROOT", "NEXTGEM_STATE_ROOT")}
    try:
        t = toy_mod.make_toy(tmp_path_factory.mktemp("toy"))
        r = toy_mod.build(t)
        assert r.returncode == 0, r.stderr[-2000:]
        t.build_stdout = r.stdout
        return t
    finally:
        for k, v in old.items():
            if v is not None:
                os.environ[k] = v


def pytest_runtest_logreport(report):
    if report.when == "call" or (report.when == "setup" and report.outcome != "passed"):
        outcome = report.outcome
        if hasattr(report, "wasxfail"):
            outcome = "xfailed" if report.skipped else "xpassed"
        elif report.when == "setup" and report.failed:
            outcome = "error"
        _OUTCOMES[report.nodeid] = outcome


def pytest_sessionfinish(session, exitstatus):
    if not _OUTCOMES or os.environ.get("NEXTGEM_NO_TEST_RECORD"):
        return
    cfg = session.config
    tests_dir = Path(__file__).resolve().parent
    args = [Path(a.split("::")[0]).resolve() for a in cfg.args]
    full = (args == [tests_dir] or args == [tests_dir.parent]) and not cfg.option.keyword and not cfg.option.markexpr
    from nextgem import frozen

    head = subprocess.run(["git", "-C", str(tests_dir.parent), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    frozen.record_test_run(dict(_OUTCOMES), full, head)
