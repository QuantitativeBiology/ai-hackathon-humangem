"""The shared Gurobi flock slot pool (DESIGN.md §11.4): n slots, the (n+1)-th waits, release frees,
and the lock is held across processes."""

import subprocess
import sys

import pytest

from nextgem import gurobi_pool as gp


def test_pool_limits_and_release(tmp_path):
    a = gp.hold_slot(2, tmp_path, "t", timeout=0)
    b = gp.hold_slot(2, tmp_path, "t", timeout=0)
    assert {a.index, b.index} == {0, 1}
    assert gp.try_slot(2, tmp_path, "t") is None
    with pytest.raises(gp.SlotTimeout):
        gp.hold_slot(2, tmp_path, "t", timeout=0.2, poll=0.05)
    a.release()
    c = gp.hold_slot(2, tmp_path, "t", timeout=0)
    assert c.index == 0
    with gp.gurobi_slot(n=3, directory=tmp_path, prefix="t", timeout=0) as d:
        assert d.index == 2 and d.held
    assert not d.held
    b.release(), c.release()


def test_slot_blocks_other_process(tmp_path):
    held = gp.hold_slot(1, tmp_path, "x", timeout=0)
    code = ("import sys; from nextgem import gurobi_pool as gp; "
            f"sys.exit(0 if gp.try_slot(1, {str(tmp_path)!r}, 'x') is None else 1)")
    assert subprocess.run([sys.executable, "-c", code]).returncode == 0  # busy while we hold it
    held.release()
    code_free = code.replace("is None else 1", "is not None else 1")
    assert subprocess.run([sys.executable, "-c", code_free]).returncode == 0


def test_env_override(monkeypatch):
    monkeypatch.setenv("NEXTGEM_GUROBI_SLOTS", "7")
    assert gp.n_slots() == 7
