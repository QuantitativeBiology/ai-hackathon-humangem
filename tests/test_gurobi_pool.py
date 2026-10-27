"""The shared Gurobi flock slot pool (DESIGN.md §11.4): one global size, the (n+1)-th waits, release frees,
the lock is held across processes and released when the holder is killed."""

import os
import signal
import subprocess
import sys
import time

import pytest

from nextgem import gurobi_pool as gp


def test_pool_limits_and_release(tmp_path, monkeypatch):
    monkeypatch.setenv("NEXTGEM_GUROBI_SLOTS", "2")
    a = gp.hold_slot(tmp_path, "t", timeout=0)
    b = gp.hold_slot(tmp_path, "t", timeout=0)
    assert {a.index, b.index} == {0, 1}
    assert gp.try_slot(tmp_path, "t") is None
    with pytest.raises(gp.SlotTimeout):
        gp.hold_slot(tmp_path, "t", timeout=0.2, poll=0.05)
    a.release()
    with gp.gurobi_slot(directory=tmp_path, prefix="t", timeout=0) as c:
        assert c.index == 0 and c.held
    assert not c.held
    b.release()


def test_size_is_global_not_per_caller(monkeypatch):
    import inspect

    monkeypatch.delenv("NEXTGEM_GUROBI_SLOTS", raising=False)
    assert gp.n_slots() == gp.POOL_SLOTS
    monkeypatch.setenv("NEXTGEM_GUROBI_SLOTS", "7")
    assert gp.n_slots() == 7
    for fn in (gp.try_slot, gp.hold_slot):
        assert "n" not in inspect.signature(fn).parameters  # no per-caller size (review NOTE 6b)


def _holder(tmp_path, prefix):
    code = ("import sys, time; from nextgem import gurobi_pool as gp; "
            f"s = gp.hold_slot({str(tmp_path)!r}, {prefix!r}, timeout=0); print('held', flush=True); time.sleep(60)")
    p = subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.PIPE, text=True,
                         env={**os.environ, "NEXTGEM_GUROBI_SLOTS": "1"})
    assert p.stdout.readline().strip() == "held"
    return p


def test_slot_blocks_other_process_and_crash_releases(tmp_path, monkeypatch):
    monkeypatch.setenv("NEXTGEM_GUROBI_SLOTS", "1")
    p = _holder(tmp_path, "x")
    try:
        assert gp.try_slot(tmp_path, "x") is None  # busy while the other process holds it
        p.send_signal(signal.SIGKILL)
        p.wait(10)
        t0 = time.monotonic()
        s = gp.hold_slot(tmp_path, "x", timeout=5, poll=0.05)  # freed by the crash
        assert s.held and time.monotonic() - t0 < 5
        s.release()
    finally:
        if p.poll() is None:
            p.kill()
