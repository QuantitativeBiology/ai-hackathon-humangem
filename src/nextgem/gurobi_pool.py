"""Shared Gurobi session pool (DESIGN.md §11.4): a ``flock`` slot pool on /dev/shm.

Every process that opens a Gurobi (WLS) environment first holds one slot. The pool size is ONE global
number, the licence session limit allotted to the project: ``POOL_SLOTS`` (8), overridable only through
the environment variable ``NEXTGEM_GUROBI_SLOTS``; callers cannot pass their own size (review NOTE 6b).
All workstreams use the same directory and prefix, so W1, W2, W4 and later stages share one pool.
A worker waits for a slot and never falls back silently to another solver. A crashed holder frees its
slot (flock state dies with the process).

Processes that must NOT use Gurobi set ``cobra.Configuration().solver = "glpk"`` before loading any
model; building a cobra model under Gurobi opens a WLS environment outside this pool (review B1).

    from nextgem.gurobi_pool import hold_slot
    slot = hold_slot()          # blocks until a slot is free; keep the object alive while solving
    ...
    slot.release()              # or let process exit release it (flock is dropped on close/exit)

    with gurobi_slot():         # context-manager form
        ...
"""

from __future__ import annotations

import fcntl
import os
import time
from contextlib import contextmanager
from pathlib import Path

DEFAULT_DIR = Path("/dev/shm")
DEFAULT_PREFIX = "nextgem_gurobi_slot"
POOL_SLOTS = 8  # global licence pool for the project (coordinator, 2026-10-27: W2 finished)


class SlotTimeout(TimeoutError):
    """No slot became free within the timeout."""


class Slot:
    def __init__(self, fh, index: int, path: Path):
        self._fh, self.index, self.path = fh, index, path

    @property
    def held(self) -> bool:
        return self._fh is not None

    def release(self) -> None:
        if self._fh is not None:
            try:
                fcntl.flock(self._fh, fcntl.LOCK_UN)
            finally:
                self._fh.close()
                self._fh = None

    def __del__(self):  # pragma: no cover - best effort
        self.release()


def n_slots() -> int:
    """The global pool size: NEXTGEM_GUROBI_SLOTS if set, else POOL_SLOTS."""
    return int(os.environ.get("NEXTGEM_GUROBI_SLOTS", POOL_SLOTS))


def try_slot(directory: Path | str = DEFAULT_DIR, prefix: str = DEFAULT_PREFIX) -> Slot | None:
    """Take a free slot without waiting; None if all slots of the global pool are held."""
    n = n_slots()
    directory = Path(directory)
    for k in range(n):
        path = directory / f"{prefix}_{k}.lock"
        fh = open(path, "a")
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            fh.close()
            continue
        return Slot(fh, k, path)
    return None


def hold_slot(directory: Path | str = DEFAULT_DIR, prefix: str = DEFAULT_PREFIX,
              timeout: float | None = None, poll: float = 5.0) -> Slot:
    """Block until a slot is free (or ``timeout`` seconds pass: SlotTimeout)."""
    t0 = time.monotonic()
    while True:
        slot = try_slot(directory, prefix)
        if slot is not None:
            return slot
        if timeout is not None and time.monotonic() - t0 >= timeout:
            raise SlotTimeout(f"no free Gurobi slot among {n_slots()} after {timeout} s")
        time.sleep(poll)


@contextmanager
def gurobi_slot(**kw):
    slot = hold_slot(**kw)
    try:
        yield slot
    finally:
        slot.release()
