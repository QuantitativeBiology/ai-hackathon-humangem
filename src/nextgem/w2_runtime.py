"""Run helpers for the W2 scripts: parent loading through the shared Gurobi pool, and provenance.

Every W2 script runs under ``nextgem.run.provenance()`` (DESIGN §11.4) through ``Manifest`` below, which appends a
row to ``results/w2/PROVENANCE.tsv`` (argv, commit, dirty flag, full SHA-256 of inputs and outputs, solver and
versions, workers, wall/CPU time, peak RSS). ``allow_dirty=True`` is passed because the W2 code is run uncommitted
(the PI commits); the row records ``dirty=True``. The earlier ``results/w2/MANIFEST.tsv`` (pre-A9 runs) is kept as
history.

Gurobi: any process that loads a parent on Gurobi first holds one slot of ``nextgem.gurobi_pool`` (once per
process); processes that only build bounds load the parents on GLPK.
"""

from __future__ import annotations

import hashlib
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from nextgem import paths

PARENTS = ("human1", "human2")


def w2_dir() -> Path:
    d = paths.results_dir() / "w2"
    d.mkdir(parents=True, exist_ok=True)
    return d


_SLOT = None


def ensure_gurobi_slot():
    """Hold one slot of the shared Gurobi pool for the life of this process (idempotent)."""
    global _SLOT
    if _SLOT is None:
        from nextgem import gurobi_pool

        _SLOT = gurobi_pool.hold_slot()
    return _SLOT


def load_parent(name: str, solver: str = "gurobi"):
    """Load a parent release (tagged) with the given solver, single-threaded (Gurobi only inside a pool slot)."""
    import cobra

    if solver == "gurobi":
        ensure_gurobi_slot()
    from nextgem.model import load_model

    cobra.Configuration().solver = solver
    m = load_model(name)
    set_threads(m)
    return m


def set_threads(model, threads: int = 1) -> None:
    if "gurobi" in model.solver.interface.__name__:
        model.solver.problem.Params.Threads = threads


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_state() -> tuple[str, bool]:
    root = paths.PROJECT_ROOT
    head = subprocess.run(["git", "-C", str(root), "rev-parse", "--short", "HEAD"], capture_output=True,
                          text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain", "--", "src", "scripts", "tests"],
                           capture_output=True, text=True).stdout.strip()
    return head, bool(dirty)


class Manifest:
    """Context manager running a W2 script under ``nextgem.run.provenance`` (row in results/w2/PROVENANCE.tsv)."""

    def __init__(self, script: str, solver: str = "gurobi", workers: int = 1, inputs=()):
        self.script, self.solver, self.workers, self.inputs = script, solver, workers, list(inputs)
        self.extra: dict = {}

    def __enter__(self):
        from nextgem.run import provenance

        self._cm = provenance(w2_dir() / "PROVENANCE.tsv", inputs=self.inputs, solver=self.solver,
                              workers=self.workers, allow_dirty=True, author="W2 engineer agent")
        self.prov = self._cm.__enter__()
        self.t0 = time.perf_counter()
        return self

    def add(self, *outs: Path) -> None:
        self.prov.output(*outs)

    def add_input(self, *ins: Path) -> None:
        self.prov.input(*ins)

    def __exit__(self, exc_type, exc, tb):
        if exc_type is None:
            self.prov.note(**self.extra)
        r = self._cm.__exit__(exc_type, exc, tb)
        if exc_type is None:
            print(f"[provenance] {self.script}: wall {time.perf_counter() - self.t0:.1f} s")
        return r


def atomic_write_tsv(df: pd.DataFrame, path: Path, **kw) -> Path:
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, **kw)
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    tmp.replace(path)
    return path
