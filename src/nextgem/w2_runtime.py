"""Run helpers for the W2 scripts: parent loading, solver threads, and a runtime manifest (DESIGN §6.7, §11.4).

``nextgem.run.provenance()`` (DESIGN §11.4) does not exist yet, and it would refuse to run on the uncommitted W2
code anyway. Until it exists, every W2 script appends one row per run to ``results/w2/MANIFEST.tsv`` with argv,
commit (and whether src/scripts/tests were dirty), solver, workers, wall and CPU time, peak RSS and output SHA-256.
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


def load_parent(name: str, solver: str = "gurobi"):
    """Load a parent release (tagged) with the given solver, single-threaded."""
    import cobra

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
    """Context manager: times a script and appends a row to results/w2/MANIFEST.tsv."""

    def __init__(self, script: str, solver: str = "gurobi", workers: int = 1):
        self.script, self.solver, self.workers = script, solver, workers
        self.outputs: list[Path] = []
        self.extra: dict = {}

    def __enter__(self):
        self.t0, self.c0 = time.perf_counter(), time.process_time()
        return self

    def add(self, *outs: Path) -> None:
        self.outputs.extend(Path(o) for o in outs)

    def __exit__(self, exc_type, exc, tb):
        if exc_type is not None:
            return False
        wall = time.perf_counter() - self.t0
        self_cpu = time.process_time() - self.c0
        child = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu = self_cpu + child.ru_utime + child.ru_stime
        rss_gb = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, child.ru_maxrss) / 1e6
        head, dirty = git_state()
        row = dict(time=time.strftime("%Y-%m-%d %H:%M:%S"), script=self.script, argv=" ".join(sys.argv[1:]),
                   commit=head, code_dirty=dirty, solver=self.solver, workers=self.workers,
                   wall_s=round(wall, 1), cpu_s=round(cpu, 1), peak_rss_gb=round(rss_gb, 2),
                   outputs=";".join(f"{o.relative_to(paths.PROJECT_ROOT)}:{sha256(o)[:16]}" for o in self.outputs
                                    if o.is_file()),
                   **self.extra)
        f = w2_dir() / "MANIFEST.tsv"
        df = pd.DataFrame([row])
        df.to_csv(f, sep="\t", index=False, mode="a", header=not f.exists())
        print(f"[manifest] {self.script}: wall {wall:.1f} s, cpu {cpu:.1f} s, rss {rss_gb:.2f} GB")
        return False


def atomic_write_tsv(df: pd.DataFrame, path: Path, **kw) -> Path:
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, **kw)
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    tmp.replace(path)
    return path
