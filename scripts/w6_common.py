"""Shared helpers for the W6 scripts (scripts/w6_*.py).

``run_main`` runs a script's ``main`` under ``nextgem.run.provenance()`` (DESIGN §11.4), appending one row to
``results/w6/MANIFEST.tsv``. Production runs refuse a dirty tree (allow_dirty=False, the §11.4 default); the
``--allow-dirty`` flag exists only for label-free pilots and is recorded in the row (``dirty=True``, note ``pilot``).

Solver rule: every script sets ``cobra.Configuration().solver = "glpk"`` before loading a model; W6 never opens a
Gurobi session (GLPK is preferred; Gurobi would go only through ``nextgem.gurobi_pool``).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT / "results" / "w6"
TMP_DIR = OUT_DIR / "tmp"  # git-ignored (results/**/tmp/)
MANIFEST = OUT_DIR / "MANIFEST.tsv"
CONFIG_DIR = OUT_DIR / "config"
AUTHOR = "W6 engineer (Claude Opus)"
MAX_WORKERS = 6
_CURRENT: list = [None]


def current():
    if _CURRENT[0] is None:
        raise RuntimeError("no provenance context: run the script through run_main()")
    return _CURRENT[0]


def atomic_write_tsv(df: pd.DataFrame, path: Path, float_format: str | None = "%.6g", header_lines=()) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w") as fh:
        for h in header_lines:
            fh.write(f"# {h}\n")
        df.to_csv(fh, sep="\t", index=False, float_format=float_format)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    current_out(path)
    return path


def atomic_write_parquet(df: pd.DataFrame, path: Path, register: bool = True) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_parquet(tmp, compression="zstd")
    os.replace(tmp, path)
    if register:
        current_out(path)
    return path


def atomic_write_json(obj, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=True, default=str) + "\n")
    os.replace(tmp, path)
    current_out(path)
    return path


def current_out(path: Path) -> None:
    if _CURRENT[0] is not None and TMP_DIR not in Path(path).resolve().parents:
        _CURRENT[0].output(path)


def run_main(main, *, inputs=(), seeds=(), solver="glpk", workers=1, allow_dirty=False, budget_cpu_h=None,
             note: dict | None = None) -> int:
    import nextgem
    from nextgem.run import provenance

    if PROJECT / "src" not in Path(nextgem.__file__).resolve().parents:
        raise RuntimeError(f"imported nextgem ({nextgem.__file__}) is not the package next to this script; in a pinned "
                           "worktree run with PYTHONPATH=<worktree>/src (DESIGN §11.4)")
    from threadpoolctl import threadpool_limits

    threadpool_limits(1)  # BLAS/OpenMP: one thread per process (also set through the environment for workers)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with provenance(MANIFEST, inputs=inputs, seeds=seeds, solver=solver, workers=workers, allow_dirty=allow_dirty,
                    author=AUTHOR) as prov:
        _CURRENT[0] = prov
        prov.note(labels="none (W6 reads no CRISPR label and no label-derived source)", budget_cpu_h=budget_cpu_h,
                  **(note or {}))
        rc = main()
        if rc:
            raise SystemExit(rc)
    return 0


def load_config(name_or_path: str):
    from nextgem.w6 import RunConfig

    p = Path(name_or_path)
    if not p.exists():
        p = CONFIG_DIR / f"{name_or_path}.json"
    return RunConfig.from_json(json.loads(p.read_text()), PROJECT), p


def run_dir(tag: str) -> Path:
    return OUT_DIR / tag


def tmp_dir(tag: str | None = None) -> Path:
    return TMP_DIR / tag if tag else TMP_DIR


class CPUWatchdog:
    """Kills the job at 2x its declared CPU budget (DESIGN §11.5): counts this process plus finished worker tasks."""

    def __init__(self, budget_cpu_h: float | None):
        import time

        self.limit_s = None if budget_cpu_h is None else 2 * budget_cpu_h * 3600
        self.t0 = time.process_time()
        self.worker_s = 0.0

    def add(self, cpu_s: float) -> None:
        import time

        self.worker_s += cpu_s
        total = (time.process_time() - self.t0) + self.worker_s
        if self.limit_s is not None and total > self.limit_s:
            raise RuntimeError(f"CPU budget exceeded: {total / 3600:.2f} CPU-h > 2 x declared budget")
