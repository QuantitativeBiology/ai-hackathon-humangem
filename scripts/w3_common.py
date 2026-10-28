"""Shared helpers for the W3 scripts (scripts/w3_*.py).

``run_main`` runs a script's ``main`` under ``nextgem.run.provenance()`` (DESIGN §11.4), appending one row to
``results/w3/MANIFEST.tsv`` (argv, commit, dirty flag, full SHA-256 of inputs and outputs, versions, seeds, solver
and its version, workers, wall and CPU time including waited-for children, peak RSS, author).

The W3 code is run uncommitted (only the PI commits), so ``allow_dirty=True`` is passed and every row records
``dirty=True`` (same convention as W2's ``w2_runtime.Manifest``).

Solver rule (brief): every script sets ``cobra.Configuration().solver`` BEFORE loading a model. GLPK is the default;
Gurobi is used only through ``nextgem.gurobi_pool`` (W3 holds at most ``W3_GUROBI_SLOTS`` = 3 slots).
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT / "results" / "w3"
TMP_DIR = OUT_DIR / "tmp"  # git-ignored (results/**/tmp/)
MANIFEST = OUT_DIR / "MANIFEST.tsv"
AUTHOR = "W3 engineer (Claude Opus)"
W3_GUROBI_SLOTS = 3
MAX_WORKERS = 8
_CURRENT: list = [None]


def current():
    if _CURRENT[0] is None:
        raise RuntimeError("no provenance context: run the script through run_main()")
    return _CURRENT[0]


def atomic_write_tsv(df: pd.DataFrame, path: Path, float_format: str | None = "%.6g") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, float_format=float_format)
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    return path


def atomic_write_parquet(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_parquet(tmp, index=False, compression="zstd")
    os.replace(tmp, path)
    return path


def atomic_write_text(text: str, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)
    return path


def frozen_label_hashes() -> dict:
    out = {}
    for ln in (PROJECT / "FROZEN.tsv").read_text().splitlines()[1:]:
        f = ln.split("\t")
        if len(f) >= 4 and f[1] == "label_file":
            out[f[0]] = f[3]
    return out


def run_main(main, *, inputs=(), seeds=(), solver="glpk", workers=1, reads_labels=False) -> int:
    from nextgem.run import provenance

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with provenance(MANIFEST, inputs=inputs, seeds=seeds, solver=solver, workers=workers, allow_dirty=True,
                    author=AUTHOR) as prov:
        _CURRENT[0] = prov
        if reads_labels:
            prov.note(label_files_frozen_sha256=frozen_label_hashes(), label_api="nextgem.labels.train()/fold(k) only")
        rc = main()
        if rc:
            raise SystemExit(rc)
    return 0


def glpk_parent(name: str):
    """Load a parent on GLPK (no Gurobi session), single-threaded, tagged with its release."""
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem.model import load_model

    return load_model(name)
