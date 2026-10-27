"""Shared helpers for the W1 scripts (scripts/w1_*.py). Label-free.

``run_main`` runs a script's ``main`` under ``nextgem.run.provenance()`` (DESIGN §11.4), which appends one
row to ``results/w1/MANIFEST.tsv`` (commit, dirty flag, full input/output SHA-256, versions, seeds, solver
and version, wall/CPU time, peak RSS).
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT / "results" / "w1"
TMP_DIR = OUT_DIR / "tmp"  # git-ignored (results/**/tmp/)
MANIFEST = OUT_DIR / "MANIFEST.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def atomic_write_parquet(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_parquet(tmp, index=False, compression="zstd")
    os.replace(tmp, path)


def atomic_write_tsv(df: pd.DataFrame, path: Path, float_format: str | None = "%.6g") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, float_format=float_format)
    os.replace(tmp, path)


DEPMAP = PROJECT / "data" / "depmap"
BASELINE = PROJECT / "results" / "baseline"
AUTHOR = "W1 analyst (Claude Opus)"
_CURRENT: list = [None]


def current():
    """The active ``nextgem.run.Provenance`` of this script run (set by :func:`run_main`)."""
    if _CURRENT[0] is None:
        raise RuntimeError("no provenance context: run the script through run_main()")
    return _CURRENT[0]


def frozen_label_hashes() -> dict:
    """SHA-256 of the label tier files as recorded in FROZEN.tsv (label files are not hashed outside
    check_frozen.py; DESIGN §11.1). Recorded in the manifest 'extra' field of label-reading scripts."""
    out = {}
    for ln in (PROJECT / "FROZEN.tsv").read_text().splitlines()[1:]:
        f = ln.split("\t")
        if len(f) >= 4 and f[1] == "label_file":
            out[f[0]] = f[3]
    return out


def run_main(main, *, inputs=(), seeds=(), solver="", workers=1, reads_labels=False) -> int:
    """Run ``main()`` under ``nextgem.run.provenance`` (DESIGN §11.4); one row in results/w1/MANIFEST.tsv.
    ``main`` registers further inputs/outputs with ``current().input(...)`` / ``current().output(...)``.
    Uncommitted code is refused unless NEXTGEM_ALLOW_DIRTY=1 (used for the pre-commit M2 runs; the row then
    records dirty=True)."""
    from nextgem.run import provenance

    allow = os.environ.get("NEXTGEM_ALLOW_DIRTY") == "1"
    with provenance(MANIFEST, inputs=inputs, seeds=seeds, solver=solver, workers=workers, allow_dirty=allow,
                    author=AUTHOR) as prov:
        _CURRENT[0] = prov
        if reads_labels:
            prov.note(label_files_frozen_sha256=frozen_label_hashes())
        rc = main()
        if rc:
            raise SystemExit(rc)
    return 0
