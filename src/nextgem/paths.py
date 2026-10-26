"""All project paths (DESIGN.md §11.1: "All data paths come from nextgem.paths").

Two environment variables redirect the roots, so that toy tests and the leakage canary never touch
the real files:

- ``NEXTGEM_DATA_ROOT``  (DESIGN §11.1): replaces ``<project>/data``.
- ``NEXTGEM_STATE_ROOT`` (M1 infrastructure): replaces ``<project>`` as the root of the *state* the
  pipeline writes and checks: ``results/`` (ledgers, splits, label reports), ``FROZEN.tsv`` and
  ``logs/``. Code (``src/``, ``scripts/``) always comes from the project checkout.

Functions, not constants, so a changed environment is honoured at call time.
"""

from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
SRC_DIR = PROJECT_ROOT / "src" / "nextgem"


def data_root() -> Path:
    return Path(os.environ.get("NEXTGEM_DATA_ROOT") or PROJECT_ROOT / "data").resolve()


def state_root() -> Path:
    return Path(os.environ.get("NEXTGEM_STATE_ROOT") or PROJECT_ROOT).resolve()


def is_redirected() -> bool:
    """True if either root is redirected (toy test, canary or sandbox run)."""
    return data_root() != (PROJECT_ROOT / "data").resolve() or state_root() != PROJECT_ROOT


# ---- data (git-ignored)
def depmap_dir() -> Path:
    return data_root() / "depmap"


def locked_dir() -> Path:
    """Validation/test labels, fingerprints, sealed outputs. Hook- and tripwire-protected."""
    return data_root() / "locked"


def labels_dir() -> Path:
    """Training labels and label-free label metadata (ids, label-derived gene lists)."""
    return data_root() / "labels"


def sealed_dir(look: str) -> Path:
    return locked_dir() / "sealed" / look


# ---- state (committed unless noted)
def results_dir() -> Path:
    return state_root() / "results"


def splits_dir() -> Path:
    return results_dir() / "splits"


def baseline_dir() -> Path:
    return results_dir() / "baseline"


def labels_results_dir() -> Path:
    return results_dir() / "labels"


def frozen_tsv() -> Path:
    return state_root() / "FROZEN.tsv"


def design_md() -> Path:
    return state_root() / "DESIGN.md"


def logs_dir() -> Path:
    """Git-ignored."""
    return state_root() / "logs"


def ledger_path(name: str) -> Path:
    if name not in ("test", "validation"):
        raise ValueError(f"unknown ledger {name!r}")
    return results_dir() / f"{name}_ledger.tsv"


def authorisations_dir() -> Path:
    return results_dir() / "authorisations"


def test_runs_tsv() -> Path:
    return results_dir() / "test_runs.tsv"
