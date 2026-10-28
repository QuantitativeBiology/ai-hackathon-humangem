"""F4 (post hoc / exploratory, DESIGN §17 A15; RESEARCH_TREE D19): line-specific constraints beyond the medium, on top
of P1. Shared definitions and pure helpers. Label-free: nothing here opens a label file.

P1 = Human1 + W3 B-loo:cofactor_pool_biomass on the line's documented medium M2(f=1) (the M6 cell H1+Bloo/M2; line ->
medium signature from results/m6/scores/line_assignment.tsv, asserted equal to W2's results/w2/controls/line_sweeps.tsv).

Two arms, fixed by D19 before any F4 build or score:

* F4a (genomic loss): per line, the pre-registered DESIGN §W4 loss rule (``nextgem.context_inputs.lost_genes``, rule
  r2-11: damaging code 2; or code 1 and (relative CN < 0.5 or log2(TPM+1) < 1); or relative CN < 0.3; lines without CN
  use mutations only, flagged) gives the lost model genes. The order-independent drop (``nextgem.context.accept_losses``)
  runs on the P1 parent with the line's M2(f=1) medium; the accepted losses are made non-functional (their GPR closure
  ``context.loss_closure`` is closed, never re-opened) and the KO sweep is re-run with the A4 KO gene-state convention
  (a lost gene is non-functional everywhere: a KO of its isozyme partner closes the shared reaction; the KO of a lost
  gene changes nothing). The generic parent has nothing closed by extraction, so the repair rule has nothing to repair.
  Lines without an accepted loss equal P1.
* F4b (E-Flux on P1): W4 family F4 with Stage A's carried configuration F4|R1|rna|real|none (E-Flux, global TPM-1
  threshold, real RNA, no genomic constraint; ``nextgem.context``), built on the P1 parent with the line's M2(f=1)
  medium. A12(b)-style pre-check first on the 30 Stage-A panel lines (nested panel LOLO, Δ̄ vs P1 >= +0.005 to proceed).
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
F4 = PROJECT / "results" / "f4"
MANIFEST = F4 / "MANIFEST.tsv"
LOSSDIR = F4 / "losses"
SWDIR = F4 / "sweeps"
TRDIR = F4 / "train"
EFDIR = F4 / "eflux_panel"
TMP = F4 / "tmp"  # git-ignored
AUTHOR = "F4 engineer agent (Claude Opus; A15 post hoc)"
MAX_WORKERS = 12  # GLPK workers (the F2 agent runs concurrently)
PARENT = "human1"
BLOO = "B-loo:cofactor_pool_biomass"
M2 = "M2(f=1)"
JOB_CPU_LIMIT_S = 1800  # A10
LOSS_RULE = "r2-11"

P1_METHOD = "M6:H1+Bloo/M2(f=1)"
F4A_METHOD = "F4a:H1+Bloo+loss/M2(f=1)"
F4B_METHOD = "F4b:H1+Bloo+EFlux_R1/M2(f=1)"
F4B_CONFIG = ("F4", "R1", "rna", "real", "none")  # Stage A's carried F4 configuration (train_selection.json)

TG_MIN_DELTA = 0.005
TG_MIN_FOLDS = 9
PANEL_MIN_DELTA = 0.005  # D19 / A12(b) mirror: F4b proceeds to all training lines only if Δ̄ vs P1 >= +0.005


# ============================================================================= pure helpers
def tg_pass(delta_lin: float, delta_lo: float, n_pos_folds: int) -> bool:
    """§9 TG: Δ̄_lin >= 0.005, CI lower bound > 0, >= 9 of 12 training folds positive."""
    return bool(np.isfinite(delta_lin) and delta_lin >= TG_MIN_DELTA and delta_lo > 0 and n_pos_folds >= TG_MIN_FOLDS)


def panel_proceeds(delta_lin: float) -> bool:
    """F4b pre-check: proceed to all training lines iff the panel Δ̄_lin vs P1 is >= +0.005 (NaN never proceeds)."""
    return bool(np.isfinite(delta_lin) and delta_lin >= PANEL_MIN_DELTA - 1e-12)


def join_genes(genes) -> str:
    return ";".join(sorted(genes))


def split_genes(s) -> list[str]:
    if s is None or (isinstance(s, float) and np.isnan(s)) or s == "":
        return []
    return [g for g in str(s).split(";") if g]


def job_key(sig: str, accepted) -> str:
    """Sweep job id of an F4a line state: (M2(f=1) signature, accepted-loss set), order-independent."""
    acc = join_genes(accepted)
    return f"{sig}|{hashlib.sha1(acc.encode()).hexdigest()[:12]}"


def sweep_mode(n_closed: int) -> str:
    """'full' if the accepted losses close >= 1 reaction (the wild-type LP changes, every gene can move); else
    'partial' (bounds unchanged: only genes whose KO closure changes in the lost state need re-solving)."""
    return "full" if n_closed > 0 else "partial"


# ============================================================================= line sets and media (label-free)
def scored_training_lines(split: pd.DataFrame, lf: pd.DataFrame) -> list[str]:
    """The W1/M6/F1 primary training set: training tier, not the rare-line pool, with an expression matrix."""
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    return sorted(l for l in elig.index if lf.loc[l, "has_expression_matrix"])


def line_signatures() -> pd.Series:
    """ModelID -> Human1 M2(f=1) W2 sweep id, from results/m6/scores/line_assignment.tsv, asserted equal to W2's
    line_sweeps.tsv mapping (the one M6 and F1 used)."""
    la = pd.read_csv(PROJECT / "results" / "m6" / "scores" / "line_assignment.tsv", sep="\t")
    a = la[(la.parent == PARENT) & (la.step == M2)].set_index("ModelID").w2_sweep_id
    ls = pd.read_csv(PROJECT / "results" / "w2" / "controls" / "line_sweeps.tsv", sep="\t")
    b = ls[(ls.parent == PARENT) & (ls.ladder_step == M2)].set_index("ModelID").sweep_id
    if a.index.duplicated().any() or not a.sort_index().equals(b.reindex(a.index).sort_index()):
        raise RuntimeError("line_assignment.tsv and line_sweeps.tsv disagree on the Human1 M2(f=1) signatures")
    return a


def signature_bounds(sigs) -> tuple[dict, dict]:
    """W2 sweep id -> exchange bounds (rebuilt from the frozen media library, as M6 and F1 did) and -> info."""
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from m6_sweeps import bounds_by_id

    return bounds_by_id({PARENT: {M2: sorted(set(sigs))}})


def variant():
    """The W3 B-loo:cofactor_pool_biomass spec exactly as swept and validated (nextgem.w3.Variant)."""
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from m6_sweeps import variant_json
    from nextgem import w3

    return w3.Variant.from_json(variant_json(PARENT, BLOO))


def variant_json() -> dict:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from m6_sweeps import variant_json as vj

    return vj(PARENT, BLOO)


def p1_parent(bounds: dict, sig: str, solver: str = "glpk"):
    """The P1 parent as a ``nextgem.context.Parent``: Human1 (loaded under ``solver``), B-loo applied permanently, the
    signature's M2(f=1) bounds as medium, biomass objective, A10 LP limits (``context.set_solver_params``)."""
    import cobra

    from nextgem import context as cx
    from nextgem.model import GEMS, load_model

    cobra.Configuration().solver = solver
    m = load_model(PARENT)
    if cx.solver_name(m) != solver:
        raise RuntimeError(f"parent loaded with {cx.solver_name(m)}, expected {solver}")
    variant().apply(m)
    return cx.make_parent(m, f"{PARENT}+{BLOO}", medium=bounds, medium_name=sig, genes_tsv=GEMS[PARENT].genes_tsv)


# ============================================================================= io helpers
def atomic_write_tsv(df: pd.DataFrame, path: Path, float_format: str | None = "%.6g") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, float_format=float_format)
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


def atomic_write_json(obj, path: Path) -> Path:
    return atomic_write_text(json.dumps(obj, indent=1, default=float) + "\n", path)


def provenance(*, inputs=(), seeds=(), solver="", workers=1, allow_dirty=False, manifest: Path | None = None):
    from nextgem.run import provenance as _prov

    F4.mkdir(parents=True, exist_ok=True)
    return _prov(manifest or MANIFEST, inputs=inputs, seeds=seeds, solver=solver, workers=workers,
                 allow_dirty=allow_dirty, author=AUTHOR)


def in_tmp(p: Path) -> bool:
    return "tmp" in Path(p).resolve().relative_to(PROJECT.resolve()).parts


def blas_single_thread() -> None:
    """Pin BLAS/OpenMP to one thread (inherited by spawned workers)."""
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[k] = "1"
