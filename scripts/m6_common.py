"""Shared definitions and pure helpers for the M6 integration factorial (DESIGN.md §10.1, §17 A11-A12;
RESEARCH_TREE D10). Label-free: nothing here opens a label file.

Cells (fixed in D10 before any factorial scoring). A cell is (parent, biomass variant, medium step):
  reference G*  = H1 / M0
  candidates    = (1) H1 + B-loo / M0, (2) H1 / M2(f=1), (3) H1 + B-loo / M2(f=1)   (one-SE order, simplest first)
  descriptive   = the H2 counterparts (H2 / M0 = generic H2, and the three H2 cells)
B-loo = W3 ``B-loo:cofactor_pool_biomass`` (spec results/w3/structure/variants_<parent>.json, applied with
``nextgem.w3.Variant``). M2(f=1) = W2 documented base medium + supplements + serum model at f = 1, gate-free, per
line as in results/w2/controls/line_sweeps.tsv (lines outside the attribution set get the integrated-pipeline medium
that W2 assigned there). M0 = Ham's (GR task inputs).

Score of a line under a cell: e = 1 − ratio_r of the KO sweep of (parent, variant, the line's medium signature),
mapped to the gene universe U with the W1 Entrez convention (absent gene -> 0; 6 decimals).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
M6 = PROJECT / "results" / "m6"
MANIFEST = M6 / "MANIFEST.tsv"
SWDIR = M6 / "sweeps"
SCDIR = M6 / "scores"
TRDIR = M6 / "train"
VALDIR = M6 / "val"
TMP = M6 / "tmp"  # git-ignored
AUTHOR = "M6 engineer agent (Claude Opus)"
MAX_WORKERS = 8

PARENTS = ("human1", "human2")
SHORT = {"human1": "H1", "human2": "H2"}
T0 = "T0"
BLOO = "B-loo:cofactor_pool_biomass"
VARIANTS = (T0, BLOO)
M0, M2 = "M0", "M2(f=1)"
STEPS = (M0, M2)

# A10 limits (DESIGN §17): 300 s per LP + one retry (nextgem.w3.solve_status), 1,800 s CPU per job (sweep);
# a timeout is a failed build.
JOB_CPU_LIMIT_S = 1800
GR_LB = 1.0  # GR task of metabolicTasks_Essential.txt: OUT biomass[e] LB 1 on Ham's inputs
FAIL_STATUSES = ("failed-timeout", "failed-nonoptimal", "failed-infeasible", "failed-wt")

TG_MIN_DELTA = 0.005  # §9 TG
TG_MIN_FOLDS = 9
FOLD_MIN_LINES = 5
CALL = 0.99  # e >= 0.99 <=> ratio_r < 0.01 (§6.1 fixed cutoff)
MEVALONATE_CORE = ("HMGCR", "MVK", "PMVK", "MVD")


@dataclass(frozen=True)
class Cell:
    key: str  # short id used in tables
    parent: str
    variant: str
    step: str
    role: str  # reference | candidate | generic_H2 | descriptive

    @property
    def method(self) -> str:
        """Method name in score files and the look spec (generic parents keep the harness names H1 / H2)."""
        if self.variant == T0 and self.step == M0:
            return SHORT[self.parent]
        return f"M6:{SHORT[self.parent]}{'+Bloo' if self.variant == BLOO else ''}/{self.step}"


def _cells() -> list[Cell]:
    out = []
    for p in PARENTS:
        s = SHORT[p]
        out.append(Cell(f"{s}", p, T0, M0, "reference" if p == "human1" else "generic_H2"))
        role = "candidate" if p == "human1" else "descriptive"
        out.append(Cell(f"{s}+Bloo/M0", p, BLOO, M0, role))
        out.append(Cell(f"{s}/M2", p, T0, M2, role))
        out.append(Cell(f"{s}+Bloo/M2", p, BLOO, M2, role))
    return out


CELLS = _cells()
CELL = {c.key: c for c in CELLS}
REFERENCE = CELL["H1"]
CANDIDATES = [CELL["H1/M2"], CELL["H1+Bloo/M0"], CELL["H1+Bloo/M2"]]  # DESIGN §9 order (medium < tasks < combination; D10 as corrected)
CANDIDATE_LABEL = {"H1/M2": "(1)", "H1+Bloo/M0": "(2)", "H1+Bloo/M2": "(3)"}


def sweep_key(variant: str, w2_sweep_id: str) -> str:
    """A sweep is (variant, W2 bounds signature); the W2 sweep id carries the parent ("human1:<signature>")."""
    return f"{variant}|{w2_sweep_id}"


def line_sweep_ids(line_sweeps: pd.DataFrame, parent: str, step: str) -> pd.Series:
    """ModelID -> W2 sweep id of the line's medium at ``step`` (M0 or M2(f=1), gate-free) for ``parent``."""
    s = line_sweeps[(line_sweeps.parent == parent) & (line_sweeps.ladder_step == step)]
    if s.ModelID.duplicated().any():
        raise ValueError(f"duplicate ModelID for {parent} {step}")
    return s.set_index("ModelID").sweep_id


def cell_sweep_keys(cell: Cell, line_sweeps: pd.DataFrame) -> pd.Series:
    """ModelID -> sweep key of ``cell``."""
    return line_sweep_ids(line_sweeps, cell.parent, cell.step).map(lambda s: sweep_key(cell.variant, s))


def cell_matrix(cell: Cell, line_sweeps: pd.DataFrame, E_sweeps: pd.DataFrame, lines) -> pd.DataFrame:
    """lines x U score matrix of ``cell`` (E_sweeps: sweep key x U e-matrix of the cell's parent)."""
    keys = cell_sweep_keys(cell, line_sweeps).reindex(list(lines))
    if keys.isna().any():
        raise KeyError(f"{int(keys.isna().sum())} lines have no {cell.step} sweep for {cell.parent}")
    missing = sorted(set(keys) - set(E_sweeps.index))
    if missing:
        raise KeyError(f"{cell.key}: sweeps missing from the e-matrix: {missing[:3]}")
    M = E_sweeps.loc[keys.to_numpy()]
    M.index = list(lines)
    return M


def e_from_ratio(ratio: pd.DataFrame, emap: pd.Series, U) -> pd.DataFrame:
    """sweeps x U e = 1 − ratio_r (absent gene -> 0), 6 decimals; NaN rows (failed sweeps) stay NaN."""
    r = ratio.reindex(columns=emap.index)
    failed = r.isna().all(axis=1)
    e = (1.0 - r).fillna(0.0)
    e.columns = emap.to_numpy()
    e = e.reindex(columns=list(U), fill_value=0.0).round(6)
    e.loc[failed] = np.nan
    return e


def gr_fail(ratio: pd.Series, wt_growth: float, lb: float = GR_LB) -> pd.Series:
    """GR-task call from a biomass KO sweep: the KO cannot reach biomass >= ``lb`` (KO growth = ratio_r x WT)."""
    return (ratio * wt_growth) < lb


def interaction(lt0: pd.Series, lt1: pd.Series, lt2: pd.Series, lt3: pd.Series) -> pd.Series:
    """Per-line 2x2 interaction (B-loo x medium) on tAP_line: (c3 − c2) − (c1 − c0), lines finite in all four."""
    df = pd.concat([lt0, lt1, lt2, lt3], axis=1, keys=["c0", "c1", "c2", "c3"]).dropna()
    return (df.c3 - df.c2) - (df.c1 - df.c0)


def tg_pass(delta_lin: float, delta_lo: float, n_pos_folds: int) -> bool:
    """§9 TG: Δ̄_lin >= 0.005, CI lower bound > 0, >= 9 of 12 training folds positive."""
    return bool(np.isfinite(delta_lin) and delta_lin >= TG_MIN_DELTA and delta_lo > 0 and n_pos_folds >= TG_MIN_FOLDS)


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
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
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


def provenance(*, inputs=(), seeds=(), solver="", workers=1, allow_dirty=False):
    """``nextgem.run.provenance`` with the M6 manifest (results/m6/MANIFEST.tsv)."""
    from nextgem.run import provenance as _prov

    M6.mkdir(parents=True, exist_ok=True)
    return _prov(MANIFEST, inputs=inputs, seeds=seeds, solver=solver, workers=workers, allow_dirty=allow_dirty,
                 author=AUTHOR)
