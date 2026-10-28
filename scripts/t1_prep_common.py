"""Shared definitions and pure helpers for the T1 preparation (scripts/t1_prep_*.py; DESIGN §6.5, §6.8, §10.1-10.3;
RESEARCH_TREE D11). Label-free: nothing here opens a label file.

P1 (D11) = Human1 + W3 B-loo:cofactor_pool_biomass on documented medium M2(f=1), gate-free, generic parent, no
context family = M6 cell (3) ``M6:H1+Bloo/M2(f=1)``; P1(H2) = the same pipeline on Human2 (``M6:H2+Bloo/M2(f=1)``).
The attribution arm (ii) "medium alone" (M0 -> best gate-free documented variant, generic parent) = M6 cell
``M6:H1/M2(f=1)``. Curation (W5) and context (M5) were not selected (DESIGN §10.1, A12c), so (iii) and (v) are
"lever not selected" and (iv) "constraints alone" = P1 itself vs G* (P1 has generic GPRs).

Outputs live in results/m6/t1/ (tracked: validation-keyed and gene-keyed files) and results/m6/tmp/t1/ (git-ignored:
test-line-keyed files, regenerable by scripts/t1_prep_build.py and pinned by SHA-256 in the P1 spec).
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as mc  # noqa: E402

PROJECT = mc.PROJECT
M6 = mc.M6
T1 = M6 / "t1"  # tracked
T1_TMP = M6 / "tmp" / "t1"  # git-ignored (test-line-keyed files)
MANIFEST = T1 / "MANIFEST.tsv"
AUTHOR = "T1-preparation engineer agent (Claude Opus)"
MAX_WORKERS = 8
SPEC_PATH = M6 / "p1_spec.json"

# ---- method names (score columns of the P1 spec)
P1_H1, P1_H2, MED_H1 = "P1[H1]", "P1[H2]", "M2[H1]"
SOURCE_COLUMN = {  # spec method -> column of results/m6/scores/scores_<tier>.large.parquet (M6 cell method names)
    "H1": "H1",
    "H2": "H2",
    P1_H1: mc.CELL["H1+Bloo/M2"].method,
    P1_H2: mc.CELL["H2+Bloo/M2"].method,
    MED_H1: mc.CELL["H1/M2"].method,
}
METHODS = ["H1", "H2", P1_H1, P1_H2, MED_H1]
ARMS = {P1_H1: "H1", P1_H2: "H2"}
ARM_PARENT = {P1_H1: "human1", P1_H2: "human2"}
ARM_CELL = {P1_H1: mc.CELL["H1+Bloo/M2"], P1_H2: mc.CELL["H2+Bloo/M2"]}
PRIMARY = P1_H1
LEVER_NOT_SELECTED = "lever not selected"
ATTRIBUTION = {
    "i": [PRIMARY, "G*"],  # P1 vs G* (fixed by the harness)
    "ii": [MED_H1, "H1"],  # medium alone: M0 -> M2(f=1) on the generic parent (H1 = G*)
    "iii": LEVER_NOT_SELECTED,  # curation alone (W5 not selected)
    "iv": [PRIMARY, "G*"],  # constraints alone = P1 with generic GPRs vs G*; P1 has generic GPRs -> P1 itself
    "v": LEVER_NOT_SELECTED,  # curation given constraints (W5 not selected)
}
FBA_CUTOFF = 0.01  # §6.1 fixed cutoff: ratio_r < 0.01
RATIO_DECIMALS = 6


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(p: Path) -> str:
    return str(Path(p).resolve().relative_to(PROJECT))


# ============================================================================= pure helpers (tested)
def fba_call_from_e(e) -> np.ndarray:
    """§6.1 FBA<0.01 call from a score e = 1 - ratio_r (6 decimals): ratio_r < 0.01. Absent genes (e = 0) -> False."""
    r = np.round(1.0 - np.asarray(e, dtype=float), RATIO_DECIMALS)
    return r < FBA_CUTOFF


def rename_scores(df: pd.DataFrame, lines) -> pd.DataFrame:
    """Subset an M6 line_gene score table to ``lines`` and the spec methods (renamed from the M6 cell names)."""
    keep = ["ModelID", "entrez"] + [SOURCE_COLUMN[m] for m in METHODS]
    out = df.loc[df.ModelID.isin(set(lines)), keep].rename(columns={v: k for k, v in SOURCE_COLUMN.items()})
    out = out.sort_values(["ModelID", "entrez"]).reset_index(drop=True)
    out["entrez"] = out["entrez"].astype("int64")
    for m in METHODS:
        out[m] = out[m].astype("float64")
    return out


def p57_pgr(task_long: pd.DataFrame, gr_task: str, U, emap: pd.Series) -> pd.DataFrame:
    """Per U gene: p57 = the KO makes >= 1 of the tasks infeasible; pgr = the GR task is infeasible.

    ``task_long``: (task, gene, feasible) over the parent's model genes; ``emap``: model gene id -> Entrez.
    Genes of U absent from the parent are False (§6.8)."""
    fails = task_long.loc[~task_long.feasible.astype(bool)]
    p57_genes = set(fails.gene)
    gr_genes = set(fails.gene[fails.task == gr_task])
    ent = emap.dropna()
    p57_e = {int(ent[g]) for g in p57_genes if g in ent.index}
    gr_e = {int(ent[g]) for g in gr_genes if g in ent.index}
    idx = pd.Index([int(g) for g in U], name="entrez")
    return pd.DataFrame({"p57": idx.isin(list(p57_e)), "pgr": idx.isin(list(gr_e))}, index=idx)


def calls_long(arm: str, gene_calls: pd.DataFrame, e_lines: pd.DataFrame) -> pd.DataFrame:
    """Long calls of one arm: (ModelID, entrez, method, p57, pgr, fba). p57/pgr are line-constant (``gene_calls``
    indexed by entrez); fba is line-specific from the arm's score matrix ``e_lines`` (lines x U)."""
    genes = list(e_lines.columns)
    gc = gene_calls.reindex(genes)
    if gc.isna().any().any():
        raise ValueError("gene calls do not cover every gene of the score matrix")
    n_l, n_g = e_lines.shape
    return pd.DataFrame({
        "ModelID": np.repeat(e_lines.index.to_numpy(), n_g),
        "entrez": np.tile(np.asarray(genes, dtype="int64"), n_l),
        "method": arm,
        "p57": np.tile(gc.p57.to_numpy(bool), n_l),
        "pgr": np.tile(gc.pgr.to_numpy(bool), n_l),
        "fba": fba_call_from_e(e_lines.to_numpy()).ravel(),
    })


def task_calls_table(task_long: pd.DataFrame, emap: pd.Series) -> pd.DataFrame:
    """(entrez int, task str, fails bool), one row per (Entrez-mapped model gene, task)."""
    ent = emap.dropna()
    t = task_long[task_long.gene.isin(ent.index)]
    out = pd.DataFrame({"entrez": t.gene.map(ent).astype("int64").to_numpy(), "task": t.task.astype(str).to_numpy(),
                        "fails": ~t.feasible.astype(bool).to_numpy()})
    if out.duplicated(["entrez", "task"]).any():
        raise ValueError("duplicate (entrez, task) rows: the Entrez map is not one-to-one")
    return out.sort_values(["task", "entrez"]).reset_index(drop=True)


def skew_decision(delta) -> dict:
    """§6.5: BCa for C1/C1b iff the training-tier distribution of line ΔtAP has |skew| > 1 (metrics rule)."""
    from scipy import stats

    from nextgem import metrics

    v = pd.Series(delta, dtype=float).dropna()
    bca = metrics.skew_triggers_bca(v)
    return {"skew": float(stats.skew(v.to_numpy())), "n_lines": int(v.size), "bca": bool(bca),
            "bootstrap_method": "bca" if bca else "percentile"}


# ============================================================================= I/O
def atomic_write_parquet(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_parquet(tmp, index=False, compression="zstd")
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    return path


def atomic_write_tsv(df: pd.DataFrame, path: Path, float_format: str | None = "%.6g") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, float_format=float_format)
    os.replace(tmp, path)
    return path


def atomic_write_json(obj, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=False, default=str) + "\n")
    os.replace(tmp, path)
    return path


def provenance(*, inputs=(), seeds=(), solver="", workers=1, allow_dirty=True):
    """``nextgem.run.provenance`` with the T1-prep manifest. Runs under 1 CPU-h may run from the working tree
    (allow_dirty, recorded as dirty=True); longer runs come from a pinned clean worktree."""
    from nextgem.run import provenance as _prov

    T1.mkdir(parents=True, exist_ok=True)
    return _prov(MANIFEST, inputs=inputs, seeds=seeds, solver=solver, workers=workers, allow_dirty=allow_dirty,
                 author=AUTHOR)
