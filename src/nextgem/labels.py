"""The single label door (DESIGN.md §5, §11.1). The only module that opens label-bearing files.

Raw sources (``<data root>/depmap``): CRISPRGeneEffect (Chronos), CRISPRGeneDependency (probability),
and, as label-derived gene lists (§11.2), CRISPRInferredCommonEssentials, AchillesCommonEssentialControls,
AchillesNonessentialControls. (ScreenGeneEffect is protected but not used.)

Files written once at M1 by ``scripts/build_labels.py`` (:func:`build`):

=====================================  ==========================================================
``<data>/labels/labels_train.parquet``  training-tier labels (long format, see ``COLUMNS``)
``<data>/locked/labels_val.parquet``    validation tier (sealed; tripwire: evaluate_val.py only)
``<data>/locked/labels_test.parquet``   test tier (sealed; tripwire: evaluate_test.py only)
``<data>/locked/fingerprints.npz``      salted row fingerprints for the content scanner (§11.1)
``<data>/labels/label_derived_lists.parquet``  DepMap common-essential/control lists (§11.2)
``<data>/labels/label_ids.json``        ModelIDs and gene header, no values (:func:`ids`)
``<state>/results/labels/dropcounts.tsv``     NaN accounting per gene and per library class (§5.1)
``<state>/results/labels/build_report.json``  counts only (ID-set check, tiers, NaN totals, kappa)
``<state>/FROZEN.tsv``                  SHA-256 of the three tier files; primary label (kappa rule)
=====================================  ==========================================================

Label variants (``config.LABELS``): ``L_prob`` = dependency > 0.5 (primary unless the kappa rule fires),
``L_ge05`` = gene effect < -0.5, ``L_ge06`` = gene effect < -0.6 (secondary, §6.8),
``L_ge10`` = gene effect < -1. A pair whose source value is NaN has label NA and is excluded for all
methods; the continuous ``gene_effect`` and ``dep_prob`` are kept for the §6.8 threshold grid and B2c.

Public API (what fitting code may receive): :func:`train`, :func:`fold`, :func:`folds`, :func:`ids`,
:func:`label_derived_list` (``purpose`` = reporting/annotation). CV goes through ``nextgem.cv.lolo``.
:func:`validation` and :func:`test` work only when the entry point is ``scripts/evaluate_val.py`` /
``scripts/evaluate_test.py`` (enforced by the tripwire, not by this module).
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from nextgem import config, paths

RAW_FILES = {"gene_effect": "CRISPRGeneEffect.csv", "dep_prob": "CRISPRGeneDependency.csv"}
LIST_FILES = {
    "depmap_inferred_common_essentials": "CRISPRInferredCommonEssentials.csv",
    "achilles_common_essential_controls": "AchillesCommonEssentialControls.csv",
    "achilles_nonessential_controls": "AchillesNonessentialControls.csv",
}
COLUMNS = ["ModelID", "entrez", "gene_effect", "dep_prob", *config.LABELS]
LABEL_NAMES = list(config.LABELS)
TIER_CATEGORY = {"train": "train", "validation": "val", "test": "test"}  # tripwire categories
GENE_COL = re.compile(r"^(?P<symbol>.+) \((?P<entrez>\d+)\)$")
MODEL_ID = re.compile(r"^ACH-\d+$")

# Fingerprint rules (scanner, §11.1a): source -> (decimals, lower, upper) on the value; cells outside
# (lower, upper) are not fingerprinted (probabilities near 0/1 are not distinctive).
FINGERPRINT_RULES = {"gene_effect": (2, -np.inf, np.inf), "dep_prob": (3, 0.05, 0.95)}


# ----------------------------------------------------------------------------- tripwire link
def _tw():
    tw = sys.modules.get("nextgem_tripwire")
    if tw is None or not tw.installed():
        raise RuntimeError(
            "nextgem label tripwire is not active in this interpreter. Use /data/tmpA/nextgem/.venv/bin/python "
            "(sitecustomize installed by scripts/install_tripwire.py); labels.py refuses to run without it."
        )
    return tw


# ----------------------------------------------------------------------------- file locations
def tier_file(tier: str) -> Path:
    if tier == "train":
        return paths.labels_dir() / "labels_train.parquet"
    if tier == "validation":
        return paths.locked_dir() / "labels_val.parquet"
    if tier == "test":
        return paths.locked_dir() / "labels_test.parquet"
    raise ValueError(tier)


def fingerprints_file() -> Path:
    return paths.locked_dir() / "fingerprints.npz"


def lists_file() -> Path:
    return paths.labels_dir() / "label_derived_lists.parquet"


def ids_file() -> Path:
    return paths.labels_dir() / "label_ids.json"


# ----------------------------------------------------------------------------- label-free inputs
def sha256_file(path: Path) -> str:
    """SHA-256 of an *unprotected* file (protected files: :func:`sha256_protected`)."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_split(verify: bool = True) -> pd.DataFrame:
    """``results/splits/lineage_split.tsv`` (label-free), checked against ``SPLIT_SHA256``."""
    f = paths.splits_dir() / "lineage_split.tsv"
    if not f.exists():
        raise FileNotFoundError(f"{f} does not exist: draw the split first (scripts/make_splits.py, DESIGN §4.1)")
    if verify:
        sha_f = paths.splits_dir() / "SPLIT_SHA256"
        recorded = sha_f.read_text().split()[0] if sha_f.exists() else ""
        if recorded != sha256_file(f):
            raise RuntimeError(f"{f} does not match {sha_f} (recorded {recorded[:12]!r}); the split is frozen")
    split = pd.read_csv(f, sep="\t", dtype={"ModelID": str})
    if not split.ModelID.is_unique:
        raise ValueError("split has duplicate ModelIDs")
    if not set(split.tier) <= set(config.TIERS):
        raise ValueError(f"unknown tiers {set(split.tier) - set(config.TIERS)}")
    return split


def gene_universe() -> pd.DataFrame:
    """U (DESIGN §3.3): Entrez IDs of H1 or H2 genes that are CRISPR columns, from the baseline maps.

    Columns: entrez (int), in_human1, in_human2. Joins use Entrez IDs only.
    """
    parts = []
    for m in ("human1", "human2"):
        f = paths.baseline_dir() / f"{m}_gene_mapping.tsv"
        df = pd.read_csv(f, sep="\t", dtype=str, usecols=["entrez", "depmap_column"])
        e = df.loc[df.depmap_column.notna(), "entrez"].dropna().astype(int).unique()
        parts.append(pd.DataFrame({"entrez": e, "model": m}))
    long = pd.concat(parts)
    u = long.assign(v=True).pivot_table(index="entrez", columns="model", values="v", aggfunc="any", fill_value=False)
    u = u.reindex(columns=["human1", "human2"], fill_value=False).astype(bool)
    return u.rename(columns={"human1": "in_human1", "human2": "in_human2"}).reset_index().sort_values("entrez", ignore_index=True)


def parse_gene_columns(columns) -> pd.DataFrame:
    """``SYMBOL (Entrez)`` headers -> DataFrame(column, symbol, entrez:int). Non-gene columns are dropped."""
    rows = [(c, m["symbol"], int(m["entrez"])) for c in columns if (m := GENE_COL.match(str(c)))]
    return pd.DataFrame(rows, columns=["column", "symbol", "entrez"])


# ----------------------------------------------------------------------------- private raw readers
def _raw_path(source: str) -> Path:
    return paths.depmap_dir() / RAW_FILES[source]


def _read_raw_header(source: str) -> list[str]:
    with _tw().scope("raw"):
        return list(pd.read_csv(_raw_path(source), nrows=0).columns)


def _read_raw_index(source: str) -> list[str]:
    with _tw().scope("raw"):
        return pd.read_csv(_raw_path(source), usecols=[0]).iloc[:, 0].astype(str).tolist()


def _read_raw(source: str, columns: list[str]) -> pd.DataFrame:
    """Rows = ModelID, columns = the requested gene columns (float64)."""
    with _tw().scope("raw"):
        head = list(pd.read_csv(_raw_path(source), nrows=0).columns)
        df = pd.read_csv(_raw_path(source), usecols=[head[0], *columns], index_col=0,
                         dtype={c: "float64" for c in columns})
    df.index = df.index.astype(str)
    df.index.name = "ModelID"
    if not df.index.is_unique:
        raise ValueError(f"{RAW_FILES[source]}: duplicate ModelIDs")
    return df


def _read_list(name: str) -> pd.DataFrame:
    """A DepMap gene list -> DataFrame(entrez, symbol); every cell matching ``SYMBOL (Entrez)`` is a member."""
    with _tw().scope("raw"):
        raw = pd.read_csv(paths.depmap_dir() / LIST_FILES[name], dtype=str)
    cells = pd.Series(raw.values.ravel()).dropna().astype(str)
    parsed = parse_gene_columns(cells)
    return parsed.drop(columns="column").drop_duplicates("entrez").sort_values("entrez", ignore_index=True)


# ----------------------------------------------------------------------------- helpers
def binarize(values: pd.Series | np.ndarray, op: str, threshold: float) -> pd.Series:
    """Label from a continuous source; NaN -> NA. Returns pandas nullable Int8."""
    v = pd.Series(values, dtype="float64")
    lab = (v > threshold) if op == ">" else (v < threshold)
    return lab.astype("Int8").mask(v.isna())


def add_labels(df: pd.DataFrame) -> pd.DataFrame:
    for name, (src, op, thr) in config.LABELS.items():
        df[name] = binarize(df[src].to_numpy(), op, thr).to_numpy()
    return df


def ge_threshold_label(df: pd.DataFrame, threshold: float) -> pd.Series:
    """§6.8 robustness grid: gene effect < threshold (NA where the gene effect is NaN)."""
    return binarize(df["gene_effect"].to_numpy(), "<", threshold)


def cohen_kappa(a, b) -> float:
    """Cohen's kappa of two binary vectors (pairs with NA in either are dropped)."""
    a = pd.Series(a, dtype="Float64")
    b = pd.Series(b, dtype="Float64")
    ok = a.notna().to_numpy() & b.notna().to_numpy()
    x = a.to_numpy(dtype=float, na_value=np.nan)[ok]
    y = b.to_numpy(dtype=float, na_value=np.nan)[ok]
    n = len(x)
    if n == 0:
        return float("nan")
    po = float(np.mean(x == y))
    p1, q1 = float(np.mean(x)), float(np.mean(y))
    pe = p1 * q1 + (1 - p1) * (1 - q1)
    return float("nan") if pe == 1 else (po - pe) / (1 - pe)


def ce_fit(df: pd.DataFrame, label: str, fraction: float = config.CE_FIT_FRACTION) -> set[int]:
    """CE_fit (§5.3): genes positive in >= ``fraction`` of the lines in ``df`` (the fitting lines)."""
    d = df[["entrez", label]].dropna()
    frac = d.groupby("entrez")[label].mean()
    return set(frac[frac >= fraction].index.astype(int))


def _parquet_bytes(df: pd.DataFrame) -> bytes:
    buf = io.BytesIO()
    df.to_parquet(buf, index=False, compression="zstd", engine="pyarrow")
    return buf.getvalue()


def _atomic_write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def fingerprint_keys(model_ids, entrez, values, source: str, salt: bytes) -> tuple[np.ndarray, np.ndarray]:
    """Salted 32-bit fingerprints of (ModelID, Entrez, rounded value) cells.

    Returns (hashes uint32, mask of cells that were fingerprintable under the source's rule).
    """
    dec, lo, hi = FINGERPRINT_RULES[source]
    v = np.asarray(values, dtype="float64")
    mask = np.isfinite(v) & (v > lo) & (v < hi)
    r = np.round(v[mask], dec) + 0.0  # + 0.0 turns -0.0 into 0.0
    mids = np.asarray(model_ids, dtype=object)[mask]
    ents = np.asarray(entrez)[mask]
    out = np.empty(len(r), dtype=np.uint32)
    fmt = f"{{}}|{{}}|{source}|{{:.{dec}f}}"
    for i, (m, e, x) in enumerate(zip(mids, ents, r)):
        out[i] = int.from_bytes(hashlib.blake2b(fmt.format(m, int(e), x).encode(), digest_size=4, key=salt).digest(), "little")
    return out, mask


# ----------------------------------------------------------------------------- build (M1, once)
def build(log=print) -> dict:
    """Build every label file. Called only by ``scripts/build_labels.py`` (the tripwire enforces it).

    Idempotent: if FROZEN.tsv already holds the tier-file hashes, the files are rebuilt in memory and
    must reproduce those hashes exactly; nothing is overwritten and the kappa rule is not re-logged.
    """
    from nextgem import frozen, ledger  # local: avoid import cycles

    tw = _tw()
    t0 = time.time()
    # review S5: redirected roots are for toy data only
    if paths.data_root() == (paths.PROJECT_ROOT / "data").resolve():
        if paths.state_root() != paths.PROJECT_ROOT:
            raise RuntimeError("the real data root needs the real state root (FROZEN.tsv, ledgers); refusing")
    elif any(tw.is_real_data(_raw_path(s)) for s in RAW_FILES):
        raise RuntimeError("NEXTGEM_DATA_ROOT is redirected but its raw label files are the real ones; refusing")
    split = load_split()
    U = gene_universe()

    # --- ID sets (DESIGN §3.2: CRISPRGeneDependency ID sets must equal CRISPRGeneEffect's)
    head = {s: _read_raw_header(s) for s in RAW_FILES}
    rows = {s: _read_raw_index(s) for s in RAW_FILES}
    gcols = {s: parse_gene_columns(head[s][1:]) for s in RAW_FILES}
    id_check = {
        "rows_equal": set(rows["gene_effect"]) == set(rows["dep_prob"]),
        "gene_columns_equal": set(gcols["gene_effect"].column) == set(gcols["dep_prob"].column),
        "rows_only_in_gene_effect": len(set(rows["gene_effect"]) - set(rows["dep_prob"])),
        "rows_only_in_dep_prob": len(set(rows["dep_prob"]) - set(rows["gene_effect"])),
        "cols_only_in_gene_effect": len(set(gcols["gene_effect"].column) - set(gcols["dep_prob"].column)),
        "cols_only_in_dep_prob": len(set(gcols["dep_prob"].column) - set(gcols["gene_effect"].column)),
    }
    if not (id_check["rows_equal"] and id_check["gene_columns_equal"]):
        raise RuntimeError(f"CRISPRGeneDependency and CRISPRGeneEffect ID sets differ (DESIGN 3.2): {id_check}")
    for s in RAW_FILES:
        if gcols[s].entrez.duplicated().any():
            raise ValueError(f"{RAW_FILES[s]}: duplicated Entrez IDs in the header")
        if not all(MODEL_ID.match(r) for r in rows[s]):
            raise ValueError(f"{RAW_FILES[s]}: first column is not a ModelID index")

    # --- Entrez join to U
    ge_cols = gcols["gene_effect"]
    in_u = ge_cols[ge_cols.entrez.isin(U.entrez)]
    u_missing = sorted(set(U.entrez) - set(ge_cols.entrez))
    if u_missing:
        raise ValueError(f"{len(u_missing)} U genes have no CRISPR column (U is defined on CRISPR columns)")
    dp_by_entrez = gcols["dep_prob"].set_index("entrez").column

    # --- lines: the split defines the cohort
    split_ids = set(split.ModelID)
    raw_ids = set(rows["gene_effect"])
    if raw_ids - split_ids:
        raise ValueError(f"{len(raw_ids - split_ids)} label rows are not in the split")
    lines_missing = sorted(split_ids - raw_ids)

    ge = _read_raw("gene_effect", in_u.column.tolist())
    ge.columns = in_u.entrez.to_numpy()
    dp_cols = [dp_by_entrez[e] for e in in_u.entrez if e in dp_by_entrez.index]
    dp = _read_raw("dep_prob", dp_cols) if dp_cols else pd.DataFrame(index=ge.index)
    dp.columns = [int(GENE_COL.match(c)["entrez"]) for c in dp.columns]
    dp = dp.reindex(index=ge.index, columns=ge.columns)

    long = pd.DataFrame({
        "ModelID": np.repeat(ge.index.to_numpy(), ge.shape[1]),
        "entrez": np.tile(ge.columns.to_numpy().astype("int32"), ge.shape[0]),
        "gene_effect": ge.to_numpy().ravel(),
        "dep_prob": dp.to_numpy().ravel(),
    })
    long = add_labels(long)
    long = long.sort_values(["ModelID", "entrez"], ignore_index=True)

    # --- NaN accounting (dropcounts.tsv)
    meta = split.set_index("ModelID")
    long_tier = long.ModelID.map(meta.tier)
    lib = long.ModelID.map(meta["library_class"]) if "library_class" in meta else pd.Series("unknown", index=long.index)
    drops = []
    for src in RAW_FILES:
        nan = long[src].isna()
        g = pd.DataFrame({"entrez": long.entrez, "nan": nan}).groupby("entrez").nan.agg(["size", "sum"])
        drops.append(pd.DataFrame({"scope": "gene", "key": g.index.astype(str), "tier": "all", "source": src,
                                   "n_pairs": g["size"].to_numpy(), "n_nan": g["sum"].to_numpy()}))
        lg = pd.DataFrame({"lib": lib.fillna("unknown"), "tier": long_tier, "nan": nan}).groupby(["lib", "tier"]).nan.agg(["size", "sum"])
        drops.append(pd.DataFrame({"scope": "library_class", "key": [k[0] for k in lg.index], "tier": [k[1] for k in lg.index],
                                   "source": src, "n_pairs": lg["size"].to_numpy(), "n_nan": lg["sum"].to_numpy()}))
        tg = pd.DataFrame({"tier": long_tier, "nan": nan}).groupby("tier").nan.agg(["size", "sum"])
        drops.append(pd.DataFrame({"scope": "tier_total", "key": "all", "tier": tg.index, "source": src,
                                   "n_pairs": tg["size"].to_numpy(), "n_nan": tg["sum"].to_numpy()}))
    drops.append(pd.DataFrame({"scope": "line_missing_from_labels", "key": lines_missing,
                               "tier": [meta.tier[m] for m in lines_missing], "source": "both", "n_pairs": 0, "n_nan": 0}))
    dropcounts = pd.concat(drops, ignore_index=True)
    dropcounts[["n_pairs", "n_nan"]] = dropcounts[["n_pairs", "n_nan"]].astype(int)

    # --- tier files (bytes first: hash, then compare with FROZEN.tsv or write)
    tier_bytes = {t: _parquet_bytes(long[long_tier.to_numpy() == t].reset_index(drop=True)[COLUMNS]) for t in config.TIERS}
    tier_sha = {t: hashlib.sha256(b).hexdigest() for t, b in tier_bytes.items()}
    frozen_rows = frozen.read_rows()
    frozen_keys = {f"labels.{t}": r for r in frozen_rows for t in config.TIERS if r["item"] == f"labels.{t}"}
    rebuild_check = None
    if frozen_keys:
        mism = [t for t in config.TIERS if frozen_keys.get(f"labels.{t}", {}).get("sha256_or_value") != tier_sha[t]]
        if mism:
            raise RuntimeError(f"labels are frozen (FROZEN.tsv) and this rebuild differs for tiers {mism}; refusing to overwrite")
        rebuild_check = "reproduced FROZEN.tsv hashes exactly"
        log("labels already frozen: rebuild reproduced every tier hash; nothing overwritten")
    else:
        for t in config.TIERS:
            _write_protected(tier_file(t), tier_bytes[t], TIER_CATEGORY[t])

    # --- label-derived lists, ids, fingerprints, hook-probe sentinel
    lists = pd.concat([_read_list(n).assign(list_name=n) for n in LIST_FILES], ignore_index=True)[["list_name", "entrez", "symbol"]]
    ids_doc = {"model_ids": sorted(raw_ids), "gene_columns": head["gene_effect"][1:],
               "U_entrez": [int(e) for e in in_u.entrez], "note": "identifiers only, no label values (labels.ids())"}
    if rebuild_check is None:
        _write_protected(lists_file(), _parquet_bytes(lists), "lists")
        _atomic_write_bytes(ids_file(), json.dumps(ids_doc, indent=0).encode())
        salt = os.urandom(16)
        fp = {"salt": np.frombuffer(salt, dtype=np.uint8)}
        # training and held-out (validation + test) cells are fingerprinted separately, so that training
        # cells cannot pad a held-out probe up to the hit thresholds (review M6)
        heldout = long_tier.isin(["validation", "test"]).to_numpy()
        for src in RAW_FILES:
            for part, sel in (("train", ~heldout), ("heldout", heldout)):
                sub = long[sel]
                h, _ = fingerprint_keys(sub.ModelID.to_numpy(), sub.entrez.to_numpy(), sub[src].to_numpy(), src, salt)
                fp[f"{src}_{part}"] = np.unique(h)
        buf = io.BytesIO()
        np.savez(buf, **fp)
        _write_protected(fingerprints_file(), buf.getvalue(), "fingerprints")
        probe = paths.locked_dir() / "HOOK_PROBE"
        if not probe.exists():
            _write_protected(probe, b"sentinel for the PreToolUse hook liveness probe (DESIGN 11.1); no data\n", "locked")

    # --- kappa rule (§5.1): training lines, pooled over all U pairs; logged before any metric
    tr = long[long_tier.to_numpy() == "train"]
    kappa = cohen_kappa(tr.L_prob, tr.L_ge05)
    ce = ce_fit(tr, "L_prob")
    tr_noce = tr[~tr.entrez.isin(ce)]
    kappa_noce = cohen_kappa(tr_noce.L_prob, tr_noce.L_ge05)
    primary = "L_ge05" if (np.isfinite(kappa) and kappa < config.KAPPA_RULE_THRESHOLD) else "L_prob"
    kappa_doc = {"rule": "kappa(L_prob, L_ge05) on training lines, pooled over U pairs; L_ge05 primary if < 0.5",
                 "kappa": round(kappa, 6), "kappa_without_ce_fit_info_only": round(kappa_noce, 6),
                 "n_pairs": int((tr.L_prob.notna() & tr.L_ge05.notna()).sum()), "primary_label": primary}
    if not any(r["item"] == "label.primary" for r in frozen_rows):
        ledger.append("validation", "non_look", look_id="kappa_rule", workstream="M1",
                      detail=json.dumps(kappa_doc, sort_keys=True), author="scripts/build_labels.py")
        frozen.upsert_rows([
            {"item": f"labels.{t}", "kind": "label_file", "path": str(tier_file(t).relative_to(paths.data_root())),
             "sha256_or_value": tier_sha[t], "note": "built by scripts/build_labels.py"} for t in config.TIERS
        ] + [{"item": "label.primary", "kind": "value", "path": "", "sha256_or_value": primary,
              "note": f"kappa rule (DESIGN 5.1); kappa={kappa_doc['kappa']}"}])

    # --- reports (counts only)
    out_dir = paths.labels_results_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    dc = out_dir / "dropcounts.tsv"
    dropcounts.to_csv(dc.with_suffix(".tmp"), sep="\t", index=False)
    os.replace(dc.with_suffix(".tmp"), dc)
    lists_u = lists[lists.entrez.isin(U.entrez)].groupby("list_name").size().to_dict()
    report = {
        "id_sets": id_check,
        "n_lines_labels": len(raw_ids), "n_lines_split": len(split_ids), "lines_missing_from_labels": len(lines_missing),
        "n_gene_columns": len(ge_cols), "n_U": int(len(U)), "n_U_in_crispr": int(len(in_u)),
        "lines_per_tier": {t: int((split.tier == t).sum()) for t in config.TIERS},
        "pairs_per_tier": {t: int((long_tier == t).sum()) for t in config.TIERS},
        "nan_pairs": {s: int(long[s].isna().sum()) for s in RAW_FILES},
        "label_derived_lists_in_U": {k: int(v) for k, v in lists_u.items()},
        "tier_sha256": tier_sha, "kappa_rule": kappa_doc, "rebuild_check": rebuild_check,
        "runtime_s": round(time.time() - t0, 1),
    }
    rp = out_dir / "build_report.json"
    rp.write_text(json.dumps(report, indent=1, sort_keys=True))
    return report


def _write_protected(path: Path, data: bytes, category: str) -> None:
    with _tw().scope(category):
        _atomic_write_bytes(path, data)


# ----------------------------------------------------------------------------- public API
def ids() -> dict:
    """ModelIDs and the gene header of CRISPRGeneEffect, without values (written by the M1 build)."""
    return json.loads(ids_file().read_text())


def _train_frame() -> pd.DataFrame:
    with _tw().scope("train"):
        df = pd.read_parquet(tier_file("train"))
    split = load_split()
    meta = split.set_index("ModelID")
    df["lineage"] = df.ModelID.map(meta.lineage)
    df["stratum"] = df.ModelID.map(meta.stratum)
    return df


def folds() -> list[str]:
    """LOLO folds: the eligible training lineages (rare lines are never a fold), sorted."""
    split = load_split()
    tr = split[(split.tier == "train") & (split.stratum != "rare")]
    return sorted(tr.lineage.unique())


def _fold_name(k) -> str:
    fl = folds()
    if isinstance(k, (int, np.integer)):
        return fl[int(k)]
    if k not in fl:
        raise KeyError(f"{k!r} is not a training fold")
    return k


def train(columns: list[str] | None = None) -> pd.DataFrame:
    """All training-tier labels (long format + lineage, stratum)."""
    df = _train_frame()
    return df if columns is None else df[["ModelID", "entrez", "lineage", *columns]]


def fold(k, inner=None, columns: list[str] | None = None) -> pd.DataFrame:
    """Training labels minus fold ``k`` (and minus inner fold ``inner`` for nested LOLO)."""
    out = {_fold_name(k)} | ({_fold_name(inner)} if inner is not None else set())
    df = train(columns)
    return df[~df.lineage.isin(out)].reset_index(drop=True)


def _heldout(k, inner=None, columns: list[str] | None = None) -> pd.DataFrame:
    """Labels of held-out fold ``k`` (or inner fold within training-minus-k). For ``nextgem.cv`` only."""
    name = _fold_name(inner) if inner is not None else _fold_name(k)
    df = train(columns)
    return df[df.lineage == name].reset_index(drop=True)


def validation() -> pd.DataFrame:
    """Validation-tier labels. The tripwire allows this only from scripts/evaluate_val.py."""
    with _tw().scope("val"):
        return pd.read_parquet(tier_file("validation"))


def test() -> pd.DataFrame:
    """Test-tier labels. The tripwire allows this only from scripts/evaluate_test.py."""
    with _tw().scope("test"):
        return pd.read_parquet(tier_file("test"))


def label_derived_list(name: str, purpose: str) -> set[int]:
    """DepMap common-essential/control lists (§11.2): reporting/annotation only, never fitting."""
    if purpose not in config.LABEL_DERIVED_PURPOSES:
        raise PermissionError(f"label-derived source {name!r} requires purpose in {config.LABEL_DERIVED_PURPOSES}")
    with _tw().scope("lists"):
        df = pd.read_parquet(lists_file())
    if name not in set(df.list_name):
        raise KeyError(name)
    return set(df.loc[df.list_name == name, "entrez"].astype(int))


def sha256_protected(path: Path | str) -> str:
    """Hash a protected file (FROZEN.tsv / ledger checks). Tripwire: check_frozen.py, verify_ledgers.py."""
    with _tw().scope("hash"):
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()


def load_fingerprints() -> dict | None:
    """Fingerprints for the content scanner (keys ``<source>_train`` and ``<source>_heldout``).
    Sealed-equivalent: with the salt they are brute-forceable, so the tripwire lets only a frozen
    scan_label_copies.py read them, and the scanner never prints match counts. None if not built."""
    f = fingerprints_file()
    if not f.exists():
        return None
    with _tw().scope("fingerprints"):
        with np.load(f) as z:
            return {"salt": z["salt"].tobytes(), **{k: z[k] for k in z.files if k != "salt"}}
