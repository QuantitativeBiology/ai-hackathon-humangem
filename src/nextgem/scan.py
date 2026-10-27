"""Content-based label-copy scanner (DESIGN.md §11.1, PI decision r2-14). Used by scripts/scan_label_copies.py.

A file is a **hit** if any of:

(a) *fingerprints*: its (ModelID, Entrez, value) cells match the salted fingerprints written by
    ``labels.build`` (gene effect rounded to 2 dp; dependency probability in (0.05, 0.95) rounded to
    3 dp; values and their negation). Held-out (validation + test) and training cells are matched
    **separately** (review M6): a hit needs >= 5 held-out matches and >= 50% of the held-out-keyed
    cells, or >= 20 training matches and >= 50% of the training-keyed cells, so training cells cannot
    pad a probe. Works for wide (either orientation) and long tables, subsets, any file name, every row.
(b) *shape*: ModelID keys, >= 50 ``SYMBOL (Entrez)`` gene columns, and either >= 1% of the values
    < -0.3 (Chronos-like) or a 0/1 matrix keyed by validation or test ModelIDs (any ModelID if no split
    exists yet); unkeyed 2-D arrays (``.npy``/``.npz``) with >= 10 rows x >= 50 columns and >= 1% < -0.3.
(c) *name*: a protected label file name outside its canonical location.
(d) *binary container* (staged / pushed files under results/ only): ``.npy .npz .pkl .pickle .joblib
    .h5 .hdf5 .zip`` unless allow-listed in ``results/scanner_allowlist.tsv`` (path, sha256).

``test_keyed`` (pre-commit, until the T2 unseal ledger row exists): any table keyed by test ModelIDs x
genes.

Formats: csv/tsv/txt/tab (also .gz/.bz2/.xz/.zip), parquet, feather/arrow, json (also compressed),
npy/npz (``allow_pickle=False``), pickles through a restricted unpickler (pandas/numpy containers only).
ModelIDs are recognised with or without the ``ACH-`` prefix (bare zero-padded 6-digit strings, matched
against the known ID set).

**Output is coarse by design (review M6):** a path and a reason class; never values, match counts or
per-column detail, so the scanner cannot be used as an oracle.
"""

from __future__ import annotations

import io
import json
import os
import pickle
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from nextgem import labels as labels_mod

COMPRESSIONS = {".gz": "gzip", ".bz2": "bz2", ".xz": "xz", ".zip": "zip"}
TEXT_TABLES = (".csv", ".tsv", ".txt", ".tab")
TABLE_BASES = TEXT_TABLES + (".parquet", ".pq", ".feather", ".arrow", ".json", ".npy", ".npz", ".pkl", ".pickle")
BINARY_CONTAINERS = (".npy", ".npz", ".pkl", ".pickle", ".joblib", ".h5", ".hdf5", ".zip")
UNSCANNABLE = (".joblib", ".h5", ".hdf5", ".rds", ".rdata")
SKIP_DIRS = {".git", ".venv", "__pycache__", ".ipynb_checkpoints", "node_modules", ".mypy_cache", ".pytest_cache"}
MAX_FULL_READ_BYTES = 200 * 1024 * 1024  # files up to this size are read completely
MAX_ROWS_LARGE = 200_000  # larger files: first rows only
FP_MIN_HELDOUT, FP_MIN_TRAIN, FP_MIN_FRAC = 5, 20, 0.5
MAX_HELDOUT_CELLS, MAX_OTHER_CELLS = 500_000, 40_000
SHAPE_MIN_GENE_COLS = 50
SHAPE_NEG_THRESHOLD = -0.3
SHAPE_NEG_FRAC = 0.01
ACH = re.compile(r"^ACH-\d{6}$", re.I)
BARE = re.compile(r"^\d{6}$")
GENE_COL = re.compile(r"^.+ \((\d+)\)$")


@dataclass
class Result:
    path: str
    hits: list[str] = field(default_factory=list)
    note: str = ""


@dataclass
class Context:
    fp: dict | None = None
    heldout_ids: set[str] | None = None
    test_ids: set[str] | None = None
    known_ids: set[str] | None = None
    container_rule: bool = False
    canonical_split: Path | None = None  # results/splits/lineage_split.tsv of the state root
    canonical_split_sha: str | None = None  # pinned by SPLIT_SHA256 (and split_meta.json, if present)


def canonical_split_pin(splits_dir: Path) -> str | None:
    """SHA-256 of the canonical split file as pinned by make_splits.py: SPLIT_SHA256, which must agree with
    split_meta.json when that exists. None if there is no completed split (or the pins disagree)."""
    import json

    f = splits_dir / "SPLIT_SHA256"
    if not f.exists():
        return None
    sha = (f.read_text().split() or [""])[0]
    if not re.fullmatch(r"[0-9a-f]{64}", sha):
        return None  # PENDING or malformed
    meta = splits_dir / "split_meta.json"
    if meta.exists():
        try:
            if json.loads(meta.read_text()).get("split_sha256") != sha:
                return None
        except ValueError:
            return None
    return sha
    allowlist: dict[str, str] = field(default_factory=dict)


# ----------------------------------------------------------------------------- reading
def _split_suffix(name: str) -> tuple[str, str | None]:
    for ext, comp in COMPRESSIONS.items():
        if name.endswith(ext) and not name.endswith(".npz"):
            inner = name[: -len(ext)]
            if inner.endswith(TABLE_BASES):
                return inner, comp
    return name, None


class _SafeUnpickler(pickle.Unpickler):
    """Only pandas/numpy containers and plain builtins (review item 7: scan pickles without running code)."""

    SAFE_BUILTINS = {"list", "dict", "tuple", "set", "frozenset", "slice", "range", "complex", "bytearray",
                     "object", "int", "float", "str", "bool", "bytes"}
    BAD = re.compile(r"eval|exec|read|load|open|system|query|call|import|compile|getattr|apply", re.I)

    def find_class(self, module, name):
        if module == "builtins" and name in self.SAFE_BUILTINS:
            return super().find_class(module, name)
        if module.split(".")[0] in ("pandas", "numpy", "collections", "datetime", "copyreg", "_codecs") \
                and not self.BAD.search(name):
            return super().find_class(module, name)
        raise pickle.UnpicklingError(f"blocked {module}.{name}")


def _frames_from_object(obj) -> list[pd.DataFrame]:
    if isinstance(obj, pd.DataFrame):
        idx = obj.index
        if isinstance(idx, pd.RangeIndex) or (idx.name is None and pd.api.types.is_integer_dtype(idx.dtype)):
            return [obj.reset_index(drop=True)]  # positional index: not a key
        return [obj.reset_index()]
    if isinstance(obj, pd.Series):
        return _frames_from_object(obj.to_frame())
    if isinstance(obj, np.ndarray):
        return [pd.DataFrame(obj)] if obj.ndim == 2 else []
    if isinstance(obj, dict):
        out = []
        try:
            out.append(pd.DataFrame(obj))
        except Exception:
            pass
        for v in list(obj.values())[:50]:
            if isinstance(v, (pd.DataFrame, pd.Series, np.ndarray, dict)):
                out += _frames_from_object(v)
        return out
    if isinstance(obj, list):
        try:
            return [pd.DataFrame(obj)]
        except Exception:
            return []
    return []


def read_frames(path: str | Path, data: bytes | None = None, max_rows: int | None = None) -> list[pd.DataFrame] | None:
    """Read a file (or its staged bytes) into DataFrames. None if unreadable/unsupported."""
    name = str(path).lower()
    base, comp = _split_suffix(name)
    src = io.BytesIO(data) if data is not None else str(path)
    try:
        if base.endswith((".parquet", ".pq")):
            import pyarrow.parquet as pq

            pf = pq.ParquetFile(src)
            batches, n = [], 0
            for b in pf.iter_batches(batch_size=50_000):
                batches.append(b.to_pandas())
                n += len(batches[-1])
                if max_rows and n >= max_rows:
                    break
            return _frames_from_object(pd.concat(batches)) if batches else []
        if base.endswith((".feather", ".arrow")):
            return _frames_from_object(pd.read_feather(src))
        if base.endswith(".npy"):
            return _frames_from_object(np.load(src, allow_pickle=False))
        if base.endswith(".npz"):
            with np.load(src, allow_pickle=False) as z:
                return [f for k in z.files for f in _frames_from_object(z[k])]
        if base.endswith((".pkl", ".pickle", ".json")):
            raw = data if data is not None else Path(path).read_bytes()
            if comp:
                import bz2
                import gzip
                import lzma
                import zipfile

                if comp == "zip":
                    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
                        raw = zf.read(zf.namelist()[0])
                else:
                    raw = {"gzip": gzip.decompress, "bz2": bz2.decompress, "xz": lzma.decompress}[comp](raw)
            if base.endswith(".json"):
                try:
                    return _frames_from_object(pd.read_json(io.BytesIO(raw)))
                except Exception:
                    return _frames_from_object(json.loads(raw))
            return _frames_from_object(_SafeUnpickler(io.BytesIO(raw)).load())
        if base.endswith(TEXT_TABLES):
            if re.search(r"\.(tsv|tab)$", base):
                sep = "\t"
            elif base.endswith(".csv"):
                sep = ","
            else:
                sep = None
            return [pd.read_csv(src, sep=sep, engine="python" if sep is None else "c", nrows=max_rows,
                                compression=comp, **({} if sep is None else {"low_memory": False}))]
    except Exception:
        return None
    return None


# ----------------------------------------------------------------------------- key detection
MIN_BARE_INT_ROWS = 10  # integer columns shorter than this are never read as bare ModelIDs


def norm_mid(v, known: set[str] | None) -> str | None:
    if isinstance(v, (int, np.integer)) and not isinstance(v, bool):
        return f"ACH-{int(v):06d}" if known and f"ACH-{int(v):06d}" in known else None
    s = str(v).strip()
    if ACH.match(s):
        return s.upper()
    if known and BARE.match(s) and f"ACH-{s}" in known:
        return f"ACH-{s}"
    return None


def _mid_array(values, known) -> np.ndarray:
    return np.array([norm_mid(v, known) for v in values], dtype=object)


def _is_mid_col(series: pd.Series, known) -> bool:
    s = series.dropna()
    if s.empty:
        return False
    if pd.api.types.is_integer_dtype(s.dtype):  # bare IDs that lost their zero padding (e.g. CSV round trip)
        # guard changelog 2026-10-27: tiny integer columns (group ids, ranks) are not ModelIDs
        return bool(known) and len(s) >= MIN_BARE_INT_ROWS and s.is_unique and float(np.mean([f"ACH-{int(v):06d}" in known for v in s])) >= 0.8
    if s.dtype != object:
        return False
    sample = s.iloc[:: max(1, len(s) // 500)]
    return np.mean([norm_mid(v, known) is not None for v in sample]) >= 0.8


def _entrez_of(name) -> int | None:
    s = str(name)
    m = GENE_COL.match(s)
    if m:
        return int(m.group(1))
    if s.isdigit() and 0 < int(s) < 10**9 and not s.startswith("0"):
        return int(s)
    return None


@dataclass
class Layout:
    kind: str  # "wide", "long" or "array"
    models: np.ndarray | None = None
    gene_cols: list = field(default_factory=list)
    entrez: np.ndarray | None = None
    symbol_entrez_cols: int = 0
    matrix: pd.DataFrame | None = None
    model_col: str | None = None
    value_cols: list = field(default_factory=list)
    frame: pd.DataFrame | None = None


def detect(df: pd.DataFrame, known: set[str] | None = None) -> list[Layout]:
    out: list[Layout] = []
    if df is None or df.empty:
        return out
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    key = next((c for c in df.columns if _is_mid_col(df[c], known)), None)
    gene_cols = [c for c in df.columns if _entrez_of(c) is not None]
    if key is not None and gene_cols:
        mat = df[gene_cols]
        if not all(pd.api.types.is_numeric_dtype(t) for t in mat.dtypes):
            mat = mat.apply(pd.to_numeric, errors="coerce")
        out.append(Layout("wide", models=_mid_array(df[key], known), gene_cols=gene_cols,
                          entrez=np.array([_entrez_of(c) for c in gene_cols]),
                          symbol_entrez_cols=sum(bool(GENE_COL.match(c)) for c in gene_cols), matrix=mat))
    mid_cols = [c for c in df.columns if norm_mid(c, known) is not None]
    if len(mid_cols) >= 10:
        gkey = next((c for c in df.columns if c not in mid_cols and df[c].dtype == object
                     and df[c].dropna().astype(str).map(lambda x: _entrez_of(x) is not None).mean() >= 0.8), None)
        if gkey is not None:
            t = df.set_index(gkey)[mid_cols].T
            t.columns = [str(c) for c in t.columns]
            gcols = [c for c in t.columns if _entrez_of(c) is not None]
            out.append(Layout("wide", models=_mid_array(mid_cols, known), gene_cols=gcols,
                              entrez=np.array([_entrez_of(c) for c in gcols]),
                              symbol_entrez_cols=sum(bool(GENE_COL.match(c)) for c in gcols),
                              matrix=t[gcols].apply(pd.to_numeric, errors="coerce")))
    if key is not None and df.shape[1] <= 64:  # long tables are narrow
        for c in df.columns:
            if c == key:
                continue
            vals = df[c].dropna()
            if vals.empty or vals.nunique() <= 1:
                continue
            parsed = vals.astype(str).map(_entrez_of)
            if parsed.notna().mean() >= 0.95:
                num = [x for x in df.columns if x not in (key, c) and pd.api.types.is_numeric_dtype(df[x])]
                ent = df[c].astype(str).map(_entrez_of)
                out.append(Layout("long", models=_mid_array(df[key], known), model_col=key, entrez=ent.to_numpy(),
                                  value_cols=num, frame=df))
                break
    if not out and df.shape[0] >= 10 and df.shape[1] >= SHAPE_MIN_GENE_COLS and all(
            pd.api.types.is_numeric_dtype(t) for t in df.dtypes):
        out.append(Layout("array", matrix=df))
    return out


# ----------------------------------------------------------------------------- checks
def _row_selection(models: np.ndarray, heldout: set[str] | None, n_cols: int) -> np.ndarray:
    """All rows keyed by held-out lines (capped) plus a strided sample of the other rows."""
    ho = np.array([m in heldout for m in models]) if heldout else np.zeros(len(models), bool)
    idx_ho = np.flatnonzero(ho)[: max(1, MAX_HELDOUT_CELLS // max(1, n_cols))]
    other = np.flatnonzero(~ho)
    k = max(1, MAX_OTHER_CELLS // max(1, n_cols))
    idx_other = other[:: max(1, len(other) // k)][:k] if len(other) else other
    return np.concatenate([idx_ho, idx_other])


def _cells(lay: Layout, heldout):
    if lay.kind == "wide":
        rows = _row_selection(lay.models, heldout, lay.matrix.shape[1])
        sub = lay.matrix.iloc[rows]
        mids = np.repeat(lay.models[rows], sub.shape[1])
        ents = np.tile(lay.entrez, len(rows))
        yield mids, ents, sub.to_numpy(dtype=float).ravel()
    elif lay.kind == "long":
        rows = _row_selection(lay.models, heldout, 1)
        f = lay.frame.iloc[rows]
        mids, ents = lay.models[rows], lay.entrez[rows]
        ok = pd.notna(ents) & pd.notna(mids)
        for c in lay.value_cols:
            yield mids[ok], ents[ok].astype(int), f[c].to_numpy(dtype=float)[ok]


def _n_match(h: np.ndarray, ref: np.ndarray | None) -> int:
    if ref is None or len(ref) == 0 or len(h) == 0:
        return 0
    pos = np.searchsorted(ref, h)
    pos[pos >= len(ref)] = 0
    return int((ref[pos] == h).sum())


def fingerprint_hit(layouts: list[Layout], ctx: Context) -> bool:
    fp = ctx.fp
    if not fp:
        return False
    heldout = ctx.heldout_ids or set()
    for lay in layouts:
        if lay.kind == "array":
            continue
        for mids, ents, vals in _cells(lay, heldout):
            ho = np.array([m in heldout for m in mids], dtype=bool)
            for src in labels_mod.RAW_FILES:
                for sign in (1.0, -1.0):
                    for part, sel, min_n in (("heldout", ho, FP_MIN_HELDOUT), ("train", ~ho, FP_MIN_TRAIN)):
                        if not sel.any():
                            continue
                        h, _ = labels_mod.fingerprint_keys(mids[sel], ents[sel], sign * vals[sel], src, fp["salt"])
                        n = _n_match(h, fp.get(f"{src}_{part}"))
                        if n >= min_n and n / max(1, len(h)) >= FP_MIN_FRAC:
                            return True
    return False


def shape_reasons(layouts: list[Layout], ctx: Context) -> list[str]:
    out = []
    for lay in layouts:
        if lay.kind == "array":
            v = lay.matrix.to_numpy(dtype=float)
            v = v[np.isfinite(v)]
            if v.size and float(np.mean(v < SHAPE_NEG_THRESHOLD)) >= SHAPE_NEG_FRAC:
                out.append("shape: unkeyed Chronos-like array")
            continue
        if lay.kind != "wide" or lay.symbol_entrez_cols < SHAPE_MIN_GENE_COLS:
            continue
        sym_cols = [c for c in lay.gene_cols if GENE_COL.match(c)]
        v = lay.matrix[sym_cols].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        if float(np.mean(v < SHAPE_NEG_THRESHOLD)) >= SHAPE_NEG_FRAC:
            out.append("shape: Chronos-like ModelID x gene matrix")
        elif set(np.unique(v)) <= {0.0, 1.0}:
            models = {m for m in lay.models if m}
            if (models & ctx.heldout_ids) if ctx.heldout_ids is not None else models:
                out.append("shape: binary matrix keyed by held-out ModelIDs")
    return out


def test_keyed(layouts: list[Layout], test_ids: set[str]) -> bool:
    for lay in layouts:
        if lay.kind == "array":
            continue
        if set(m for m in lay.models if m) & test_ids and (lay.kind == "long" or len(lay.gene_cols) >= SHAPE_MIN_GENE_COLS):
            return True
    return False


# ----------------------------------------------------------------------------- per file
def canonical_label_location(path: Path) -> bool:
    """True for the label stores themselves (raw originals, data/labels, data/locked)."""
    from nextgem import paths

    real = path.resolve()
    roots = {paths.data_root(), (paths.PROJECT_ROOT / "data").resolve()}
    for r in roots:
        for d in (r / "depmap", r / "labels", r / "locked"):
            if real == d or str(real).startswith(str(d) + os.sep):
                return True
    return False


def scan_file(path: Path, ctx: Context, data: bytes | None = None, rel: str | None = None) -> Result:
    import hashlib
    import sys

    tw = sys.modules.get("nextgem_tripwire")
    res = Result(str(path))
    name = path.name.lower()
    if tw is not None and tw.classify(str(path)) is not None:
        if not canonical_label_location(path):
            res.hits.append("name: protected label file name outside its canonical location")
        else:
            res.note = "canonical label store (skipped)"
        return res
    # the canonical split file (ModelID, tier, covariates; no labels) is exempt only at its own location and
    # only if its content matches the make_splits.py pin; any copy or modification is scanned as usual
    if ctx.canonical_split_sha and ctx.canonical_split is not None and path.resolve() == ctx.canonical_split.resolve():
        content = data if data is not None else path.read_bytes()
        if hashlib.sha256(content).hexdigest() == ctx.canonical_split_sha:
            res.note = "canonical split (pinned by SPLIT_SHA256; skipped)"
            return res
    if ctx.container_rule and rel and rel.startswith("results/") and name.endswith(BINARY_CONTAINERS):
        content = data if data is not None else path.read_bytes()
        if ctx.allowlist.get(rel) != hashlib.sha256(content).hexdigest():
            res.hits.append("binary container under results/ (not allow-listed)")
    base, _ = _split_suffix(name)
    if base.endswith(UNSCANNABLE):
        res.note = "unscannable binary container"
        return res
    if not base.endswith(TABLE_BASES):
        res.note = "not a table"
        return res
    size = len(data) if data is not None else path.stat().st_size
    frames = read_frames(path, data, None if size <= MAX_FULL_READ_BYTES else MAX_ROWS_LARGE)
    if frames is None:
        res.note = "unreadable"
        return res
    lays = [lay for f in frames for lay in detect(f, ctx.known_ids)]
    if not lays:
        res.note = "no ModelID keys"
        return res
    if fingerprint_hit(lays, ctx):
        res.hits.append("fingerprint: matches label values")
    res.hits += sorted(set(shape_reasons(lays, ctx)))
    if ctx.test_ids and test_keyed(lays, ctx.test_ids):
        res.hits.append("test-keyed: table keyed by test ModelIDs x genes")
    return res


def iter_files(roots: list[Path], exclude: list[Path], registered: set[Path]):
    excl = [e.resolve() for e in exclude]
    for root in roots:
        root = Path(root)
        if root.is_file():
            yield root
            continue
        for dirpath, dirnames, filenames in os.walk(root, onerror=lambda e: None):
            dp = Path(dirpath)
            keep = []
            for d in dirnames:
                full = (dp / d).resolve()
                if d in SKIP_DIRS or any(full == e or str(full).startswith(str(e) + os.sep) for e in excl):
                    continue
                keep.append(d)
            dirnames[:] = keep
            for fn in filenames:
                p = dp / fn
                try:
                    if not p.is_file() or p.resolve() in registered:
                        continue
                except OSError:
                    continue
                yield p
