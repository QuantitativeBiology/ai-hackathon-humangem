"""Label tripwire (DESIGN.md §11.1, PI decision r2-15).

Loaded for every interpreter launch of ``.venv`` by ``sitecustomize.py`` (installed by
``scripts/install_tripwire.py``), as the standalone module ``nextgem_tripwire``: it does not import
the ``nextgem`` package and uses only the standard library, so it costs ~1 ms per launch.

What it guards
--------------
Protected paths ("categories"):

- ``raw``          the raw label files, by basename: CRISPRGeneEffect, CRISPRGeneDependency,
                   ScreenGeneEffect, CRISPRInferredCommonEssentials, AchillesCommonEssentialControls,
                   AchillesNonessentialControls (any directory; hard links to the real files by inode);
- ``train``/``val``/``test``  ``labels_train*``, ``labels_val*``, ``labels_test*`` (any directory);
- ``lists``        ``label_derived_lists*`` (label-derived gene lists, §11.2);
- ``fingerprints``/``sealed``/``locked``  anything under ``<data root>/locked`` (and always under the
                   real ``<project>/data/locked``, whatever ``NEXTGEM_DATA_ROOT`` says).

Hooks: ``sys.addaudithook`` for ``open`` (builtins/io/os, hence pandas' readers), ``os.rename``,
``os.link``, ``os.symlink``, ``os.remove``, ``shutil.copyfile`` and process launches
(``subprocess.Popen``, ``os.system``, ``os.exec``, ``os.posix_spawn``, ``os.spawn``) whose arguments
name a protected path; wrappers for the pyarrow readers and writers, which raise no audit event
(``pyarrow.parquet.read_table``/``ParquetFile``/``ParquetDataset``/``write_table``,
``pyarrow.dataset.dataset``, ``pyarrow.csv.read_csv``/``open_csv``, ``pyarrow.feather``,
``pyarrow.json.read_json``, ``pyarrow.memory_map``/``input_stream``), installed at import time.

Rule
----
An access is allowed only if **both** hold (DESIGN §11.1):

1. a frame of ``src/nextgem/labels.py`` is on the stack and has opened a *scope* naming the category
   (``labels.py`` opens scopes; nothing else should); and
2. the scope is a public, fold-scoped training API scope (read of ``train`` or ``lists``), **or** the
   entry point (``__main__.__file__``, which must be ``<checkout>/scripts/<name>`` of the same
   checkout as that ``labels.py``) is permitted for that category and mode (``ENTRY_PERMS``).

Writes of raw-named files are denied only inside the data roots (originals are read-only); elsewhere a
raw-named write is a toy fixture and its later *read* is still guarded.

Every allowed or denied event on a protected path is logged to ``<state root>/logs/label_access.log``.

Residual risk (documented in DESIGN §11.1): interpreters outside ``.venv``, ``python -S``/``-I``, a
``sitecustomize.py`` earlier on ``sys.path``, native libraries that open files themselves (other than
the wrapped pyarrow entry points), and deliberate circumvention. The PreToolUse hook and the content
scanner (``scripts/scan_label_copies.py``) are the other two layers.
"""

from __future__ import annotations

import contextlib
import importlib.abc
import importlib.machinery
import os
import re
import sys
import threading
import time

VERSION = "1"

_HERE = os.path.dirname(os.path.realpath(__file__))  # <checkout>/src/nextgem
PROJECT_ROOT = os.path.dirname(os.path.dirname(_HERE))
REAL_DATA_ROOT = os.path.join(PROJECT_ROOT, "data")

RAW_RE = re.compile(
    r"CRISPRGeneEffect|CRISPRGeneDependency|ScreenGeneEffect|CRISPRInferredCommonEssentials|"
    r"AchillesCommonEssentialControls|AchillesNonessentialControls"
)
TIER_RE = re.compile(r"labels_(train|val|test)")
LISTS_RE = re.compile(r"label_derived_lists")
LOCKED_SEG_RE = re.compile(r"(^|/)data/locked(/|$)")
# For process arguments: any mention of a protected name or a locked directory.
ARG_RE = re.compile(RAW_RE.pattern + r"|labels_(train|val|test)|label_derived_lists|data/locked")

_TIER_CAT = {"train": "train", "val": "val", "test": "test"}

# Entry points (basename of __main__.__file__ under <checkout>/scripts/) -> permissions.
# "hash" = read-only hashing of any protected file via labels.sha256_protected().
ENTRY_PERMS: dict[str, dict[str, set[str]]] = {
    "build_labels.py": {
        "read": {"raw", "train", "val", "test", "lists", "fingerprints", "locked"},
        "write": {"train", "val", "test", "lists", "fingerprints", "locked"},
    },
    "evaluate_val.py": {"read": {"val"}, "write": set()},
    "evaluate_test.py": {"read": {"val", "test", "sealed"}, "write": {"sealed"}},
    "check_frozen.py": {"read": {"hash"}, "write": set()},
    "verify_ledgers.py": {"read": {"hash"}, "write": set()},
    "scan_label_copies.py": {"read": {"fingerprints"}, "write": set()},
}
PUBLIC_READ = {"train", "lists"}  # fold-scoped training API (labels.train(), labels.fold(k), lists)

_SKIP_SUFFIXES = (".py", ".pyc", ".so", ".pth", ".pyi", ".typed")

_state = threading.local()
_installed = False
_raw_inodes: set[tuple[int, int]] | None = None
_realpath_cache: dict[str, str] = {}
_labels_files: dict[str, bool] = {}


class TripwireDenied(PermissionError):
    """Raised for a denied access to a protected label path."""


# ----------------------------------------------------------------------------- classification
def _data_roots() -> list[str]:
    roots = [os.path.realpath(REAL_DATA_ROOT)]
    env = os.environ.get("NEXTGEM_DATA_ROOT")
    if env:
        roots.append(os.path.realpath(env))
    return roots


def _locked_dirs() -> list[str]:
    return [os.path.join(r, "locked") for r in _data_roots()]


def _state_root() -> str:
    return os.path.realpath(os.environ.get("NEXTGEM_STATE_ROOT") or PROJECT_ROOT)


def _load_raw_inodes() -> dict[tuple[int, int], tuple[str, bool]]:
    """(dev, ino) -> (category, is_real) for the raw label files of every data root and for the files of
    the *real* data/locked and data/labels stores (hard links to them are caught by inode)."""
    global _raw_inodes
    if _raw_inodes is None:
        inodes: dict[tuple[int, int], tuple[str, bool]] = {}
        real_root = os.path.realpath(REAL_DATA_ROOT)
        for root in _data_roots():
            is_real = root == real_root
            d = os.path.join(root, "depmap")
            try:
                names = os.listdir(d)
            except OSError:
                names = []
            for n in names:
                if RAW_RE.search(n):
                    try:
                        st = os.stat(os.path.join(d, n))
                        inodes[(st.st_dev, st.st_ino)] = ("raw", is_real)
                    except OSError:
                        pass
        for sub in ("locked", "labels"):
            base = os.path.join(real_root, sub)
            for dirpath, _dirs, files in os.walk(base):
                for n in files:
                    full = os.path.join(dirpath, n)
                    cat = _locked_category(full.replace(os.sep, "/")) if sub == "locked" else (
                        _TIER_CAT[m.group(1)] if (m := TIER_RE.search(n)) else ("lists" if LISTS_RE.search(n) else None))
                    if cat is None:
                        continue
                    try:
                        st = os.stat(full)
                        inodes[(st.st_dev, st.st_ino)] = (cat, True)
                    except OSError:
                        pass
        _raw_inodes = inodes
    return _raw_inodes


def is_real_data(path) -> bool:
    """True if ``path`` is (or links to) a file of the real project data root."""
    p = _as_str(path)
    if not p:
        return False
    if _under(os.path.realpath(p), os.path.realpath(REAL_DATA_ROOT)):
        return True
    try:
        st = os.stat(p)
    except OSError:
        return False
    hit = _load_raw_inodes().get((st.st_dev, st.st_ino))
    return bool(hit and hit[1])


def _as_str(path) -> str | None:
    if isinstance(path, int):
        return None
    try:
        return os.fsdecode(os.fspath(path))
    except TypeError:
        return None


def _under(path: str, root: str) -> bool:
    return path == root or path.startswith(root + os.sep)


def classify(path) -> str | None:
    """Category of a path, or None if it is not protected."""
    p = _as_str(path)
    if not p:
        return None
    if p.endswith(_SKIP_SUFFIXES):  # imports: name-based check only, no realpath/stat
        if not (RAW_RE.search(p) or TIER_RE.search(p) or LISTS_RE.search(p) or LOCKED_SEG_RE.search(p)):
            return None
    real = os.path.realpath(p)
    for cand in (p, real):
        norm = cand.replace(os.sep, "/")
        if LOCKED_SEG_RE.search(norm):
            return _locked_category(norm)
    for ld in _locked_dirs():
        if _under(real, ld):
            return _locked_category(real.replace(os.sep, "/"))
    for cand in (os.path.basename(p), os.path.basename(real)):
        if RAW_RE.search(cand):
            return "raw"
        m = TIER_RE.search(cand)
        if m:
            return _TIER_CAT[m.group(1)]
        if LISTS_RE.search(cand):
            return "lists"
    if not p.endswith(_SKIP_SUFFIXES):
        inodes = _load_raw_inodes()
        if inodes:
            try:
                st = os.stat(real)
            except OSError:
                return None
            hit = inodes.get((st.st_dev, st.st_ino))
            if hit:
                return hit[0]
    return None


def _locked_category(norm: str) -> str:
    base = norm.rsplit("/", 1)[-1]
    if "/locked/sealed" in norm:
        return "sealed"
    m = TIER_RE.search(base)
    if m:
        return _TIER_CAT[m.group(1)]
    if base.startswith("fingerprints"):
        return "fingerprints"
    return "locked"


def _in_data_root(path: str) -> bool:
    real = os.path.realpath(path)
    return any(_under(real, r) for r in _data_roots())


# ----------------------------------------------------------------------------- caller context
def _is_labels_file(filename: str) -> bool:
    hit = _labels_files.get(filename)
    if hit is None:
        real = os.path.realpath(filename)
        parts = real.split(os.sep)
        hit = len(parts) >= 3 and parts[-1] == "labels.py" and parts[-2] == "nextgem" and parts[-3] == "src"
        _labels_files[filename] = hit
    return hit


def _labels_frame_file() -> str | None:
    f = sys._getframe(1)
    while f is not None:
        fn = f.f_code.co_filename
        if fn.endswith("labels.py") and _is_labels_file(fn):
            return os.path.realpath(fn)
        f = f.f_back
    return None


def entry_point(labels_file: str | None = None) -> str | None:
    """Basename of the entry script if it lives in ``<checkout>/scripts`` of the labels.py checkout."""
    main = sys.modules.get("__main__")
    mf = getattr(main, "__file__", None)
    if not mf:
        return None
    real = os.path.realpath(mf)
    checkout = os.path.dirname(os.path.dirname(os.path.dirname(labels_file))) if labels_file else PROJECT_ROOT
    if os.path.dirname(real) != os.path.join(checkout, "scripts"):
        return None
    return os.path.basename(real)


def _scope() -> set[str] | None:
    stack = getattr(_state, "scopes", None)
    return stack[-1] if stack else None


@contextlib.contextmanager
def scope(*categories: str):
    """Opened by labels.py around each protected access. Names the categories it intends to touch."""
    stack = getattr(_state, "scopes", None)
    if stack is None:
        stack = _state.scopes = []
    stack.append(set(categories))
    try:
        yield
    finally:
        stack.pop()


def policy(category: str, write: bool, entry: str | None, sc: set[str]) -> tuple[bool, str, bool]:
    """The entry-point policy table alone. Returns (allowed, reason, needs_trust)."""
    perms = ENTRY_PERMS.get(entry or "", {"read": set(), "write": set()})
    if "hash" in sc:
        if write:
            return False, "hash scope is read-only", False
        if "hash" in perms["read"]:
            return True, "hash", True
        return False, f"entry point {entry!r} may not hash protected files", False
    if category not in sc:
        return False, f"scope {sorted(sc)} does not cover {category}", False
    if not write and category in PUBLIC_READ:
        return True, "fold-scoped training API", False
    if category in perms["write" if write else "read"]:
        return True, f"entry point {entry}", True
    return False, f"entry point {entry!r} may not {'write' if write else 'read'} {category}", False


# ----------------------------------------------------------------------------- code trust (review M4)
_sha_cache: dict[str, tuple[float, int, str]] = {}
_frozen_cache: dict[str, object] = {}


def _sha256(path: str) -> str | None:
    import hashlib

    try:
        st = os.stat(path)
    except OSError:
        return None
    hit = _sha_cache.get(path)
    if hit and hit[0] == st.st_mtime and hit[1] == st.st_size:
        return hit[2]
    h = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return None
    _sha_cache[path] = (st.st_mtime, st.st_size, h.hexdigest())
    return h.hexdigest()


def _parse_frozen(text: str) -> dict[str, str]:
    rows = {}
    for ln in text.splitlines()[1:]:
        f = ln.split("\t")
        if len(f) >= 4 and f[1] == "file":
            rows[f[2]] = f[3]
    return rows


def frozen_code_rows(committed: bool) -> dict[str, str]:
    """path -> SHA-256 of ``file`` rows of the *project's* FROZEN.tsv (``HEAD:`` if committed, else the
    working tree). The project is the checkout this tripwire was loaded from (fixed by sitecustomize)."""
    key = "head" if committed else "work"
    if committed:
        if key not in _frozen_cache:
            import subprocess

            r = subprocess.run(["git", "-C", PROJECT_ROOT, "show", "HEAD:FROZEN.tsv"], capture_output=True, text=True)
            _frozen_cache[key] = _parse_frozen(r.stdout) if r.returncode == 0 else {}
        return _frozen_cache[key]  # type: ignore[return-value]
    f = os.path.join(PROJECT_ROOT, "FROZEN.tsv")
    try:
        with open(f) as fh:
            return _parse_frozen(fh.read())
    except OSError:
        return {}


def trusted(labels_file: str, entry: str | None, real_data: bool) -> tuple[bool, str]:
    """labels.py and the entry script must both carry the SHA-256 recorded for them in FROZEN.tsv:
    the committed HEAD:FROZEN.tsv for real data; for toy data (no real file involved) the working-tree
    FROZEN.tsv is accepted too."""
    committed = frozen_code_rows(True)
    work = {} if real_data else frozen_code_rows(False)

    def ok(rel: str, path: str) -> bool:
        sha = _sha256(path)
        return sha is not None and sha in {committed.get(rel), work.get(rel)} - {None}

    if not ok("src/nextgem/labels.py", labels_file):
        return False, "labels.py differs from its FROZEN.tsv SHA-256 (or is not frozen)"
    checkout = os.path.dirname(os.path.dirname(os.path.dirname(labels_file)))
    if not entry or not ok(f"scripts/{entry}", os.path.join(checkout, "scripts", entry)):
        return False, f"entry point {entry!r} differs from its FROZEN.tsv SHA-256 (or is not frozen)"
    return True, "trusted"


def decide(category: str, write: bool, path=None) -> tuple[bool, str, str | None, set[str] | None]:
    labels_file = _labels_frame_file()
    sc = _scope()
    entry = entry_point(labels_file) if labels_file else entry_point()
    if labels_file is None:
        return False, "no labels.py frame on the stack", entry, sc
    if sc is None:
        return False, "labels.py opened no scope", entry, sc
    ok, reason, needs_trust = policy(category, write, entry, sc)
    if ok and needs_trust:
        ok, why = trusted(labels_file, entry, is_real_data(path) if path is not None else True)
        if not ok:
            reason = why
    return ok, reason, entry, sc


# ----------------------------------------------------------------------------- logging
def _log(decision: str, event: str, category: str, write: bool, path: str, entry, sc, reason: str) -> None:
    try:
        d = os.path.join(_state_root(), "logs")
        os.makedirs(d, exist_ok=True)
        line = "\t".join([
            time.strftime("%Y-%m-%dT%H:%M:%S"), str(os.getpid()), decision, event, category,
            "w" if write else "r", path, str(entry), ",".join(sorted(sc)) if sc else "-", reason,
        ]) + "\n"
        fd = os.open(os.path.join(d, "label_access.log"), os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
        try:
            os.write(fd, line.encode())
        finally:
            os.close(fd)
    except Exception:
        pass


# ----------------------------------------------------------------------------- checks
def check(path, write: bool = False, event: str = "manual") -> None:
    """Raise TripwireDenied if this access to ``path`` is not allowed. No-op for unprotected paths."""
    if getattr(_state, "busy", False):
        return
    _state.busy = True
    try:
        category = classify(path)
        if category is None:
            return
        p = _as_str(path) or ""
        if category == "raw" and write and not _in_data_root(p):
            _log("allow", event, category, write, p, None, None, "raw-named write outside data roots (toy fixture)")
            return
        if event == "os.remove" and not _in_data_root(p):
            _log("allow", event, category, write, p, None, None, "deletion of a copy outside the data roots")
            return
        ok, reason, entry, sc = decide(category, write, p)
        _log("allow" if ok else "DENY", event, category, write, p, entry, sc, reason)
        if not ok:
            raise TripwireDenied(
                f"nextgem tripwire: {'write' if write else 'read'} of protected label path ({category}) denied: "
                f"{reason}. Labels are read only through src/nextgem/labels.py (DESIGN.md §11.1). Path: {p}"
            )
    finally:
        _state.busy = False


_GLOB_TARGETS = ["depmap", "locked", "labels", "sealed", "CRISPRGeneEffect.csv", "CRISPRGeneDependency.csv",
                 "ScreenGeneEffect.csv", "CRISPRInferredCommonEssentials.csv", "AchillesCommonEssentialControls.csv",
                 "AchillesNonessentialControls.csv", "labels_train.parquet", "labels_val.parquet",
                 "labels_test.parquet", "label_derived_lists.parquet", "fingerprints.npz"]


def _glob_reaches_labels(tok: str) -> bool:
    """A path-like token whose glob segment could expand to a label directory or file name."""
    if "/" not in tok or not re.search(r"[*?\[]", tok):
        return False
    import fnmatch

    segs = tok.split("/")
    for i, seg in enumerate(segs):
        if not re.search(r"[*?\[]", seg) or not any(fnmatch.fnmatchcase(t, seg) for t in _GLOB_TARGETS):
            continue
        literal = re.sub(r"\[[^\]]*\]|[*?]", "", seg)
        if len(literal) >= 2 or (i > 0 and segs[i - 1] in ("data", "depmap", "locked", "labels")):
            return True
    return False


def check_args(values, event: str) -> None:
    """Deny a process launch whose arguments name a protected path."""
    if getattr(_state, "busy", False):
        return
    for v in values:
        s = _as_str(v) if not isinstance(v, (list, tuple)) else None
        if isinstance(v, (list, tuple)):
            check_args(v, event)
            continue
        if not s:
            continue
        hit = ARG_RE.search(s) or any(ld in s for ld in _locked_dirs()) or any(  # path globs that could reach labels
            _glob_reaches_labels(tok) for tok in re.split(r"[\s;|&'\"]+", s))
        if hit:
            _state.busy = True
            try:
                _log("DENY", event, "process", False, s[:300], entry_point(), None, "process argument names a protected path")
            finally:
                _state.busy = False
            raise TripwireDenied(
                f"nextgem tripwire: process launch naming a protected label path denied ({event}). "
                "Labels are read only through src/nextgem/labels.py (DESIGN.md §11.1)."
            )


_WRITE_FLAGS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND


def _is_write(mode, flags) -> bool:
    if isinstance(mode, str):
        return any(c in mode for c in "wax+")
    if isinstance(flags, int):
        return bool(flags & _WRITE_FLAGS)
    return False


def _audit(event: str, args) -> None:
    if event == "open":
        path = args[0]
        if path is None or isinstance(path, int):
            return
        check(path, _is_write(args[1] if len(args) > 1 else None, args[2] if len(args) > 2 else None), "open")
    elif event in ("os.rename", "os.link", "os.symlink"):
        check(args[0], True, event)
        check(args[1], True, event)
    elif event == "os.remove":
        check(args[0], True, event)
    elif event == "shutil.copyfile":
        check(args[0], False, event)
        check(args[1], True, event)
    elif event == "subprocess.Popen":
        check_args([args[0], args[1], args[2]], event)
    elif event == "os.system":
        check_args([args[0]], event)
    elif event in ("os.exec", "os.posix_spawn", "os.spawn"):
        check_args(list(args[:3]), event)


_WATCHED = {"open", "os.rename", "os.link", "os.symlink", "os.remove", "shutil.copyfile",
            "subprocess.Popen", "os.system", "os.exec", "os.posix_spawn", "os.spawn"}


def _hook(event: str, args) -> None:
    if event not in _WATCHED:
        return
    try:
        _audit(event, args)
    except TripwireDenied:
        raise
    except Exception:  # never break unrelated I/O because of a guard bug
        return


# ----------------------------------------------------------------------------- pyarrow wrappers
_READ_FUNCS = {
    "pyarrow.parquet.core": ["read_table", "read_pandas", "read_metadata", "read_schema"],
    "pyarrow.parquet": ["read_table", "read_pandas", "read_metadata", "read_schema"],
    "pyarrow.dataset": ["dataset", "parquet_dataset"],
    "pyarrow.csv": ["read_csv", "open_csv"],
    "pyarrow.feather": ["read_table", "read_feather"],
    "pyarrow.json": ["read_json"],
    "pyarrow.orc": ["read_table"],
    "pyarrow": ["memory_map", "input_stream"],
}
_WRITE_FUNCS = {
    "pyarrow.parquet.core": ["write_table", "write_to_dataset"],
    "pyarrow.parquet": ["write_table", "write_to_dataset"],
    "pyarrow.dataset": ["write_dataset"],
    "pyarrow.csv": ["write_csv"],
    "pyarrow.feather": ["write_feather"],
    "pyarrow.orc": ["write_table"],
    "pyarrow": ["output_stream"],
}
_CLASS_INITS = {
    "pyarrow.parquet.core": [("ParquetFile", False), ("ParquetDataset", False), ("ParquetWriter", True)],
    "pyarrow.orc": [("ORCFile", False)],
}
_TARGETS = set(_READ_FUNCS) | set(_WRITE_FUNCS) | set(_CLASS_INITS)


_PATH_KWARGS = ("source", "where", "path", "path_or_paths", "base_dir", "root_path", "output_file",
                "dest", "sink", "input_file", "file")


def _check_call_args(args, kwargs, write: bool, via: str) -> None:
    """Check the path-like arguments: the first two positionals and the known path keywords."""
    cands = list(args[:2]) + [kwargs[k] for k in _PATH_KWARGS if k in kwargs]
    for v in cands:
        if isinstance(v, (str, bytes, os.PathLike)):
            check(v, write, via)
        elif isinstance(v, (list, tuple)) and v and all(isinstance(x, (str, bytes, os.PathLike)) for x in v):
            for x in v:
                check(x, write, via)


def _wrap_func(fn, write: bool, via: str):
    if getattr(fn, "_nextgem_wrapped", False):
        return fn

    def wrapper(*args, **kwargs):
        _check_call_args(args, kwargs, write, via)
        return fn(*args, **kwargs)

    wrapper._nextgem_wrapped = True
    wrapper.__wrapped__ = fn
    wrapper.__name__ = getattr(fn, "__name__", "wrapped")
    wrapper.__doc__ = getattr(fn, "__doc__", None)
    return wrapper


def _patch_module(mod) -> None:
    name = mod.__name__
    for fname in _READ_FUNCS.get(name, []):
        if hasattr(mod, fname):
            setattr(mod, fname, _wrap_func(getattr(mod, fname), False, f"{name}.{fname}"))
    for fname in _WRITE_FUNCS.get(name, []):
        if hasattr(mod, fname):
            setattr(mod, fname, _wrap_func(getattr(mod, fname), True, f"{name}.{fname}"))
    for cname, write in _CLASS_INITS.get(name, []):
        cls = getattr(mod, cname, None)
        if cls is not None and not getattr(cls.__init__, "_nextgem_wrapped", False):
            try:
                cls.__init__ = _wrap_func(cls.__init__, write, f"{name}.{cname}")
            except (TypeError, AttributeError):
                pass


class _PatchingLoader(importlib.abc.Loader):
    def __init__(self, inner):
        self._inner = inner

    def create_module(self, spec):
        return self._inner.create_module(spec)

    def exec_module(self, module):
        self._inner.exec_module(module)
        try:
            _patch_module(module)
        except Exception:
            pass

    def __getattr__(self, item):
        return getattr(self._inner, item)


class _PyarrowFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path, target=None):
        if fullname not in _TARGETS:
            return None
        spec = importlib.machinery.PathFinder.find_spec(fullname, path)
        if spec is None or spec.loader is None:
            return None
        spec.loader = _PatchingLoader(spec.loader)
        return spec


# ----------------------------------------------------------------------------- install
def install() -> None:
    global _installed
    if _installed:
        return
    sys.addaudithook(_hook)
    sys.meta_path.insert(0, _PyarrowFinder())
    for name in _TARGETS:  # already imported (unusual this early): patch now
        if name in sys.modules:
            _patch_module(sys.modules[name])
    _installed = True


def installed() -> bool:
    return _installed
