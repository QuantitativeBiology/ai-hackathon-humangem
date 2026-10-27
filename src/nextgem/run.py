"""Provenance wrapper for result-producing scripts (DESIGN.md §11.4, §6.7). Minimal implementation.

Usage::

    from nextgem.run import provenance
    with provenance("results/w1/MANIFEST.tsv", inputs=[...], seeds=[11], solver="gurobi", workers=8,
                    allow_dirty=False) as prov:
        ...                       # do the work
        prov.input(path)          # (optional) register inputs known only inside the body
        prov.output(path)         # register every output file
        prov.note(key=value)      # optional extra fields (JSON in the 'extra' column)

On exit (only if the body succeeded) one row is appended to the manifest with: timestamp, script, argv,
git commit, dirty flag (``git status --porcelain -- src scripts tests``), full SHA-256 of every input and
output, seeds, solver and its version, python / cobra / optlang / gurobipy versions, workers, wall time, CPU
time (this process plus waited-for children), peak RSS and author.

``allow_dirty=False`` (DESIGN §11.4 default) refuses to start when ``src/``, ``scripts/`` or ``tests/`` are
dirty. Inputs that are protected label files are not hashed here (that needs ``labels.sha256_protected`` from
a permitted entry point); pass their FROZEN.tsv hash through ``prov.note`` instead.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

from nextgem import paths

MANIFEST_COLUMNS = ["timestamp", "script", "argv", "commit", "dirty", "inputs", "outputs", "seeds", "solver",
                    "solver_version", "versions", "workers", "wall_s", "cpu_s", "peak_rss_mb", "author", "extra"]


class DirtyTreeError(RuntimeError):
    pass


def sha256_file(path: Path | str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit(root: Path | None = None) -> str:
    root = Path(root or paths.PROJECT_ROOT)
    r = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else ""


def git_dirty(root: Path | None = None, scopes=("src", "scripts", "tests")) -> bool:
    root = Path(root or paths.PROJECT_ROOT)
    r = subprocess.run(["git", "-C", str(root), "status", "--porcelain", "--", *scopes], capture_output=True, text=True)
    return bool(r.stdout.strip())


def versions() -> dict:
    out = {"python": platform.python_version()}
    for mod in ("cobra", "optlang", "numpy", "pandas", "sklearn"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:  # noqa: BLE001 - optional
            out[mod] = None
    try:
        import gurobipy

        out["gurobipy"] = ".".join(map(str, gurobipy.gurobi.version()))
    except Exception:  # noqa: BLE001 - optional
        out["gurobipy"] = None
    return out


def _rel(p: Path) -> str:
    p = Path(p).resolve()
    for root in (paths.state_root(), paths.PROJECT_ROOT):
        try:
            return str(p.relative_to(root))
        except ValueError:
            continue
    return str(p)


@dataclass
class Provenance:
    manifest: Path
    inputs: list
    seeds: list
    solver: str
    workers: int
    author: str
    outputs: list = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    def output(self, *paths_):
        self.outputs.extend(Path(p) for p in paths_)

    def input(self, *paths_):
        self.inputs.extend(Path(p) for p in paths_)

    def note(self, **kw):
        self.extra.update(kw)


def _solver_version(solver: str, vers: dict) -> str:
    if solver.lower().startswith("gurobi"):
        return vers.get("gurobipy") or ""
    if solver.lower().startswith("glpk"):
        try:
            import swiglpk

            return str(swiglpk.glp_version())
        except Exception:  # noqa: BLE001
            return ""
    return ""


@contextlib.contextmanager
def provenance(manifest, *, inputs=(), seeds=(), solver: str = "", workers: int = 1, allow_dirty: bool = False,
               author: str = "", argv=None):
    """Context manager writing one manifest row after a successful body (see module docstring)."""
    if not allow_dirty and git_dirty():
        raise DirtyTreeError("src/, scripts/ or tests/ have uncommitted changes; refusing to run (DESIGN §11.4)")
    t0, c0 = time.time(), time.process_time()
    ch0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    prov = Provenance(Path(manifest), [Path(p) for p in inputs], list(seeds), solver, workers, author)
    yield prov
    ch1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (time.process_time() - c0) + (ch1.ru_utime - ch0.ru_utime) + (ch1.ru_stime - ch0.ru_stime)
    rss = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, ch1.ru_maxrss) / 1024
    vers = versions()
    argv = list(argv if argv is not None else sys.argv)
    row = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "script": Path(argv[0]).name if argv else "",
        "argv": " ".join(argv[1:]), "commit": git_commit(), "dirty": git_dirty(),
        "inputs": json.dumps({_rel(p): sha256_file(p) for p in prov.inputs}, sort_keys=True),
        "outputs": json.dumps({_rel(p): sha256_file(p) for p in prov.outputs}, sort_keys=True),
        "seeds": json.dumps(prov.seeds), "solver": solver, "solver_version": _solver_version(solver, vers),
        "versions": json.dumps(vers, sort_keys=True), "workers": workers, "wall_s": round(time.time() - t0, 1),
        "cpu_s": round(cpu, 1), "peak_rss_mb": round(rss, 1), "author": author,
        "extra": json.dumps(prov.extra, sort_keys=True, default=str),
    }
    m = prov.manifest
    m.parent.mkdir(parents=True, exist_ok=True)
    new = not m.exists()
    with open(m, "a") as fh:
        if new:
            fh.write("\t".join(MANIFEST_COLUMNS) + "\n")
        fh.write("\t".join(str(row[c]).replace("\t", " ").replace("\n", " ") for c in MANIFEST_COLUMNS) + "\n")
        fh.flush()
        os.fsync(fh.fileno())
    prov.row = row
