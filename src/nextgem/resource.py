"""W8 resource: SBML L3V1 FBC v2 export and validation gates (DESIGN.md section 9 W8, section 12 row W8).

Label-free. Nothing here reads labels. Callers set ``cobra.Configuration().solver`` before loading models.

* ``export_model`` writes ``{ModelID}_{parent}-{version}_{pipelineID}.xml.gz`` (deterministic gzip, mtime 0)
  with the build metadata recorded in the SBML model notes.
* ``validate_model`` runs the per-model gates: libSBML (no new errors vs the parent), cobra round-trip
  (identical reaction / metabolite / gene sets and bounds), WT growth, 20 random KOs, viability tasks.
* ``load_checked`` is the only sanctioned way to read a resource file; it rejects corrupted files.
"""

from __future__ import annotations

import gzip
import json
import math
import tempfile
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import cobra
import libsbml
import numpy as np

from nextgem.config import SEEDS
from nextgem.run import sha256_file

REL_TOL = 1e-6
N_RANDOM_KO = 20
ABS_FLOOR = 1e-9  # LP noise floor: lethal KOs return 0 or ~1e-12 depending on the solve path
NOTE_KEYS = ("release_status", "model_id", "parent", "parent_version", "pipeline_id", "medium", "medium_source_tier",
             "git_commit", "input_shas", "repaired_reactions", "n_repaired_reactions", "build_time_s")
MANIFEST_COLUMNS = ["file", "model_id", "parent", "version", "pipeline_id", "sha256", "size_bytes", "build_time_s",
                    "n_reactions", "n_metabolites", "n_genes", "n_closed_reactions", "medium", "medium_source_tier",
                    "benchmarked", "flags", "validation_status", "validation_detail", "release_status"]


class ResourceError(RuntimeError):
    """A resource file that fails to load or validate."""


@dataclass
class Meta:
    """Context metadata written to the SBML notes and to the manifest."""

    model_id: str  # ModelID (ACH-...) or "generic"
    parent: str  # "human1" | "human2"
    parent_version: str  # e.g. "2.1.1"
    pipeline_id: str  # e.g. "F0R1" or "none" for a generic parent
    medium: str = "M0"
    medium_source_tier: str = "generic (Ham's, GR task)"
    git_commit: str = ""
    input_shas: dict = field(default_factory=dict)
    repaired_reactions: list = field(default_factory=list)
    n_repaired_reactions: int | None = None  # when the list itself was not stored by the build
    build_time_s: float = float("nan")
    benchmarked: bool = False
    flags: str = ""
    release_status: str = "release"  # "demo, not a release" for demos

    def stem(self) -> str:
        return f"{self.model_id}_{self.parent}-{self.parent_version}_{self.pipeline_id}"


def set_notes(model: cobra.Model, meta: Meta) -> None:
    n_rep = meta.n_repaired_reactions if meta.n_repaired_reactions is not None else len(meta.repaired_reactions)
    model.notes.update({
        "release_status": meta.release_status, "model_id": meta.model_id, "parent": meta.parent,
        "parent_version": meta.parent_version, "pipeline_id": meta.pipeline_id, "medium": meta.medium,
        "medium_source_tier": meta.medium_source_tier, "git_commit": meta.git_commit,
        "input_shas": json.dumps(meta.input_shas, sort_keys=True),
        "repaired_reactions": ";".join(sorted(meta.repaired_reactions)) if meta.repaired_reactions else "not_recorded",
        "n_repaired_reactions": str(n_rep), "build_time_s": f"{meta.build_time_s:.1f}",
    })


def read_notes(model: cobra.Model) -> dict:
    import html

    return {k: html.unescape(str(v)) for k, v in model.notes.items()}


def write_sbml_gz(model: cobra.Model, path: Path) -> None:
    """SBML L3V1 FBC v2 (cobra default), gzipped deterministically."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "m.xml"
        cobra.io.write_sbml_model(model, str(tmp))
        raw = tmp.read_bytes()
    with open(path, "wb") as fh, gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0, compresslevel=6) as gz:
        gz.write(raw)


def export_model(model: cobra.Model, meta: Meta, out_dir: Path) -> Path:
    """Write ``model`` (with ``meta`` in its notes) under ``out_dir``; returns the file path."""
    old_id, old_notes = model.id, dict(model.notes)
    try:
        model.id = f"{meta.model_id}_{meta.parent}_{meta.pipeline_id}".replace("-", "_")
        set_notes(model, meta)
        path = Path(out_dir) / f"{meta.stem()}.xml.gz"
        write_sbml_gz(model, path)
    finally:
        model.id = old_id
        model.notes.clear()
        model.notes.update(old_notes)
    return path


# ----------------------------------------------------------------------------- reading
def _read_bytes(path: Path) -> bytes:
    path = Path(path)
    try:
        data = path.read_bytes()
        if path.suffix == ".gz" or data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        raise ResourceError(f"{path.name}: unreadable or corrupt gzip ({e})") from e
    return data


def libsbml_errors(path: Path) -> Counter:
    """Counter of libSBML error ids (severity >= error, parse plus consistency check) for a file."""
    doc = libsbml.readSBMLFromString(_read_bytes(path).decode("utf-8", errors="replace"))
    errs = Counter()
    for i in range(doc.getNumErrors()):
        e = doc.getError(i)
        if e.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
            errs[e.getErrorId()] += 1
    if not errs:  # consistency checks are only meaningful on a parsable document
        doc.checkConsistency()
        for i in range(doc.getNumErrors()):
            e = doc.getError(i)
            if e.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
                errs[e.getErrorId()] += 1
    return errs


def new_errors(child: Counter, parent: Counter) -> dict:
    return {k: v - parent.get(k, 0) for k, v in child.items() if v > parent.get(k, 0)}


def load_checked(path: Path) -> cobra.Model:
    """Read a resource file; raise ``ResourceError`` if it is corrupt (bad gzip, malformed XML, fatal libSBML)."""
    data = _read_bytes(path)
    doc = libsbml.readSBMLFromString(data.decode("utf-8", errors="replace"))
    fatal = [doc.getError(i).getMessage()[:120] for i in range(doc.getNumErrors())
             if doc.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_FATAL
             or doc.getError(i).getErrorId() < 10000 and doc.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR]
    if doc.getModel() is None or fatal:
        raise ResourceError(f"{Path(path).name}: not a valid SBML document ({fatal[:1]})")
    try:
        return cobra.io.read_sbml_model(data.decode("utf-8"))
    except Exception as e:  # noqa: BLE001 - any cobra/libsbml failure is a rejection
        raise ResourceError(f"{Path(path).name}: cobra could not read the model ({e})") from e


# ----------------------------------------------------------------------------- gates
def roundtrip_diff(a: cobra.Model, b: cobra.Model) -> dict:
    """Differences between two models: id sets and bounds. Empty dict = identical."""
    out = {}
    for what in ("reactions", "metabolites", "genes"):
        ia, ib = {x.id for x in getattr(a, what)}, {x.id for x in getattr(b, what)}
        if ia != ib:
            out[f"{what}_only_a"], out[f"{what}_only_b"] = len(ia - ib), len(ib - ia)
    if "reactions_only_a" not in out and "reactions_only_b" not in out:
        bad = [r.id for r in a.reactions
               if not (math.isclose(r.lower_bound, b.reactions.get_by_id(r.id).lower_bound, rel_tol=1e-12, abs_tol=1e-12)
                       and math.isclose(r.upper_bound, b.reactions.get_by_id(r.id).upper_bound, rel_tol=1e-12,
                                        abs_tol=1e-12))]
        if bad:
            out["bound_mismatch"] = len(bad)
            out["bound_mismatch_examples"] = bad[:5]
        gbad = [r.id for r in a.reactions if r.gene_reaction_rule != b.reactions.get_by_id(r.id).gene_reaction_rule]
        if gbad:
            out["gpr_mismatch"] = len(gbad)
    return out


def _growth(model: cobra.Model) -> float:
    v = model.slim_optimize(error_value=math.nan)
    return 0.0 if model.solver.status == "infeasible" else v


def _rel_equal(x: float, y: float, tol: float = REL_TOL) -> bool:
    if math.isnan(x) or math.isnan(y):
        return math.isnan(x) and math.isnan(y)
    return abs(x - y) <= tol * max(abs(x), abs(y), 1e-12) or abs(x - y) <= ABS_FLOOR


def random_ko_compare(a: cobra.Model, b: cobra.Model, n: int = N_RANDOM_KO, seed: int = SEEDS[0]) -> dict:
    genes = sorted(g.id for g in a.genes)
    rng = np.random.default_rng(seed)
    pick = [genes[i] for i in rng.choice(len(genes), size=min(n, len(genes)), replace=False)]
    bad = []
    for g in pick:
        va, vb = [], []
        for m, out in ((a, va), (b, vb)):
            with m:
                m.genes.get_by_id(g).knock_out()
                out.append(_growth(m))
        if not _rel_equal(va[0], vb[0]):
            bad.append((g, va[0], vb[0]))
    return {"n": len(pick), "genes": pick, "mismatches": bad}


def viability_tasks(model: cobra.Model, tasks_file: Path, types=("GR", "ER")) -> dict:
    """Pass/fail of the viability tasks (raven task semantics, boundaries closed) on a model as it is."""
    from raven_toolbox.tasks.check import task_name_maps
    from raven_toolbox.tasks.tasklist import parse_task_list

    from nextgem import context as cx
    from nextgem.knockouts import solver_name  # noqa: F401 - (kept for solver-name symmetry in logs)

    tasks = [t for t in parse_task_list(str(tasks_file)) if t.id in types]
    maps = task_name_maps(model)
    out = {}
    for t in tasks:
        saved = None
        try:
            with model:
                saved, err = cx._apply_task(model, t, maps)
                if err:
                    out[t.id] = False
                    continue
                model.objective = {}
                model.slim_optimize(error_value=math.nan)
                out[t.id] = (model.solver.status == "optimal") != t.should_fail
        finally:
            cx._restore_mets(model, saved)
    return out


@dataclass
class Validation:
    status: str  # "pass" | "fail"
    gates: dict
    detail: str = ""


def validate_model(source: cobra.Model, path: Path, parent_xml: Path | None, tasks_file: Path | None,
                   parent_errors: Counter | None = None, n_ko: int = N_RANDOM_KO) -> Validation:
    """Run all per-model gates for ``path`` against the in-memory ``source`` model it was exported from.

    The task gate runs on the re-read model ``b`` (what a user gets). ``tasks_file=None`` skips it and
    reports the gate as "skipped" (status fail for a resource build, allowed for toy tests).
    """
    gates: dict = {}
    try:
        b = load_checked(path)
    except ResourceError as e:
        return Validation("fail", {"load": f"rejected: {e}"}, str(e))
    if parent_errors is None:
        parent_errors = libsbml_errors(parent_xml) if parent_xml else Counter()
    ne = new_errors(libsbml_errors(path), parent_errors)
    gates["libsbml_no_new_errors"] = (not ne, ne)
    d = roundtrip_diff(source, b)
    gates["roundtrip_identical"] = (not d, d)
    wa, wb = _growth(source), _growth(b)
    gates["wt_growth_match"] = (_rel_equal(wa, wb), (wa, wb))
    ko = random_ko_compare(source, b, n_ko)
    gates["random_ko_match"] = (not ko["mismatches"], {"n": ko["n"], "mismatches": ko["mismatches"][:3]})
    if tasks_file is None:
        gates["viability_tasks"] = ("skipped", {})
    else:
        t = viability_tasks(b, tasks_file)
        gates["viability_tasks"] = (all(t.values()) and bool(t), t)
    failed = [k for k, (ok, _) in gates.items() if ok is not True]
    return Validation("pass" if not failed else "fail", gates, ",".join(failed))


def manifest_row(path: Path, meta: Meta, model: cobra.Model, val: Validation, n_closed: int = 0) -> dict:
    return {
        "file": str(path), "model_id": meta.model_id, "parent": meta.parent, "version": meta.parent_version,
        "pipeline_id": meta.pipeline_id, "sha256": sha256_file(path), "size_bytes": Path(path).stat().st_size,
        "build_time_s": f"{meta.build_time_s:.1f}", "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites), "n_genes": len(model.genes), "n_closed_reactions": n_closed,
        "medium": meta.medium, "medium_source_tier": meta.medium_source_tier,
        "benchmarked": str(bool(meta.benchmarked)).lower(), "flags": meta.flags,
        "validation_status": val.status, "validation_detail": val.detail, "release_status": meta.release_status,
    }


def write_manifest(rows: list[dict], path: Path) -> None:
    import csv

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def line_model_from_mask(parent: cobra.Model, hexmask: str) -> tuple[cobra.Model, int]:
    """Reconstruct a Stage-A per-line model from the stored ``open_mask`` (hex of packbits over the parent's
    reaction order; 1 = open): closed reactions get bounds (0, 0). Valid for families whose line model is a
    closed set only (F0, F2, F3, F3b, F5), not E-Flux (F4), whose scaled bounds are not in the mask."""
    bits = np.unpackbits(np.frombuffer(bytes.fromhex(hexmask), dtype=np.uint8))[: len(parent.reactions)].astype(bool)
    if len(bits) != len(parent.reactions):
        raise ResourceError("open_mask length does not match the parent reaction count")
    m = parent.copy()
    n = 0
    for r, open_ in zip(m.reactions, bits):
        if not open_:
            r.bounds = (0.0, 0.0)
            n += 1
    return m, n


def now() -> float:
    return time.perf_counter()
