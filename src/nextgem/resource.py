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
             "git_commit", "input_shas", "repaired_reactions", "n_repaired_reactions", "build_time_s", "flags")
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
    extra: dict = field(default_factory=dict)  # further notes (str values kept, others JSON-encoded)
    repaired_note: str = "not_recorded"  # what an empty repaired_reactions list means for this pipeline

    def stem(self) -> str:
        return f"{self.model_id}_{self.parent}-{self.parent_version}_{self.pipeline_id}"


def set_notes(model: cobra.Model, meta: Meta) -> None:
    n_rep = meta.n_repaired_reactions if meta.n_repaired_reactions is not None else len(meta.repaired_reactions)
    reserved = set(NOTE_KEYS) | {"version"}
    clash = reserved & set(meta.extra)
    if clash:
        raise ValueError(f"extra notes may not override {sorted(clash)}")
    model.notes.update({
        "release_status": meta.release_status, "model_id": meta.model_id, "parent": meta.parent,
        "parent_version": meta.parent_version, "pipeline_id": meta.pipeline_id, "medium": meta.medium,
        "medium_source_tier": meta.medium_source_tier, "git_commit": meta.git_commit,
        "input_shas": json.dumps(meta.input_shas, sort_keys=True),
        "repaired_reactions": ";".join(sorted(meta.repaired_reactions)) if meta.repaired_reactions else meta.repaired_note,
        "n_repaired_reactions": str(n_rep), "build_time_s": f"{meta.build_time_s:.1f}", "flags": meta.flags,
    })
    model.notes.update({k: (v if isinstance(v, str) else json.dumps(v, sort_keys=True)) for k, v in meta.extra.items()})


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
    """Biomass optimum; 0 if infeasible; NaN if the LP stays non-optimal after the A10 retry (DESIGN §17 A10:
    300 s per LP and one retry, via ``nextgem.w3.solve_status``)."""
    from nextgem.w3 import SolverInfraError, solve_status

    try:
        st = solve_status(model)
    except SolverInfraError:
        return math.nan
    return 0.0 if st == "infeasible" else model.solver.objective.value


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
    keys = task_keys(tasks)
    out = {}
    for key, t in zip(keys, tasks):
        saved = None
        try:
            with model:
                saved, err = cx._apply_task(model, t, maps)
                if err:
                    out[key] = False
                    continue
                model.objective = {}
                model.slim_optimize(error_value=math.nan)
                out[key] = (model.solver.status == "optimal") != t.should_fail
        finally:
            cx._restore_mets(model, saved)
    return out


def task_keys(tasks) -> list[str]:
    """Unique result keys for a task list. Raven task ids are categories (9 tasks share id "ER"), so a dict keyed by
    ``t.id`` keeps only the last task of each category; repeated ids get a 1-based suffix (ER1..ER9), unique ids
    (GR) stay as they are."""
    n = Counter(t.id for t in tasks)
    seen: Counter = Counter()
    out = []
    for t in tasks:
        seen[t.id] += 1
        out.append(t.id if n[t.id] == 1 else f"{t.id}{seen[t.id]}")
    return out


@dataclass
class Validation:
    status: str  # "pass" | "fail"
    gates: dict
    detail: str = ""


def validate_model(source: cobra.Model, path: Path, parent_xml: Path | None, tasks_file: Path | None,
                   parent_errors: Counter | None = None, n_ko: int = N_RANDOM_KO,
                   extra_gates: dict | None = None) -> Validation:
    """Run all per-model gates for ``path`` against the in-memory ``source`` model it was exported from.

    The task gate runs on the re-read model ``b`` (what a user gets). ``tasks_file=None`` skips it and
    reports the gate as "skipped" (status fail for a resource build, allowed for toy tests).
    ``extra_gates``: name -> callable(source, b) returning (ok: bool, detail), run after the five gates.
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
    for name, fn in (extra_gates or {}).items():
        ok, detail = fn(source, b)
        gates[name] = (bool(ok), detail)
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


# ----------------------------------------------------------------------------- P1 resource (RESEARCH_TREE D11)
# P1 = generic parent + W3 B-loo:cofactor_pool_biomass on the W2 documented medium M2(f=1), gate-free, no context
# family. P1 models differ across lines only by the medium signature (the W2 bounds hash), so the resource stores one
# SBML per (parent, medium signature) and a line manifest mapping every line to its model file.
PIPELINE_P1 = "P1"
P1_DESCRIPTION = ("P1 (RESEARCH_TREE D11): generic Human-GEM parent + W3 B-loo:cofactor_pool_biomass (the cofactor/"
                  "vitamin pool species removed from the biomass reaction MAR13082) on the W2 documented medium "
                  "M2(f=1) (base medium + supplements + serum model at f = 1), gate-free; no context family, no "
                  "pruning, no repaired reactions; objective MAR13082")
P1_FLAG_ORDER = ("W3-D1", "W3-MEV", "W2-SHUF", "W2-MOPEN")
LICENCE = ("CC BY 4.0, inherited from Human-GEM. Medium bounds come from the W2 documented-media library "
           "(vendor sheets; no HPLM/Plasmax composition is used).")
FALLBACK_STRATA = ("imputed", "fallback_sheet_missing", "unmappable")
REF_WT_REL_TOL = 1e-5  # the M6 sweep tables store WT growth at 6 significant digits
P1_MODEL_COLUMNS = MANIFEST_COLUMNS + [
    "w2_sweep_id", "medium_signature", "medium_keys", "serum", "n_uptake_open", "n_lines", "line_strata",
    "wt_growth", "m6_wt_growth", "mevalonate_core_lethal",
    "gate_libsbml", "gate_roundtrip", "gate_wt_growth", "gate_random_ko", "gate_viability", "gate_reload_m6",
    "job_status", "cpu_s"]
GATE_COLUMNS = {"gate_libsbml": "libsbml_no_new_errors", "gate_roundtrip": "roundtrip_identical",
                "gate_wt_growth": "wt_growth_match", "gate_random_ko": "random_ko_match",
                "gate_viability": "viability_tasks", "gate_reload_m6": "reload_vs_m6"}
LINE_COLUMNS = ["ModelID", "parent", "version", "pipeline_id", "model_file", "model_sha256", "w2_sweep_id",
                "medium_signature", "medium_key", "medium_text_used", "medium_source", "medium_stratum",
                "medium_fallback", "serum_free", "serum_pct_documented", "lineage", "cancer", "screened", "flags",
                "validation_status", "release_status"]


def p1_flags(mevalonate_core_lethal: bool) -> str:
    """Flags carried by a P1 model (D11): W3-D1, W2-SHUF and W2-MOPEN always; W3-MEV where the model's own KO
    sweep shows the mevalonate core (HMGCR, MVK, PMVK, MVD) non-lethal (every serum medium; A11b)."""
    flags = set(P1_FLAG_ORDER) - ({"W3-MEV"} if mevalonate_core_lethal else set())
    return ";".join(f for f in P1_FLAG_ORDER if f in flags)


def build_medium_variant_model(parent: cobra.Model, bounds: dict, variant) -> cobra.Model:
    """A copy of ``parent`` with uptake bounds ``bounds`` (cobra medium: every other uptake closed, secretion
    unchanged), the biomass variant applied permanently and the biomass objective; the same steps, in the same
    order, as ``nextgem.w3.ko_sweep`` applies them inside its context."""
    from nextgem.w3 import BIOMASS

    m = parent.copy()
    m.medium = dict(bounds)
    variant.apply(m)
    m.objective = BIOMASS
    return m


def n_closed(model: cobra.Model) -> int:
    """Reactions with bounds (0, 0) (kept in the file, never deleted)."""
    return sum(1 for r in model.reactions if r.lower_bound == 0 and r.upper_bound == 0)


def ko_ratio_vector(model: cobra.Model, genes: list[str]) -> tuple[float, np.ndarray]:
    """WT growth and KO/WT growth ratios for ``genes`` (``Gene.knock_out()``; 6 decimals, clipped to [0, 1]: the
    W2/W3/M6 sweep convention). Every gene is solved explicitly (no shortcut)."""
    from nextgem.knockouts import RATIO_DECIMALS, WT_GROWTH_TOL

    wt = _growth(model)
    out = np.full(len(genes), math.nan)
    if not wt > WT_GROWTH_TOL:
        return wt, out
    for i, g in enumerate(genes):
        with model:
            model.genes.get_by_id(g).knock_out()
            out[i] = _growth(model) / wt
    return wt, np.clip(np.round(out, RATIO_DECIMALS), 0.0, 1.0) + 0.0


def reload_sample(ref_ratio: dict, n_partial: int = 25, n_one: int = 15, seed: int = SEEDS[0]) -> list[str]:
    """Deterministic gene sample for the reload check: up to ``n_partial`` genes whose reference KO ratio is < 1
    and ``n_one`` genes at ratio 1 (sorted ids, ``seed``)."""
    rng = np.random.default_rng(seed)
    ids = sorted(g for g, v in ref_ratio.items() if not (isinstance(v, float) and math.isnan(v)))
    low = [g for g in ids if ref_ratio[g] < 1]
    one = [g for g in ids if ref_ratio[g] >= 1]
    pick = [low[i] for i in sorted(rng.choice(len(low), size=min(n_partial, len(low)), replace=False))] if low else []
    pick += [one[i] for i in sorted(rng.choice(len(one), size=min(n_one, len(one)), replace=False))] if one else []
    return pick


def reload_check(source: cobra.Model, b: cobra.Model, ref_ratio: dict, ref_wt: float, genes: list[str]
                 ) -> tuple[bool, dict]:
    """Reloaded model ``b`` vs in-memory ``source`` vs the M6 sweep (``ref_ratio``: gene -> ratio_r; ``ref_wt``):
    WT growth (source = b within REL_TOL; = M6 within REF_WT_REL_TOL, the stored precision) and the KO ratio
    vector on ``genes`` (all three equal within 1e-6, i.e. at the sweeps' 6-decimal rounding)."""
    wa, ra = ko_ratio_vector(source, genes)
    wb, rb = ko_ratio_vector(b, genes)
    ref = np.array([ref_ratio[g] for g in genes], dtype=float)
    ok_wt_ab = _rel_equal(wa, wb)
    ok_wt_ref = (not math.isnan(wb)) and abs(wb - ref_wt) <= REF_WT_REL_TOL * max(abs(ref_wt), 1e-12)
    d_ab = np.abs(ra - rb)
    d_ref = np.abs(rb - ref)
    ok_ab = bool(np.all(d_ab <= 1e-6))
    ok_ref = bool(np.all(d_ref <= 1e-6))
    bad = [(g, float(x), float(y), float(z)) for g, x, y, z in zip(genes, ra, rb, ref)
           if not (abs(x - y) <= 1e-6 and abs(y - z) <= 1e-6)]
    detail = {"wt_source": wa, "wt_reloaded": wb, "wt_m6": ref_wt, "n_genes": len(genes),
              "n_partial": int((ref < 1).sum()), "max_abs_diff_source_reloaded": float(np.nanmax(d_ab)) if len(genes) else 0.0,
              "max_abs_diff_reloaded_m6": float(np.nanmax(d_ref)) if len(genes) else 0.0, "mismatches": bad[:5]}
    return ok_wt_ab and ok_wt_ref and ok_ab and ok_ref, detail


def gate_cells(val: Validation) -> dict:
    """pass / FAIL / skipped per manifest gate column."""
    out = {}
    for col, key in GATE_COLUMNS.items():
        ok = val.gates.get(key, (None,))[0]
        out[col] = "pass" if ok is True else ("skipped" if ok == "skipped" else ("FAIL" if key in val.gates else "not_run"))
    return out


def line_manifest(lines, models):
    """One row per (line, parent): the line's own medium metadata + its model file (``models``: one row per
    (parent, w2_sweep_id) with ``file``, ``sha256``, ``version``, ``pipeline_id``, ``medium_signature``, ``flags``,
    ``validation_status``, ``release_status``). Every line must map to exactly one model."""
    import pandas as pd

    need = ["ModelID", "parent", "w2_sweep_id"]
    if lines.duplicated(["ModelID", "parent"]).any():
        raise ResourceError("a (line, parent) pair occurs twice in the line table")
    mcols = ["parent", "w2_sweep_id", "file", "sha256", "version", "pipeline_id", "medium_signature", "flags",
             "validation_status", "release_status"]
    if models.duplicated(["parent", "w2_sweep_id"]).any():
        raise ResourceError("two model files for one (parent, medium signature)")
    m = lines.merge(models[mcols], on=["parent", "w2_sweep_id"], how="left", validate="many_to_one")
    miss = m.file.isna()
    if miss.any():
        raise ResourceError(f"{int(miss.sum())} (line, parent) pairs have no model file, e.g. "
                            f"{m.loc[miss, need].head(3).to_dict('records')}")
    m = m.rename(columns={"file": "model_file", "sha256": "model_sha256"})
    m["medium_fallback"] = m.medium_stratum.isin(FALLBACK_STRATA).map(
        {True: "classic RPMI + 10% FBS (W2 integrated-pipeline fallback)", False: ""})
    for c in LINE_COLUMNS:
        if c not in m:
            m[c] = ""
    # Human2 rows first (DESIGN W8 release order), then by ModelID
    return m[LINE_COLUMNS].sort_values(["parent", "ModelID"], ascending=[False, True]).reset_index(drop=True)


LINE_NOTE_FIELDS = ("medium_key", "medium_text_used", "medium_source", "medium_stratum", "medium_fallback",
                    "lineage", "cancer")


def load_line(resource_dir: Path, model_id: str, parent: str, manifest: str = "manifest.tsv") -> cobra.Model:
    """The P1 model of one line (read from its shared per-medium file, sha256-checked) with the line's own
    metadata in its notes and id."""
    import pandas as pd

    resource_dir = Path(resource_dir)
    man = pd.read_csv(resource_dir / manifest, sep="\t", dtype=str, keep_default_na=False)
    row = man[(man.ModelID == model_id) & (man.parent == parent)]
    if len(row) != 1:
        raise ResourceError(f"{model_id}/{parent}: {len(row)} manifest rows")
    row = row.iloc[0]
    path = resource_dir / row.model_file
    if sha256_file(path) != row.model_sha256:
        raise ResourceError(f"{path.name}: sha256 differs from the manifest")
    m = load_checked(path)
    notes = read_notes(m)
    notes.update({"shared_model_file": row.model_file, "shared_model_sha256": row.model_sha256,
                  "shared_model_id": notes.get("model_id", ""), "model_id": model_id,
                  **{f"line_{k}": str(row[k]) for k in LINE_NOTE_FIELDS}})
    m.notes.clear()
    m.notes.update(notes)
    m.id = f"{model_id}_{parent}_{row.pipeline_id}".replace("-", "_")
    return m


def materialise_line(resource_dir: Path, model_id: str, parent: str, out_dir: Path) -> Path:
    """Write a per-line SBML ``{ModelID}_{parent}-{version}_{pipeline}.xml.gz`` (DESIGN W8 naming) on demand."""
    import pandas as pd

    m = load_line(resource_dir, model_id, parent)
    man = pd.read_csv(Path(resource_dir) / "manifest.tsv", sep="\t", dtype=str, keep_default_na=False)
    row = man[(man.ModelID == model_id) & (man.parent == parent)].iloc[0]
    path = Path(out_dir) / f"{model_id}_{parent}-{row.version}_{row.pipeline_id}.xml.gz"
    write_sbml_gz(m, path)
    return path
