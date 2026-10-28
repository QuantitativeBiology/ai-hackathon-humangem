"""W8: build, export and validate the P1 resource (DESIGN §9 W8; RESEARCH_TREE D11). Label-free.

P1 = generic Human-GEM parent + W3 B-loo:cofactor_pool_biomass on the W2 documented medium M2(f=1), gate-free, no
context family. P1 models differ across lines only by their medium signature (the W2 bounds hash; 20 per parent),
so the resource stores ONE gzipped SBML per (parent, medium signature) plus a line manifest:

  resource/v{X.Y.Z}/medium-{signature}_{parent}-{version}_P1.xml.gz   one model per (parent, medium signature)
  resource/models.tsv     one row per model file: sha256, size, counts, flags, the 5 gates + reload gate
  resource/manifest.tsv   one row per (line, parent): the line's own medium metadata and its model file
  results/w8/p1/gates.jsonl, summary.json, MANIFEST row (results/w8/MANIFEST.tsv)

A per-line SBML (DESIGN naming {ModelID}_{parent}-{version}_P1.xml.gz) is materialised on demand by
scripts/w8_line_sbml.py; it equals the shared file plus the line's notes.

Inputs (all label-free): results/m6/scores/line_assignment.tsv (line -> W2 sweep id at M2(f=1)), the W2 media
library and floor tables (bounds rebuilt by ``w2_engine.build_sweeps`` through ``m6_sweeps.bounds_by_id``; the
rebuilt bounds must hash to the W2 signature), the W3 variant spec, the M6 sweep catalog / ratio_r (reference for the
reload gate), results/m6/scores/w3_mev.tsv (W3-MEV flag), Model.csv (lineage / cancer columns only).

Gates per model (resource.validate_model): libSBML no new errors vs the raw parent; cobra round-trip identical
sets/bounds/GPRs; WT growth; 20 random KOs; viability tasks (GR + 9 ER of the parent's essential-task file); plus
``reload_vs_m6``: the reloaded model's WT growth and KO ratio vector on a 40-gene sample equal the in-memory model and
the M6 sweep of that (parent, variant, signature).

Solver and limits (A10): GLPK; 300 s per LP + one retry (``w3.solve_status``); 1,800 s CPU per model job
(``ITIMER_PROF``, checked between Python calls). A timeout is a failed build, recorded, never re-run with relaxed
limits.

Usage:
  .venv/bin/python scripts/w8_build_p1.py --pilot 1 --workers 2          # 1 model per parent, results/w8/tmp/pilot
  PYTHONPATH=$PWD/src <venv>/python scripts/w8_build_p1.py --workers 8 --release 1.0.0   # pinned clean worktree
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import multiprocessing
import os
import signal
import sys
import time
from collections import Counter
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as mc  # noqa: E402

PROJECT = Path(__file__).resolve().parents[1]
W8 = PROJECT / "results" / "w8"
MANIFEST = W8 / "MANIFEST.tsv"
AUTHOR = "W8 resource engineer agent (Claude Opus)"
MAX_WORKERS = 8
JOB_CPU_LIMIT_S = 1800
PARENT_ORDER = ("human2", "human1")  # DESIGN W8: Human2-parent models are released first
_W: dict = {}


class JobCPUTimeout(Exception):
    pass


def _on_prof(signum, frame):
    raise JobCPUTimeout()


def _init(parent: str):
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem.model import load_model

    _W.update(parent=parent, model=load_model(parent))
    signal.signal(signal.SIGPROF, _on_prof)


def _parent_errors(parent: str) -> dict:
    from nextgem import resource as rs
    from nextgem.model import GEMS

    return dict(rs.libsbml_errors(GEMS[parent].sbml))


def _job(job: dict) -> dict:
    from nextgem import media, resource as rs, w3
    from nextgem.model import GEMS

    assert job["parent"] == _W["parent"]
    t0, c0 = time.time(), time.process_time()
    out = {k: job[k] for k in ("parent", "w2_sweep_id", "signature")}
    out.update(status="ok", row=None, gates=None)
    signal.setitimer(signal.ITIMER_PROF, float(job["cpu_limit"]))
    try:
        if media.bounds_signature(job["bounds"]) != job["signature"]:
            raise RuntimeError(f"{job['w2_sweep_id']}: rebuilt bounds do not hash to the W2 signature")
        variant = w3.Variant.from_json(job["variant_json"])
        m = rs.build_medium_variant_model(_W["model"], job["bounds"], variant)
        meta = rs.Meta(**job["meta"])
        meta.build_time_s = time.process_time() - c0  # notes: CPU-s to build the model (copy, medium, variant)
        path = rs.export_model(m, meta, Path(job["out_dir"]))
        genes = rs.reload_sample(job["ref_ratio"])
        extra = {"reload_vs_m6": lambda a, b: rs.reload_check(a, b, job["ref_ratio"], job["ref_wt"], genes)}
        val = rs.validate_model(m, path, None, GEMS[job["parent"]].essential_tasks,
                                parent_errors=Counter(job["parent_errors"]), extra_gates=extra)
        row = rs.manifest_row(path, meta, m, val, rs.n_closed(m))
        row["file"] = str(Path(path).relative_to(Path(job["resource_root"])))
        rl = val.gates.get("reload_vs_m6", (None, {}))[1]
        row.update(job["model_info"], wt_growth=rl.get("wt_source", float("nan")), m6_wt_growth=job["ref_wt"],
                   **rs.gate_cells(val))
        out.update(row=row, gates={k: [v[0], v[1]] for k, v in val.gates.items()})
    except JobCPUTimeout:
        out["status"] = "failed-timeout"
    except Exception as e:  # noqa: BLE001 - any build/export failure is a failed build, recorded
        out["status"] = f"failed-error: {type(e).__name__}: {str(e)[:300]}"
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0.0)
    out.update(cpu_s=time.process_time() - c0, wall_s=time.time() - t0)
    return out


# ------------------------------------------------------------------------------------------------- main side
def line_table() -> pd.DataFrame:
    """(line, parent) -> W2 sweep id at M2(f=1) + the line's medium metadata (label-free)."""
    from nextgem import paths

    la = pd.read_csv(mc.SCDIR / "line_assignment.tsv", sep="\t")
    la = la[la.step == mc.M2][["ModelID", "parent", "w2_sweep_id", "medium_key", "stratum"]]
    lm = pd.read_csv(PROJECT / "results" / "w2" / "media_library" / "line_media.tsv", sep="\t",
                     usecols=["ModelID", "lineage", "medium_text_used", "medium_source", "serum_free", "serum_pct",
                              "medium_key", "stratum"])
    chk = la.merge(lm[["ModelID", "medium_key", "stratum"]], on="ModelID", suffixes=("", "_w2"))
    if len(chk) != len(la) or (chk.medium_key != chk.medium_key_w2).any() or (chk.stratum != chk.stratum_w2).any():
        raise RuntimeError("line_assignment.tsv and line_media.tsv disagree on medium key / stratum")
    model_csv = paths.data_root() / "depmap" / "Model.csv"
    mm = pd.read_csv(model_csv, usecols=["ModelID", "OncotreeLineage", "OncotreePrimaryDisease"]).set_index("ModelID")
    # DESIGN W8 resource scope rule (scripts/design_numbers.py): non-cancer = Fibroblast/Normal lineage or a
    # "Non-Cancerous" primary disease
    noncancer = (mm.OncotreeLineage.isin(["Fibroblast", "Normal"])
                 | mm.OncotreePrimaryDisease.astype(str).str.contains("Non-Cancerous", case=False))
    t = la.merge(lm.drop(columns=["medium_key", "stratum"]), on="ModelID", how="left")
    t = t.rename(columns={"stratum": "medium_stratum", "serum_pct": "serum_pct_documented"})
    t["cancer"] = t.ModelID.map(~noncancer).map({True: "true", False: "false"}).fillna("unknown")
    t["screened"] = "true"  # line_assignment covers the 1,361 CRISPR-screened lines (W2/M6 cohort)
    return t, model_csv


def reference(parent: str, sweep_ids: list[str], m6_sweeps) -> tuple[dict, dict]:
    """M6 reference per W2 sweep id: ratio_r (gene -> value) and WT growth of the B-loo sweep."""
    cat = pd.read_csv(mc.SWDIR / "catalog.tsv", sep="\t")
    cat = cat[(cat.parent == parent) & (cat.variant == mc.BLOO) & (cat.step == mc.M2)].set_index("w2_sweep_id")
    reused = m6_sweeps.reused_rows()[parent]
    new = pd.read_parquet(mc.SWDIR / f"ratio_r_{parent}.parquet").set_index("sweep")
    ratios, wts = {}, {}
    for sid in sweep_ids:
        k = mc.sweep_key(mc.BLOO, sid)
        row = new.loc[k] if k in new.index else reused.loc[k]
        if cat.loc[sid, "status"] != "ok":
            raise RuntimeError(f"{k}: M6 sweep status {cat.loc[sid, 'status']}")
        ratios[sid] = {g: float(v) for g, v in row.items()}
        wts[sid] = float(cat.loc[sid, "wt_growth"])
    return ratios, wts


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--parents", nargs="+", default=list(PARENT_ORDER), choices=PARENT_ORDER)
    ap.add_argument("--release", default="1.0.0", help="resource version X.Y.Z (directory resource/vX.Y.Z)")
    ap.add_argument("--release-status", default="release candidate (P1 per D11; built before T1)")
    ap.add_argument("--pilot", type=int, default=0, help="first N signatures per parent (most lines first); "
                    "outputs under results/w8/tmp/pilot, tagged 'pilot, not a release'")
    ap.add_argument("--allow-dirty", action="store_true", help="pilots only")
    a = ap.parse_args(argv)
    if a.workers > MAX_WORKERS:
        raise SystemExit(f"at most {MAX_WORKERS} worker processes")
    if a.allow_dirty and not a.pilot:
        raise SystemExit("--allow-dirty is for pilots only (DESIGN §11.4: runs > 1 CPU-h come from a pinned clean "
                         "worktree)")
    import nextgem

    if not a.pilot and PROJECT / "src" not in Path(nextgem.__file__).resolve().parents:
        raise SystemExit(f"imported nextgem ({nextgem.__file__}) is not the package next to this script; run "
                         "with PYTHONPATH=<worktree>/src (pinned-worktree guard, §11.4)")
    from nextgem.model import GEMS
    from nextgem.run import provenance

    shared, per_parent = input_sets(a.parents)
    inputs = shared + [f for p in a.parents for f in per_parent[p]]
    with provenance(MANIFEST, inputs=inputs, seeds=[], solver="glpk", workers=a.workers,
                    allow_dirty=a.allow_dirty, author=AUTHOR) as prov:
        rc = run(a, prov, shared, per_parent)
    return rc


def input_sets(parents) -> tuple[list[Path], dict]:
    """Shared inputs (line assignment, media library, floor, M6 catalog/MEV, Model.csv) and per-parent inputs
    (parent SBML, task file, W3 variant spec, the reference sweeps)."""
    from nextgem import paths
    from nextgem.model import GEMS

    w2 = PROJECT / "results" / "w2"
    shared = [mc.SCDIR / "line_assignment.tsv", mc.SCDIR / "w3_mev.tsv", mc.SWDIR / "catalog.tsv",
              w2 / "media_library" / "line_media.tsv", w2 / "media_library" / "compositions.tsv",
              w2 / "media_library" / "exchange_map.tsv", w2 / "controls" / "line_sweeps.tsv",
              w2 / "floor" / "floor.tsv", w2 / "floor" / "requirements.tsv", w2 / "floor" / "floor_mult.tsv",
              paths.data_root() / "depmap" / "Model.csv"]
    per = {p: [GEMS[p].sbml, GEMS[p].essential_tasks, PROJECT / "results" / "w3" / "structure" / f"variants_{p}.json",
               mc.SWDIR / f"ratio_r_{p}.parquet", PROJECT / "results" / "w3" / "sweeps" / f"ratio_r_{p}.parquet",
               w2 / "ko" / f"ratio_r_{p}.parquet"] for p in parents}
    return shared, per


def _rel(p: Path) -> str:
    """Project-relative path for the notes (data/ may be a symlink into the main tree in a pinned worktree)."""
    from nextgem import paths

    p = Path(p)
    for base in (PROJECT, PROJECT.resolve(), paths.data_root().parent):
        try:
            return str(p.relative_to(base))
        except ValueError:
            continue
    return str(p)


def run(a, prov, shared, per_parent) -> int:
    import m6_sweeps

    from nextgem import resource as rs
    from nextgem.model import GEMS
    from nextgem.run import git_commit, git_dirty, sha256_file

    t_start = time.time()
    pilot = bool(a.pilot)
    res_root = (W8 / "tmp" / "pilot") if pilot else (PROJECT / "resource")
    out_dir = res_root / f"v{a.release}"
    release_status = "pilot, not a release" if pilot else a.release_status
    lines, _model_csv = line_table()
    lines = lines[lines.parent.isin(a.parents)]
    mev = pd.read_csv(mc.SCDIR / "w3_mev.tsv", sep="\t").set_index("key")
    commit, dirty = git_commit(), git_dirty()
    shared_shas = {_rel(f): sha256_file(f) for f in shared}

    need = {p: {mc.M2: sorted(lines[lines.parent == p].w2_sweep_id.unique())} for p in a.parents}
    bounds, info = m6_sweeps.bounds_by_id(need)
    jobs = []
    for p in a.parents:
        sub = lines[lines.parent == p]
        counts = sub.w2_sweep_id.value_counts()
        ids = sorted(need[p][mc.M2], key=lambda s: (-counts[s], s))
        if a.pilot:
            ids = ids[: a.pilot]
        ratios, wts = reference(p, ids, m6_sweeps)
        vj = m6_sweeps.variant_json(p, mc.BLOO)
        ver = GEMS[p].version.lstrip("v")
        p_shas = {**shared_shas, **{_rel(f): sha256_file(f) for f in per_parent[p]}}
        for sid in ids:
            s = sub[sub.w2_sweep_id == sid]
            sig = sid.split(":", 1)[1]
            strata = s.medium_stratum.value_counts().sort_index().to_dict()
            keys = sorted(set(s.medium_key))
            core_lethal = bool(mev.loc[mc.sweep_key(mc.BLOO, sid), "core_lethal"])
            flags = rs.p1_flags(core_lethal)
            tier = ";".join(f"{k}:{v}" for k, v in strata.items())
            extra = {
                "pipeline": rs.P1_DESCRIPTION, "human_gem_release": f"Human-GEM {GEMS[p].version}",
                "medium_step": "M2(f=1), gate-free", "medium_signature": sig, "w2_sweep_id": sid,
                "medium_keys": ";".join(keys), "serum": str(bool(info[sid]["serum"])).lower(),
                "n_uptake_open": str(len(bounds[sid])), "w3_variant": json.dumps(vj, sort_keys=True),
                "n_lines": str(len(s)), "line_strata": tier, "git_dirty": str(bool(dirty)).lower(),
                "licence": rs.LICENCE, "lines_manifest": "manifest.tsv (one row per line; this file is shared by "
                "every line with this medium signature)",
                "mevalonate_core_lethal": str(core_lethal).lower(),
            }
            meta = dict(model_id=f"medium-{sig}", parent=p, parent_version=ver, pipeline_id=rs.PIPELINE_P1,
                        medium=f"M2(f=1) {';'.join(keys)}", medium_source_tier=tier, git_commit=commit,
                        input_shas=p_shas, repaired_reactions=[], n_repaired_reactions=0,
                        repaired_note="none (generic parent; no context family, nothing pruned or re-opened)",
                        benchmarked=True, flags=flags, release_status=release_status, extra=extra)
            model_info = {"w2_sweep_id": sid, "medium_signature": sig, "medium_keys": ";".join(keys),
                          "serum": str(bool(info[sid]["serum"])).lower(), "n_uptake_open": len(bounds[sid]),
                          "n_lines": len(s), "line_strata": tier, "mevalonate_core_lethal": str(core_lethal).lower()}
            jobs.append(dict(parent=p, w2_sweep_id=sid, signature=sig, bounds=bounds[sid], variant_json=vj,
                             meta=meta, model_info=model_info, ref_ratio=ratios[sid], ref_wt=wts[sid],
                             out_dir=str(out_dir), resource_root=str(res_root), cpu_limit=JOB_CPU_LIMIT_S))
    per = ", ".join(f"{p}: {sum(j['parent'] == p for j in jobs)}" for p in a.parents)
    print(f"{len(jobs)} model jobs ({per}); {len(lines)} (line, parent) rows; out {out_dir}", flush=True)

    ctx = multiprocessing.get_context("spawn")
    with cf.ProcessPoolExecutor(max_workers=min(2, len(a.parents)), mp_context=ctx) as pool:
        perr = dict(zip(a.parents, pool.map(_parent_errors, a.parents)))
    print("parent libSBML errors:", perr, f"{time.time() - t_start:.0f} s", flush=True)
    res = []
    for p in a.parents:
        pj = [j for j in jobs if j["parent"] == p]
        for j in pj:
            j["parent_errors"] = perr[p]
        with cf.ProcessPoolExecutor(max_workers=min(a.workers, len(pj)), mp_context=ctx, initializer=_init,
                                    initargs=(p,)) as pool:
            for n, r in enumerate(pool.map(_job, pj, chunksize=1), 1):
                res.append(r)
                row = r["row"] or {}
                print(f"  [{p}] {n}/{len(pj)} {r['w2_sweep_id']}: {r['status']} {row.get('validation_status', '-')} "
                      f"{row.get('validation_detail', '')} {r['cpu_s']:.0f} s CPU, {time.time() - t_start:.0f} s",
                      flush=True)

    # ---- tables
    rows = []
    for r in res:
        if r["row"] is not None:
            row = dict(r["row"])
        else:  # failed build: no file
            j = next(j for j in jobs if j["parent"] == r["parent"] and j["w2_sweep_id"] == r["w2_sweep_id"])
            row = {c: "" for c in rs.P1_MODEL_COLUMNS}
            row.update(j["model_info"], parent=r["parent"], pipeline_id=rs.PIPELINE_P1,
                       version=j["meta"]["parent_version"], validation_status="fail",
                       validation_detail=r["status"], release_status=release_status)
        row.update(job_status=r["status"], cpu_s=round(r["cpu_s"], 1))
        rows.append(row)
    models = pd.DataFrame(rows)[rs.P1_MODEL_COLUMNS]
    models["_o"] = models.parent.map({p: i for i, p in enumerate(PARENT_ORDER)})
    models = models.sort_values(["_o", "w2_sweep_id"]).drop(columns="_o")
    built = models[models.file != ""]
    sel = lines[lines.w2_sweep_id.isin(set(built.w2_sweep_id))]
    man = rs.line_manifest(sel, built)
    outs = [mc.atomic_write_tsv(models, res_root / "models.tsv", float_format=None),
            mc.atomic_write_tsv(man, res_root / "manifest.tsv", float_format=None)]
    gdir = (W8 / "tmp" / "pilot") if pilot else (W8 / "p1")
    gl = "\n".join(json.dumps({"parent": r["parent"], "w2_sweep_id": r["w2_sweep_id"], "status": r["status"],
                               "gates": r["gates"]}, default=str, sort_keys=True) for r in res) + "\n"
    outs.append(mc.atomic_write_text(gl, gdir / "gates.jsonl"))
    files = [res_root / f for f in built.file]
    disk = sum(f.stat().st_size for f in files)
    summ = {
        "release": a.release, "release_status": release_status, "pilot": a.pilot, "parents": a.parents,
        "n_models": int(len(models)), "n_models_built": int(len(built)),
        "n_models_pass": int((models.validation_status == "pass").sum()),
        "per_parent": {p: {"models": int((models.parent == p).sum()),
                           "pass": int(((models.parent == p) & (models.validation_status == "pass")).sum()),
                           "lines": int((man.parent == p).sum())} for p in a.parents},
        "gate_pass": {c: int((models[c] == "pass").sum()) for c in rs.GATE_COLUMNS},
        "job_status": models.job_status.value_counts().to_dict(),
        "n_line_rows": int(len(man)), "n_lines": int(man.ModelID.nunique()),
        "model_bytes": int(disk), "job_cpu_s": round(float(sum(r["cpu_s"] for r in res)), 1),
        "wall_s": round(time.time() - t_start, 1), "workers": a.workers, "solver": "glpk",
        "a10_limits": {"lp_s": 300, "lp_retries": 1, "job_cpu_s": JOB_CPU_LIMIT_S},
        "commit": commit, "dirty": bool(dirty),
    }
    outs.append(mc.atomic_write_text(json.dumps(summ, indent=1) + "\n", gdir / "summary.json"))
    outs.append(mc.atomic_write_text(resource_readme(summ, models, man), res_root / "README.md"))
    prov.output(*outs, *files)
    prov.note(**{k: v for k, v in summ.items() if k not in ("a10_limits", "per_parent", "gate_pass")})
    print(json.dumps(summ, indent=1), flush=True)
    rc = rs.gate_exit_code(models)
    print(f"exit code {rc} ({summ['n_models_pass']}/{summ['n_models']} models pass every gate)", flush=True)
    return rc


def resource_readme(summ: dict, models: pd.DataFrame, man: pd.DataFrame) -> str:
    """resource/README.md: layout, how to get a line's model, licence (DESIGN W8 data SF11)."""
    from nextgem import resource as rs

    per = "\n".join(f"| {p} | {v['models']} | {v['pass']} | {v['lines']} |" for p, v in summ["per_parent"].items())
    scr = man.drop_duplicates("ModelID").screened.astype(str).str.lower().value_counts().to_dict()
    fails = models[models.validation_status != "pass"]
    fail_txt = "\n".join(f"- `{r.file}` ({r.medium_keys}; {r.n_lines} line(s)): failed `{r.validation_detail}`; "
                         "see results/w8/README.md" for r in fails.itertuples()) or "- none"
    size_mb = sum(int(x) for x in models.size_bytes if str(x)) / 1e6
    return f"""# nextgem P1 resource v{summ['release']} ({summ['release_status']})

{rs.P1_DESCRIPTION}.

P1 models differ across lines only by the line's medium signature, so this resource holds **one SBML model per
(parent, medium signature)** and a manifest that maps every line to its model:

- `v{summ['release']}/medium-<signature>_<parent>-<Human-GEM version>_P1.xml.gz`: SBML L3V1 FBC v2, gzipped. Notes
  record the pipeline, Human-GEM release, medium signature/keys/tiers, W3 variant, git commit, input SHA-256s and
  flags. Closed reactions are kept with bounds (0, 0).
- `models.tsv`: one row per model file (sha256, size, counts, flags, validation gates).
- `manifest.tsv`: one row per (line, parent): the line's own medium key, text, source and stratum (tier),
  fallback flag (imputed / sheet missing / unmappable lines get classic RPMI + 10% FBS), lineage, and its model file.

| Parent | Models | Gates passed | Lines |
|---|---|---|---|
{per}

Lines: {scr.get('true', 0)} CRISPR-screened lines (W2/M6 medium assignment; `screened = true`) and
{scr.get('false', 0)} unscreened DepMap 26Q3 models in the DESIGN W8 scope (default expression profile, cancer lineage;
medium by W2's rule from the Model.csv onboarded text, else the flagged classic RPMI + 10% FBS fallback;
`screened = false`). Model files: {len(models)}, {size_mb:.1f} MB.

Gates per model (models.tsv): libSBML no new errors vs the parent, cobra round-trip, WT growth, 20 random KOs,
viability tasks (GR + 9 ER), and the reloaded model's WT growth and 40-gene KO vector vs the reference KO sweep.
Models failing a gate:
{fail_txt}

A per-line file (`<ModelID>_<parent>-<version>_P1.xml.gz`, the shared model plus the line's notes) is written on
demand: `python scripts/w8_line_sbml.py --lines ACH-000001 --parent human2`; in Python,
`nextgem.resource.load_line("resource", "ACH-000001", "human2")`.

Flags (RESEARCH_TREE D11): W3-D1, W2-SHUF, W2-MOPEN on every model; W3-MEV on every serum medium (mevalonate core
non-lethal). See results/w8/README.md.

Licence: {rs.LICENCE} The resource contains no per-protein values and no DepMap screen data.
"""


if __name__ == "__main__":
    raise SystemExit(main())
