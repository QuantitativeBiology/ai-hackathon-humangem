"""W8: extend the P1 resource to the DepMap 26Q3 models outside the screened cohort (DESIGN W8 scope). Label-free.

DESIGN W8 scope: every 26Q3 model with a default expression profile and a cancer lineage (scripts/design_numbers.py
rule: not Fibroblast/Normal lineage, primary disease not "Non-Cancerous"). The P1 v1.0.0 build covered the 1,361
CRISPR-screened lines (the W2/M6 cohort); this script adds the in-scope models that are not screened.

Steps
  --step assign  W2's frozen assignment rule (``media.assign_line_media``, called on the target lines; its screen
                 branch cannot fire for unscreened lines, so a line gets its Model.csv onboarded text, else the
                 imputed classic RPMI + 10% FBS, with W2's strata/flags), medium key (``media.medium_key``), and
                 the M2(f=1) gate-free bounds per parent (``media.uptake_bounds``) -> W2 signature. Every key that
                 also occurs in the screened cohort must reproduce that cohort's signature (asserted).
                 -> results/w8/unscreened/{line_media_unscreened.tsv, signatures.tsv}
  --step sweeps  for each NEW (parent, signature): the B-loo KO reference sweep with the M6 engine
                 (``w3.ko_sweep``, A10 limits) -> results/w8/unscreened/{ratio_r_<parent>.parquet, sweeps.tsv}
  --step build   export + the 6 gates for the new models (``w8_build_p1._job``), then append to
                 resource/models.tsv and rebuild resource/manifest.tsv (screened rows unchanged) and README.

Existing model files are never rewritten (their notes count only the screened lines that share them; the manifest
is the authority on which lines use a file).

Usage: .venv/bin/python scripts/w8_extend_unscreened.py --step assign|sweeps|build|readme [--workers 8] [--allow-dirty]
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
from pathlib import Path
from unittest import mock

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as mc  # noqa: E402
import w8_build_p1 as wb  # noqa: E402

PROJECT = wb.PROJECT
OUT = PROJECT / "results" / "w8" / "unscreened"
RES = PROJECT / "resource"
EXPR_FILE = "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
_W: dict = {}


# ------------------------------------------------------------------------------------------------ assign
def design_scope() -> tuple[list[str], pd.DataFrame]:
    """DESIGN W8 scope (design_numbers.py rule): default-expression models in Model.csv with a cancer lineage."""
    from nextgem import paths

    dm = paths.depmap_dir()
    ex = pd.read_csv(dm / EXPR_FILE, usecols=["ModelID", "IsDefaultEntryForModel"])
    E = set(ex.loc[ex.IsDefaultEntryForModel == "Yes", "ModelID"])
    model = pd.read_csv(dm / "Model.csv", usecols=["ModelID", "OncotreeLineage", "OncotreePrimaryDisease"])
    model = model.set_index("ModelID")
    ma = model.loc[sorted(E & set(model.index))]
    noncancer = ma.OncotreeLineage.isin(["Fibroblast", "Normal"]) | ma.OncotreePrimaryDisease.astype(str).str.contains(
        "Non-Cancerous", case=False)
    return sorted(ma.index[~noncancer]), model


def assign(screened_sig: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """W2 assignment + medium key + M2(f=1) signature per parent for the in-scope unscreened models."""
    from nextgem import media

    scope, model = design_scope()
    screened = set(media.screened_lines())
    target = [m for m in scope if m not in screened]
    with mock.patch.object(media, "screened_lines", lambda: target):
        df = media.assign_line_media(set(media.SHEET_FILES))
    assert list(df.ModelID) == target and set(df.tier) <= {"onboarded", "imputed"}, set(df.tier)
    sheets = media.load_sheets()
    on_disk = set(media.SHEET_FILES)
    keys, meds = {}, {}
    for t in sorted(set(df.medium_text_used)):
        p = media.parse_formulation(t)
        if not (p.parts and all(b in on_disk for b, _ in p.parts)):
            raise RuntimeError(f"medium text used is not buildable from the sheets on disk: {t!r}")
        med = media.build_medium(p, sheets)
        k = media.medium_key(p, med, sheets)
        med.key = k
        keys[t] = k
        meds.setdefault(k, med)
    df["medium_key"] = df.medium_text_used.map(keys)
    w2 = PROJECT / "results" / "w2" / "floor"
    pcs = media.floor_constraints(w2 / "floor.tsv", w2 / "requirements.tsv", w2 / "floor_mult.tsv")
    sig_rows = []
    for p in wb.PARENT_ORDER:
        known = screened_sig[screened_sig.parent == p].groupby("medium_key").w2_sweep_id.agg(set)
        known_ids = set(screened_sig[screened_sig.parent == p].w2_sweep_id)
        for k, med in sorted(meds.items()):
            b = media.uptake_bounds(None, "M2", med, pcs[p], f=1.0, gate_open=True)
            sid = f"{p}:{media.bounds_signature(b)}"
            if k in known.index and known[k] != {sid}:
                raise RuntimeError(f"{p} {k}: rebuilt signature {sid} differs from the screened cohort's {known[k]}")
            sig_rows.append(dict(parent=p, medium_key=k, w2_sweep_id=sid, new=sid not in known_ids,
                                 serum=bool(med.serum), n_uptake_open=len(b), n_lines=int((df.medium_key == k).sum()),
                                 bounds_json=json.dumps(b, sort_keys=True)))
    sigs = pd.DataFrame(sig_rows)
    out = df[["ModelID", "lineage", "tier", "text", "base_class", "stratum", "medium_text_used", "medium_source",
              "serum_free", "serum_pct", "medium_key"]].copy()
    out["cancer"] = "true"
    out["screened"] = "false"
    return out, sigs


def screened_signatures() -> pd.DataFrame:
    la = pd.read_csv(mc.SCDIR / "line_assignment.tsv", sep="\t")
    return la[la.step == mc.M2][["ModelID", "parent", "w2_sweep_id", "medium_key"]]


# ------------------------------------------------------------------------------------------------ sweeps
def _init_sweep(parent: str):
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import w3
    from nextgem.model import load_model

    m = load_model(parent)
    _W.update(parent=parent, model=m, genes=sorted(g.id for g in m.genes), closures=w3.gene_closures(m))
    signal.signal(signal.SIGPROF, wb._on_prof)


def _sweep_job(job: dict) -> dict:
    from nextgem import w3
    from nextgem.knockouts import WT_GROWTH_TOL

    t0, c0 = time.time(), time.process_time()
    status, r = "ok", None
    signal.setitimer(signal.ITIMER_PROF, float(mc.JOB_CPU_LIMIT_S))
    try:
        r = w3.ko_sweep(_W["model"], _W["genes"], _W["closures"], job["bounds"], w3.Variant.from_json(job["variant_json"]))
    except wb.JobCPUTimeout:
        status = "failed-timeout"
    except w3.SolverInfraError as e:
        status = "failed-timeout" if "time" in str(e).lower() else "failed-nonoptimal"
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0.0)
    out = {k: job[k] for k in ("parent", "w2_sweep_id", "medium_key")}
    if r is not None:
        if r["wt_status"] != "optimal":
            status = "failed-infeasible"
        elif not r["wt_growth"] > WT_GROWTH_TOL:
            status = "failed-wt"
        out.update({k: r[k] for k in ("wt_growth", "wt_status", "n_solved", "n_genes_lp", "n_infeasible",
                                      "shortcut_verify_max_dev")})
    out.update(status=status, ratio_r=None if r is None else r["ratio_r"], genes=_W["genes"],
               cpu_s=time.process_time() - c0, wall_s=time.time() - t0)
    return out


# ------------------------------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", required=True, choices=("assign", "sweeps", "build", "readme"))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--release", default="1.0.0")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    if a.workers > wb.MAX_WORKERS:
        raise SystemExit(f"at most {wb.MAX_WORKERS} worker processes")
    import nextgem
    from nextgem import paths
    from nextgem.model import GEMS
    from nextgem.run import provenance

    if not a.allow_dirty and PROJECT / "src" not in Path(nextgem.__file__).resolve().parents:
        raise SystemExit(f"imported nextgem ({nextgem.__file__}) is not the package next to this script; run "
                         "with PYTHONPATH=<worktree>/src (pinned-worktree guard, §11.4)")
    shared, per = wb.input_sets(wb.PARENT_ORDER)
    inputs = shared + [paths.depmap_dir() / EXPR_FILE, paths.depmap_dir() / "ModelCondition.csv",
                       paths.depmap_dir() / "ScreenSequenceMap.csv", paths.depmap_dir() / "CRISPRScreenMap.csv"]
    if a.step != "assign":
        inputs += [OUT / "signatures.tsv", OUT / "line_media_unscreened.tsv"]
        inputs += [f for p in wb.PARENT_ORDER for f in per[p]]
    if a.step == "build":
        inputs += [OUT / "sweeps.tsv"] + [OUT / f"ratio_r_{p}.parquet" for p in wb.PARENT_ORDER
                                          if (OUT / f"ratio_r_{p}.parquet").exists()]
        inputs += [RES / "models.tsv", RES / "manifest.tsv"]
    with provenance(wb.MANIFEST, inputs=inputs, seeds=[], solver="glpk" if a.step != "assign" else "",
                    workers=a.workers if a.step != "assign" else 1, allow_dirty=a.allow_dirty,
                    author=wb.AUTHOR) as prov:
        rc = {"assign": run_assign, "sweeps": run_sweeps, "build": run_build, "readme": run_readme}[a.step](
            a, prov, GEMS, per)
    return rc


def run_assign(a, prov, GEMS, per) -> int:
    t0 = time.time()
    lines, sigs = assign(screened_signatures())
    outs = [mc.atomic_write_tsv(lines, OUT / "line_media_unscreened.tsv", float_format=None),
            mc.atomic_write_tsv(sigs, OUT / "signatures.tsv", float_format=None)]
    summ = {"n_lines": int(len(lines)), "tier": lines.tier.value_counts().to_dict(),
            "stratum": lines.stratum.value_counts().to_dict(), "n_keys": int(lines.medium_key.nunique()),
            "new_signatures": {p: int(sigs[(sigs.parent == p) & sigs.new].w2_sweep_id.nunique())
                               for p in wb.PARENT_ORDER},
            "lines_on_new_signatures": {p: int(sigs[(sigs.parent == p) & sigs.new].n_lines.sum())
                                        for p in wb.PARENT_ORDER}, "wall_s": round(time.time() - t0, 1)}
    outs.append(mc.atomic_write_text(json.dumps(summ, indent=1) + "\n", OUT / "assign_summary.json"))
    prov.output(*outs)
    prov.note(step="assign", **{k: v for k, v in summ.items() if k in ("n_lines", "n_keys", "new_signatures")})
    print(json.dumps(summ, indent=1))
    return 0


def _new_sigs() -> pd.DataFrame:
    s = pd.read_csv(OUT / "signatures.tsv", sep="\t")
    return s[s.new].drop_duplicates(["parent", "w2_sweep_id"])


def run_sweeps(a, prov, GEMS, per) -> int:
    import m6_sweeps
    import numpy as np

    t0 = time.time()
    new = _new_sigs()
    jobs = [dict(parent=r.parent, w2_sweep_id=r.w2_sweep_id, medium_key=r.medium_key, bounds=json.loads(r.bounds_json),
                 variant_json=m6_sweeps.variant_json(r.parent, mc.BLOO)) for r in new.itertuples()]
    print(f"{len(jobs)} new B-loo reference sweeps", flush=True)
    res = []
    ctx = multiprocessing.get_context("spawn")
    for p in wb.PARENT_ORDER:
        pj = [j for j in jobs if j["parent"] == p]
        if not pj:
            continue
        with cf.ProcessPoolExecutor(max_workers=min(a.workers, len(pj)), mp_context=ctx, initializer=_init_sweep,
                                    initargs=(p,)) as pool:
            for r in pool.map(_sweep_job, pj, chunksize=1):
                res.append(r)
                print(f"  [{p}] {r['w2_sweep_id']} {r['medium_key']}: {r['status']} WT {r.get('wt_growth', float('nan')):.6g} "
                      f"{r['cpu_s']:.0f} s", flush=True)
    meta = pd.DataFrame([{k: v for k, v in r.items() if k not in ("ratio_r", "genes")} for r in res])
    meta["key"] = [mc.sweep_key(mc.BLOO, s) for s in meta.w2_sweep_id]
    outs = [mc.atomic_write_tsv(meta, OUT / "sweeps.tsv")]  # %.6g: the M6 convention for WT growth
    for p in wb.PARENT_ORDER:
        rr = [r for r in res if r["parent"] == p and r["status"] == "ok"]
        if rr:
            R = pd.DataFrame(np.vstack([r["ratio_r"] for r in rr]), columns=rr[0]["genes"],
                             index=pd.Index([mc.sweep_key(mc.BLOO, r["w2_sweep_id"]) for r in rr], name="sweep"))
            outs.append(mc.atomic_write_parquet(R.reset_index(), OUT / f"ratio_r_{p}.parquet"))
    prov.output(*outs)
    summ = {"n_sweeps": len(res), "status": meta.status.value_counts().to_dict() if len(meta) else {},
            "job_cpu_s": round(float(meta.cpu_s.sum()), 1) if len(meta) else 0.0, "wall_s": round(time.time() - t0, 1)}
    prov.note(step="sweeps", **summ)
    print(json.dumps(summ, indent=1))
    return 0 if all(r["status"] == "ok" for r in res) else 3


def run_build(a, prov, GEMS, per) -> int:
    from collections import Counter  # noqa: F401

    from nextgem import resource as rs
    from nextgem.run import git_commit, git_dirty, sha256_file

    t0 = time.time()
    lines = pd.read_csv(OUT / "line_media_unscreened.tsv", sep="\t", dtype=str, keep_default_na=False)
    sigs = pd.read_csv(OUT / "signatures.tsv", sep="\t")
    new = _new_sigs()
    models = pd.read_csv(RES / "models.tsv", sep="\t", dtype=str, keep_default_na=False)
    old_man = pd.read_csv(RES / "manifest.tsv", sep="\t", dtype=str, keep_default_na=False)
    models, old_man = drop_previous_extension(models, old_man, new)
    release_status = models.release_status.iloc[0]
    commit, dirty = git_commit(), git_dirty()
    shared, _ = wb.input_sets(wb.PARENT_ORDER)
    extra_inputs = [OUT / "line_media_unscreened.tsv", OUT / "signatures.tsv", OUT / "sweeps.tsv"]
    shared_shas = {wb._rel(f): sha256_file(f) for f in shared + extra_inputs}
    sw = pd.read_csv(OUT / "sweeps.tsv", sep="\t").set_index("w2_sweep_id") if len(new) else None
    jobs = []
    for r in new.itertuples():
        p, sid = r.parent, r.w2_sweep_id
        if sw.loc[sid, "status"] != "ok":
            raise RuntimeError(f"{sid}: reference sweep status {sw.loc[sid, 'status']}")
        ref = pd.read_parquet(OUT / f"ratio_r_{p}.parquet").set_index("sweep").loc[mc.sweep_key(mc.BLOO, sid)]
        ref_ratio = {g: float(v) for g, v in ref.items()}
        sig = sid.split(":", 1)[1]
        s = lines[lines.medium_key.isin(set(sigs[(sigs.parent == p) & (sigs.w2_sweep_id == sid)].medium_key))]
        strata = s.stratum.value_counts().sort_index().to_dict()
        tier = ";".join(f"{k}:{v}" for k, v in strata.items())
        keys = sorted(set(s.medium_key))
        core_lethal = mevalonate_core_lethal(p, ref_ratio)
        vj = __import__("m6_sweeps").variant_json(p, mc.BLOO)
        bounds = json.loads(r.bounds_json)
        p_shas = {**shared_shas, **{wb._rel(f): sha256_file(f) for f in per[p]},
                  wb._rel(OUT / f"ratio_r_{p}.parquet"): sha256_file(OUT / f"ratio_r_{p}.parquet")}
        ver = GEMS[p].version.lstrip("v")
        extra = {"pipeline": rs.P1_DESCRIPTION, "human_gem_release": f"Human-GEM {GEMS[p].version}",
                 "medium_step": "M2(f=1), gate-free", "medium_signature": sig, "w2_sweep_id": sid,
                 "medium_keys": ";".join(keys), "serum": str(bool(r.serum)).lower(), "n_uptake_open": str(len(bounds)),
                 "w3_variant": json.dumps(vj, sort_keys=True), "n_lines": str(len(s)), "line_strata": tier,
                 "git_dirty": str(bool(dirty)).lower(), "licence": rs.LICENCE,
                 "lines_manifest": "manifest.tsv (one row per line; unscreened-cohort medium signature, W8 extension)",
                 "mevalonate_core_lethal": str(core_lethal).lower(),
                 "reference_sweep": "results/w8/unscreened (W8 extension; same engine as M6: w3.ko_sweep)"}
        meta = dict(model_id=f"medium-{sig}", parent=p, parent_version=ver, pipeline_id=rs.PIPELINE_P1,
                    medium=f"M2(f=1) {';'.join(keys)}", medium_source_tier=tier, git_commit=commit,
                    input_shas=p_shas, repaired_reactions=[], n_repaired_reactions=0,
                    repaired_note="none (generic parent; no context family, nothing pruned or re-opened)",
                    benchmarked=False, flags=rs.p1_flags(core_lethal), release_status=release_status, extra=extra)
        info = {"w2_sweep_id": sid, "medium_signature": sig, "medium_keys": ";".join(keys),
                "serum": str(bool(r.serum)).lower(), "n_uptake_open": len(bounds), "n_lines": len(s),
                "line_strata": tier, "mevalonate_core_lethal": str(core_lethal).lower()}
        jobs.append(dict(parent=p, w2_sweep_id=sid, signature=sig, bounds=bounds, variant_json=vj, meta=meta,
                         model_info=info, ref_ratio=ref_ratio, ref_wt=float(sw.loc[sid, "wt_growth"]),
                         out_dir=str(RES / f"v{a.release}"), resource_root=str(RES), cpu_limit=wb.JOB_CPU_LIMIT_S))
    res = []
    if jobs:
        from nextgem import resource as _rs  # noqa: F401
        ctx = multiprocessing.get_context("spawn")
        with cf.ProcessPoolExecutor(max_workers=2, mp_context=ctx) as pool:
            perr = dict(zip(wb.PARENT_ORDER, pool.map(wb._parent_errors, wb.PARENT_ORDER)))
        for p in wb.PARENT_ORDER:
            pj = [j for j in jobs if j["parent"] == p]
            for j in pj:
                j["parent_errors"] = perr[p]
            if not pj:
                continue
            with cf.ProcessPoolExecutor(max_workers=min(a.workers, len(pj)), mp_context=ctx, initializer=wb._init,
                                        initargs=(p,)) as pool:
                for r in pool.map(wb._job, pj, chunksize=1):
                    res.append(r)
                    row = r["row"] or {}
                    print(f"  [{p}] {r['w2_sweep_id']}: {r['status']} {row.get('validation_status', '-')} "
                          f"{row.get('validation_detail', '')} {r['cpu_s']:.0f} s", flush=True)
    rows = []
    for r in res:
        if r["row"] is not None:
            row = dict(r["row"])
        else:
            j = next(j for j in jobs if j["w2_sweep_id"] == r["w2_sweep_id"])
            row = {c: "" for c in rs.P1_MODEL_COLUMNS}
            row.update(j["model_info"], parent=r["parent"], pipeline_id=rs.PIPELINE_P1,
                       version=j["meta"]["parent_version"], validation_status="fail", validation_detail=r["status"],
                       release_status=release_status)
        row.update(job_status=r["status"], cpu_s=round(r["cpu_s"], 1))
        rows.append({k: str(v) for k, v in row.items()})
    all_models = pd.concat([models, pd.DataFrame(rows, columns=rs.P1_MODEL_COLUMNS)], ignore_index=True)
    built = all_models[all_models.file != ""]
    # unscreened lines: one row per parent, mapped by the parent's signature of the line's medium key
    sk = sigs[["parent", "medium_key", "w2_sweep_id"]]
    ul = lines.merge(sk, on="medium_key", how="left", validate="many_to_many")
    ul = ul.rename(columns={"stratum": "medium_stratum", "serum_pct": "serum_pct_documented"})
    new_man = rs.line_manifest(ul[["ModelID", "parent", "w2_sweep_id", "medium_key", "medium_text_used",
                                   "medium_source", "medium_stratum", "serum_free", "serum_pct_documented",
                                   "lineage", "cancer", "screened"]], built)
    man = pd.concat([old_man, new_man.astype(str)], ignore_index=True)
    if man.duplicated(["ModelID", "parent"]).any():
        raise RuntimeError("a line occurs in both cohorts")
    man["_o"] = man.parent.map({p: i for i, p in enumerate(wb.PARENT_ORDER)})
    man = man.sort_values(["_o", "ModelID"]).drop(columns="_o").reset_index(drop=True)
    all_models["_o"] = all_models.parent.map({p: i for i, p in enumerate(wb.PARENT_ORDER)})
    all_models = all_models.sort_values(["_o", "w2_sweep_id"]).drop(columns="_o")
    outs = [mc.atomic_write_tsv(all_models, RES / "models.tsv", float_format=None),
            mc.atomic_write_tsv(man, RES / "manifest.tsv", float_format=None)]
    gl = "\n".join(json.dumps({"parent": r["parent"], "w2_sweep_id": r["w2_sweep_id"], "status": r["status"],
                               "gates": r["gates"]}, default=str, sort_keys=True) for r in res) + "\n"
    outs.append(mc.atomic_write_text(gl, OUT / "gates.jsonl"))
    files = [RES / f for f in pd.DataFrame(rows, columns=rs.P1_MODEL_COLUMNS).file if f]
    summ = {"n_new_models": len(rows), "n_new_pass": sum(r["validation_status"] == "pass" for r in rows),
            "gate_pass": {c: sum(r.get(c) == "pass" for r in rows) for c in rs.GATE_COLUMNS},
            "n_unscreened_lines": int(lines.ModelID.nunique()), "n_line_rows_total": int(len(man)),
            "n_lines_total": int(man.ModelID.nunique()), "n_models_total": int(len(all_models)),
            "new_model_bytes": int(sum(f.stat().st_size for f in files)),
            "resource_model_bytes": int(sum((RES / f).stat().st_size for f in built.file)),
            "job_cpu_s": round(float(sum(r["cpu_s"] for r in res)), 1), "wall_s": round(time.time() - t0, 1),
            "commit": commit, "dirty": bool(dirty)}
    outs.append(mc.atomic_write_text(json.dumps(summ, indent=1) + "\n", OUT / "build_summary.json"))
    full = json.loads((PROJECT / "results" / "w8" / "p1" / "summary.json").read_text())
    per_parent = {p: {"models": int((all_models.parent == p).sum()),
                      "pass": int(((all_models.parent == p) & (all_models.validation_status == "pass")).sum()),
                      "lines": int((man.parent == p).sum())} for p in wb.PARENT_ORDER}
    outs.append(mc.atomic_write_text(wb.resource_readme(dict(full, per_parent=per_parent), all_models, man),
                                     RES / "README.md"))
    prov.output(*outs, *files)
    prov.note(step="build", **{k: v for k, v in summ.items() if k != "gate_pass"})
    print(json.dumps(summ, indent=1))
    new_rows = pd.DataFrame(rows, columns=rs.P1_MODEL_COLUMNS)
    rc = rs.gate_exit_code(new_rows) if len(new_rows) else 0
    print(f"exit code {rc} ({summ['n_new_pass']}/{summ['n_new_models']} new models pass every gate)", flush=True)
    return rc


def run_readme(a, prov, GEMS, per) -> int:
    """Regenerate resource/README.md from the current models.tsv / manifest.tsv (no model is touched)."""
    models = pd.read_csv(RES / "models.tsv", sep="\t", dtype=str, keep_default_na=False)
    man = pd.read_csv(RES / "manifest.tsv", sep="\t", dtype=str, keep_default_na=False)
    full = json.loads((PROJECT / "results" / "w8" / "p1" / "summary.json").read_text())
    per_parent = {p: {"models": int((models.parent == p).sum()),
                      "pass": int(((models.parent == p) & (models.validation_status == "pass")).sum()),
                      "lines": int((man.parent == p).sum())} for p in wb.PARENT_ORDER}
    out = mc.atomic_write_text(wb.resource_readme(dict(full, per_parent=per_parent), models, man), RES / "README.md")
    prov.output(out)
    return 0


def drop_previous_extension(models: pd.DataFrame, man: pd.DataFrame, new: pd.DataFrame):
    """Make --step build idempotent: remove the model rows of the signatures being (re)built and every earlier
    unscreened manifest row (screened = false); the screened cohort and its 40 models are kept as they are."""
    key = set(zip(new.parent, new.w2_sweep_id))
    keep = [(p, s) not in key for p, s in zip(models.parent, models.w2_sweep_id)]
    return models[keep].reset_index(drop=True), man[man.screened.astype(str).str.lower() != "false"].reset_index(drop=True)


def mevalonate_core_lethal(parent: str, ref_ratio: dict) -> bool:
    """W3-MEV rule of results/m6/scores/w3_mev.tsv: HMGCR, MVK, PMVK, MVD all lethal (ratio_r < 0.01)."""
    from nextgem.model import GEMS

    g = pd.read_csv(GEMS[parent].genes_tsv, sep="\t")
    sym = dict(zip(g.genes, g.geneSymbols))
    ids = [k for k, v in sym.items() if v in mc.MEVALONATE_CORE]
    if len(ids) != len(mc.MEVALONATE_CORE):
        raise RuntimeError(f"{parent}: mevalonate core genes not found once each")
    return all(ref_ratio[i] < 0.01 for i in ids)


if __name__ == "__main__":
    raise SystemExit(main())
