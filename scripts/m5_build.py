"""M5 / amendment A12(b): label-free builds for the panel pre-check of the context-model axis (DESIGN.md §17 A12 b).

On the 30 Stage-A panel lines (results/w4/stage_a_panel.tsv), with the parent = generic Human1 + the validated W3
variant B-loo:cofactor_pool_biomass (spec: results/w3/structure/variants_human1.json, applied with nextgem.w3.Variant):

  (1) generic   the parent itself (one line-constant KO sweep; checked against W3's sweep row M0|B-loo:cofactor_pool_biomass)
  (2) F2        GIMME, R1, RNA, real expression, no genomic constraint  (Stage A's F2 configuration carried to M6)
  (3) F4        E-Flux, R1, RNA, real expression, no genomic constraint (Stage A's full-panel-best F4 configuration)

Everything else is Stage A's (scripts/w4_stage_a.py, amendment A4): medium M0 (Ham's, from the parent's own GR task),
the W4 builders in nextgem.context / nextgem.context_run (extraction -> repair rule (F2) -> WT growth -> KO sweep of
all parent genes with the A4 KO gene-state convention -> 57-task profile -> open mask), declared solvers (GLPK for F2
and F4), one fresh spawned process per build, 300 s LP limit with one retry (context.LP_TIME_LIMIT_S), a per-job CPU
kill at 1,800 s (failed-timeout, A10). The only change from Stage A is the parent (Human1 + B-loo biomass instead of
Human2), as A12(b) specifies; R1 thresholds are global (TPM 1), so no Human2-fitted gene statistic is used.

No label is read. Outputs (default results/m5/builds/): builds.tsv, gates.tsv, ko.parquet (build_id, gene, ratio_r,
status, solver), build_summary.json (repair and gate statistics, generic check), MANIFEST rows in results/m5/MANIFEST.tsv.

Usage:
  .venv/bin/python scripts/m5_build.py --dry-run
  .venv/bin/python scripts/m5_build.py --pilot --lines 2 --out results/m5/tmp/pilot --workers 3     # cost pilot
  .venv/bin/python scripts/m5_build.py --workers 8                                                   # pinned worktree
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import psutil

from nextgem import context as cx
from nextgem import context_inputs as ci
from nextgem import context_run as cr
from nextgem import paths
from nextgem import w4_analysis as wa
from nextgem.config import SEEDS

PARENT = "human1"
VARIANT_NAME = "B-loo:cofactor_pool_biomass"
PARENT_TAG = f"{PARENT}+{VARIANT_NAME}"
MEDIUM = "M0"
W3_STRUCT = paths.results_dir() / "w3" / "structure"
VARIANTS_FILE = W3_STRUCT / "variants_human1.json"
W3_SWEEP = paths.results_dir() / "w3" / "sweeps" / "ratio_r_human1.parquet"
W3_SWEEP_ROW = f"{MEDIUM}|{VARIANT_NAME}"
W4 = paths.results_dir() / "w4"
M5 = paths.results_dir() / "m5"
OUT = M5 / "builds"
MANIFEST = M5 / "MANIFEST.tsv"
CONFIGS = (  # (family, rule, modality, variant, genomic): Stage A's carried configurations (train_selection.json)
    ("F2", "R1", "rna", "real", ("none",)),
    ("F4", "R1", "rna", "real", ("none",)),
)
JOB_CPU_LIMIT_S = 1800.0  # A10 per-job CPU kill
MAX_WORKERS = 8
MIN_FREE_GB = 10.0
AUTHOR = "M5 engineer agent"
CODE_FILES = ["src/nextgem/context.py", "src/nextgem/context_inputs.py", "src/nextgem/context_run.py",
              "src/nextgem/knockouts.py", "src/nextgem/model.py", "src/nextgem/w3.py", "scripts/m5_build.py"]


# ----------------------------------------------------------------------------- parent (worker side)
def load_variant(path: Path = VARIANTS_FILE, name: str = VARIANT_NAME):
    """The W3 variant spec exactly as swept and validated (results/w3/structure/variants_human1.json)."""
    from nextgem import w3

    spec = json.loads(Path(path).read_text())
    hits = [v for group in spec.values() for v in group if v["name"] == name]
    if len(hits) != 1:
        raise KeyError(f"{name}: {len(hits)} specs in {path}")
    return w3.Variant.from_json(hits[0])


def parent_model(solver: str = "glpk"):
    """Human1 under ``solver`` (set before the model is built: no WLS session for GLPK) with the B-loo biomass applied
    permanently (outside any ``with model:`` block)."""
    import cobra

    from nextgem.model import load_model

    cobra.Configuration().solver = solver
    m = load_model(PARENT)
    load_variant().apply(m)
    return m


def init_parent(solver: str = "glpk") -> None:
    """Worker init: the Stage-A ``context_run.init_worker`` with the Human1 + B-loo parent (fills ``context_run._W``,
    which ``context_run.worker_run`` / ``run_job`` use)."""
    from nextgem.model import GEMS

    t0, c0 = time.time(), time.process_time()
    m = parent_model(solver)
    if cx.solver_name(m) != solver:
        raise RuntimeError(f"parent loaded with {cx.solver_name(m)}, expected {solver}")
    parent = cx.make_parent(m, PARENT_TAG, medium=None, medium_name=MEDIUM, genes_tsv=GEMS[PARENT].genes_tsv)
    cr._W.update(parent=parent, solver=solver, load_wall_s=time.time() - t0, load_cpu_s=time.process_time() - c0,
                 fastcc_cache=None)


def generic_build(parent, line: str = "generic", genes=None, tasks: bool = True) -> dict:
    """The parent itself as a line model (nothing closed, no repair): KO sweep with the same engine, task profile."""
    c0 = time.process_time()
    lm = cx.LineModel(line=line, family="F0", closed=set())
    lm.info["_s_r"] = pd.Series(dtype="float64")
    cx.finish(parent, lm)
    row = {"build_id": "generic|none", "line": line, "family": "generic", "rule": "-", "modality": "-",
           "variant": "-", "genomic": "none", **cx.summary_row(lm)}
    row["family"] = "generic"
    ko = cx.ko_sweep(parent, lm, genes)
    row.update({k: lm.info[k] for k in ("sweep_cpu_s", "sweep_n_solved", "sweep_n_nonoptimal", "sweep_n_lethal")})
    if tasks:
        prof = cx.task_profile(parent, lm)
        row["tasks_passed"] = int(sum(prof.values()))
        row["task_profile"] = "".join("1" if v else "0" for v in prof.values())
    row["wt_generic"] = parent.wt_generic
    row["solver"] = cx.solver_name(parent.model)
    row["peak_rss_gb"] = cr.peak_rss_gb()
    row["job_cpu_s"] = time.process_time() - c0
    ko = ko[["gene", "ratio_r", "status", "solver"]].astype({"ratio_r": "float64", "status": "category", "solver": "category"})
    return {"rows": [row], "ko": {"none": ko}, "error": None}


def _child(conn, solver: str, job):
    """One build in a fresh spawned process (Stage-A protocol: 'started' once the parent is loaded, then the result)."""
    try:
        init_parent(solver)
        conn.send(("started", None))
        if job.family == "generic":
            res = generic_build(cr._W["parent"], genes=job.genes, tasks=job.tasks)
        else:
            res = cr.worker_run(job, None)  # asserts the declared solver (A4 c), then extract/repair/sweep/tasks
        res["parent_wt_generic"] = cr._W["parent"].wt_generic
        conn.send(("ok", res))
    except BaseException:  # noqa: BLE001 - reported to the parent, which aborts the run
        conn.send(("error", traceback.format_exc()))
    finally:
        conn.close()


class Runner(cr.Runner):
    """Stage-A runner (lanes, watchdog, failed-timeout) whose children load the Human1 + B-loo parent."""

    def start(self, key, lane, job, loss_line=None, payload=None):
        parent_conn, child_conn = self.ctx.Pipe(duplex=False)
        p = self.ctx.Process(target=_child, args=(child_conn, lane, job), daemon=True)
        p.start()
        child_conn.close()
        limit = float(job.meta.get("cpu_limit_s", self.job_timeout))
        self.running[key] = (p, parent_conn, lane, time.time(), payload, False, limit, 0.0)


# ----------------------------------------------------------------------------- enumeration and inputs (parent side)
def enumerate_jobs(lines) -> list[dict]:
    specs = [dict(build_id="generic", line="generic", family="generic", rule="-", modality="-", variant="-",
                  genomic=("none",))]
    for fam, rule, mod, var, gen in CONFIGS:
        for line in lines:
            specs.append(dict(build_id="|".join(["M5", line, fam, rule, mod, var]), line=line, family=fam, rule=rule,
                              modality=mod, variant=var, genomic=gen))
    return specs


class Inputs:
    def __init__(self, lines):
        import cobra

        from nextgem.model import GEMS

        cobra.Configuration().solver = "glpk"  # structure only: no WLS session in the parent
        m = cobra.io.read_sbml_model(str(GEMS[PARENT].sbml))
        self.genes = ci.gene_table(m, GEMS[PARENT].genes_tsv)
        del m
        st = json.loads((W4 / "threshold_stats.json").read_text())
        gm = pd.read_parquet(W4 / "threshold_gene_means.parquet").set_index("gene")
        # R1 is global (log2(TPM+1) = 1); gene means are carried only to construct the Stage-A object
        self.rna_stats = ci.ThresholdStats("rna", st["rna"]["p25"], st["rna"]["p75"], gm.rna_mean.dropna(),
                                           st["rna"]["floor_value"], st["rna"]["n_lines"])
        self.rna = ci.read_expression(self.genes, lines)

    def scores(self, spec) -> pd.Series:
        assert spec["rule"] == "R1" and spec["modality"] == "rna" and spec["variant"] == "real"
        return ci.signed_scores(spec["line"], spec["rule"], spec["modality"], self.genes, self.rna, self.rna_stats)


def input_files() -> list[Path]:
    from nextgem.model import GEMS

    dm = paths.depmap_dir()
    return [VARIANTS_FILE, W3_SWEEP, W4 / "stage_a_panel.tsv", W4 / "threshold_stats.json",
            W4 / "threshold_gene_means.parquet", paths.splits_dir() / "lineage_split.tsv", dm / ci.EXPR_FILE,
            GEMS[PARENT].sbml, GEMS[PARENT].genes_tsv, GEMS[PARENT].essential_tasks,
            W4 / "stage_a_analysis" / "train_selection.json"]


def code_hash(params: dict) -> str:
    h = hashlib.sha256()
    for f in CODE_FILES:
        h.update(f.encode())
        h.update((paths.PROJECT_ROOT / f).read_bytes())
    h.update(VARIANTS_FILE.read_bytes())
    h.update(json.dumps(params, sort_keys=True).encode())
    return h.hexdigest()[:16]


def atomic_pickle(obj, path: Path) -> None:
    tmp = path.with_suffix(".part")
    with open(tmp, "wb") as fh:
        pickle.dump(obj, fh)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


# ----------------------------------------------------------------------------- label-free summaries
def viability_idx() -> list[int]:
    from raven_toolbox.tasks.tasklist import parse_task_list

    from nextgem.model import GEMS

    return wa.viability_task_indices(parse_task_list(str(GEMS[PARENT].essential_tasks)))


def task_lp_statuses(s: str) -> dict:
    """'desc=3;desc2=optimal...' -> counts of non-integer outcomes (an integer = reactions re-opened by that LP)."""
    out = {}
    for item in (s or "").split(";"):
        if "=" not in item:
            continue
        v = item.rsplit("=", 1)[1]
        key = "solved" if v.lstrip("-").isdigit() else v
        out[key] = out.get(key, 0) + 1
    return out


def summarise(builds: pd.DataFrame, gates: pd.DataFrame, generic_check: dict, timeouts: list) -> dict:
    out = {"parent": PARENT_TAG, "medium": MEDIUM, "configs": [list(c[:4]) + [list(c[4])] for c in CONFIGS],
           "generic_check": generic_check, "timeouts": timeouts, "families": {}}
    for fam, g in builds[builds.family != "generic"].groupby("family"):
        gg = gates[gates.family == fam]
        tl = {}
        for s in g.get("viability_task_lps", pd.Series(dtype=str)).fillna(""):
            for k, v in task_lp_statuses(s).items():
                tl[k] = tl.get(k, 0) + v
        failed = g.failed.fillna(False).astype(bool)
        out["families"][fam] = {
            "n_builds": int(len(g)), "n_failed": int(failed.sum()),
            "n_failed_timeout": int((g.get("error", pd.Series(dtype=str)) == "failed-timeout").sum()),
            "repair_growth_lp": g.get("repair_growth_lp", pd.Series(dtype=str)).fillna("NA").value_counts().to_dict(),
            "repair_task_lps": tl,
            "n_reopened": {"min": float(g.n_reopened.min()), "median": float(g.n_reopened.median()),
                           "max": float(g.n_reopened.max()), "n_builds_with_reopen": int((g.n_reopened > 0).sum())},
            "repair_frac_max": float(g.repair_frac_reactions.max()),
            "wt_ratio_generic": {"min": float(g.wt_ratio_generic.min()), "median": float(g.wt_ratio_generic.median()),
                                 "max": float(g.wt_ratio_generic.max())},
            "wt_growth": {"min": float(g.wt_growth.min()), "max": float(g.wt_growth.max())},
            "sweep_n_lethal": {"min": float(g.sweep_n_lethal.min()), "median": float(g.sweep_n_lethal.median()),
                               "max": float(g.sweep_n_lethal.max())},
            "sweep_n_nonoptimal_total": int(g.sweep_n_nonoptimal.fillna(0).sum()),
            "gate_frac": {k: float(gg[k].mean()) for k in ("gate_growth", "gate_viability", "gate_repair", "gate_all")},
            "config_gate_pass": bool(gg.gate_all.mean() >= wa.CONFIG_MIN_PASS_FRAC),
            "job_cpu_s": {"mean": float(g.job_cpu_s.mean()), "max": float(g.job_cpu_s.max())},
            "peak_rss_gb_max": float(g.peak_rss_gb.max()),
        }
    return out


def generic_vs_w3(ko: pd.DataFrame) -> dict:
    """The generic H1 + B-loo sweep must equal W3's validated sweep row (same variant, M0)."""
    w3s = pd.read_parquet(W3_SWEEP).set_index("sweep").loc[W3_SWEEP_ROW]
    mine = ko.set_index("gene").ratio_r
    common = mine.index.intersection(w3s.index)
    d = (mine.reindex(common) - w3s.reindex(common).astype(float)).abs()
    return {"w3_row": W3_SWEEP_ROW, "n_genes_m5": int(mine.size), "n_genes_w3": int(w3s.size),
            "n_common": int(common.size), "n_mismatch_gt_1e-6": int((d > 1e-6).sum()),
            "max_abs_diff": float(d.max()) if len(d) else None,
            "n_lethal_m5": int((mine < 0.01).sum()), "n_lethal_w3": int((w3s.astype(float) < 0.01).sum())}


# ----------------------------------------------------------------------------- main
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=MAX_WORKERS, help=f"GLPK worker processes (<= {MAX_WORKERS})")
    ap.add_argument("--job-timeout", type=float, default=JOB_CPU_LIMIT_S, help="per-job CPU seconds (A10: 1800)")
    ap.add_argument("--lines", type=int, help="pilot: only the first N panel lines")
    ap.add_argument("--out", help="output directory (default results/m5/builds)")
    ap.add_argument("--pilot", action="store_true", help="cost pilot: allows a dirty tree and --lines; needs --out")
    ap.add_argument("--dry-run", action="store_true")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    from nextgem.run import provenance

    a = parse_args(argv)
    if a.workers > MAX_WORKERS or a.workers < 1:
        sys.exit(f"--workers must be in 1..{MAX_WORKERS}")
    if a.job_timeout != JOB_CPU_LIMIT_S and not a.pilot:
        sys.exit("the per-job CPU limit is fixed at 1,800 s (A10)")
    if a.lines and not a.pilot:
        sys.exit("--lines needs --pilot")
    if a.pilot and not a.out and not a.dry_run:
        sys.exit("--pilot needs --out (a pilot never writes results/m5/builds)")
    out = OUT if not a.out else (Path(a.out) if a.out.startswith("/") else paths.PROJECT_ROOT / a.out)
    if a.pilot and "tmp" not in out.relative_to(paths.PROJECT_ROOT).parts:
        sys.exit("a pilot writes under a git-ignored tmp/ directory")
    if Path(__file__).resolve().parents[1] != paths.PROJECT_ROOT:  # pinned-worktree guard (§11.4)
        sys.exit(f"nextgem imported from {paths.PROJECT_ROOT}, script from {Path(__file__).resolve().parents[1]}")
    panel = pd.read_csv(W4 / "stage_a_panel.tsv", sep="\t")
    if a.lines:
        panel = panel.head(a.lines)
    split = pd.read_csv(paths.splits_dir() / "lineage_split.tsv", sep="\t").set_index("ModelID")
    assert (split.loc[panel.ModelID, "tier"] == "train").all(), "the pre-check uses training panel lines only"
    sel = json.loads((W4 / "stage_a_analysis" / "train_selection.json").read_text())
    assert sel["stage_b_configs"] == ["F2|R1|rna|real|none"], "Stage A's carried F2 configuration changed"
    assert sel["families"]["F4"]["vs_H1"]["full_panel_best"] == "F4|R1|rna|real|none", "Stage A's F4 configuration"
    specs = enumerate_jobs(list(panel.ModelID))
    print(f"{len(specs)} builds (1 generic + {len(specs) - 1} line builds: {[c[0] for c in CONFIGS]} x "
          f"{len(panel)} lines); parent {PARENT_TAG} on {MEDIUM}; solver glpk; workers {a.workers}")
    if a.dry_run:
        return 0
    params = {"lines": list(panel.ModelID), "configs": [list(map(str, c[:4])) for c in CONFIGS],
              "job_timeout": a.job_timeout, "variant": VARIANT_NAME}
    chash = code_hash(params)
    jobs_dir = out / "tmp" / f"jobs_{chash}" if not a.pilot else out / f"jobs_{chash}"
    jobs_dir.mkdir(parents=True, exist_ok=True)
    manifest = MANIFEST if not a.pilot else out / "MANIFEST.tsv"
    with provenance(manifest, inputs=input_files(), seeds=[SEEDS[0]], solver="glpk", workers=a.workers,
                    allow_dirty=a.pilot, author=AUTHOR) as prov:
        prov.note(step="m5_build (A12 b, label-free)", code_hash=chash, pilot=a.pilot, parent=PARENT_TAG,
                  medium=MEDIUM, job_cpu_limit_s=a.job_timeout, lp_time_limit_s=cx.LP_TIME_LIMIT_S,
                  n_specs=len(specs))
        return _run(a, out, panel, specs, jobs_dir, prov)


def _run(a, out, panel, specs, jobs_dir, prov) -> int:
    t0 = time.time()
    I = Inputs(list(panel.ModelID))
    runner = Runner({"glpk": a.workers}, a.job_timeout)
    state = {"n_done": 0, "timeouts": [], "cpu": 0.0}

    def guard():
        free = psutil.virtual_memory().available / 1e9
        if free < MIN_FREE_GB:
            raise SystemExit(f"abort: available memory {free:.1f} GB < {MIN_FREE_GB} GB")

    def job_path(bid):
        return jobs_dir / (hashlib.sha1(bid.encode()).hexdigest()[:16] + ".pkl")

    todo = [s for s in specs if not job_path(s["build_id"]).exists()]
    print(f"{len(specs) - len(todo)} builds checkpointed; {len(todo)} to run", flush=True)
    items = []
    for s in todo:
        scores = pd.Series(dtype="float64") if s["family"] == "generic" else I.scores(s)
        job = cr.Job(build_id=s["build_id"], line=s["line"], family=s["family"], rule=s["rule"], modality=s["modality"],
                     variant=s["variant"], scores=scores, genomic=s["genomic"], accepted=(), tasks=True,
                     meta={"row": "M5", "solver": "glpk", "cpu_limit_s": a.job_timeout})
        items.append((s["build_id"], "glpk", job, None, s))

    def on_build(key, status, res, spec, wall):
        if status == "error":
            raise SystemExit(f"worker error in {key} (run aborted, §11.4):\n{res}")
        if status == "timeout":  # A10: failed-timeout, recorded, never re-run with relaxed limits
            res = {"rows": [{"build_id": f"{key}|{g}", "line": spec["line"], "family": spec["family"],
                             "rule": spec["rule"], "modality": spec["modality"], "variant": spec["variant"],
                             "genomic": g, "row": "M5", "failed": True, "error": "failed-timeout",
                             "job_wall_s": wall, "solver": "glpk"} for g in spec["genomic"]],
                   "ko": {}, "error": "failed-timeout"}
            state["timeouts"].append(key)
        else:
            state["cpu"] += sum(r.get("job_cpu_s", 0) or 0 for r in res["rows"])
        atomic_pickle(res, job_path(key))
        state["n_done"] += 1
        print(f"[{time.time() - t0:6.0f}s] {state['n_done']}/{len(todo)} {key} {status} wall {wall:.0f}s", flush=True)

    runner.run_all(items, on_build, guard)

    rows, kos, parent_wt = [], [], set()
    for s in specs:
        res = pickle.loads(job_path(s["build_id"]).read_bytes())
        if "parent_wt_generic" in res:
            parent_wt.add(round(res["parent_wt_generic"], 9))
        rows += res["rows"]
        for gvar, ko in res["ko"].items():
            ko = ko.copy()
            ko.insert(0, "build_id", f"{s['build_id']}|{gvar}" if s["family"] != "generic" else "generic|none")
            kos.append(ko)
    assert len(parent_wt) == 1, f"parent WT differs between processes: {parent_wt}"
    builds = pd.DataFrame(rows)
    builds.loc[builds.family == "generic", "build_id"] = "generic|none"
    lin = panel.set_index("ModelID").lineage
    builds.insert(2, "lineage", builds.line.map(lin).fillna("-"))
    gates = wa.build_gates(builds, viability_idx())
    gcols = ["build_id", "line", "lineage", "family", "rule", "genomic", "failed", "gate_growth", "gate_viability",
             "gate_repair", "gate_all", "viability_passed", "tasks_passed", "wt_growth", "wt_ratio_generic",
             "n_reopened", "repair_frac_reactions", "repair_growth_lp", "sweep_n_lethal", "job_cpu_s"]
    gates = gates[[c for c in gcols if c in gates.columns]]
    ko = pd.concat(kos, ignore_index=True)
    ko["build_id"] = ko.build_id.astype("category")
    gen_check = generic_vs_w3(ko[ko.build_id == "generic|none"])
    summary = summarise(builds, gates, gen_check, state["timeouts"])
    summary["parent_wt_generic"] = parent_wt.pop()
    out.mkdir(parents=True, exist_ok=True)
    builds.to_csv(out / "builds.tsv", sep="\t", index=False)
    gates.to_csv(out / "gates.tsv", sep="\t", index=False)
    ko.to_parquet(out / "ko.parquet", index=False)
    cpu_total = float(builds.job_cpu_s.fillna(0).sum())
    summary.update({"build_cpu_h": round(cpu_total / 3600, 4), "wall_s": round(time.time() - t0),
                    "n_line_builds": int((builds.family != "generic").sum())})
    (out / "build_summary.json").write_text(json.dumps(summary, indent=1, default=str))
    prov.output(out / "builds.tsv", out / "gates.tsv", out / "ko.parquet", out / "build_summary.json")
    prov.note(build_cpu_h=summary["build_cpu_h"], timeouts=state["timeouts"], wall_s=summary["wall_s"])
    print(json.dumps({k: summary[k] for k in ("parent_wt_generic", "generic_check", "timeouts", "build_cpu_h",
                                              "wall_s")}, indent=1, default=str))
    for fam, f in summary["families"].items():
        print(fam, json.dumps(f, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
