"""F4b step 1 (label-free): E-Flux on the P1 parent, 30-panel pre-check builds (A15 post hoc; RESEARCH_TREE D19).

Mirrors scripts/m5_build.py (A12(b)) with two changes fixed by D19: the parent is P1's, i.e. Human1 + W3
B-loo:cofactor_pool_biomass with the LINE's documented medium M2(f=1) (bounds of the line's W2 signature, rebuilt as M6
did; line -> signature from results/m6/scores/line_assignment.tsv), instead of M0; and the reference is P1.

Builds on the 30 Stage-A panel lines (results/w4/stage_a_panel.tsv; training tier):
  F4b      E-Flux, Stage A's carried configuration F4|R1|rna|real|none (global TPM-1 threshold, real RNA, no genomic
           constraint), the W4 builders (``context.extract`` -> ``context.finish`` (E-Flux: no repair) -> WT growth ->
           ``context.ko_sweep`` with the A4 KO convention -> 57-task profile), one fresh spawned process per build.
  generic  the P1 parent itself on every distinct panel signature (same W4 engine; checked bitwise-or-<=1e-6 against
           the pinned P1 sweep, as M5 checked its generic build against W3).
GLPK (declared solver for F4, A4 c); A10: 300 s per LP + one retry (``context.LP_TIME_LIMIT_S``), CPU kill at 1,800 s
per job (failed-timeout, never re-run with relaxed limits); <= 12 workers.

Label-blind gates (``w4_analysis.build_gates``): growth, viability (GR + 9 ER tasks), repair <= 5% (n/a for E-Flux).
Known units artefact (results/w4, results/m5): E-Flux WT growth is ~1-2% of the parent, so the GR task (biomass >= 1)
fails whenever growth < 1; the gate statistics are reported, as in M5.

No label is read. Outputs (default results/f4/eflux_panel/builds/): builds.tsv, gates.tsv, ko.parquet, build_summary.json.

Usage:
  .venv/bin/python scripts/f4_eflux_build.py --dry-run
  .venv/bin/python scripts/f4_eflux_build.py --pilot --lines 2 --out results/f4/tmp/eflux_pilot --workers 4
  .venv/bin/python scripts/f4_eflux_build.py --workers 12            # pinned worktree (> 1 CPU-h)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import psutil

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f4_common as fc  # noqa: E402

fc.blas_single_thread()

from nextgem import context as cx  # noqa: E402
from nextgem import context_run as cr  # noqa: E402
from nextgem import w4_analysis as wa  # noqa: E402

W4 = fc.PROJECT / "results" / "w4"
OUT = fc.EFDIR / "builds"
MIN_FREE_GB = 10.0


def _child(conn, job, sig: str, bounds: dict):
    """One build in a fresh spawned process: the P1 parent on ``sig``, then the build."""
    import m5_build

    try:
        t0, c0 = time.time(), time.process_time()
        parent = fc.p1_parent(bounds, sig)
        cr._W.update(parent=parent, solver="glpk", load_wall_s=time.time() - t0, load_cpu_s=time.process_time() - c0,
                     fastcc_cache=None)
        conn.send(("started", None))
        if job.family == "generic":
            res = m5_build.generic_build(parent, line=job.line, genes=job.genes, tasks=job.tasks)
            res["rows"][0]["build_id"] = job.build_id
            res["rows"][0]["w2_sweep_id"] = sig
        else:
            res = cr.worker_run(job, None)  # asserts the declared solver (A4 c); extract / finish / sweep / tasks
        res["parent_wt_generic"] = parent.wt_generic
        res["load_cpu_s"] = cr._W["load_cpu_s"]
        conn.send(("ok", res))
    except BaseException:  # noqa: BLE001 - reported to the parent, which aborts the run
        conn.send(("error", traceback.format_exc()))
    finally:
        conn.close()


class Runner(cr.Runner):
    def __init__(self, workers, job_timeout, bounds: dict):
        super().__init__(workers, job_timeout)
        self.bounds = bounds

    def start(self, key, lane, job, loss_line=None, payload=None):
        sig = job.meta["w2_sweep_id"]
        parent_conn, child_conn = self.ctx.Pipe(duplex=False)
        p = self.ctx.Process(target=_child, args=(child_conn, job, sig, self.bounds[sig]), daemon=True)
        p.start()
        child_conn.close()
        limit = float(job.meta.get("cpu_limit_s", self.job_timeout))
        self.running[key] = (p, parent_conn, lane, time.time(), payload, False, limit, 0.0)


def build_id(line: str) -> str:
    return "|".join(["F4b", line, *fc.F4B_CONFIG[:4]])


def enumerate_jobs(panel: pd.DataFrame, sig: pd.Series) -> list[dict]:
    specs = [dict(build_id=f"generic|{s}", line=f"generic|{s}", family="generic", w2_sweep_id=s)
             for s in sorted(set(sig.reindex(panel.ModelID)))]
    for line in panel.ModelID:
        specs.append(dict(build_id=build_id(line), line=line, family="F4", w2_sweep_id=sig[line]))
    return specs


def generic_check(ko: pd.DataFrame, builds: pd.DataFrame) -> pd.DataFrame:
    """Every W4-engine generic build equals the pinned P1 sweep of its signature (max |Δratio_r|; bitwise flag)."""
    from f4_sweeps import p1_rows

    g = builds[builds.family == "generic"]
    P1 = p1_rows(sorted(g.w2_sweep_id))
    rows = []
    for r in g.itertuples():
        mine = ko[ko.build_id.astype(str) == r.build_id].set_index("gene").ratio_r.reindex(P1.columns).astype(float)
        ref = P1.loc[r.w2_sweep_id]
        d = np.abs(mine.to_numpy() - ref.to_numpy())
        rows.append({"w2_sweep_id": r.w2_sweep_id, "bitwise_equal": bool(np.array_equal(mine.to_numpy(), ref.to_numpy())),
                     "max_abs_diff": float(np.nanmax(d)), "n_diff_gt_1e-6": int((d > 1e-6).sum()),
                     "n_lethal_w4_engine": int((mine < 0.01).sum()), "n_lethal_p1": int((ref < 0.01).sum()),
                     "pass": bool(np.nanmax(d) <= 1e-6)})
    return pd.DataFrame(rows)


def summarise(builds: pd.DataFrame, gates: pd.DataFrame, gchk: pd.DataFrame, ko: pd.DataFrame, timeouts) -> dict:
    f = builds[builds.family == "F4"]
    gg = gates[gates.family == "F4"]
    failed = f.failed.fillna(False).astype(bool)
    # lethal set vs the line's P1 (W4-engine generic build of its signature)
    lethal = {b: set(g.gene[g.ratio_r < 0.01]) for b, g in ko.groupby("build_id", observed=True)}
    same = [lethal.get(r.build_id) is not None and lethal.get(r.build_id) == lethal.get(f"generic|{r.w2_sweep_id}|none")
            for r in f.itertuples()]
    gr, er = gr_er_idx()
    prof = f.task_profile.fillna("")
    gr_fail = prof.map(lambda s: len(s) > gr and s[gr] == "0")
    er_fail = prof.map(lambda s: sum(1 for i in er if len(s) > i and s[i] == "0"))
    return {
        "parent": f"{fc.PARENT}+{fc.BLOO} on the line's {fc.M2}", "config": "|".join(fc.F4B_CONFIG),
        "n_builds": int(len(f)), "n_failed": int(failed.sum()), "timeouts": list(timeouts),
        "wt_ratio_parent": {"min": float(f.wt_ratio_generic.min()), "median": float(f.wt_ratio_generic.median()),
                            "max": float(f.wt_ratio_generic.max())},
        "wt_growth": {"min": float(f.wt_growth.min()), "median": float(f.wt_growth.median()), "max": float(f.wt_growth.max())},
        "n_wt_growth_lt_1": int((f.wt_growth < 1).sum()),
        "n_gr_task_fail": int(gr_fail.sum()), "n_builds_with_er_fail": int((er_fail > 0).sum()),
        "gate_frac": {k: float(gg[k].mean()) for k in ("gate_growth", "gate_viability", "gate_repair", "gate_all")},
        "n_gate_all_pass": int(gg.gate_all.sum()),
        "config_gate_pass": bool(gg.gate_all.mean() >= wa.CONFIG_MIN_PASS_FRAC),
        "sweep_n_lethal": {"min": float(f.sweep_n_lethal.min()), "median": float(f.sweep_n_lethal.median()),
                           "max": float(f.sweep_n_lethal.max())},
        "n_lethal_set_equal_to_P1": int(sum(bool(x) for x in same)),
        "sweep_n_nonoptimal_total": int(f.sweep_n_nonoptimal.fillna(0).sum()),
        "tasks_passed": {"min": float(f.tasks_passed.min()), "median": float(f.tasks_passed.median()),
                         "max": float(f.tasks_passed.max())},
        "generic_check": {"n": int(len(gchk)), "n_pass": int(gchk["pass"].sum()),
                          "n_bitwise": int(gchk.bitwise_equal.sum()), "max_abs_diff": float(gchk.max_abs_diff.max())},
        "job_cpu_s": {"mean": float(f.job_cpu_s.mean()), "max": float(f.job_cpu_s.max())},
        "peak_rss_gb_max": float(builds.peak_rss_gb.max()),
    }


def m5_viability_idx() -> list[int]:
    import m5_build

    return m5_build.viability_idx()


def gr_er_idx() -> tuple[int, list[int]]:
    """Profile positions of the GR task and of the 9 ER tasks (Human1 essential task list)."""
    from raven_toolbox.tasks.tasklist import parse_task_list

    from nextgem.model import GEMS

    t = parse_task_list(str(GEMS[fc.PARENT].essential_tasks))
    vi = m5_viability_idx()
    gr = [i for i in vi if t[i].id == "GR"]
    assert len(gr) == 1
    return gr[0], [i for i in vi if t[i].id == "ER"]


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=fc.MAX_WORKERS)
    ap.add_argument("--lines", type=int, help="pilot: only the first N panel lines")
    ap.add_argument("--out", help="output directory (default results/f4/eflux_panel/builds)")
    ap.add_argument("--pilot", action="store_true", help="cost pilot: dirty tree and --lines allowed; needs --out under tmp/")
    ap.add_argument("--dry-run", action="store_true")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    from nextgem import paths
    from nextgem.config import SEEDS

    a = parse_args(argv)
    if not 1 <= a.workers <= fc.MAX_WORKERS:
        sys.exit(f"--workers must be in 1..{fc.MAX_WORKERS}")
    if a.lines and not a.pilot:
        sys.exit("--lines needs --pilot")
    out = OUT if not a.out else (Path(a.out) if a.out.startswith("/") else fc.PROJECT / a.out)
    if a.pilot and not a.dry_run and not fc.in_tmp(out):
        sys.exit("a pilot writes under a git-ignored tmp/ directory")
    if Path(__file__).resolve().parents[1] != paths.PROJECT_ROOT:  # pinned-worktree guard (§11.4)
        sys.exit(f"nextgem imported from {paths.PROJECT_ROOT}, script from {Path(__file__).resolve().parents[1]}")
    panel = pd.read_csv(W4 / "stage_a_panel.tsv", sep="\t")
    if a.lines:
        panel = panel.head(a.lines)
    split = pd.read_csv(paths.splits_dir() / "lineage_split.tsv", sep="\t").set_index("ModelID")
    assert (split.loc[panel.ModelID, "tier"] == "train").all(), "the pre-check uses training panel lines only"
    sel = json.loads((W4 / "stage_a_analysis" / "train_selection.json").read_text())
    assert sel["families"]["F4"]["vs_H1"]["full_panel_best"] == "|".join(fc.F4B_CONFIG), "Stage A's F4 configuration"
    sig = fc.line_signatures()
    specs = enumerate_jobs(panel, sig)
    n_gen = sum(s["family"] == "generic" for s in specs)
    print(f"{len(specs)} builds ({n_gen} generic P1 signatures + {len(specs) - n_gen} F4b E-Flux lines); parent "
          f"{fc.PARENT}+{fc.BLOO} on each line's {fc.M2}; GLPK; workers {a.workers}")
    if a.dry_run:
        return 0
    from nextgem import context_inputs as ci
    from nextgem.model import GEMS

    inputs = [W4 / "stage_a_panel.tsv", W4 / "threshold_stats.json", W4 / "threshold_gene_means.parquet",
              W4 / "stage_a_analysis" / "train_selection.json", paths.splits_dir() / "lineage_split.tsv",
              paths.depmap_dir() / ci.EXPR_FILE, GEMS[fc.PARENT].sbml, GEMS[fc.PARENT].genes_tsv,
              GEMS[fc.PARENT].essential_tasks, fc.PROJECT / "results" / "w3" / "structure" / "variants_human1.json",
              fc.PROJECT / "results" / "m6" / "scores" / "line_assignment.tsv",
              fc.PROJECT / "results" / "w2" / "controls" / "line_sweeps.tsv",
              fc.PROJECT / "results" / "m6" / "sweeps" / "ratio_r_human1.parquet",
              fc.PROJECT / "results" / "w3" / "sweeps" / "ratio_r_human1.parquet"]
    manifest = fc.MANIFEST if not a.pilot else out / "MANIFEST.tsv"
    with fc.provenance(inputs=inputs, seeds=[SEEDS[0]], solver="glpk", workers=a.workers, allow_dirty=a.pilot,
                       manifest=manifest) as prov:
        prov.note(step="f4_eflux_build (F4b panel pre-check builds, label-free, D19)", pilot=a.pilot,
                  job_cpu_limit_s=fc.JOB_CPU_LIMIT_S, lp_time_limit_s=cx.LP_TIME_LIMIT_S, n_specs=len(specs))
        return _run(a, out, panel, specs, prov)


def _run(a, out, panel, specs, prov) -> int:
    import m5_build

    t0 = time.time()
    I = m5_build.Inputs(list(panel.ModelID))
    bounds, info = fc.signature_bounds({s["w2_sweep_id"] for s in specs})
    runner = Runner({"glpk": a.workers}, float(fc.JOB_CPU_LIMIT_S), bounds)
    jobs_dir = (out / "tmp" if not a.pilot else out) / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)
    state = {"n_done": 0, "timeouts": []}

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
        if s["family"] == "generic":
            job = cr.Job(build_id=s["build_id"], line=s["line"], family="generic", rule="-", modality="-", variant="-",
                         scores=pd.Series(dtype="float64"), tasks=True,
                         meta={"row": "F4b", "solver": "glpk", "cpu_limit_s": fc.JOB_CPU_LIMIT_S,
                               "w2_sweep_id": s["w2_sweep_id"]})
        else:
            fam, rule, mod, var, gen = fc.F4B_CONFIG
            sc = I.scores(dict(line=s["line"], rule=rule, modality=mod, variant=var))
            job = cr.Job(build_id=s["build_id"], line=s["line"], family=fam, rule=rule, modality=mod, variant=var,
                         scores=sc, genomic=(gen,), accepted=(), tasks=True,
                         meta={"row": "F4b", "solver": "glpk", "cpu_limit_s": fc.JOB_CPU_LIMIT_S,
                               "w2_sweep_id": s["w2_sweep_id"]})
        items.append((s["build_id"], "glpk", job, None, s))

    def on_build(key, status, res, spec, wall):
        if status == "error":
            raise SystemExit(f"worker error in {key} (run aborted, §11.4):\n{res}")
        if status == "timeout":  # A10
            res = {"rows": [{"build_id": f"{key}|none", "line": spec["line"], "family": spec["family"],
                             "w2_sweep_id": spec["w2_sweep_id"], "failed": True, "error": "failed-timeout",
                             "job_wall_s": wall, "solver": "glpk"}], "ko": {}, "error": "failed-timeout"}
            state["timeouts"].append(key)
        m5_build.atomic_pickle(res, job_path(key))
        state["n_done"] += 1
        print(f"[{time.time() - t0:6.0f}s] {state['n_done']}/{len(todo)} {key} {status} wall {wall:.0f}s", flush=True)

    runner.run_all(items, on_build, guard)

    rows, kos, load_cpu = [], [], 0.0
    for s in specs:
        res = pickle.loads(job_path(s["build_id"]).read_bytes())
        load_cpu += float(res.get("load_cpu_s", 0.0))
        for r in res["rows"]:
            r.setdefault("w2_sweep_id", s["w2_sweep_id"])
            if s["family"] == "generic":
                r["build_id"] = f"{s['build_id']}|none"
                r["line"] = s["line"]
            rows.append(r)
        for gvar, ko in res["ko"].items():
            ko = ko.copy()
            ko.insert(0, "build_id", f"{s['build_id']}|{gvar}")
            kos.append(ko)
    builds = pd.DataFrame(rows)
    lin = panel.set_index("ModelID").lineage
    builds.insert(2, "lineage", builds.line.map(lin).fillna("-"))
    gates = wa.build_gates(builds, m5_viability_idx())
    gcols = ["build_id", "line", "lineage", "family", "w2_sweep_id", "failed", "gate_growth", "gate_viability",
             "gate_repair", "gate_all", "viability_passed", "tasks_passed", "wt_growth", "wt_ratio_generic",
             "n_reopened", "repair_frac_reactions", "sweep_n_lethal", "job_cpu_s"]
    gates = gates[[c for c in gcols if c in gates.columns]]
    ko = pd.concat(kos, ignore_index=True)
    gchk = generic_check(ko, builds)
    summary = summarise(builds, gates, gchk, ko, state["timeouts"])
    cpu_total = float(builds.job_cpu_s.fillna(0).sum()) + load_cpu
    summary.update({"build_cpu_h": round(cpu_total / 3600, 4), "wall_s": round(time.time() - t0),
                    "medium_keys": {s: info[s]["medium_keys"] for s in sorted(gchk.w2_sweep_id)}})
    out.mkdir(parents=True, exist_ok=True)
    fc.atomic_write_tsv(builds, out / "builds.tsv")
    fc.atomic_write_tsv(gates, out / "gates.tsv")
    fc.atomic_write_tsv(gchk, out / "generic_check.tsv")
    ko["build_id"] = ko.build_id.astype("category")
    fc.atomic_write_parquet(ko, out / "ko.parquet")
    fc.atomic_write_json(summary, out / "build_summary.json")
    prov.output(out / "builds.tsv", out / "gates.tsv", out / "generic_check.tsv", out / "ko.parquet",
                out / "build_summary.json")
    prov.note(build_cpu_h=summary["build_cpu_h"], timeouts=state["timeouts"], wall_s=summary["wall_s"])
    print(json.dumps(summary, indent=1, default=str))
    print(gchk.to_string(index=False))
    return 0 if bool(gchk["pass"].all()) else 3


if __name__ == "__main__":
    sys.exit(main())
