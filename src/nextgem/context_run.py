"""Worker-side runner for W4 line builds (pilots and the Stage-A sweep).

One worker process = one loaded parent (one solver environment; one Gurobi WLS session when the
solver is Gurobi). Jobs carry the line's signed gene scores (computed in the parent process from
label-free inputs) and, for genomic variants, the accepted-loss set; the worker extracts, applies
the genomic constraint post-extraction (A4: gene KO + repair + sweep, never a rebuild), repairs,
sweeps all genes and records the 57-task capability profile.
"""

from __future__ import annotations

import os
import resource
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from nextgem import context as cx

_W: dict = {}


@dataclass
class Job:
    build_id: str
    line: str
    family: str
    rule: str
    modality: str
    variant: str  # "real" | "perm" | "noise<k>" | "gimme_seed<k>"
    scores: pd.Series
    genomic: tuple = ("none",)  # ("none",), ("none", "loss") or ("loss",)
    accepted: tuple = ()
    gimme_noise_seed: int | None = None
    tasks: bool = True
    keep_ko: bool = True
    meta: dict = field(default_factory=dict)


def init_worker(parent_name: str = "human2", solver: str = "gurobi", fastcc_cache=None, medium=None,
                medium_name: str = "M0") -> None:
    from nextgem.model import load_model

    t0, c0 = time.time(), time.process_time()
    m = load_model(parent_name)
    if solver != "gurobi":
        m.solver = solver
    parent = cx.make_parent(m, parent_name, medium=medium, medium_name=medium_name)
    if fastcc_cache is not None and os.path.exists(fastcc_cache):
        cx.consistent_reactions(parent, fastcc_cache)
    _W.update(parent=parent, solver=solver, load_wall_s=time.time() - t0, load_cpu_s=time.process_time() - c0,
              fastcc_cache=fastcc_cache)


def peak_rss_gb() -> float:
    """Peak RSS of this process's address space (VmHWM, GB). Not ru_maxrss: on Linux that survives
    fork + exec, so spawned workers would report their parent's peak."""
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) / 1e6
    except OSError:
        pass
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def open_mask(parent, lm) -> str:
    closed = lm.final_closed
    bits = np.array([rid not in closed for rid in parent.rxn_ids], dtype=bool)
    return np.packbits(bits).tobytes().hex()


def run_job(job: Job) -> dict:
    """Returns {"rows": [summary dicts], "ko": {variant_key: DataFrame}, "error": str|None}."""
    parent = _W["parent"]
    if job.family == "F3" and parent.consistent is None:
        cx.consistent_reactions(parent, _W.get("fastcc_cache"))
    out = {"rows": [], "ko": {}, "error": None}
    c0, t0 = time.process_time(), time.time()
    try:
        s = job.scores.reindex(parent.genes)
        base = cx.extract(parent, job.family, s, gimme_noise_seed=job.gimme_noise_seed)
    except cx.BuildFailed as e:
        out["error"] = f"extract: {e}"
        out["rows"].append({**_ids(job, "none"), "failed": True, "error": str(e), "solver": _W["solver"],
                            "peak_rss_gb": peak_rss_gb(), "job_cpu_s": time.process_time() - c0})
        return out
    extract_cpu = base.info["extract_cpu_s"]
    for gvar in job.genomic:
        lm = cx.LineModel(line=job.line, family=job.family, closed=set(base.closed), eflux=dict(base.eflux),
                          off=set(base.off), lin=dict(base.lin), cap_max=base.cap_max, info=dict(base.info))
        c1 = time.process_time()
        try:
            cx.finish(parent, lm, set(job.accepted) if gvar == "loss" else None)
            row = {**_ids(job, gvar), **cx.summary_row(lm)}
            if lm.info["failed"]:
                row["error"] = "no growth after repair"
            else:
                ko = cx.ko_sweep(parent, lm)
                row.update({k: lm.info[k] for k in ("sweep_cpu_s", "sweep_n_solved", "sweep_n_nonoptimal", "sweep_n_lethal")})
                if job.keep_ko:
                    out["ko"][gvar] = ko[["gene", "ratio_r", "status", "solver"]].astype(
                        {"ratio_r": "float32", "status": "category", "solver": "category"})  # §11.4: solver per KO row
            if job.tasks:
                tc = time.process_time()
                prof = cx.task_profile(parent, lm)
                row["tasks_passed"] = int(sum(prof.values()))
                row["task_profile"] = "".join("1" if v else "0" for v in prof.values())
                row["tasks_cpu_s"] = time.process_time() - tc
            row["open_mask"] = open_mask(parent, lm)
            row["viability_task_lps"] = ";".join(f"{k}={v}" for k, v in lm.info.get("repair_task_lps", {}).items())
        except cx.BuildFailed as e:
            row = {**_ids(job, gvar), "failed": True, "error": str(e)}
        row["extract_cpu_s"] = extract_cpu if gvar == job.genomic[0] else 0.0  # extraction counted once
        row["variant_cpu_s"] = time.process_time() - c1
        row["solver"] = _W["solver"]
        row["peak_rss_gb"] = peak_rss_gb()
        row["pid"] = os.getpid()
        out["rows"].append(row)
    out["rows"][-1]["job_cpu_s"] = time.process_time() - c0
    out["rows"][-1]["job_wall_s"] = time.time() - t0
    return out


def _ids(job: Job, gvar: str) -> dict:
    return {"build_id": f"{job.build_id}|{gvar}", "line": job.line, "family": job.family, "rule": job.rule,
            "modality": job.modality, "variant": job.variant, "genomic": gvar, **job.meta}
