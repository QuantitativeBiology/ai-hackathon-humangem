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
import traceback
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
    genes: tuple | None = None  # KO sweep subset (tests); None = all parent genes
    meta: dict = field(default_factory=dict)


def init_worker(parent_name: str = "human2", solver: str = "gurobi", fastcc_cache=None, medium=None,
                medium_name: str = "M0") -> None:
    """Load the parent directly under ``solver``: cobra's default solver is set BEFORE the model is built,
    so a GLPK worker never creates a gurobipy model / WLS environment (review B1). A Gurobi worker must
    already hold a slot of ``nextgem.gurobi_pool``."""
    import cobra

    from nextgem.model import load_model

    t0, c0 = time.time(), time.process_time()
    cobra.Configuration().solver = solver
    m = load_model(parent_name)
    if cx.solver_name(m) != solver:
        raise RuntimeError(f"parent loaded with {cx.solver_name(m)}, expected {solver}")
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
                ko = cx.ko_sweep(parent, lm, job.genes)
                row.update({k: lm.info[k] for k in ("sweep_cpu_s", "sweep_n_solved", "sweep_n_nonoptimal", "sweep_n_lethal")})
                if job.keep_ko:
                    out["ko"][gvar] = ko[["gene", "ratio_r", "status", "solver"]].astype(
                        {"ratio_r": "float64", "status": "category", "solver": "category"})  # §11.4: solver per KO row
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


# ============================================================================= one-job-per-process runner
def _cpu(pid) -> float:
    try:
        import psutil

        t = psutil.Process(pid).cpu_times()
        return t.user + t.system
    except Exception:  # noqa: BLE001 - process already gone
        return 0.0


def _child(conn, solver: str, job: Job, loss_line, fastcc_cache):
    """Runs in a fresh spawned process: set up the parent under ``solver``, run ONE job, send the result."""
    try:
        if solver == "gurobi":
            from nextgem import gurobi_pool

            _W["slot"] = gurobi_pool.hold_slot()  # held until this process exits
        init_worker("human2", solver, fastcc_cache)
        conn.send(("started", None))  # the watchdog clock starts now (slot wait and model load excluded)
        conn.send(("ok", worker_run(job, loss_line)))
    except BaseException:  # noqa: BLE001 - reported to the parent, which aborts the run
        conn.send(("error", traceback.format_exc()))
    finally:
        conn.close()


def worker_run(job: Job, loss_line: str | None):
    want = job.meta.get("solver") or (cx.LOSS_RULE_SOLVER if loss_line is not None else cx.DECLARED_SOLVER[job.family])
    have = cx.solver_name(_W["parent"].model)
    if have != want:  # always asserted, also under --pilot (review NOTE 6a)
        raise RuntimeError(f"{job.build_id}: parent solver {have} != declared {want} (amendment A4 c)")
    if loss_line is not None:  # phase 0: loss rule for one line
        lost, expr = job.meta.pop("_lost"), job.meta.pop("_expr")
        c0 = time.process_time()
        res = cx.accept_losses(_W["parent"], lost, expr)
        res["cpu_s"] = time.process_time() - c0
        res["peak_rss_gb"] = peak_rss_gb()
        return {"loss": res}
    return run_job(job)


class Runner:
    """Per-solver lanes of one-job processes (spawned, one job each), with a watchdog per job.

    A job is killed and reported as ``timeout`` when, after its solver is ready, its process has used more than
    its CPU limit (``job.meta["cpu_limit_s"]``, default ``job_timeout``) or more than twice that in wall time.
    """

    def __init__(self, workers: dict, job_timeout: float, fastcc_cache=None):
        import multiprocessing as mp

        self.ctx = mp.get_context("spawn")
        self.workers, self.job_timeout, self.fastcc = workers, job_timeout, fastcc_cache
        self.running = {}  # key -> (proc, conn, lane, t0, payload, started, cpu_limit_s, cpu_at_start)

    def n_running(self, lane):
        return sum(1 for v in self.running.values() if v[2] == lane)

    def start(self, key, lane, job, loss_line=None, payload=None):
        parent_conn, child_conn = self.ctx.Pipe(duplex=False)
        p = self.ctx.Process(target=_child, args=(child_conn, lane, job, loss_line, self.fastcc), daemon=True)
        p.start()
        child_conn.close()
        limit = float(job.meta.get("cpu_limit_s", self.job_timeout))
        self.running[key] = (p, parent_conn, lane, time.time(), payload, False, limit, 0.0)

    def poll(self, timeout=1.0):
        """Yield (key, status, result, payload) for finished or timed-out jobs."""
        from multiprocessing.connection import wait as mpwait

        conns = {v[1]: k for k, v in self.running.items()}
        ready = mpwait(list(conns), timeout=timeout) if conns else []
        done = []
        for c in ready:
            k = conns[c]
            p, conn, lane, t0, payload, started, limit, cpu0 = self.running[k]
            try:
                status, res = conn.recv()
            except EOFError:
                status, res = "error", f"worker died (exit code {p.exitcode})"
            if status == "started":
                self.running[k] = (p, conn, lane, time.time(), payload, True, limit, _cpu(p.pid))
                continue
            self.running.pop(k)
            p.join(30)
            done.append((k, status, res, payload, time.time() - t0))
        now = time.time()
        for k, (p, conn, lane, t0, payload, started, limit, cpu0) in list(self.running.items()):
            if started and (_cpu(p.pid) - cpu0 > limit or now - t0 > 2 * limit):
                p.kill()
                p.join(30)
                self.running.pop(k)
                done.append((k, "timeout", None, payload, now - t0))
        return done

    def kill_all(self):
        for p, *_ in self.running.values():
            p.kill()
        self.running.clear()

    def run_all(self, items, on_done, guard=None):
        """items: list of (key, lane, job, loss_line, payload); calls on_done(key, status, res, payload, wall)."""
        queues = {lane: [it for it in items if it[1] == lane] for lane in self.workers}
        try:
            while any(queues.values()) or self.running:
                for lane, q in queues.items():
                    while q and self.n_running(lane) < self.workers[lane]:
                        key, ln, job, loss_line, payload = q.pop(0)
                        self.start(key, ln, job, loss_line, payload)
                for args in self.poll():
                    on_done(*args)
                if guard:
                    guard()
        except BaseException:
            self.kill_all()
            raise
