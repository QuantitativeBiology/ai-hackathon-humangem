"""W4 Stage A: the frozen LP-only build table on the 30-line training panel (DESIGN.md §9 W4).

Rows (label-free; Human2 parent; medium M0 until W2 fixes per-line media):
  A1  F0, F2, F3, F4 x {R1, R3} x {real, permuted} RNA                         panel (30)      480 builds
  A2  F0, F2, F3, F4 x R2 x real RNA                                              panel (30)      120
  A3  F0, F2, F3, F4 x {R1, R2, R3} x {protein, rna_protein}                      panel w/ ProCan 24 x n_p
  A4  + genomic constraint post-extraction (accepted losses -> KO + repair + sweep) on every real build of A1-A3
  A5  F5 LCP, fixed label-free 24-point grid (amendment A4 d): w_r = 1, w_p in {0, 0.5, 1}, w_c in {0, 1},
      tau in {-1.5, -1, -0.5, 0}; gene off iff w_r z_RNA + w_p z_prot + w_c loss < tau (z against training
      lines; loss = -1 for accepted losses; z_prot = 0 without proteomics). Grid points that give identical
      scores on a line are built once; a5_grid_map.tsv maps all 720 (line, point) pairs to their build.
  A6  F3b CORDA, scope by the M3 pilot band (--corda-band): full = {R1, R3} x {real, perm} on 30 lines;
      fallback = R1 x {real, perm} on 15 lines; drop = none (the drop is reported)
  A7  stability: F2 GIMME (objective perturbation, seeds SEEDS[0:3]) and F3 FASTCORE (expression noise
      N(0, --noise-sd) on log2 values, seeds SEEDS[0:3]) on 20 panel lines, R1 real RNA

Every build: extraction -> (genomic) -> repair (F0/F2/F3/F3b) -> WT growth -> KO sweep of all genes ->
57-task capability profile -> open-reaction bitmask. Outputs (results/w4/stage_a/):
  builds.tsv              one row per build x genomic variant (cost, repair size, growth, tasks, open_mask)
  ko_<row>.parquet        build_id, gene, ratio_r (float32), status
  accepted_losses.json    per panel line: lost / contradicted / dropped / accepted genes
  MANIFEST.tsv            argv, commit, dirty, solver pool, workers, wall, CPU, peak RSS, outputs
Checkpoints in results/w4/stage_a/tmp/ (git-ignored); a rerun resumes.

Solvers (amendment A4 c): one declared solver per family, context.DECLARED_SOLVER (GLPK: F0, F2, F4, F5;
Gurobi: F3, F3b); the loss rule (phase 0) runs on GLPK. Any other solver is refused unless --pilot.
Guards: refuses a dirty src/scripts/tests tree (§11.4) unless --pilot; --pilot also requires --out; declared
CPU budget (--budget-cpu-h) with abort at 2x (§11.5); abort if available memory < 10 GB (§11.4); every
Gurobi worker holds a slot of the shared pool nextgem.gurobi_pool (--gurobi-slots, default 4).
Run long jobs (> 1 CPU-h) from a pinned worktree: /data/tmpA/nextgem_runs/<commit> (§11.4).

Usage:
  .venv/bin/python scripts/w4_stage_a.py --dry-run
  .venv/bin/python scripts/w4_stage_a.py --gurobi-workers 4 --glpk-workers 20 --corda-band drop --budget-cpu-h 100
  .venv/bin/python scripts/w4_stage_a.py --pilot --out results/w4/stage_a_smoke_a5 --rows A5 --lines 2 \
      --gurobi-workers 0 --glpk-workers 2                                          # smoke test
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import pickle
import resource
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np
import pandas as pd
import psutil

from nextgem import context_inputs as ci
from nextgem import context_run as cr
from nextgem import gurobi_pool, paths
from nextgem.context import DECLARED_SOLVER
from nextgem.config import SEEDS

W4 = paths.results_dir() / "w4"
OUT = W4 / "stage_a"  # --out overrides (smoke tests never write here)
TMP = OUT / "tmp"
FASTCC = W4 / "fastcc_human2_M0.json"
LP_FAMILIES = ("F0", "F2", "F3", "F4")
LOSS_RULE_SOLVER = "glpk"  # phase 0 (family-independent accepted losses); fixed for all reruns
MIN_FREE_GB = 10.0


# ----------------------------------------------------------------------------- enumeration
def enumerate_jobs(panel: pd.DataFrame, corda_band: str, rows=None, a7_lines=20) -> list[dict]:
    lines = list(panel.ModelID)
    prot_lines = list(panel[panel.has_procan.astype(bool)].ModelID)
    specs = []

    def add(row, line, fam, rule, mod, var, genomic, **kw):
        specs.append(dict(row=row, line=line, family=fam, rule=rule, modality=mod, variant=var, genomic=genomic, **kw))

    for f in LP_FAMILIES:
        for rule in ("R1", "R3"):
            for line in lines:
                add("A1", line, f, rule, "rna", "real", ("none", "loss"))  # A4 re-sweep rides on the real build
                add("A1", line, f, rule, "rna", "perm", ("none",))
        for line in lines:
            add("A2", line, f, "R2", "rna", "real", ("none", "loss"))
        for rule in ("R1", "R2", "R3"):
            for mod in ("protein", "rna_protein"):
                for line in prot_lines:
                    add("A3", line, f, rule, mod, "real", ("none", "loss"))
    if corda_band == "full":
        for rule in ("R1", "R3"):
            for var in ("real", "perm"):
                for line in lines:
                    add("A6", line, "F3b", rule, "rna", var, ("none",))
    elif corda_band == "fallback":
        sub = [lines[i] for i in range(0, len(lines), 2)][:15]  # every other panel line: 1-2 per lineage
        for var in ("real", "perm"):
            for line in sub:
                add("A6", line, "F3b", "R1", "rna", var, ("none",))
    for line in lines:
        mod = "rna_protein" if line in prot_lines else "rna"
        for (w_r, w_p, w_c, tau) in ci.LCP_GRID:
            add("A5", line, "F5", "lcp", mod, f"wr{w_r:g}_wp{w_p:g}_wc{w_c:g}_tau{tau:g}", ("none",),
                lcp=(w_r, w_p, w_c, tau))
    for k, seed in enumerate(SEEDS[0:3]):
        for line in lines[:a7_lines]:
            add("A7", line, "F2", "R1", "rna", f"gimme_seed{seed}", ("none",), gimme_noise_seed=seed)
            add("A7", line, "F3", "R1", "rna", f"noise_seed{seed}", ("none",), noise_seed=seed)
    if rows:
        specs = [s for s in specs if s["row"] in rows]
    for s in specs:
        s["build_id"] = "|".join([s["row"], s["line"], s["family"], s["rule"], s["modality"], s["variant"]])
    # longest first: CORDA, then FASTCORE, GIMME, F0, F4
    order = {"F3b": 0, "F3": 1, "F2": 2, "F0": 3, "F5": 4, "F4": 5}
    return sorted(specs, key=lambda s: (order[s["family"]], s["build_id"]))


# ----------------------------------------------------------------------------- worker side
def worker_init(solver: str, n_slots: int):
    if solver == "gurobi":
        cr._W["slot"] = gurobi_pool.hold_slot(n_slots)  # held for the life of the worker process
    cr.init_worker("human2", solver, str(FASTCC))


def worker_run(job: cr.Job, accepted_detail_line: str | None):
    want = LOSS_RULE_SOLVER if accepted_detail_line is not None else DECLARED_SOLVER[job.family]
    if cr._W["solver"] != want and not job.meta.get("pilot"):
        raise RuntimeError(f"{job.build_id}: worker solver {cr._W['solver']} != declared {want} (amendment A4 c)")
    if accepted_detail_line is not None:  # phase 0: loss rule for one line
        from nextgem import context as cx

        lost, expr = job.meta.pop("_lost"), job.meta.pop("_expr")
        c0 = time.process_time()
        res = cx.accept_losses(cr._W["parent"], lost, expr)
        res["cpu_s"] = time.process_time() - c0
        res["peak_rss_gb"] = cr.peak_rss_gb()
        return {"loss": res}
    return cr.run_job(job)


# ----------------------------------------------------------------------------- parent side
def git_state():
    root = paths.PROJECT_ROOT
    head = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain", "--", "src", "scripts", "tests"],
                           capture_output=True, text=True).stdout.strip()
    return head, bool(dirty)


class Inputs:
    def __init__(self, panel):
        import cobra

        m = cobra.io.read_sbml_model(str(paths.data_root() / "Human-GEM-main/model/Human-GEM.xml"))
        self.genes = ci.gene_table(m)
        del m
        st = json.loads((W4 / "threshold_stats.json").read_text())
        gm = pd.read_parquet(W4 / "threshold_gene_means.parquet").set_index("gene")
        self.rna_stats = ci.ThresholdStats("rna", st["rna"]["p25"], st["rna"]["p75"], gm.rna_mean.dropna(),
                                           st["rna"]["floor_value"], st["rna"]["n_lines"])
        self.prot_stats = ci.ThresholdStats("protein", st["protein"]["p25"], st["protein"]["p75"],
                                            gm.protein_mean.dropna(), st["protein"]["floor_value"], st["protein"]["n_lines"])
        lines = list(panel.ModelID)
        self.rna = ci.read_expression(self.genes, lines)
        self.rna_perm = ci.permute_across_lines(self.rna.loc[lines], SEEDS[0])
        self.cn = ci.read_cn(self.genes, lines)
        self.dmg = ci.read_damaging(self.genes, lines)
        self.prot = ci.read_procan(self.genes, lines)
        self.zstats = ci.ZStats(gm.rna_mean, gm.rna_sd, gm.protein_mean_floor, gm.protein_sd_floor,
                                st["protein"]["floor_value"])

    def scores(self, spec, noise_sd, accepted=None):
        if spec["family"] == "F5":
            return ci.lcp_scores(spec["line"], *spec["lcp"], self.genes, self.rna, self.prot, self.zstats,
                                 accepted[spec["line"]]["accepted"])
        rna = self.rna_perm if spec["variant"] == "perm" else self.rna
        if spec.get("noise_seed") is not None:
            rng = np.random.default_rng(spec["noise_seed"])
            rna = rna.copy()
            rna.loc[spec["line"]] = rna.loc[spec["line"]] + rng.normal(0, noise_sd, rna.shape[1])
        return ci.signed_scores(spec["line"], spec["rule"], spec["modality"], self.genes, rna, self.rna_stats,
                                self.prot, self.prot_stats)


def job_path(build_id: str):
    return TMP / "jobs" / (hashlib.sha1(build_id.encode()).hexdigest()[:16] + ".pkl")


def atomic_pickle(obj, path):
    tmp = path.with_suffix(".part")
    with open(tmp, "wb") as fh:
        pickle.dump(obj, fh)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rows", nargs="*", help="subset of A1 A2 A3 A6 A7 (A4 rides on A1-A3)")
    ap.add_argument("--corda-band", choices=["full", "fallback", "drop"], default="drop")
    ap.add_argument("--gurobi-workers", type=int, default=4)
    ap.add_argument("--gurobi-slots", type=int, default=4)
    ap.add_argument("--glpk-workers", type=int, default=20)
    ap.add_argument("--noise-sd", type=float, default=0.25, help="A7 FASTCORE expression noise (log2 units)")
    ap.add_argument("--budget-cpu-h", type=float, default=100.0)
    ap.add_argument("--limit", type=int, help="pilot: only the first N jobs per family")
    ap.add_argument("--families", nargs="*", help="pilot: subset of families")
    ap.add_argument("--lines", type=int, help="pilot: only the first N panel lines")
    ap.add_argument("--solver-override", nargs="*", default=[], help="pilot: FAMILY=SOLVER (or loss=SOLVER)")
    ap.add_argument("--out", help="output directory (default results/w4/stage_a)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--pilot", action="store_true", help="smoke/pilot mode: allows a dirty tree, subsets and "
                    "solver overrides; requires --out")
    ap.add_argument("--no-tasks", action="store_true")
    a = ap.parse_args()
    global OUT, TMP
    if a.out:
        OUT = paths.PROJECT_ROOT / a.out if not a.out.startswith("/") else type(OUT)(a.out)
        TMP = OUT / "tmp"
    if (a.limit or a.families or a.lines or a.solver_override) and not a.pilot:
        sys.exit("--limit/--families/--lines/--solver-override need --pilot")
    if a.pilot and not a.out and not a.dry_run:
        sys.exit("--pilot needs --out (never write a pilot run into results/w4/stage_a)")
    solver_of = dict(DECLARED_SOLVER, loss=LOSS_RULE_SOLVER)
    for item in a.solver_override:
        fam, _, sol = item.partition("=")
        if fam not in solver_of or sol not in ("glpk", "gurobi"):
            sys.exit(f"bad --solver-override {item!r}")
        solver_of[fam] = sol

    panel = pd.read_csv(W4 / "stage_a_panel.tsv", sep="\t")
    if a.lines:
        panel = panel.head(a.lines)
    split = pd.read_csv(paths.splits_dir() / "lineage_split.tsv", sep="\t").set_index("ModelID")
    assert (split.loc[panel.ModelID, "tier"] == "train").all(), "Stage A uses training lines only"
    specs = enumerate_jobs(panel, a.corda_band, a.rows)
    if a.families:
        specs = [s for s in specs if s["family"] in a.families]
    if a.limit:
        seen = {}
        specs = [s for s in specs if seen.setdefault(s["family"], []).append(1) or len(seen[s["family"]]) <= a.limit]
    n_by_row = pd.Series([s["row"] for s in specs]).value_counts().sort_index().to_dict()
    n_resweeps = sum(1 for s in specs if "loss" in s["genomic"])
    print(f"{len(specs)} builds {n_by_row} + {n_resweeps} A4 genomic re-sweeps (A5: grid points before "
          f"deduplication; A6 band {a.corda_band})")
    print("solvers:", {f: solver_of[f] for f in sorted({s['family'] for s in specs})}, "loss rule:", solver_of["loss"])
    if a.dry_run:
        return
    head, dirty = git_state()
    if dirty and not a.pilot:
        sys.exit("src/scripts/tests are dirty: commit first (DESIGN §11.4), or --pilot for a smoke test")
    workers = {"gurobi": a.gurobi_workers, "glpk": a.glpk_workers}
    needed = {solver_of[s["family"]] for s in specs} | {solver_of["loss"]}
    missing = [k for k in needed if workers[k] < 1]
    if missing:
        sys.exit(f"declared solvers need workers: {missing} (set --{missing[0]}-workers)")
    (TMP / "jobs").mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    t0, c0 = time.time(), time.process_time()
    I = Inputs(panel)
    ctx = mp.get_context("spawn")
    pools = {k: ProcessPoolExecutor(workers[k], ctx, initializer=worker_init, initargs=(k, a.gurobi_slots))
             for k in sorted(needed)}
    cpu_used = 0.0
    budget_s = a.budget_cpu_h * 3600

    def guard():
        free = psutil.virtual_memory().available / 1e9
        if free < MIN_FREE_GB:
            raise SystemExit(f"abort: available memory {free:.1f} GB < {MIN_FREE_GB} GB")
        if cpu_used > 2 * budget_s:
            raise SystemExit(f"abort: CPU {cpu_used / 3600:.1f} h > 2 x budget {a.budget_cpu_h} h")

    # ---- phase 0: accepted losses per panel line (family-independent; on the parent + medium)
    acc_path = OUT / "accepted_losses.json"
    if acc_path.exists():
        accepted = json.loads(acc_path.read_text())
    else:
        futs = {}
        for line in panel.ModelID:
            lost, flags = ci.lost_genes(line, I.dmg, I.cn, I.rna)
            j = cr.Job(build_id=f"loss|{line}", line=line, family="-", rule="-", modality="-", variant="-",
                       scores=pd.Series(dtype=float), meta={"_lost": lost, "_expr": I.rna.loc[line], "pilot": a.pilot})
            futs[pools[solver_of["loss"]].submit(worker_run, j, line)] = (line, flags)
        accepted = {}
        for f, (line, flags) in futs.items():
            r = f.result()["loss"]
            cpu_used += r["cpu_s"]
            accepted[line] = {**{k: r[k] for k in ("lost", "contradicted", "dropped_joint", "accepted")}, **flags}
        acc_path.write_text(json.dumps(accepted, indent=0))

    # ---- A5: score every grid point, build each distinct score vector once per line, map all points
    a5 = [s for s in specs if s["row"] == "A5"]
    if a5:
        canon, gmap, keep = {}, [], set()
        for s in a5:
            sc = I.scores(s, a.noise_sd, accepted)
            h = hashlib.sha1(np.round(sc.reindex(I.genes.index).to_numpy(), 9).tobytes()).hexdigest()
            key = (s["line"], h)
            if key not in canon:
                canon[key] = s["build_id"]
                keep.add(s["build_id"])
                s["_scores"] = sc
            gmap.append({"line": s["line"], "w_r": s["lcp"][0], "w_p": s["lcp"][1], "w_c": s["lcp"][2],
                         "tau": s["lcp"][3], "build_id": canon[key], "n_off": int((sc < 0).sum())})
        pd.DataFrame(gmap).to_csv(OUT / "a5_grid_map.tsv", sep="\t", index=False)
        specs = [s for s in specs if s["row"] != "A5" or s["build_id"] in keep]
        print(f"A5: {len(a5)} grid points -> {len(keep)} distinct builds")

    # ---- phase 1: builds
    todo = [s for s in specs if not job_path(s["build_id"]).exists()]
    print(f"{len(specs) - len(todo)} builds already checkpointed; {len(todo)} to run")
    # one declared solver per family (amendment A4 c); one queue and one in-flight cap per pool
    def pool_of(fam):
        return solver_of[fam]

    queues = {k: [s for s in todo if pool_of(s["family"]) == k] for k in pools}
    caps = {"gurobi": 2 * a.gurobi_workers, "glpk": 2 * a.glpk_workers}
    running = {}

    def submit(spec):
        job = cr.Job(build_id=spec["build_id"], line=spec["line"], family=spec["family"], rule=spec["rule"],
                     modality=spec["modality"], variant=spec["variant"],
                     scores=spec["_scores"] if "_scores" in spec else I.scores(spec, a.noise_sd, accepted),
                     genomic=spec["genomic"], accepted=tuple(accepted[spec["line"]]["accepted"]),
                     gimme_noise_seed=spec.get("gimme_noise_seed"), tasks=not a.no_tasks,
                     meta={"row": spec["row"], "pilot": a.pilot})
        k = pool_of(spec["family"])
        running[pools[k].submit(worker_run, job, None)] = (k, spec)

    n_done = 0
    while any(queues.values()) or running:
        for k, q in queues.items():
            while q and sum(1 for v in running.values() if v[0] == k) < caps[k]:
                submit(q.pop(0))
        done, _ = wait(list(running), return_when=FIRST_COMPLETED)
        for f in done:
            _, spec = running.pop(f)
            res = f.result()
            cpu_used += sum(r.get("variant_cpu_s", 0) + r.get("extract_cpu_s", 0) for r in res["rows"])
            atomic_pickle(res, job_path(spec["build_id"]))
            n_done += 1
            if n_done % 20 == 0:
                print(f"[{time.time() - t0:7.0f}s] {n_done}/{len(todo)} builds; CPU {cpu_used / 3600:.2f} h", flush=True)
        guard()
    for p in pools.values():
        p.shutdown()

    # ---- consolidate
    rows, kos = [], {}
    for s in specs:
        res = pickle.loads(job_path(s["build_id"]).read_bytes())
        rows += res["rows"]
        for gvar, ko in res["ko"].items():
            ko = ko.copy()
            ko.insert(0, "build_id", f"{s['build_id']}|{gvar}")
            kos.setdefault(s["row"] if gvar == "none" else "A4", []).append(ko)
    pd.DataFrame(rows).to_csv(OUT / "builds.tsv", sep="\t", index=False)
    for row, parts in kos.items():
        df = pd.concat(parts, ignore_index=True)
        df["build_id"] = df.build_id.astype("category")
        df.to_parquet(OUT / f"ko_{row}.parquet", index=False)
    wall = time.time() - t0
    man = {"argv": " ".join(sys.argv), "commit": head, "dirty": dirty, "parent": "human2", "medium": "M0",
           "seeds": json.dumps(SEEDS[:3]), "gurobi_workers": a.gurobi_workers, "glpk_workers": a.glpk_workers,
           "solvers": json.dumps({f: solver_of[f] for f in sorted({s["family"] for s in specs})} | {"loss": solver_of["loss"]}),
           "pilot": a.pilot, "n_builds": len(specs), "wall_s": round(wall),
           "worker_cpu_h": round(cpu_used / 3600, 3), "parent_cpu_s": round(time.process_time() - c0),
           "parent_peak_rss_gb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6,
           "worker_peak_rss_gb_max": max((r.get("peak_rss_gb") or 0) for r in rows) if rows else None,
           "outputs": ",".join(sorted(p.name for p in OUT.glob("*.parquet")))}
    mpath = OUT / "MANIFEST.tsv"
    pd.DataFrame([man]).to_csv(mpath, sep="\t", index=False, mode="a", header=not mpath.exists())
    print(json.dumps(man, indent=1))


if __name__ == "__main__":
    main()
