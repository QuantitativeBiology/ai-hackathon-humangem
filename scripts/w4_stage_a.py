"""W4 Stage A: the frozen LP-only build table on the 30-line training panel (DESIGN.md §9 W4, amendment A4).

Rows (label-free; Human2 parent; medium M0 per A4 f):
  A1  F0, F2, F3, F4 x {R1, R3} x {real, permuted} RNA                         panel (30)      480 builds
  A2  F0, F2, F3, F4 x R2 x real RNA                                              panel (30)      120
  A3  F0, F2, F3, F4 x {R1, R2, R3} x {protein, rna_protein}                      panel w/ ProCan 24 x n_p
  A4  + genomic constraint post-extraction (accepted losses -> KO + repair + sweep) on every real build of A1-A3
  A5  F5 LCP, fixed label-free 24-point grid (A4 d): w_r = 1, w_p in {0, 0.5, 1}, w_c in {0, 1},
      tau in {-1.5, -1, -0.5, 0}; gene off iff w_r z_RNA + w_p z_prot + w_c loss < tau (z against training
      lines; loss = -1 for accepted losses; z_prot = 0 without proteomics). Grid points that give identical
      scores on a line are built once; a5_grid_map.tsv maps all (line, point) pairs to their build.
  A6' F3b CORDA (amendment A7, replaces A6): {R1, R3} x real RNA on 10 representative panel lines (one per lineage,
      closest to the lineage median log2(TPM+1) profile; results/w4/a6prime_lines.tsv), 20 builds on Gurobi; CPU
      limit 5 h per build (3x the pilot cost; failed-timeout above). Tagged "descriptive, not eligible for P1".
  A7  stability, 20 panel lines = the first 2 panel lines of each of the 10 panel lineages: F2 GIMME (cost noise
      1e-3 U(0,1)) and F3 FASTCORE (N(0, --noise-sd) on log2 values), seeds SEEDS[0:3]; noise drawn
      independently per line (seed sequence [seed, panel position])

Every build: extraction -> (genomic) -> repair (F0/F2/F3/F5) -> WT growth -> KO sweep of all genes ->
57-task capability profile -> open-reaction bitmask.

Execution (review B2): every build runs in its OWN fresh process (spawned, one job, then exit), so a build never
depends on what a worker solved before (solver warm starts, variable order). Every LP has a time limit (context
LP_TIME_LIMIT_S); a job that exceeds --job-timeout wall seconds is killed and recorded as failed-timeout; any
other worker exception aborts the run with no consolidated output (§11.4).

Solvers (A4 c): context.DECLARED_SOLVER (GLPK: F0, F2, F4, F5; Gurobi: F3, F3b); loss rule on GLPK. cobra's solver
is set before any model load, so GLPK processes never open a WLS session (review B1). Every Gurobi process holds a
slot of nextgem.gurobi_pool (one global size). Every job asserts that its parent model runs the declared solver
(after --solver-override, which needs --pilot).

Provenance (review B3): main() runs under nextgem.run.provenance(); the manifest row (results/w4/stage_a/MANIFEST.tsv)
hashes every input (omics files, Human-GEM, split, panel, threshold statistics, FASTCC cache) and every output,
and refuses a dirty src/scripts/tests tree unless --pilot. Checkpoints live in tmp/jobs_<code hash>/ (hash of
the W4 code files and the parameters), so a resume after a code change can never mix versions.

Outputs (results/w4/stage_a/): builds.tsv, ko_<row>.parquet (build_id, gene, ratio_r float64, status, solver),
accepted_losses.json, a5_grid_map.tsv, MANIFEST.tsv. Run long jobs from a pinned worktree (§11.4); the script
refuses to run if the imported nextgem package is not the one next to this script.

Usage:
  .venv/bin/python scripts/w4_stage_a.py --dry-run
  .venv/bin/python scripts/w4_stage_a.py --gurobi-workers 8 --glpk-workers 20 --budget-cpu-h 100
  .venv/bin/python scripts/w4_stage_a.py --pilot --out results/w4/stage_a_smoke --lines 2 \
      --gurobi-workers 2 --glpk-workers 6                                           # smoke test
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psutil

from nextgem import context_inputs as ci
from nextgem import context_run as cr
from nextgem import gurobi_pool, paths
from nextgem.config import SEEDS
from nextgem.context import DECLARED_SOLVER, LOSS_RULE_SOLVER
from nextgem.context_run import Runner, worker_run  # noqa: F401  (worker_run re-exported for tests)

W4 = paths.results_dir() / "w4"
OUT = W4 / "stage_a"  # --out overrides (pilot runs never write here)
FASTCC = W4 / "fastcc_human2_M0.json"
LP_FAMILIES = ("F0", "F2", "F3", "F4")
MIN_FREE_GB = 10.0
AUTHOR = "W4 engineer agent"
CODE_FILES = ["src/nextgem/context.py", "src/nextgem/context_inputs.py", "src/nextgem/context_run.py",
              "src/nextgem/gurobi_pool.py", "src/nextgem/knockouts.py", "src/nextgem/model.py", "scripts/w4_stage_a.py"]


# ----------------------------------------------------------------------------- enumeration
def a7_lines(panel: pd.DataFrame, per_lineage: int = 2) -> list[str]:
    """The first ``per_lineage`` panel rows of every panel lineage (panel order), i.e. 20 lines from 10 lineages."""
    return list(panel.groupby("lineage", sort=False).head(per_lineage).ModelID)


CORDA_CPU_LIMIT_S = 5 * 3600.0  # amendment A7: 3x the pilot cost per model


def family_tag(fam: str) -> str:
    return {"F3b": ci.A6P_TAG, "F5": "label-calibrated (F5)"}.get(fam, "")


def enumerate_jobs(panel: pd.DataFrame, a6p_lines=(), rows=None) -> list[dict]:
    lines = list(panel.ModelID)
    pos = {m: i for i, m in enumerate(lines)}
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
    for line in lines:
        mod = "rna_protein" if line in prot_lines else "rna"
        for (w_r, w_p, w_c, tau) in ci.LCP_GRID:
            add("A5", line, "F5", "lcp", mod, f"wr{w_r:g}_wp{w_p:g}_wc{w_c:g}_tau{tau:g}", ("none",),
                lcp=(w_r, w_p, w_c, tau))
    for line in a6p_lines:
        for rule in ("R1", "R3"):
            add("A6p", line, "F3b", rule, "rna", "real", ("none",), cpu_limit_s=CORDA_CPU_LIMIT_S)
    for seed in SEEDS[0:3]:
        for line in a7_lines(panel):
            per_line = [seed, pos[line]]  # independent per line, deterministic (review NOTE 7c)
            add("A7", line, "F2", "R1", "rna", f"gimme_seed{seed}", ("none",), gimme_noise_seed=per_line)
            add("A7", line, "F3", "R1", "rna", f"noise_seed{seed}", ("none",), noise_seed=per_line)
    if rows:
        specs = [s for s in specs if s["row"] in rows]
    for s in specs:
        s["build_id"] = "|".join([s["row"], s["line"], s["family"], s["rule"], s["modality"], s["variant"]])
    order = {"F3b": 0, "F3": 1, "F2": 2, "F0": 3, "F5": 4, "F4": 5}  # longest first
    return sorted(specs, key=lambda s: (order[s["family"]], s["build_id"]))


# ----------------------------------------------------------------------------- parent side
class Inputs:
    def __init__(self, panel):
        import cobra

        cobra.Configuration().solver = "glpk"  # structure only: no WLS session in the parent (review B1)
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


def a6p_subset(panel: pd.DataFrame, split: pd.DataFrame) -> pd.DataFrame:
    """Amendment A7 rule on a (pilot) subset panel; for the full panel this equals results/w4/a6prime_lines.tsv."""
    import cobra

    cobra.Configuration().solver = "glpk"
    genes = ci.gene_table(cobra.io.read_sbml_model(str(paths.data_root() / "Human-GEM-main/model/Human-GEM.xml")))
    sp = split.reset_index()
    tr = sp[(sp.tier == "train") & sp.has_expression.astype(bool) & sp.lineage.isin(set(panel.lineage))]
    out = ci.a6prime_lines(panel, ci.read_expression(genes, tr.ModelID), sp)
    out["status"] = ci.A6P_TAG
    return out


def input_files() -> list[Path]:
    D, dm = paths.data_root(), paths.depmap_dir()
    return [W4 / "stage_a_panel.tsv", W4 / "threshold_stats.json", W4 / "threshold_gene_means.parquet", FASTCC,
            paths.splits_dir() / "lineage_split.tsv", dm / ci.EXPR_FILE, dm / ci.CN_FILE, dm / ci.DAMAGING_FILE,
            dm / "Model.csv", D / ci.PROCAN_FILE, D / "Human-GEM-main/model/Human-GEM.xml",
            D / "Human-GEM-main/model/genes.tsv", D / "Human-GEM-main/data/metabolicTasks/metabolicTasks_Essential.txt"]


def code_hash(params: dict) -> str:
    h = hashlib.sha256()
    for f in CODE_FILES:
        h.update(f.encode())
        h.update((paths.PROJECT_ROOT / f).read_bytes())
    h.update(json.dumps(params, sort_keys=True).encode())
    return h.hexdigest()[:16]


def atomic_pickle(obj, path):
    tmp = path.with_suffix(".part")
    with open(tmp, "wb") as fh:
        pickle.dump(obj, fh)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rows", nargs="*", help="subset of A1 A2 A3 A5 A6p A7 (A4 rides on A1-A3)")
    ap.add_argument("--max-corda", type=int, help="pilot: at most N CORDA (A6') builds")
    ap.add_argument("--gurobi-workers", type=int, default=8, help="<= the global Gurobi pool size")
    ap.add_argument("--glpk-workers", type=int, default=20)
    ap.add_argument("--job-timeout", type=float, default=1800.0, help="wall seconds per build process")
    ap.add_argument("--noise-sd", type=float, default=0.25, help="A7 FASTCORE expression noise (log2 units)")
    ap.add_argument("--budget-cpu-h", type=float, default=100.0)
    ap.add_argument("--limit", type=int, help="pilot: only the first N jobs per family")
    ap.add_argument("--families", nargs="*", help="pilot: subset of families")
    ap.add_argument("--lines", type=int, help="pilot: only the first N panel lines")
    ap.add_argument("--solver-override", nargs="*", default=[], help="pilot: FAMILY=SOLVER (or loss=SOLVER)")
    ap.add_argument("--out", help="output directory (default results/w4/stage_a)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--pilot", action="store_true", help="smoke/pilot mode: allows a dirty tree and subsets; "
                    "requires --out")
    ap.add_argument("--no-tasks", action="store_true")
    return ap.parse_args(argv)


def main(argv=None):
    from nextgem.run import provenance

    a = parse_args(argv)
    global OUT
    if a.out:
        OUT = Path(a.out) if a.out.startswith("/") else paths.PROJECT_ROOT / a.out
    if (a.limit or a.families or a.lines or a.solver_override or a.max_corda) and not a.pilot:
        sys.exit("--limit/--families/--lines/--solver-override/--max-corda need --pilot")
    if a.pilot and not a.out and not a.dry_run:
        sys.exit("--pilot needs --out (never write a pilot run into results/w4/stage_a)")
    if a.gurobi_workers > gurobi_pool.n_slots():
        sys.exit(f"--gurobi-workers {a.gurobi_workers} > global Gurobi pool {gurobi_pool.n_slots()}")
    if Path(__file__).resolve().parents[1] != paths.PROJECT_ROOT:  # pinned-worktree guard (§11.4)
        sys.exit(f"nextgem imported from {paths.PROJECT_ROOT}, script from {Path(__file__).resolve().parents[1]}")
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
    a6p_file = W4 / "a6prime_lines.tsv"
    if a.lines:  # pilot subset: the same A7 rule applied to the subset panel
        a6p = a6p_subset(panel, split)
    else:
        a6p = pd.read_csv(a6p_file, sep="\t")
    a6p_lines = list(a6p.ModelID)
    specs = enumerate_jobs(panel, a6p_lines, a.rows)
    if a.max_corda is not None:
        k = [0]
        specs = [s for s in specs if s["family"] != "F3b" or (k.__setitem__(0, k[0] + 1) or k[0] <= a.max_corda)]
    if a.families:
        specs = [s for s in specs if s["family"] in a.families]
    if a.limit:
        seen = {}
        specs = [s for s in specs if seen.setdefault(s["family"], []).append(1) or len(seen[s["family"]]) <= a.limit]
    n_by_row = pd.Series([s["row"] for s in specs]).value_counts().sort_index().to_dict()
    n_resweeps = sum(1 for s in specs if "loss" in s["genomic"])
    print(f"{len(specs)} builds {n_by_row} + {n_resweeps} A4 genomic re-sweeps (A5: grid points before "
          f"deduplication; A6' lines {a6p_lines})")
    print("solvers:", {f: solver_of[f] for f in sorted({s['family'] for s in specs})}, "loss rule:", solver_of["loss"])
    if a.dry_run:
        return 0
    workers = {"gurobi": a.gurobi_workers, "glpk": a.glpk_workers}
    needed = sorted({solver_of[s["family"]] for s in specs} | {solver_of["loss"]})
    missing = [k for k in needed if workers[k] < 1]
    if missing:
        sys.exit(f"declared solvers need workers: {missing}")
    workers = {k: workers[k] for k in needed}

    params = {"rows": sorted(n_by_row), "lines": list(panel.ModelID), "solvers": solver_of, "noise_sd": a.noise_sd,
              "tasks": not a.no_tasks, "a6p_lines": a6p_lines, "corda_cpu_limit_s": CORDA_CPU_LIMIT_S}
    chash = code_hash(params)
    jobs_dir = OUT / "tmp" / f"jobs_{chash}"
    stale = [p.name for p in (OUT / "tmp").glob("jobs_*") if p.name != jobs_dir.name] if (OUT / "tmp").exists() else []
    if stale:
        print(f"note: checkpoints of other code versions are ignored: {stale}")
    jobs_dir.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    a6p.to_csv(OUT / "a6prime_lines.tsv", sep="\t", index=False)
    with provenance(OUT / "MANIFEST.tsv", inputs=input_files() + ([] if a.lines else [a6p_file]), seeds=SEEDS[:3],
                    solver=",".join(f"{k}:{v}" for k, v in sorted(solver_of.items())),
                    workers=sum(workers.values()), allow_dirty=a.pilot, author=AUTHOR) as prov:
        prov.note(code_hash=chash, pilot=a.pilot, workers=workers, job_timeout_s=a.job_timeout,
                  n_specs=len(specs), parent="human2", medium="M0")
        prov.output(OUT / "a6prime_lines.tsv")
        rc = _run(a, panel, specs, solver_of, workers, jobs_dir, prov)
    return rc


def _run(a, panel, specs, solver_of, workers, jobs_dir, prov):
    t0 = time.time()
    I = Inputs(panel)
    runner = Runner(workers, a.job_timeout, str(FASTCC))
    state = {"cpu": 0.0, "n_done": 0, "timeouts": []}
    budget_s = a.budget_cpu_h * 3600

    def guard():
        free = psutil.virtual_memory().available / 1e9
        if free < MIN_FREE_GB:
            raise SystemExit(f"abort: available memory {free:.1f} GB < {MIN_FREE_GB} GB")
        if state["cpu"] > 2 * budget_s:
            raise SystemExit(f"abort: CPU {state['cpu'] / 3600:.1f} h > 2 x budget {a.budget_cpu_h} h")

    # ---- phase 0: accepted losses per panel line (family-independent; parent + medium; one process per line)
    acc_path = jobs_dir / "accepted_losses.json"
    if acc_path.exists():
        accepted = json.loads(acc_path.read_text())
    else:
        accepted, items = {}, []
        for line in panel.ModelID:
            lost, flags = ci.lost_genes(line, I.dmg, I.cn, I.rna)
            j = cr.Job(build_id=f"loss|{line}", line=line, family="-", rule="-", modality="-", variant="-",
                       scores=pd.Series(dtype=float),
                       meta={"_lost": lost, "_expr": I.rna.loc[line], "solver": solver_of["loss"]})
            items.append((f"loss|{line}", solver_of["loss"], j, line, (line, flags)))

        def on_loss(key, status, res, payload, wall):
            if status != "ok":
                raise SystemExit(f"loss rule failed for {key}: {status} {res}")
            r = res["loss"]
            state["cpu"] += r["cpu_s"]
            line, flags = payload
            accepted[line] = {**{k: r[k] for k in ("lost", "contradicted", "dropped_joint", "accepted")}, **flags}

        runner.run_all(items, on_loss)
        accepted = {line: accepted[line] for line in panel.ModelID}
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
        prov.output(OUT / "a5_grid_map.tsv")
        specs = [s for s in specs if s["row"] != "A5" or s["build_id"] in keep]
        print(f"A5: {len(a5)} grid points -> {len(keep)} distinct builds")

    # ---- phase 1: builds, one fresh process each
    def job_path(build_id):
        return jobs_dir / (hashlib.sha1(build_id.encode()).hexdigest()[:16] + ".pkl")

    todo = [s for s in specs if not job_path(s["build_id"]).exists()]
    print(f"{len(specs) - len(todo)} builds already checkpointed; {len(todo)} to run", flush=True)
    items = []
    for spec in todo:
        lane = solver_of[spec["family"]]
        job = cr.Job(build_id=spec["build_id"], line=spec["line"], family=spec["family"], rule=spec["rule"],
                     modality=spec["modality"], variant=spec["variant"],
                     scores=spec["_scores"] if "_scores" in spec else I.scores(spec, a.noise_sd, accepted),
                     genomic=spec["genomic"], accepted=tuple(accepted[spec["line"]]["accepted"]),
                     gimme_noise_seed=spec.get("gimme_noise_seed"), tasks=not a.no_tasks,
                     meta={"row": spec["row"], "solver": lane,
                           "cpu_limit_s": spec.get("cpu_limit_s", a.job_timeout)})
        items.append((spec["build_id"], lane, job, None, spec))

    def on_build(key, status, res, spec, wall):
        if status == "error":
            raise SystemExit(f"worker error in {key} (run aborted, §11.4):\n{res}")
        if status == "timeout":  # failed-timeout build: recorded, run continues
            res = {"rows": [{"build_id": f"{key}|{g}", "line": spec["line"], "family": spec["family"],
                             "rule": spec["rule"], "modality": spec["modality"], "variant": spec["variant"],
                             "genomic": g, "row": spec["row"], "failed": True, "error": "failed-timeout",
                             "job_wall_s": wall, "solver": solver_of[spec["family"]]} for g in spec["genomic"]],
                   "ko": {}, "error": "failed-timeout"}
            state["timeouts"].append(key)
        else:
            state["cpu"] += sum(r.get("variant_cpu_s", 0) + r.get("extract_cpu_s", 0) for r in res["rows"])
        atomic_pickle(res, job_path(key))
        state["n_done"] += 1
        if state["n_done"] % 20 == 0:
            print(f"[{time.time() - t0:7.0f}s] {state['n_done']}/{len(todo)} builds; build CPU "
                  f"{state['cpu'] / 3600:.2f} h", flush=True)

    runner.run_all(items, on_build, guard)

    # ---- consolidate
    rows, kos = [], {}
    for s in specs:
        res = pickle.loads(job_path(s["build_id"]).read_bytes())
        for r in res["rows"]:
            r["status_tag"] = family_tag(s["family"])
            r["p1_eligible"] = s["family"] != "F3b"  # amendment A7: CORDA is descriptive only
        rows += res["rows"]
        for gvar, ko in res["ko"].items():
            ko = ko.copy()
            ko.insert(0, "build_id", f"{s['build_id']}|{gvar}")
            if s["family"] == "F3b":
                ko["status_tag"] = pd.Categorical([ci.A6P_TAG] * len(ko))
            kos.setdefault(s["row"] if gvar == "none" else "A4", []).append(ko)
    pd.DataFrame(rows).to_csv(OUT / "builds.tsv", sep="\t", index=False)
    (OUT / "accepted_losses.json").write_text(json.dumps(accepted, indent=0))
    prov.output(OUT / "builds.tsv", OUT / "accepted_losses.json")
    for row, parts in kos.items():
        df = pd.concat(parts, ignore_index=True)
        df["build_id"] = df.build_id.astype("category")
        df.to_parquet(OUT / f"ko_{row}.parquet", index=False)
        prov.output(OUT / f"ko_{row}.parquet")
    prov.note(n_builds=len(specs), build_cpu_h=round(state["cpu"] / 3600, 3), timeouts=state["timeouts"],
              worker_peak_rss_gb_max=max((r.get("peak_rss_gb") or 0) for r in rows) if rows else None,
              wall_s=round(time.time() - t0))
    print(f"done: {len(specs)} builds, {len(state['timeouts'])} timeouts, build CPU {state['cpu'] / 3600:.2f} h, "
          f"wall {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
