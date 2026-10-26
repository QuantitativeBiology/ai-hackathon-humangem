"""Draw the lineage split exactly once (DESIGN.md §4.1, as amended by A1). Label-free.

Production (run once, by the PI, after review and amendment A1):

    .venv/bin/python scripts/make_splits.py --seed S --freeze-commit 61b0beb --code-commit <C>

Refuses (exit 2, nothing written) if:
- production mode is given ``--reference-commit`` other than the freeze commit, ``--repo`` other than the
  project, or redirected roots (``NEXTGEM_DATA_ROOT``/``NEXTGEM_STATE_ROOT``) (review M1);
- S is not a positive integer, or any attempt seed S+a (a = 0..99) lies in the retired range
  20261026..20261125 (S in 20260927..20261125); re-checked per attempt inside the loop (V-2);
- git (V-1 as amended by A1): the freeze commit is not an ancestor of C, or C not of HEAD;
  ``git diff --quiet C HEAD -- src scripts tests DESIGN.md FROZEN.tsv`` fails; those paths are dirty or
  have untracked files; when C differs from the freeze commit, DESIGN.md changed since the freeze
  anywhere except added rows of the §17 amendment table (the placeholder row may be replaced), or
  HEAD:FROZEN.tsv lacks the current SHA-256 of make_splits.py, splits.py and cohort.py;
- ``results/splits/SEED_SOURCE.md`` is not committed (unmodified) or does not contain S and the freeze SHA;
- ``results/splits/SPLIT_SHA256`` already exists (the split is never redrawn);
- pins (review M3): an input file's SHA-256 differs from ``data/checksums.tsv`` at the freeze commit;
  the ``split.*`` value rows of FROZEN.tsv (HEAD in production) differ from what is computed now (input
  manifest, line->lineage map, per-lineage expression/proteomics counts, draw parameters);
  per-lineage expression/proteomics counts differ from ``coverage_by_lineage`` in the frozen JSON;
- cross-checks: strata; the 20 reference draws; (production only) the acceptance statistics for seeds
  1-1000; ``splits.draw``/``checks`` against the functions extracted verbatim from
  ``scripts/design_numbers.py`` at the freeze commit (seeds 1-1000, and S..S+99); and the accepted
  attempt of ``splits.run_attempts`` against a loop run on the frozen functions.

Writes ``results/splits/``: ``lineage_split.tsv``, ``split_summary.tsv``, ``attempt_log.tsv``,
``split_meta.json`` and ``SPLIT_SHA256``. Test-ledger rows ``split_started`` (committed before the draw;
the SPLIT_SHA256 marker is created only after that commit, review S2) and ``split_completed``.

Sandbox (tests and review only): ``--sandbox DIR`` writes everything under DIR, keeps the ledger
uncommitted, reports the git checks without enforcing them (``--enforce-git`` enforces them against
``--repo``), and requires S >= 500000 (review M2). Simulated seeds 38..137 are never evaluated outside
production: the 1-1000 checks skip them in sandbox mode.

``--compute-pins`` writes the ``split.*`` pin rows into the working-tree FROZEN.tsv (no draw).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CHECK_PATHS = ["src", "scripts", "tests", "DESIGN.md", "FROZEN.tsv"]
PINNED_CODE = ["scripts/make_splits.py", "src/nextgem/splits.py", "src/nextgem/cohort.py"]
SEED_SOURCE_REL = "results/splits/SEED_SOURCE.md"
SANDBOX_MIN_SEED = 500_000
EXCLUDED_SIM = range(38, 138)  # never evaluated outside production (coordinator / review M2)
PLACEHOLDER_ROW = re.compile(r"^\|\s*—\s*\|\s*—\s*\|\s*—\s*\|\s*—\s*\|\s*—\s*\|\s*$")


class Refused(SystemExit):
    def __init__(self, msg: str):
        print(f"REFUSED: {msg}", file=sys.stderr)
        super().__init__(2)


def sh(repo: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canon_sha(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def resolve(repo: Path, rev: str) -> str:
    r = sh(repo, "rev-parse", "--verify", f"{rev}^{{commit}}")
    if r.returncode != 0:
        raise Refused(f"{rev!r} is not a commit in {repo}")
    return r.stdout.strip()


# ----------------------------------------------------------------------------- git checks (V-1 / A1)
def design_additions_failures(repo: Path, freeze: str, head: str) -> list[str]:
    """DESIGN.md may change after the freeze only by added rows of the §17 amendment table (A1)."""
    diff = sh(repo, "diff", "-U0", freeze, head, "--", "DESIGN.md").stdout
    if not diff.strip():
        return []
    text = sh(repo, "show", f"{head}:DESIGN.md").stdout.splitlines()
    start = next((i for i, ln in enumerate(text) if re.match(r"###\s*Amendment log", ln)), None)
    span = set()
    if start is not None:
        i = start + 1
        while i < len(text) and (text[i].startswith("|") or not text[i].strip()):
            if text[i].startswith("|"):
                span.add(i + 1)  # 1-based
            elif span:
                break
            i += 1
    fails, new_ln = [], 0
    for ln in diff.splitlines():
        m = re.match(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,\d+)? @@", ln)
        if m:
            new_ln = int(m.group(1))
            continue
        if ln.startswith(("+++", "---")):
            continue
        if ln.startswith("-"):
            if not PLACEHOLDER_ROW.match(ln[1:]):
                fails.append("DESIGN.md lost or changed a line since the freeze commit (only §17 amendment rows may be added)")
        elif ln.startswith("+"):
            if new_ln not in span or not ln[1:].startswith("|"):
                fails.append(f"DESIGN.md line {new_ln} was added outside the §17 amendment table")
            new_ln += 1
    return sorted(set(fails))


def git_checks(repo: Path, freeze: str, code: str, seed: int) -> list[str]:
    fails = []
    head = resolve(repo, "HEAD")
    if sh(repo, "merge-base", "--is-ancestor", freeze, head).returncode != 0:
        fails.append(f"freeze commit {freeze[:10]} is not an ancestor of HEAD {head[:10]}")
    if sh(repo, "merge-base", "--is-ancestor", freeze, code).returncode != 0:
        fails.append(f"code commit {code[:10]} does not descend from the freeze commit {freeze[:10]}")
    if sh(repo, "merge-base", "--is-ancestor", code, head).returncode != 0:
        fails.append(f"code commit {code[:10]} is not an ancestor of HEAD")
    if sh(repo, "diff", "--quiet", code, head, "--", *CHECK_PATHS).returncode != 0:
        changed = sh(repo, "diff", "--name-only", code, head, "--", *CHECK_PATHS).stdout.split()
        fails.append(f"code/DESIGN/FROZEN differ between {code[:10]} and HEAD: {changed[:10]}")
    dirty = sh(repo, "status", "--porcelain", "--untracked-files=all", "--", *CHECK_PATHS).stdout.strip()
    if dirty:
        fails.append(f"working tree has uncommitted or untracked changes under {CHECK_PATHS}: {dirty.splitlines()[:10]}")
    if code != freeze:
        frozen_sha = frozen_rows(repo, head, "file")
        for rel in PINNED_CODE:
            p = repo / rel
            if frozen_sha.get(rel) is None or not p.exists() or frozen_sha[rel] != sha256(p):
                fails.append(f"FROZEN.tsv at HEAD does not hold the current SHA-256 of {rel}")
        fails += design_additions_failures(repo, freeze, head)
    if sh(repo, "cat-file", "-e", f"HEAD:{SEED_SOURCE_REL}").returncode != 0:
        fails.append(f"{SEED_SOURCE_REL} is not committed")
    else:
        if sh(repo, "diff", "--quiet", "HEAD", "--", SEED_SOURCE_REL).returncode != 0:
            fails.append(f"{SEED_SOURCE_REL} has uncommitted changes")
        text = sh(repo, "show", f"HEAD:{SEED_SOURCE_REL}").stdout
        if not re.search(rf"(?<![0-9]){seed}(?![0-9])", text):
            fails.append(f"{SEED_SOURCE_REL} does not contain the seed {seed}")
        if freeze[:7] not in text:
            fails.append(f"{SEED_SOURCE_REL} does not name the freeze commit {freeze[:7]}")
    return fails


def frozen_rows(repo: Path, rev: str | None, kind: str) -> dict[str, str]:
    """{path (file rows) or item (value rows): sha/value} from FROZEN.tsv at ``rev`` (None = working tree)."""
    if rev is None:
        f = repo / "FROZEN.tsv"
        text = f.read_text() if f.exists() else ""
    else:
        r = sh(repo, "show", f"{rev}:FROZEN.tsv")
        text = r.stdout if r.returncode == 0 else ""
    out = {}
    for ln in text.splitlines()[1:]:
        f = ln.split("\t")
        if len(f) >= 4 and f[1] == kind:
            out[f[2] if kind == "file" else f[0]] = f[3]
    return out


# ----------------------------------------------------------------------------- pins (M3)
def compute_pins(cohort, data_root: Path) -> tuple[dict[str, str], dict[str, str]]:
    """(pin value rows, per-input SHA-256). Label-free; no draw."""
    from nextgem import config, splits
    from nextgem.cohort import INPUT_FILES

    inputs = {rel: sha256(data_root / rel) for rel in INPUT_FILES}
    params = {"MAX_SPLIT_ATTEMPTS": config.MAX_SPLIT_ATTEMPTS, "RETIRED_SEED_RANGE": list(config.RETIRED_SEED_RANGE),
              "MIN_LINES_ELIGIBLE": splits.MIN_LINES_ELIGIBLE, "STRATA": splits.STRATA, "ALLOC_TEST": splits.ALLOC_TEST,
              "ALLOC_VAL": splits.ALLOC_VAL, "rng": "numpy.random.default_rng(S + attempt)"}
    pins = {
        "split.inputs_sha256": canon_sha(inputs),
        "split.line_lineage_sha256": canon_sha({k: str(v) for k, v in cohort.lineage.items()}),
        "split.expr_by_lineage_sha256": canon_sha(cohort.by_lineage(cohort.expr)),
        "split.prot_by_lineage_sha256": canon_sha(cohort.by_lineage(cohort.proteomics_any)),
        "split.draw_params_sha256": canon_sha(params),
    }
    return pins, inputs


def check_pins(pins, inputs, freeze: str, recorded: dict[str, str], ref: dict) -> list[str]:
    fails = []
    ck = sh(PROJECT_ROOT, "show", f"{freeze}:data/checksums.tsv").stdout
    frozen_ck = {ln.split("\t")[1]: ln.split("\t")[0] for ln in ck.splitlines()[1:] if ln.count("\t") == 1}
    for rel, s in inputs.items():
        if rel in frozen_ck and frozen_ck[rel] != s:
            fails.append(f"input {rel} differs from data/checksums.tsv at the freeze commit")
    for k, v in pins.items():
        if recorded.get(k) != v:
            fails.append(f"pin {k} missing from FROZEN.tsv or different from the recorded value")
    return fails


def check_coverage(cohort, ref: dict) -> list[str]:
    cov = ref["coverage_by_lineage"]
    e, p = cohort.by_lineage(cohort.expr), cohort.by_lineage(cohort.proteomics_any)
    fails = []
    if {L: v["expr"] for L, v in cov.items()} != e:
        fails.append("per-lineage expression counts differ from coverage_by_lineage in the frozen JSON")
    if {L: v["proteomics_any"] for L, v in cov.items()} != p:
        fails.append("per-lineage proteomics counts differ from coverage_by_lineage in the frozen JSON")
    return fails


# ----------------------------------------------------------------------------- cross-checks
def frozen_attempt_loop(seed: int, fdraw, fchecks):
    """The attempt loop of §4.1 step 5 run on the frozen draw/checks functions."""
    from nextgem import config

    for attempt in range(config.MAX_SPLIT_ATTEMPTS):
        t, v = fdraw(seed + attempt)
        if all(fchecks(t, v).values()):
            return attempt, list(t), list(v)
    return None


def crosscheck(ref_commit: str, ref: dict, cohort, seed: int, production: bool, log) -> dict:
    from nextgem import splits

    sim = ref["split_simulation"]
    counts = cohort.lineage_counts()
    if counts != ref["cohort"]["lineage_counts"]:
        raise Refused("screened-line lineage counts differ from design_numbers.json")
    _, _, strata = splits.eligible_and_strata(counts)
    if strata != sim["strata"]:
        raise Refused("strata differ from design_numbers.json")
    expr_by_lin = cohort.by_lineage(cohort.expr)
    prot_by_lin = cohort.by_lineage(cohort.proteomics_any)
    for sd, d in sim["reference_draws_seeds_1_to_20"].items():
        t, v = splits.draw(int(sd), strata)
        if t != d["test"] or v != d["val"]:
            raise Refused(f"reference draw for seed {sd} not reproduced")
    sim_seeds = list(range(1, 1001)) if production else [s for s in range(1, 1001) if s not in EXCLUDED_SIM]
    if production:
        stats = splits.acceptance_stats(sim_seeds, strata, expr_by_lin, prot_by_lin)
        for k, v in stats.items():
            rv = sim[k]
            same = (v == rv) if not isinstance(v, list) else (len(v) == len(rv) and all(abs(a - b) < 1e-9 for a, b in zip(v, rv)))
            if not same:
                raise Refused(f"acceptance statistic {k} differs from design_numbers.json")
    fdraw, fchecks = splits.frozen_reference_functions(PROJECT_ROOT, ref_commit, strata, expr_by_lin, prot_by_lin)
    for sd in sim_seeds + list(range(seed, seed + 100)):
        t, v = splits.draw(sd, strata)
        if (t, v) != tuple(map(list, fdraw(sd))) or splits.checks(t, v, expr_by_lin, prot_by_lin) != fchecks(t, v):
            raise Refused("splits.draw/checks differ from the frozen design_numbers.py functions")
    log("cross-check passed: strata, 20 reference draws, "
        + ("acceptance stats (1..1000), " if production else "")
        + "frozen-function equality")
    return {"strata": strata, "expr_by_lin": expr_by_lin, "prot_by_lin": prot_by_lin, "fdraw": fdraw, "fchecks": fchecks}


# ----------------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, help="user-supplied split seed S")
    ap.add_argument("--freeze-commit", required=True, help="DESIGN freeze commit (61b0beb)")
    ap.add_argument("--code-commit", help="commit holding the reviewed M1 code (default: the freeze commit)")
    ap.add_argument("--reference-commit", help="sandbox only: reference commit for design_numbers.{py,json}")
    ap.add_argument("--sandbox", type=Path, help="write everything under this directory (tests/review only)")
    ap.add_argument("--repo", type=Path, help="sandbox only: repository for the git checks")
    ap.add_argument("--enforce-git", action="store_true", help="sandbox: enforce the git checks")
    ap.add_argument("--compute-pins", action="store_true", help="write the split.* pin rows to FROZEN.tsv (no draw)")
    a = ap.parse_args(argv)
    t0 = time.time()

    from nextgem import splits

    production = not a.sandbox and not a.compute_pins
    if production or a.compute_pins:  # review M1
        if a.repo is not None and a.repo.resolve() != PROJECT_ROOT:
            raise Refused("--repo is allowed in sandbox mode only")
        if os.environ.get("NEXTGEM_DATA_ROOT") or os.environ.get("NEXTGEM_STATE_ROOT"):
            raise Refused("NEXTGEM_DATA_ROOT/NEXTGEM_STATE_ROOT overrides are allowed in sandbox mode only")
    if production or a.compute_pins:
        ref_commit = resolve(PROJECT_ROOT, a.freeze_commit)
        if a.reference_commit and resolve(PROJECT_ROOT, a.reference_commit) != ref_commit:
            raise Refused("--reference-commit must be the freeze commit in production mode")
    else:  # sandbox: the git checks may run on a toy repo; the reference stays in the project
        ref_commit = resolve(PROJECT_ROOT, a.reference_commit or a.freeze_commit)
    ref = json.loads(sh(PROJECT_ROOT, "show", f"{ref_commit}:results/design/design_numbers.json").stdout or "null")
    if not ref:
        raise Refused(f"cannot read results/design/design_numbers.json at {ref_commit}")

    from nextgem import cohort as cohort_mod, ledger, paths

    if a.compute_pins:
        cohort = cohort_mod.load()
        fails = check_coverage(cohort, ref)
        if fails:
            raise Refused("; ".join(fails))
        pins, _ = compute_pins(cohort, paths.data_root())
        from nextgem import frozen

        frozen.upsert_rows([{"item": k, "kind": "value", "sha256_or_value": v,
                             "note": "split input pin (make_splits.py --compute-pins; review M3)"} for k, v in pins.items()])
        print(f"wrote {len(pins)} split.* pin rows to {paths.frozen_tsv()}")
        return 0

    if a.seed is None:
        raise Refused("--seed is required")
    try:
        splits.check_seed(a.seed)
    except splits.SeedRefused as e:
        raise Refused(str(e))
    if a.sandbox:
        if a.seed < SANDBOX_MIN_SEED:  # review M2
            raise Refused(f"sandbox seeds must be >= {SANDBOX_MIN_SEED}")
        sb = a.sandbox.resolve()
        if sb == PROJECT_ROOT or str(sb).startswith(str(PROJECT_ROOT) + os.sep):
            raise Refused("the sandbox must be outside the project tree")
        os.environ["NEXTGEM_STATE_ROOT"] = str(sb)
        sb.mkdir(parents=True, exist_ok=True)

    out = paths.splits_dir()
    if (out / "SPLIT_SHA256").exists():
        raise Refused(f"{out / 'SPLIT_SHA256'} exists: the split is drawn once and never redrawn")

    repo = (a.repo or PROJECT_ROOT).resolve()
    freeze = resolve(repo, a.freeze_commit)
    code = resolve(repo, a.code_commit) if a.code_commit else freeze
    fails = git_checks(repo, freeze, code, a.seed)
    enforce = production or a.enforce_git
    if fails and enforce:
        raise Refused("; ".join(fails))

    log = print
    cohort = cohort_mod.load()
    pin_fails = check_coverage(cohort, ref)
    pins, inputs = compute_pins(cohort, paths.PROJECT_ROOT / "data")
    recorded = frozen_rows(PROJECT_ROOT, "HEAD" if production else None, "value")
    pin_fails += check_pins(pins, inputs, ref_commit, recorded, ref)
    if pin_fails:
        raise Refused("; ".join(pin_fails))
    xc = crosscheck(ref_commit, ref, cohort, a.seed, production, log)

    script_shas = {rel: sha256(PROJECT_ROOT / rel) for rel in PINNED_CODE}
    seed_source = {"path": SEED_SOURCE_REL, "sha256": None, "mode": "sandbox (simulated seed)" if a.sandbox else "production"}
    if production:
        seed_source["sha256"] = sha256(repo / SEED_SOURCE_REL)
    started = {"seed": a.seed, "freeze_commit": freeze, "code_commit": code, "reference_commit": ref_commit,
               "script_sha256": script_shas, "pins": pins, "seed_source": seed_source, "sandbox": bool(a.sandbox)}
    # ledger row first (committed in production), then the exclusive marker (review S2)
    ledger.append("test", "split_started", look_id="split", detail=json.dumps(started, sort_keys=True),
                  author="scripts/make_splits.py", commit=production)
    out.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(out / "SPLIT_SHA256", os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        raise Refused("SPLIT_SHA256 appeared concurrently")
    os.write(fd, b"PENDING draw in progress\n")
    os.close(fd)

    test, val, attempts = splits.run_attempts(a.seed, xc["strata"], xc["expr_by_lin"], xc["prot_by_lin"])
    ref_loop = frozen_attempt_loop(a.seed, xc["fdraw"], xc["fchecks"])
    if ref_loop != (len(attempts) - 1, list(test), list(val)):
        raise Refused("the attempt loop differs from the loop run on the frozen functions (SPLIT_SHA256 left PENDING)")

    lt = cohort_mod.line_table(cohort)
    stratum_of = {L: s for s, ls in xc["strata"].items() for L in ls}
    lt["stratum"] = lt.lineage.map(stratum_of).fillna("rare")
    lt["tier"] = "train"
    lt.loc[lt.lineage.isin(val), "tier"] = "validation"
    lt.loc[lt.lineage.isin(test), "tier"] = "test"
    lt = lt[["ModelID", "tier", "stratum", "lineage"] + [c for c in lt.columns if c not in ("ModelID", "tier", "stratum", "lineage")]]
    lt = lt.sort_values("ModelID", ignore_index=True)
    assert lt.ModelID.is_unique and len(lt) == len(cohort.screened)
    assert not (set(test) & set(val)) and lt.groupby("lineage").tier.nunique().max() == 1
    assert set(lt.loc[lt.stratum == "rare", "tier"]) == {"train"}

    def atomic_tsv(df, name):
        tmp = out / (name + ".tmp")
        df.to_csv(tmp, sep="\t", index=False)
        os.replace(tmp, out / name)

    atomic_tsv(lt, "lineage_split.tsv")
    summ = lt.groupby(["tier", "stratum", "lineage"]).agg(
        n_lines=("ModelID", "size"), n_expression=("has_expression", "sum"), n_cn=("has_cn", "sum"),
        n_proteomics_any=("has_proteomics_any", "sum"), n_nusinow=("has_nusinow", "sum"), n_procan=("has_procan", "sum"),
        n_organoid=("model_type", lambda s: int((s == "Organoid").sum()))).reset_index()
    atomic_tsv(summ, "split_summary.tsv")
    import pandas as pd

    atomic_tsv(pd.DataFrame(attempts), "attempt_log.tsv")
    split_sha = sha256(out / "lineage_split.tsv")
    meta = {**started, "n_attempts": len(attempts), "n_redraws": len(attempts) - 1,
            "rejected_attempts": [r["attempt"] for r in attempts if not r["accepted"]],
            "git_checks": {"enforced": enforce, "failures": fails},
            "split_sha256": split_sha, "n_lines": int(len(lt)),
            "lines_per_tier": {k: int(v) for k, v in lt.tier.value_counts().items()},
            "lineages_per_tier": {k: int(v) for k, v in lt.groupby("tier").lineage.nunique().items()},
            "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "runtime_s": round(time.time() - t0, 1)}
    (out / "split_meta.json").write_text(json.dumps(meta, indent=1, sort_keys=True))
    tmp = out / "SPLIT_SHA256.tmp"
    tmp.write_text(f"{split_sha}  lineage_split.tsv\n")
    os.replace(tmp, out / "SPLIT_SHA256")
    ledger.append("test", "split_completed", look_id="split", output_sha=split_sha,
                  detail=json.dumps({"n_redraws": meta["n_redraws"], "lines_per_tier": meta["lines_per_tier"]}, sort_keys=True),
                  author="scripts/make_splits.py", commit=production)
    print(json.dumps({k: meta[k] for k in ("split_sha256", "n_lines", "lines_per_tier", "lineages_per_tier", "n_redraws",
                                           "runtime_s")}, indent=1))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
