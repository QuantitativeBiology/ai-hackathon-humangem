"""W6 end to end in one command: inputs (if missing) -> simulate -> analyze.

    .venv/bin/python scripts/w6_run_all.py --config provisional      # the provisional (pre-P1) run
    .venv/bin/python scripts/w6_run_all.py --config p1               # once P1 is frozen (results/w6/config/p1.json)

Runs above 1 CPU-h must come from a pinned clean worktree (DESIGN §11.4): /data/tmpA/nextgem_runs/<commit>; the
provenance wrapper refuses a dirty src/, scripts/ or tests/ unless --allow-dirty (pilots only) is given.
Each step is a separate process with its own MANIFEST row (results/w6/MANIFEST.tsv).
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable


def run(args: list[str]) -> None:
    env = dict(os.environ)
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[v] = "1"
    print("+", " ".join(args), flush=True)
    r = subprocess.run([PY, *args], env=env)
    if r.returncode:
        raise SystemExit(r.returncode)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="provisional")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--refresh-inputs", action="store_true")
    ap.add_argument("--pilot", action="store_true", help="2 medium keys, 6 targets, a few units per analysis")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    dirty = ["--allow-dirty"] if a.allow_dirty else []
    root = HERE.parent
    if a.refresh_inputs or not (root / "results/w6/tmp/inputs/responses.parquet").exists():
        run([str(HERE / "w6_inputs.py"), *dirty])
    sim = [str(HERE / "w6_simulate.py"), "--config", a.config, "--workers", str(a.workers), *dirty]
    ana = [str(HERE / "w6_analyze.py"), "--config", a.config, *dirty]
    if a.pilot:
        sim += ["--pilot", "2", "--pilot-targets", "6"]
        ana += ["--pilot"]
    run(sim)
    run(ana)
    return 0


if __name__ == "__main__":
    sys.exit(main())
