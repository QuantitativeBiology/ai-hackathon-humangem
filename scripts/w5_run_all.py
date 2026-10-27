"""Run the whole W5 pipeline in order with the final code and check determinism (every output SHA-256 before vs after).

Order: w5_enumerate (H2, H1) -> w5_curate -> w5_apply -> w5_factorial (H1) -> w5_score_train -> w5_val_dryrun.
Each step writes its own MANIFEST.tsv row (nextgem.run.provenance). Writes results/w5/rerun_sha_check.tsv and
results/w5/runtime.tsv (from the MANIFEST rows of this run).
Usage: OMP_NUM_THREADS=1 .venv/bin/python scripts/w5_run_all.py [--workers 6]
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "w5"
PY = str(ROOT / ".venv" / "bin" / "python")


def shas() -> dict:
    out = {}
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and "logs" not in p.parts and p.name not in ("MANIFEST.tsv", "rerun_sha_check.tsv", "runtime.tsv", "README.md") \
                and not p.name.endswith(".tmp"):
            out[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    w = str(a.workers)
    steps = [["scripts/w5_enumerate.py", "--parent", "human2", "--workers", w],
             ["scripts/w5_enumerate.py", "--parent", "human1", "--workers", w],
             ["scripts/w5_curate.py"],
             ["scripts/w5_apply.py", "--workers", w],
             ["scripts/w5_factorial.py", "--parent", "human1"],
             ["scripts/w5_score_train.py", "--workers", w],
             ["scripts/w5_val_dryrun.py", "--methods", "H1_TYMS", "H1_W5a", "H1_W5a_pid30", "H1_W5b"]]
    before = shas()
    t0 = time.strftime("%Y-%m-%dT%H:%M:%S")
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    (OUT / "logs").mkdir(exist_ok=True)
    for s in steps:
        log = OUT / "logs" / ("run_all_" + "_".join(Path(x).stem if x.endswith(".py") else x for x in s[:3] if not x.startswith("--")) + ".log")
        r = subprocess.run([PY, *s], cwd=ROOT, env=env, stdout=open(log, "w"), stderr=subprocess.STDOUT)
        print(f"{' '.join(s)} -> rc {r.returncode}", flush=True)
        if r.returncode:
            return r.returncode
    after = shas()
    rows = [{"file": f, "sha_before": before.get(f, ""), "sha_after": after.get(f, ""),
             "identical": before.get(f) == after.get(f)} for f in sorted(set(before) | set(after))]
    pd.DataFrame(rows).to_csv(OUT / "rerun_sha_check.tsv", sep="\t", index=False)
    m = pd.read_csv(OUT / "MANIFEST.tsv", sep="\t")
    m = m[m.timestamp >= t0]
    rt = m[["script", "argv", "workers", "wall_s", "cpu_s", "peak_rss_mb", "solver"]].copy()
    tot = {"script": "TOTAL", "argv": "", "workers": "", "wall_s": rt.wall_s.sum(), "cpu_s": rt.cpu_s.sum(),
           "peak_rss_mb": rt.peak_rss_mb.max(), "solver": "glpk"}
    rt = pd.concat([rt, pd.DataFrame([tot])], ignore_index=True)
    rt.to_csv(OUT / "runtime.tsv", sep="\t", index=False)
    n_diff = sum(not r["identical"] for r in rows)
    print(f"{len(rows)} outputs, {n_diff} changed; total wall {tot['wall_s']:.0f} s, CPU {tot['cpu_s']:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
