"""Run the whole W1 (M2, training tier) pipeline in order, time each step, and report which outputs changed
against the SHA-256 values of the outputs already on disk (a determinism check when re-run).

Steps: w1_inputs -> w1_gene_panels -> w1_backgrounds -> w1_p57 (human1, human2; Gurobi, 8 sessions) -> w1_b6
-> w1_score_train -> w1_binary -> w1_verify -> w1_verify --task-equality (Gurobi, 8 sessions).
The validation look (scripts/evaluate_val.py) is NOT part of this driver: it runs once, after review and freeze.

Usage: NEXTGEM_ALLOW_DIRTY=1 .venv/bin/python scripts/w1_run_all.py [--skip-gurobi]
(NEXTGEM_ALLOW_DIRTY=1 only while the W1 code is uncommitted; every step's MANIFEST row then records dirty=True.)
Output: results/w1/runtime.tsv (step, wall_s, exit code), results/w1/rerun_sha_check.tsv
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w1_common import OUT_DIR, PROJECT, atomic_write_tsv, current, run_main, sha256  # noqa: E402

PY = str(PROJECT / ".venv" / "bin" / "python")
STEPS = [
    ("inputs", ["scripts/w1_inputs.py"], {"OMP_NUM_THREADS": "2"}, False),
    ("gene_panels", ["scripts/w1_gene_panels.py"], {}, False),
    ("backgrounds", ["scripts/w1_backgrounds.py"], {}, False),
    ("p57_human1", ["scripts/w1_p57.py", "--model", "human1", "--workers", "7"], {"OMP_NUM_THREADS": "1"}, True),
    ("p57_human2", ["scripts/w1_p57.py", "--model", "human2", "--workers", "7"], {"OMP_NUM_THREADS": "1"}, True),
    ("b6", ["scripts/w1_b6.py"], {"OMP_NUM_THREADS": "8"}, False),
    ("score_train", ["scripts/w1_score_train.py", "--workers", "14"], {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}, False),
    ("binary", ["scripts/w1_binary.py"], {"OMP_NUM_THREADS": "2"}, False),
    ("verify", ["scripts/w1_verify.py"], {"OMP_NUM_THREADS": "2"}, False),
    ("verify_tasks", ["scripts/w1_verify.py", "--task-equality"], {"OMP_NUM_THREADS": "1"}, True),
]
TRACKED = ["gene_features.tsv", "line_features.tsv", "background_groups.tsv", "p57_human1.parquet", "p57_human2.parquet",
           "train_grid.tsv", "train_per_line.tsv", "train_per_gene.tsv", "train_per_lineage.tsv", "train_addressable.tsv",
           "sesoi_c1b.json", "train_binary.tsv", "train_binary_grid.tsv", "train_binary_paired.tsv", "verification.tsv",
           "verification_tasks.tsv", "tmp/b6_oof.parquet"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-gurobi", action="store_true")
    a = ap.parse_args(argv)
    before = {f: sha256(OUT_DIR / f) for f in TRACKED if (OUT_DIR / f).exists()}
    panels = PROJECT / "data" / "controls" / "gene_panels.tsv"
    before["data/controls/gene_panels.tsv"] = sha256(panels) if panels.exists() else ""
    rows = []
    logs = OUT_DIR / "tmp" / "logs"  # git-ignored (Gurobi licence banners)
    logs.mkdir(parents=True, exist_ok=True)
    for name, args, env, gurobi in STEPS:
        if gurobi and a.skip_gurobi:
            continue
        t0 = time.time()
        with open(logs / f"run_all_{name}.log", "w") as fh:
            r = subprocess.run([PY, *args], cwd=PROJECT, env={**os.environ, **env}, stdout=fh, stderr=subprocess.STDOUT)
        rows.append({"step": name, "wall_s": round(time.time() - t0, 1), "exit": r.returncode, "gurobi": gurobi})
        print(rows[-1], flush=True)
        if r.returncode != 0:
            break
    atomic_write_tsv(pd.DataFrame(rows), OUT_DIR / "runtime.tsv")
    after = {f: sha256(OUT_DIR / f) for f in TRACKED if (OUT_DIR / f).exists()}
    after["data/controls/gene_panels.tsv"] = sha256(panels) if panels.exists() else ""
    chk = pd.DataFrame([{"file": f, "before": before.get(f, "")[:16], "after": after.get(f, "")[:16],
                         "identical": before.get(f) == after.get(f)} for f in sorted(set(before) | set(after))])
    atomic_write_tsv(chk, OUT_DIR / "rerun_sha_check.tsv", float_format=None)
    print(chk.to_string(index=False))
    current().output(OUT_DIR / "runtime.tsv", OUT_DIR / "rerun_sha_check.tsv")
    current().note(steps=rows)
    return 0 if all(r["exit"] == 0 for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(run_main(main))
