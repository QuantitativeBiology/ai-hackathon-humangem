"""W3 runtime table (DESIGN §6.7, §13): one row per successful script run in results/w3/MANIFEST.tsv, plus totals.

Rows superseded by a later run of the same script and argv are kept but flagged ``superseded``; totals are given
for all runs (spent budget) and for the final runs only. Label-free.
Usage: .venv/bin/python scripts/w3_runtime.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w3_common import MANIFEST, OUT_DIR, atomic_write_tsv  # noqa: E402


def main() -> int:
    m = pd.read_csv(MANIFEST, sep="\t")
    m["argv"] = m.argv.fillna("")
    m["superseded"] = m.duplicated(["script", "argv"], keep="last")
    m["cpu_h"] = m.cpu_s / 3600
    m["wall_min"] = m.wall_s / 60
    gur = m.solver.fillna("").str.contains("gurobi")
    out = m[["timestamp", "script", "argv", "solver", "workers", "wall_min", "cpu_h", "peak_rss_mb", "superseded"]]
    atomic_write_tsv(out, OUT_DIR / "runtime.tsv")
    tot = {"all_runs_cpu_h": round(float(m.cpu_h.sum()), 2), "final_runs_cpu_h": round(float(m.cpu_h[~m.superseded].sum()), 2),
           "gurobi_cpu_h": round(float(m.cpu_h[gur].sum()), 2), "all_runs_wall_h": round(float(m.wall_s.sum() / 3600), 2),
           "n_runs": int(len(m))}
    (OUT_DIR / "runtime_totals.json").write_text(json.dumps(tot, indent=1))
    print(out.to_string(index=False))
    print(json.dumps(tot))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
