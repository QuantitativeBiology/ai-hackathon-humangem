"""Check that two runs of scripts/run_baseline_ko.py agree bitwise on ratio_r (and on genes/statuses).

Usage: .venv/bin/python scripts/compare_ko_runs.py DIR_A DIR_B [--models human1 human2]
Exit code 1 if any model differs. Raw ``ratio`` differences (solver noise) are reported, not tested.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("dir_a", type=Path)
    parser.add_argument("dir_b", type=Path)
    parser.add_argument("--models", nargs="+", default=["human1", "human2"])
    args = parser.parse_args()

    ok = True
    for name in args.models:
        a = pd.read_parquet(args.dir_a / f"{name}_ko.parquet")
        b = pd.read_parquet(args.dir_b / f"{name}_ko.parquet")
        same = (
            a.gene.equals(b.gene)
            and a.status.equals(b.status)
            and np.array_equal(a.ratio_r.to_numpy().view(np.int64), b.ratio_r.to_numpy().view(np.int64))
        )
        ok &= same
        print(
            f"{name}: ratio_r bitwise identical={same}; raw ratio bitwise-equal rows "
            f"{int((a.ratio.to_numpy() == b.ratio.to_numpy()).sum())}/{len(a)}, "
            f"max |raw diff|={(a.ratio - b.ratio).abs().max():.2g}"
        )
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
