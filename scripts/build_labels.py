"""Build the frozen label files once at M1 (DESIGN.md §5, §11.1). Entry point permitted by the tripwire.

Requires the split (results/splits/lineage_split.tsv + SPLIT_SHA256) and the baseline gene maps.
Writes the tier files (training under data/labels/, validation and test sealed under data/locked/),
fingerprints, label-derived lists, ids, results/labels/{dropcounts.tsv,build_report.json}; logs the
kappa rule in the validation ledger (committed) and writes the tier hashes and the primary label to
FROZEN.tsv. Re-running verifies that the rebuild reproduces the frozen hashes and overwrites nothing.

Prints counts and hashes only, never label values.

Run: /data/tmpA/nextgem/.venv/bin/python scripts/build_labels.py   (~2-4 min, ~6 GB RAM)
     [--data-root DIR] [--state-root DIR]   (toy tests; default: the project)
"""

from __future__ import annotations

import argparse
import json
import os
import sys


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", help="override NEXTGEM_DATA_ROOT (toy tests)")
    ap.add_argument("--state-root", help="override NEXTGEM_STATE_ROOT (toy tests)")
    a = ap.parse_args(argv)
    if a.data_root:
        os.environ["NEXTGEM_DATA_ROOT"] = a.data_root
    if a.state_root:
        os.environ["NEXTGEM_STATE_ROOT"] = a.state_root
    from nextgem import labels

    report = labels.build()
    summary = {k: report[k] for k in ("n_lines_labels", "n_lines_split", "lines_missing_from_labels", "n_U",
                                      "n_U_in_crispr", "lines_per_tier", "pairs_per_tier", "id_sets", "tier_sha256",
                                      "rebuild_check", "runtime_s")}
    summary["primary_label"] = report["kappa_rule"]["primary_label"]
    print(json.dumps(summary, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
