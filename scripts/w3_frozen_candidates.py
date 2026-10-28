"""Candidate FROZEN.tsv rows for the W3 frozen lists (DESIGN W3: B-loo components V-13, B-add metabolite IDs V-14,
pseudo-task definitions; plus the variant specs the sweeps used). The PI adds them; this script edits nothing.

Output: results/w3/frozen_candidates.tsv (FROZEN.tsv columns: item, kind, path, sha256_or_value, frozen_on,
amendment, note). Label-free.
Usage: .venv/bin/python scripts/w3_frozen_candidates.py
"""

from __future__ import annotations

import hashlib
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w3_common import OUT_DIR, PROJECT, atomic_write_tsv  # noqa: E402

FILES = [
    ("results/w3/structure/bloo_components.tsv", "W3 B-loo components per parent (V-13); incl. metpool:dihydroorotate "
     "added from the label-free M3_W5 review finding before W3 scoring"),
    ("results/w3/structure/badd_metabolites.tsv", "W3 B-add metabolite IDs per parent (V-14); skips logged "
     "(heme A absent in both; Fe-S cluster absent in Human1)"),
    ("results/w3/structure/variants_human1.json", "W3 variant specs used by the KO sweeps (Human1)"),
    ("results/w3/structure/variants_human2.json", "W3 variant specs used by the KO sweeps (Human2)"),
    ("results/w3/structure/demand_feasibility_log.tsv", "W3 D-BS / D-ER demand feasibility log (all kept)"),
    ("results/w3/tasks/pseudo_human1_tasks.tsv", "W3 pseudo-task definitions (swaps) for SEEDS[0:5], Human1 = G*"),
]


def main() -> int:
    day = time.strftime("%Y-%m-%d")
    rows = []
    for rel, note in FILES:
        p = PROJECT / rel
        if not p.exists():
            continue
        rows.append({"item": rel, "kind": "file", "path": rel,
                     "sha256_or_value": hashlib.sha256(p.read_bytes()).hexdigest(), "frozen_on": day, "amendment": "",
                     "note": note})
    sc = OUT_DIR / "structure" / "frozen_candidates_structure.tsv"
    if sc.exists():
        rows += pd.read_csv(sc, sep="\t", keep_default_na=False).to_dict("records")
    atomic_write_tsv(pd.DataFrame(rows), OUT_DIR / "frozen_candidates.tsv", float_format=None)
    print(pd.DataFrame(rows)[["item", "kind", "sha256_or_value"]].to_string(index=False, max_colwidth=70))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
