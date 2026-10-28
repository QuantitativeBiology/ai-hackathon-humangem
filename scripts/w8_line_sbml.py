"""W8: write the per-line P1 SBML of one or more lines on demand (DESIGN W8 naming
``{ModelID}_{parent}-{version}_P1.xml.gz``). Label-free.

The P1 resource stores one model per (parent, medium signature) (scripts/w8_build_p1.py); a per-line file is that
shared model (sha256-checked against resource/manifest.tsv) with the line's own metadata added to its notes
(model_id, line_medium_key, line_medium_stratum, line_medium_fallback, shared_model_file, ...). Nothing else changes.

Usage:
  .venv/bin/python scripts/w8_line_sbml.py --lines ACH-000001 [ACH-...] [--parent human2] [--out DIR]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"

from nextgem import paths, resource as rs  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lines", nargs="+", required=True)
    ap.add_argument("--parent", default="human2", choices=("human2", "human1"))
    ap.add_argument("--resource", default=str(paths.state_root() / "resource"))
    ap.add_argument("--out", default=".")
    a = ap.parse_args(argv)
    for line in a.lines:
        p = rs.materialise_line(Path(a.resource), line, a.parent, Path(a.out))
        print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
