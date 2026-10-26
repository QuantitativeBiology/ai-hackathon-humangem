"""Check FROZEN.tsv (DESIGN.md §11.3). Runs as the git pre-commit hook (``--staged``) and at the start of
every evaluate script. Entry point permitted by the tripwire to *hash* label files.

  .venv/bin/python scripts/check_frozen.py [--staged] [--no-test-runs]
  .venv/bin/python scripts/check_frozen.py --add PATH [PATH ...] [--note TEXT]
        add file rows (or refresh rows that are not yet committed in HEAD:FROZEN.tsv)
  .venv/bin/python scripts/check_frozen.py --add PATH --amendment ID
        change a committed row: ID must be in the DESIGN.md §17 amendment log

Exit 1 with one line per problem; exit 0 if clean.
"""

from __future__ import annotations

import argparse
import sys


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--staged", action="store_true", help="check the staged (index) versions (pre-commit)")
    ap.add_argument("--add", nargs="+", help="add/refresh file rows (paths relative to the project)")
    ap.add_argument("--amendment", default="", help="amendment id for changing a committed row")
    ap.add_argument("--note", default="")
    ap.add_argument("--no-test-runs", action="store_true", help="skip the test-run check (record-only commits)")
    a = ap.parse_args(argv)

    from nextgem import frozen, labels, paths
    from nextgem.ledger import git

    if a.add:
        head = git(paths.state_root(), "show", "HEAD:FROZEN.tsv", check=False)
        committed = {r["item"]: r for r in frozen.parse(head.stdout)} if head.returncode == 0 else {}
        rows = []
        for rel in a.add:
            row = frozen.file_row(rel, a.note)
            if rel in committed and committed[rel]["sha256_or_value"] != row["sha256_or_value"]:
                if not a.amendment:
                    print(f"{rel} is frozen in HEAD:FROZEN.tsv; changing it needs --amendment <id> (DESIGN §17)", file=sys.stderr)
                    return 1
                row["amendment"] = a.amendment
            rows.append(row)
        frozen.upsert_rows(rows)
        print(f"FROZEN.tsv: {len(rows)} row(s) added/refreshed")

    problems = frozen.check(staged=a.staged, hash_protected=labels.sha256_protected, test_runs=not a.no_test_runs)
    for p in problems:
        print(f"FROZEN: {p}")
    if not problems:
        n = len(frozen.read_rows())
        print(f"FROZEN.tsv OK ({n} rows)", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
