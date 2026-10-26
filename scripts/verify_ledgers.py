"""Verify the test and validation ledgers (DESIGN.md §11.3): hash chain intact; additions only across git
history; look caps; every sealed output under data/locked/sealed/<look>/ has its SHA-256 in an ``output``
row of that look; looks left open (started, not completed) are reported. Entry point permitted by the
tripwire to *hash* sealed outputs.

  .venv/bin/python scripts/verify_ledgers.py      # exit 1 on any problem
"""

from __future__ import annotations

import sys


def main() -> int:
    from nextgem import labels, ledger, paths

    problems, notes = [], []
    for name in ("test", "validation"):
        problems += ledger.verify(name)
        rows = ledger.read(name)
        for lk_id, lk in ledger.looks(rows).items():
            if lk["started"] is not None and not lk["completed"]:
                notes.append(f"{name}: look {lk_id} started but not completed (counts as a look)")
        notes.append(f"{name}: {len(rows)} rows, {len(ledger.looks(rows))} looks")
    sealed_root = paths.locked_dir() / "sealed"
    if sealed_root.exists():
        test_rows = ledger.read("test")
        outs = {}
        for r in test_rows:
            if r["event"] == "output":
                outs.setdefault(r["look_id"], set()).add(r["output_sha"])
        for look_dir in sorted(p for p in sealed_root.iterdir() if p.is_dir()):
            for f in sorted(p for p in look_dir.rglob("*") if p.is_file()):
                sha = labels.sha256_protected(f)
                if sha not in outs.get(look_dir.name, set()):
                    problems.append(f"sealed output {f.relative_to(sealed_root)} has no output row in look {look_dir.name}")
    for n in notes:
        print(n, file=sys.stderr)
    for p in problems:
        print(f"LEDGER: {p}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
