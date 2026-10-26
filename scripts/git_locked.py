"""Run a git command under the ledger commit lock (DESIGN.md §11.3: ledger appends commit under a lock
shared with the PI's commits). Use it for PI commits while jobs may append ledger rows:

  .venv/bin/python scripts/git_locked.py commit -m "message" -- path ...
"""

from __future__ import annotations

import subprocess
import sys


def main() -> int:
    from nextgem import ledger, paths

    with ledger.commit_lock(paths.state_root()):
        return subprocess.run(["git", "-C", str(paths.state_root()), *sys.argv[1:]]).returncode


if __name__ == "__main__":
    sys.exit(main())
