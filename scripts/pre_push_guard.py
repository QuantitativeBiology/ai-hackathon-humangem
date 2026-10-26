"""Git pre-push guard (DESIGN.md §11.1 "before every push", §11.3). Installed as .git/hooks/pre-push.

1. The content scanner on the blobs of every commit being pushed (not the working tree, so a label
   copy committed and later deleted is still caught; review S4).
2. scripts/check_frozen.py (including the test-run downgrade check).
3. The full test suite (``pytest tests``), recorded in results/test_runs.tsv. The pre-commit hook stays
   fast and never runs it.

Exit 1 blocks the push.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent


def main() -> int:
    from nextgem import paths

    spec = sys.stdin.read()
    problems = []
    r = subprocess.run([sys.executable, str(SCRIPTS / "scan_label_copies.py"), "--push-range"], input=spec,
                       capture_output=True, text=True)
    if r.returncode != 0:
        problems += [f"scanner: {ln}" for ln in r.stdout.splitlines() if ln.strip()] or [f"scanner failed: {r.stderr[-400:]}"]
    r = subprocess.run([sys.executable, str(SCRIPTS / "check_frozen.py")], capture_output=True, text=True)
    if r.returncode != 0:
        problems += r.stdout.splitlines() or [f"check_frozen failed: {r.stderr[-400:]}"]
    if not os.environ.get("NEXTGEM_PREPUSH_SKIP_TESTS"):  # used only by the hook's own toy test
        root = paths.state_root()
        r = subprocess.run([sys.executable, "-m", "pytest", "tests", "-q", "-p", "no:cacheprovider"],
                           cwd=root, capture_output=True, text=True)
        if r.returncode != 0:
            problems.append("test suite failed: " + (r.stdout.strip().splitlines() or ["?"])[-1])
    for p in problems:
        print(f"pre-push: {p}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
