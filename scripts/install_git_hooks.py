"""Install the nextgem git hooks (DESIGN.md §11.1, §11.3): pre-commit -> scripts/pre_commit_guard.py,
pre-push -> scripts/pre_push_guard.py (scanner on pushed blobs, check_frozen, full test suite). Both run with the project venv interpreter (so the
tripwire is active) and NEXTGEM_STATE_ROOT = the repository being committed.

  .venv/bin/python scripts/install_git_hooks.py [--repo PATH] [--force]

Refuses to overwrite an existing hook that it did not write, unless --force.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
MARK = "# nextgem guard hook"


def hook_text(script: str, args: str = "") -> str:
    return f"""#!/bin/sh
{MARK} (DESIGN.md 11.1/11.3); written by scripts/install_git_hooks.py
ROOT="$(git rev-parse --show-toplevel)" || exit 1
NEXTGEM_STATE_ROOT="$ROOT" exec "{PROJECT}/.venv/bin/python" "{PROJECT}/scripts/{script}" {args}
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", type=Path, default=PROJECT)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)
    gd = subprocess.run(["git", "-C", str(a.repo), "rev-parse", "--git-path", "hooks"], capture_output=True, text=True, check=True).stdout.strip()
    hooks = Path(gd) if Path(gd).is_absolute() else a.repo / gd
    hooks.mkdir(parents=True, exist_ok=True)
    for name, text in (("pre-commit", hook_text("pre_commit_guard.py")),
                       ("pre-push", hook_text("pre_push_guard.py"))):
        f = hooks / name
        if f.exists() and MARK not in f.read_text() and not a.force:
            print(f"{f} exists and was not written by this script; use --force", file=sys.stderr)
            return 1
        f.write_text(text)
        os.chmod(f, 0o755)
        print(f"installed {f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
