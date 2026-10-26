"""Git pre-commit guard (DESIGN.md §11.1 "Pre-commit hook", §11.3). Installed as .git/hooks/pre-commit by
scripts/install_git_hooks.py. Rejects the commit (exit 1) if any of:

1. a staged file is under data/ (except data/README.md, data/checksums.tsv) or papers/;
2. a staged file is larger than 50 MB;
3. the content scanner flags a staged file (fingerprints, shape, protected names) or, until the T2
   unseal ledger row exists, a staged table is keyed by test ModelIDs x genes
   (scripts/scan_label_copies.py --staged);
4. scripts/check_frozen.py --staged fails (frozen files, amendments; test-run downgrades only for commits
   touching src/, scripts/, tests/ or FROZEN.tsv);
5. the latest recorded test run has failures and this commit edits a failing test file together with
   src/ or scripts/ ("a failing test is never edited in the same commit as the code fix").

It is fast (two short subprocesses) and never runs the test suite; the full suite runs in pre-push
(scripts/pre_push_guard.py).
"""

from __future__ import annotations

import csv
import io
import subprocess
import sys
from pathlib import Path

MAX_BYTES = 50 * 1024 * 1024
ALLOWED_DATA = {"data/README.md", "data/checksums.tsv"}
SCRIPTS = Path(__file__).resolve().parent


def main() -> int:
    from nextgem import paths

    root = paths.state_root()
    names = [n for n in subprocess.run(["git", "-C", str(root), "diff", "--cached", "--name-only", "--diff-filter=ACMR", "-z"],
                                       capture_output=True, check=True).stdout.decode().split("\0") if n]
    problems = []
    for n in names:
        if (n.startswith("data/") and n not in ALLOWED_DATA) or n.startswith("papers/"):
            problems.append(f"{n}: files under data/ or papers/ are never committed")
        size = subprocess.run(["git", "-C", str(root), "cat-file", "-s", f":{n}"], capture_output=True, text=True)
        if size.returncode == 0 and int(size.stdout.strip() or 0) > MAX_BYTES:
            problems.append(f"{n}: larger than 50 MB")

    if names:
        r = subprocess.run([sys.executable, str(SCRIPTS / "scan_label_copies.py"), "--staged"], capture_output=True, text=True)
        if r.returncode != 0:
            problems += [f"scanner: {ln}" for ln in r.stdout.splitlines() if ln.strip()] or [f"scanner failed: {r.stderr.strip()[-400:]}"]
    # the test-run check applies only to commits that touch code, tests or the manifest (review S1);
    # record-only commits (notes, results, ledger rows) are never blocked by a failing test
    code_commit = any(n.startswith(("src/", "scripts/", "tests/")) or n == "FROZEN.tsv" for n in names)
    r = subprocess.run([sys.executable, str(SCRIPTS / "check_frozen.py"), "--staged",
                        *([] if code_commit else ["--no-test-runs"])], capture_output=True, text=True)
    if r.returncode != 0:
        problems += [ln for ln in r.stdout.splitlines() if ln.strip()] or [f"check_frozen failed: {r.stderr.strip()[-400:]}"]

    tr = paths.test_runs_tsv()
    if tr.exists():
        rows = list(csv.DictReader(io.StringIO(tr.read_text()), delimiter="\t"))
        if rows and int(rows[-1].get("n_failed") or 0) > 0:
            of = root / rows[-1]["outcomes_file"]
            failing_files = set()
            if of.exists():
                for ln in of.read_text().splitlines()[1:]:
                    nodeid, _, outcome = ln.partition("\t")
                    if outcome in ("failed", "error"):
                        failing_files.add(nodeid.split("::")[0])
            touched_tests = failing_files & set(names)
            touched_code = [n for n in names if n.startswith(("src/", "scripts/"))]
            if touched_tests and touched_code:
                problems.append(f"failing test file(s) {sorted(touched_tests)} edited in the same commit as code "
                                f"{touched_code[:5]} (DESIGN §11.3: never in the same commit)")

    for p in problems:
        print(f"pre-commit: {p}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
