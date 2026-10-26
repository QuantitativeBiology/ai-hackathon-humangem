"""Content-based label-copy scanner (DESIGN.md §11.1, PI decision r2-14). Entry point permitted by the
tripwire (when frozen) to read the label fingerprints (``data/locked/fingerprints.npz``) - nothing else.

Default scope: the project tree (results/, notes/, data/ outside data/locked and outside the registered
raw downloads in data/checksums.tsv, ...) plus the session scratchpads ``/tmp/claude-*``.
Hit rules: see ``src/nextgem/scan.py``. Output is coarse: ``HIT <path> <reason class>`` only - never
values or match counts (review M6). Exit 1 if any hit.

Usage:
  .venv/bin/python scripts/scan_label_copies.py                  # milestone scan
  .venv/bin/python scripts/scan_label_copies.py --staged         # pre-commit (staged blobs; + test-keyed rule)
  .venv/bin/python scripts/scan_label_copies.py --tracked        # all tracked files (working tree)
  .venv/bin/python scripts/scan_label_copies.py --push-range     # pre-push: blobs of the pushed commits (stdin)
  .venv/bin/python scripts/scan_label_copies.py --roots DIR ...  # explicit roots / files
  [--delete-scratch-hits]  delete hits under /tmp/claude-* and log the deletion (logs/label_copy_deletions.tsv)
"""

from __future__ import annotations

import argparse
import glob
import subprocess
import sys
import time
from pathlib import Path

ZERO = "0" * 40


def _git(root: Path, *args, binary=False):
    r = subprocess.run(["git", "-C", str(root), *args], capture_output=True, check=False)
    return r.stdout if binary else r.stdout.decode(errors="replace")


def push_blobs(root: Path, spec: str) -> list[tuple[str, bytes]]:
    """(path, blob bytes) for every file added/modified in the commits being pushed (pre-push stdin)."""
    seen, out = set(), []
    for line in spec.splitlines():
        parts = line.split()
        if len(parts) != 4 or parts[1] == ZERO:
            continue
        local, remote = parts[1], parts[3]
        rng = [local, "--not", "--remotes"] if remote == ZERO else [f"{remote}..{local}"]
        for c in _git(root, "rev-list", *rng).split():
            raw = _git(root, "diff-tree", "-r", "--no-commit-id", "--root", "--diff-filter=ACMR", "-z", c)
            f = raw.split("\0")
            for meta, path in zip(f[0::2], f[1::2]):
                blob = meta.split()[3] if len(meta.split()) >= 4 else None
                if blob and (path, blob) not in seen:
                    seen.add((path, blob))
                    out.append((path, _git(root, "cat-file", "-p", blob, binary=True)))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--roots", nargs="*", type=Path, help="directories or files to scan (default: project + scratch)")
    ap.add_argument("--staged", action="store_true", help="scan staged files (pre-commit)")
    ap.add_argument("--tracked", action="store_true", help="scan all git-tracked files")
    ap.add_argument("--push-range", action="store_true", help="scan blobs of the pushed commits (pre-push stdin)")
    ap.add_argument("--test-keyed", action="store_true", help="also flag tables keyed by test ModelIDs x genes")
    ap.add_argument("--delete-scratch-hits", action="store_true")
    ap.add_argument("--scratch-glob", default="/tmp/claude-*")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    t0 = time.time()

    import pandas as pd

    from nextgem import labels, paths, scan
    from nextgem.ledger import read as read_ledger

    root = paths.state_root()
    ctx = scan.Context(fp=labels.load_fingerprints())
    split_f = paths.splits_dir() / "lineage_split.tsv"
    if split_f.exists():
        sp = pd.read_csv(split_f, sep="\t", usecols=["ModelID", "tier"])
        ctx.heldout_ids = set(sp.ModelID[sp.tier.isin(["validation", "test"])])
        ctx.test_ids = set(sp.ModelID[sp.tier == "test"])
        ctx.known_ids = set(sp.ModelID)
    elif labels.ids_file().exists():
        ctx.known_ids = set(labels.ids()["model_ids"])
    unsealed = any(r["event"] == "unseal" for r in read_ledger("test")) if paths.ledger_path("test").exists() else False
    if not ((a.test_keyed or a.staged or a.push_range) and not unsealed):
        ctx.test_ids = None
    ctx.container_rule = a.staged or a.tracked or a.push_range
    ctx.canonical_split = paths.splits_dir() / "lineage_split.tsv"
    ctx.canonical_split_sha = scan.canonical_split_pin(paths.splits_dir())
    al = paths.results_dir() / "scanner_allowlist.tsv"
    if al.exists():
        for ln in al.read_text().splitlines()[1:]:
            p = ln.split("\t")
            if len(p) >= 2:
                ctx.allowlist[p[0]] = p[1]

    registered = set()
    ck = paths.PROJECT_ROOT / "data" / "checksums.tsv"
    if ck.exists():
        for ln in ck.read_text().splitlines()[1:]:
            parts = ln.split("\t")
            if len(parts) == 2:
                registered.add((paths.PROJECT_ROOT / "data" / parts[1]).resolve())

    scratch_roots = [Path(p) for p in glob.glob(a.scratch_glob)]
    jobs: list[tuple[Path, bytes | None, str | None]] = []
    if a.push_range:
        jobs = [(root / p, b, p) for p, b in push_blobs(root, sys.stdin.read())]
    elif a.staged or a.tracked:
        if a.staged:
            names = _git(root, "diff", "--cached", "--name-only", "--diff-filter=ACMR", "-z").split("\0")
        else:
            names = _git(root, "ls-files", "-z").split("\0")
        for n in filter(None, names):
            data = _git(root, "cat-file", "-p", f":{n}", binary=True) if a.staged else None
            jobs.append((root / n, data, n))
    else:
        roots = a.roots if a.roots else [root, *scratch_roots]
        exclude = [paths.locked_dir(), paths.PROJECT_ROOT / "data" / "locked", paths.PROJECT_ROOT / ".venv",
                   paths.PROJECT_ROOT / "papers"]
        jobs = [(p, None, None) for p in scan.iter_files(roots, exclude, registered)]

    hits, n_files = [], 0
    for p, data, rel in jobs:
        n_files += 1
        try:
            res = scan.scan_file(p, ctx, data=data, rel=rel)
        except Exception as e:  # never let one odd file stop the scan
            res = scan.Result(str(p), note=f"error {type(e).__name__}")
        if res.hits:
            hits.append(res)
            print(f"HIT\t{p}\t{'; '.join(res.hits)}")
        elif a.verbose:
            print(f"ok\t{p}\t{res.note}")

    if a.delete_scratch_hits:
        log = paths.logs_dir() / "label_copy_deletions.tsv"
        log.parent.mkdir(parents=True, exist_ok=True)
        for res in hits:
            p = Path(res.path).resolve()
            if any(str(p).startswith(str(s.resolve()) + "/") for s in scratch_roots):
                p.unlink(missing_ok=True)
                with open(log, "a") as fh:
                    fh.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')}\t{p}\t{'; '.join(res.hits)}\n")
                print(f"DELETED\t{p}\t(add this deletion to notes/label_exposure.md)")

    print(f"scanned {n_files} files, {len(hits)} hit(s); fingerprints {'loaded' if ctx.fp else 'not built yet'}; "
          f"test-keyed rule {'on' if ctx.test_ids else 'off'}; {time.time() - t0:.1f} s", file=sys.stderr)
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
