"""FROZEN.tsv manifest and its checks (DESIGN.md §11.3, PI decision 18).

``FROZEN.tsv`` (state root) columns: ``item kind path sha256_or_value frozen_on amendment note``.

- ``kind=file``: ``path`` relative to the state root; ``sha256_or_value`` = SHA-256 of its content.
- ``kind=label_file``: ``path`` relative to the data root (label tier files); hashed through
  ``labels.sha256_protected`` (tripwire: only check_frozen.py / verify_ledgers.py may hash them).
- ``kind=value``: a frozen decision, e.g. ``label.primary`` (the kappa rule outcome).

Rules checked by :func:`check` (``scripts/check_frozen.py``; git pre-commit hook; start of every
evaluate script):

1. every ``file``/``label_file`` row matches the current content (``staged=True``: the index blob);
2. against ``HEAD:FROZEN.tsv``, no row may disappear, and a row whose path/hash/value changed must name
   an ``amendment`` id that appears in the DESIGN.md §17 amendment log (an amendment row committed
   before first use); new rows are free (first freeze);
3. ``results/test_runs.tsv``: a test that passed in the previous full run may not disappear or be
   downgraded in the latest full run.
"""

from __future__ import annotations

import csv
import hashlib
import io
import re
import subprocess
import time
from pathlib import Path

from nextgem import paths
from nextgem.ledger import git

FIELDS = ["item", "kind", "path", "sha256_or_value", "frozen_on", "amendment", "note"]
KINDS = {"file", "label_file", "value"}


def parse(text: str) -> list[dict]:
    if not text.strip():
        return []
    rdr = csv.DictReader(io.StringIO(text), delimiter="\t")
    if rdr.fieldnames != FIELDS:
        raise ValueError(f"FROZEN.tsv header must be {FIELDS}")
    return [{k: (r[k] or "") for k in FIELDS} for r in rdr]


def read_rows(path: Path | None = None) -> list[dict]:
    p = path or paths.frozen_tsv()
    return parse(p.read_text()) if p.exists() else []


def write_rows(rows: list[dict], path: Path | None = None) -> None:
    p = path or paths.frozen_tsv()
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in FIELDS})
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(buf.getvalue())
    tmp.replace(p)


def upsert_rows(new: list[dict]) -> None:
    """Add rows, or replace rows of the same item (a replaced *committed* row then needs an amendment,
    which :func:`check` enforces)."""
    rows = read_rows()
    by = {r["item"]: i for i, r in enumerate(rows)}
    for r in new:
        r = {**{k: "" for k in FIELDS}, **r}
        r["frozen_on"] = r["frozen_on"] or time.strftime("%Y-%m-%d")
        if r["kind"] not in KINDS:
            raise ValueError(f"bad kind {r['kind']!r}")
        if r["item"] in by:
            rows[by[r["item"]]] = r
        else:
            by[r["item"]] = len(rows)
            rows.append(r)
    write_rows(rows)


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def file_row(rel: str, note: str = "") -> dict:
    p = paths.state_root() / rel
    return {"item": rel, "kind": "file", "path": rel, "sha256_or_value": sha256_bytes(p.read_bytes()), "note": note}


def _staged_bytes(rel: str) -> bytes | None:
    """Exact bytes of the index (staged) version of ``rel``; None if it is not in the index."""
    out = subprocess.run(["git", "-C", str(paths.state_root()), "cat-file", "-p", f":{rel}"], capture_output=True)
    return out.stdout if out.returncode == 0 else None


def _git_ignored(rel: str) -> bool:
    out = subprocess.run(["git", "-C", str(paths.state_root()), "check-ignore", "-q", "--", rel], capture_output=True)
    return out.returncode == 0


def _content_sha(rel: str, staged: bool) -> str | None:
    p = paths.state_root() / rel
    if staged and not _git_ignored(rel):  # what the commit will contain: the index only
        b = _staged_bytes(rel)
        return sha256_bytes(b) if b is not None else None
    # A5: git-ignored inputs (data/, results/**/tmp/) are never staged; they are checked in the working tree.
    return sha256_bytes(p.read_bytes()) if p.exists() else None


AMENDMENT_ID = re.compile(r"\bA\d+\b")


def amendment_table(design_text: str) -> list[str]:
    """Data rows of the DESIGN.md §17 'Amendment log' table (header and separator rows excluded)."""
    m = re.search(r"###\s*Amendment log[^\n]*\n(.*?)(?=\n#{1,3}\s|\Z)", design_text, re.S)
    if not m:
        return []
    rows = [ln for ln in m.group(1).splitlines() if ln.strip().startswith("|")]
    return [r for r in rows[2:] if not re.fullmatch(r"\|[\s|:-]*\|", r.strip())]  # skip header + separator


def amendment_ids(design_text: str) -> set[str]:
    """Explicit amendment IDs (``A1``, ``A2``, ...) in the data rows of the §17 amendment table (review M5)."""
    return {i for row in amendment_table(design_text) for i in AMENDMENT_ID.findall(row)}


def _design_text(staged: bool) -> str:
    if staged:
        b = _staged_bytes("DESIGN.md")
        if b is not None:
            return b.decode()
    p = paths.design_md()
    return p.read_text() if p.exists() else ""


def check(staged: bool = False, hash_protected=None, test_runs: bool = True) -> list[str]:
    """Return a list of problems (empty = OK). ``hash_protected(path)`` hashes label files.
    ``test_runs=False`` skips the test-run check (record-only commits, review S1)."""
    root = paths.state_root()
    problems: list[str] = []
    if staged:  # the manifest being committed (none in the index yet -> nothing frozen to check)
        b = _staged_bytes("FROZEN.tsv")
        text = b.decode() if b is not None else ""
    else:
        text = paths.frozen_tsv().read_text() if paths.frozen_tsv().exists() else ""
    try:
        rows = parse(text)
    except ValueError as e:
        return [str(e)]
    seen = set()
    for r in rows:
        if r["item"] in seen:
            problems.append(f"duplicate item {r['item']}")
        seen.add(r["item"])
        if r["kind"] not in KINDS:
            problems.append(f"{r['item']}: bad kind {r['kind']!r}")
        elif r["kind"] == "file":
            sha = _content_sha(r["path"], staged)
            if sha is None:
                problems.append(f"{r['item']}: frozen file {r['path']} is missing")
            elif sha != r["sha256_or_value"]:
                problems.append(f"{r['item']}: {r['path']} differs from its frozen SHA-256 (amend via DESIGN §17 + FROZEN.tsv)")
        elif r["kind"] == "label_file":
            p = paths.data_root() / r["path"]
            if not p.exists():
                problems.append(f"{r['item']}: frozen label file {r['path']} is missing")
            elif hash_protected is None:
                problems.append(f"{r['item']}: cannot hash label files outside check_frozen.py")
            elif hash_protected(p) != r["sha256_or_value"]:
                problems.append(f"{r['item']}: label file {r['path']} differs from its frozen SHA-256")

    head = git(root, "show", "HEAD:FROZEN.tsv", check=False)
    if head.returncode == 0:
        try:
            head_rows = {r["item"]: r for r in parse(head.stdout)}
        except ValueError:
            head_rows = {}
        cur = {r["item"]: r for r in rows}
        ids = amendment_ids(_design_text(staged))
        for item, hr in head_rows.items():
            if item not in cur:
                problems.append(f"{item}: row removed from FROZEN.tsv (frozen rows are never removed)")
                continue
            cr = cur[item]
            if (cr["path"], cr["sha256_or_value"], cr["kind"]) != (hr["path"], hr["sha256_or_value"], hr["kind"]):
                if not cr["amendment"]:
                    problems.append(f"{item}: frozen row changed without an amendment id")
                elif cr["amendment"] not in ids:
                    problems.append(f"{item}: amendment {cr['amendment']!r} is not in the DESIGN.md §17 amendment log")
    if test_runs:
        problems += check_test_runs()
    return problems


# ----------------------------------------------------------------------------- test-run ledger
TEST_RUN_FIELDS = ["timestamp", "commit", "full_run", "n_passed", "n_failed", "n_xfailed", "n_skipped", "n_error",
                   "testids_sha256", "outcomes_file"]


def record_test_run(outcomes: dict[str, str], full_run: bool, commit: str) -> Path:
    """Append a row to results/test_runs.tsv and write the per-test outcomes file it points to."""
    lines = "".join(f"{k}\t{v}\n" for k, v in sorted(outcomes.items()))
    sha_ids = sha256_bytes("".join(f"{k}\n" for k in sorted(outcomes)).encode())
    sha_out = sha256_bytes(lines.encode())
    d = paths.results_dir() / "test_runs"
    d.mkdir(parents=True, exist_ok=True)
    of = d / f"{sha_out[:16]}.tsv"
    of.write_text("nodeid\toutcome\n" + lines)
    counts = {k: sum(1 for v in outcomes.values() if v == k) for k in ("passed", "failed", "xfailed", "skipped", "error")}
    f = paths.test_runs_tsv()
    new = not f.exists()
    with open(f, "a") as fh:
        if new:
            fh.write("\t".join(TEST_RUN_FIELDS) + "\n")
        fh.write("\t".join([time.strftime("%Y-%m-%dT%H:%M:%S"), commit, str(full_run), *(str(counts[k]) for k in
                                                                                          ("passed", "failed", "xfailed", "skipped", "error")),
                            sha_ids, str(of.relative_to(paths.state_root()))]) + "\n")
    return of


def _outcomes(rel: str) -> dict[str, str]:
    p = paths.state_root() / rel
    if not p.exists():
        return {}
    rows = [ln.split("\t") for ln in p.read_text().splitlines()[1:] if ln]
    return {r[0]: r[1] for r in rows}


def retired_tests() -> set[str]:
    """Tests retired by an amendment: results/test_retirements.tsv rows (nodeid, amendment id in §17)."""
    f = paths.results_dir() / "test_retirements.tsv"
    if not f.exists():
        return set()
    ids = amendment_ids(_design_text(False))
    out = set()
    for ln in f.read_text().splitlines()[1:]:
        parts = ln.split("\t")
        if len(parts) >= 2 and parts[1].strip() in ids:
            out.add(parts[0])
    return out


def check_test_runs() -> list[str]:
    """Every test that passed in *any* earlier full run must pass in the latest full run, unless it was
    retired by a logged amendment (review S1: a downgrade is not forgotten after one more run)."""
    f = paths.test_runs_tsv()
    if not f.exists():
        return []
    rows = list(csv.DictReader(io.StringIO(f.read_text()), delimiter="\t"))
    full = [r for r in rows if r.get("full_run") == "True"]
    if len(full) < 2:
        return []
    baseline: set[str] = set()
    for r in full[:-1]:
        baseline |= {t for t, o in _outcomes(r["outcomes_file"]).items() if o == "passed"}
    last = _outcomes(full[-1]["outcomes_file"])
    problems = []
    for t in sorted(baseline - retired_tests()):
        if last.get(t) != "passed":
            problems.append(f"test {t} passed in an earlier full run but is {last.get(t, 'missing')} now "
                            "(known-answer tests may not disappear or be downgraded)")
    return problems
