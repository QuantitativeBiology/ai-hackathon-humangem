"""Hash-chained, file-locked, committed-on-append ledgers with look caps (DESIGN.md §11.3, §10.2, §4.1).

``results/test_ledger.tsv`` and ``results/validation_ledger.tsv`` (under the state root). Each row holds
the SHA-256 of the previous row (``prev_sha``) and of itself (``row_sha``). :func:`append` is the only
non-PI commit path: it appends under ``.git/nextgem_commit.lock`` (shared with ``scripts/git_locked.py``,
which the PI uses for his own commits) and runs ``git commit --only <ledger>``. If the commit fails,
the row is rolled back and :class:`LedgerError` is raised, so a look whose ``started`` row is not
committed cannot proceed (AI r2 S3).

Events
------
``started`` / ``output`` / ``completed``  a look (validation: any workstream; test: T1, T2, ...)
``non_look``                         e.g. the kappa rule (§5.1)
``split_started`` / ``split_completed``  the one split draw (test ledger; never twice)
``unseal``                           the T2 unseal (test ledger; requires a completed T2 look)

Caps (in code): validation <= 40 looks (raised only by a *committed*
``results/authorisations/val_cap_raise_<n>.md``, n = the new cap), of which <= 20 for W7. A look id
that already exists is refused; a crashed look (``started`` without ``completed``) resumes only with
the same config hash and still counts as a look.
"""

from __future__ import annotations

import contextlib
import fcntl
import hashlib
import json
import os
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

from nextgem import config, paths

COLUMNS = ["seq", "timestamp", "ledger", "event", "look_id", "workstream", "config_sha", "output_sha",
           "detail", "author", "head_commit", "prev_sha", "row_sha"]
EVENTS = {"started", "output", "completed", "non_look", "split_started", "split_completed", "unseal"}
LOCK_NAME = "nextgem_commit.lock"


class LedgerError(RuntimeError):
    pass


# ----------------------------------------------------------------------------- git helpers
def git(root: Path, *args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True, check=check)


def git_dir(root: Path) -> Path:
    out = git(root, "rev-parse", "--git-dir").stdout.strip()
    p = Path(out)
    return p if p.is_absolute() else (root / p)


def head_commit(root: Path) -> str:
    r = git(root, "rev-parse", "HEAD", check=False)
    return r.stdout.strip() if r.returncode == 0 else ""


def is_committed_clean(root: Path, rel: str) -> bool:
    """File is tracked in HEAD and unmodified in the working tree."""
    if git(root, "cat-file", "-e", f"HEAD:{rel}", check=False).returncode != 0:
        return False
    return git(root, "diff", "--quiet", "HEAD", "--", rel, check=False).returncode == 0


@contextlib.contextmanager
def commit_lock(root: Path | None = None, timeout_s: float = 120.0):
    root = Path(root or paths.state_root())
    lock = git_dir(root) / LOCK_NAME
    fd = os.open(lock, os.O_CREAT | os.O_RDWR, 0o644)
    t0 = time.time()
    try:
        while True:
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.time() - t0 > timeout_s:
                    raise LedgerError(f"could not acquire {lock} within {timeout_s} s")
                time.sleep(0.2)
        yield
    finally:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)


# ----------------------------------------------------------------------------- rows and chain
def genesis(name: str) -> str:
    return hashlib.sha256(f"nextgem {name} ledger genesis".encode()).hexdigest()


def row_hash(row: dict) -> str:
    return hashlib.sha256("\t".join(str(row[c]) for c in COLUMNS[:-1]).encode()).hexdigest()


def _clean(s) -> str:
    return re.sub(r"[\t\r\n]+", " ", str(s if s is not None else ""))


def parse(text: str) -> list[dict]:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    if lines[0].split("\t") != COLUMNS:
        raise LedgerError("ledger header does not match the expected columns")
    rows = []
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) != len(COLUMNS):
            raise LedgerError(f"malformed ledger row: {ln[:80]!r}")
        rows.append(dict(zip(COLUMNS, f)))
    return rows


def read(name: str) -> list[dict]:
    p = paths.ledger_path(name)
    return parse(p.read_text()) if p.exists() else []


def verify_chain(name: str, rows: list[dict]) -> list[str]:
    problems, prev = [], genesis(name)
    for i, r in enumerate(rows):
        if r["seq"] != str(i + 1):
            problems.append(f"{name}: row {i + 1} has seq {r['seq']}")
        if r["ledger"] != name:
            problems.append(f"{name}: row {i + 1} names ledger {r['ledger']!r}")
        if r["prev_sha"] != prev:
            problems.append(f"{name}: row {i + 1} prev_sha does not match row {i}")
        if row_hash(r) != r["row_sha"]:
            problems.append(f"{name}: row {i + 1} row_sha does not match its content")
        prev = r["row_sha"]
    return problems


# ----------------------------------------------------------------------------- look bookkeeping
def looks(rows: list[dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for r in rows:
        if r["event"] in ("started", "output", "completed"):
            lk = out.setdefault(r["look_id"], {"started": None, "outputs": [], "completed": False, "workstream": "", "config_sha": ""})
            if r["event"] == "started":
                lk["started"] = r
                lk["workstream"] = r["workstream"]
                lk["config_sha"] = r["config_sha"]
            elif r["event"] == "output":
                lk["outputs"].append(r["output_sha"])
            else:
                lk["completed"] = True
    return out


def validation_cap(root: Path | None = None) -> int:
    root = Path(root or paths.state_root())
    cap = config.VALIDATION_LOOK_CAP
    d = paths.authorisations_dir()
    if d.exists():
        for f in d.glob("val_cap_raise_*.md"):
            m = re.fullmatch(r"val_cap_raise_(\d+)\.md", f.name)
            if m and is_committed_clean(root, str(f.relative_to(root))):
                cap = max(cap, int(m.group(1)))
    return cap


def _validate(name: str, event: str, look_id: str, workstream: str, config_sha: str, output_sha: str,
              rows: list[dict], root: Path) -> dict | None:
    """Raise LedgerError if the event is not allowed. Returns an existing row to resume, or None."""
    if event not in EVENTS:
        raise LedgerError(f"unknown event {event!r}")
    lk = looks(rows)
    if event == "started":
        if not look_id or not config_sha:
            raise LedgerError("a look needs look_id and config_sha")
        if name == "validation" and not workstream:
            raise LedgerError("a validation look needs its workstream")
        if look_id in lk:
            prev = lk[look_id]
            if prev["completed"]:
                raise LedgerError(f"look {look_id!r} already exists and is completed; looks are never repeated")
            if prev["config_sha"] != config_sha:
                raise LedgerError(f"look {look_id!r} crashed earlier with config {prev['config_sha'][:12]}; "
                                  "it resumes only with the same config hash")
            return prev["started"]  # resume: still one look
        if name == "validation":
            n = len(lk)
            cap = validation_cap(root)
            if n >= cap:
                raise LedgerError(f"validation look cap reached ({n}/{cap}); raising it needs a committed "
                                  f"results/authorisations/val_cap_raise_<n>.md")
            if workstream == "W7" and sum(1 for v in lk.values() if v["workstream"] == "W7") >= config.VALIDATION_LOOK_CAP_W7:
                raise LedgerError(f"W7 validation look cap reached ({config.VALIDATION_LOOK_CAP_W7})")
    elif event in ("output", "completed"):
        if look_id not in lk or lk[look_id]["started"] is None:
            raise LedgerError(f"look {look_id!r} has no started row")
        if lk[look_id]["completed"]:
            raise LedgerError(f"look {look_id!r} is already completed")
        if event == "output" and not re.fullmatch(r"[0-9a-f]{64}", output_sha or ""):
            raise LedgerError("an output row needs a SHA-256 output_sha")
    elif event in ("split_started", "split_completed", "unseal"):
        if name != "test":
            raise LedgerError(f"{event} rows belong to the test ledger")
        evs = [r["event"] for r in rows]
        if event == "split_started" and "split_started" in evs:
            raise LedgerError("the split was already drawn (split_started exists); it is never redrawn")
        if event == "split_completed" and ("split_started" not in evs or "split_completed" in evs):
            raise LedgerError("split_completed needs exactly one prior split_started")
        if event == "unseal" and not any(k.startswith("T2") and v["completed"] for k, v in lk.items()):
            raise LedgerError("unseal needs a completed T2 look")
    return None


# ----------------------------------------------------------------------------- append
def append(name: str, event: str, *, look_id: str = "", workstream: str = "", config_sha: str = "",
           output_sha: str = "", detail: str = "", author: str = "", commit: bool = True) -> dict:
    """Append one row and commit it (``git commit --only``). Returns the row (or the resumed one)."""
    root = paths.state_root()
    path = paths.ledger_path(name)
    rel = str(path.relative_to(root))
    with commit_lock(root) if commit else contextlib.nullcontext():
        old = path.read_bytes() if path.exists() else None
        rows = parse(old.decode()) if old else []
        problems = verify_chain(name, rows)
        if problems:
            raise LedgerError("ledger chain broken; refusing to append: " + "; ".join(problems))
        resumed = _validate(name, event, look_id, workstream, config_sha, output_sha, rows, root)
        if resumed is not None:
            return resumed
        row = {"seq": str(len(rows) + 1), "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "ledger": name,
               "event": event, "look_id": _clean(look_id), "workstream": _clean(workstream),
               "config_sha": _clean(config_sha), "output_sha": _clean(output_sha), "detail": _clean(detail),
               "author": _clean(author), "head_commit": head_commit(root) if commit else "",
               "prev_sha": rows[-1]["row_sha"] if rows else genesis(name)}
        row["row_sha"] = row_hash(row)
        text = (old.decode() if old else "\t".join(COLUMNS) + "\n") + "\t".join(row[c] for c in COLUMNS) + "\n"
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".tmp")
        tmp.write_text(text)
        os.replace(tmp, path)
        if commit:
            try:
                git(root, "add", "--", rel)
                git(root, "commit", "--only", "-q", "-m", f"ledger({name}): {event} {look_id}".strip(), "--", rel)
            except subprocess.CalledProcessError as e:
                if old is None:
                    path.unlink(missing_ok=True)
                    git(root, "rm", "--cached", "-q", "--ignore-unmatch", "--", rel, check=False)
                else:
                    path.write_bytes(old)
                    git(root, "add", "--", rel, check=False)
                raise LedgerError(f"ledger commit failed; row rolled back: {(e.stderr or e.stdout or '').strip()[:400]}") from e
    return row


@dataclass
class Look:
    name: str
    look_id: str

    def output(self, sha256: str, detail: str = "", author: str = "") -> dict:
        return append(self.name, "output", look_id=self.look_id, output_sha=sha256, detail=detail, author=author)


@contextlib.contextmanager
def look(name: str, look_id: str, *, workstream: str, config_sha: str, detail: str = "", author: str = ""):
    """``started`` is committed before the body runs; ``completed`` only if the body succeeds."""
    append(name, "started", look_id=look_id, workstream=workstream, config_sha=config_sha, detail=detail, author=author)
    lk = Look(name, look_id)
    yield lk
    append(name, "completed", look_id=look_id, author=author)


# ----------------------------------------------------------------------------- verification
def verify_history(name: str) -> list[str]:
    """Every committed version of the ledger is a prefix of the next, and of the working tree."""
    root = paths.state_root()
    path = paths.ledger_path(name)
    if not path.exists():
        return []
    rel = str(path.relative_to(root))
    commits = git(root, "log", "--format=%H", "--reverse", "--", rel, check=False).stdout.split()
    problems, prev_text, prev_c = [], "", None
    for c in commits:
        r = git(root, "show", f"{c}:{rel}", check=False)
        text = r.stdout if r.returncode == 0 else ""
        if not text.startswith(prev_text):
            problems.append(f"{name}: commit {c[:10]} rewrites rows committed in {prev_c[:10] if prev_c else '?'} (additions only)")
        prev_text, prev_c = text, c
    if not path.read_text().startswith(prev_text):
        problems.append(f"{name}: working tree rewrites committed rows (additions only)")
    return problems


def verify(name: str) -> list[str]:
    rows = read(name)
    problems = verify_chain(name, rows) + verify_history(name)
    lk = looks(rows)
    if name == "validation":
        cap = validation_cap()
        if len(lk) > cap:
            problems.append(f"validation looks {len(lk)} exceed the cap {cap}")
        w7 = sum(1 for v in lk.values() if v["workstream"] == "W7")
        if w7 > config.VALIDATION_LOOK_CAP_W7:
            problems.append(f"W7 validation looks {w7} exceed {config.VALIDATION_LOOK_CAP_W7}")
    return problems


def detail_json(row: dict) -> dict:
    try:
        return json.loads(row["detail"])
    except (json.JSONDecodeError, KeyError):
        return {}
