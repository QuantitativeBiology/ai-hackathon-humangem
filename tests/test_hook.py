"""§12 "Agent guards": the PreToolUse hook (.claude/hooks/guard_labels.py) denies the liveness probe and
direct label reads; it lets ordinary reads through. (Event payloads go through stdin.)"""

from __future__ import annotations

import json
import subprocess

import pytest

import toy

HOOK = toy.PROJECT / ".claude" / "hooks" / "guard_labels.py"
RAW = "data/depmap/" + toy.RAW_GE


def _run(event):
    r = subprocess.run(["python3", str(HOOK)], input=json.dumps(event), capture_output=True, text=True)
    assert r.returncode == 0
    if not r.stdout.strip():
        return "allow"
    return json.loads(r.stdout)["hookSpecificOutput"]["permissionDecision"]


@pytest.mark.parametrize("event", [
    {"tool_name": "Read", "tool_input": {"file_path": str(toy.PROJECT / "data" / "locked" / "HOOK_PROBE")}},
    {"tool_name": "Read", "tool_input": {"file_path": str(toy.PROJECT / RAW)}},
    {"tool_name": "Bash", "tool_input": {"command": f"head -c 20 {RAW}"}},
    {"tool_name": "Bash", "tool_input": {"command": "ls && cat " + str(toy.PROJECT / "data" / "labels" / "labels_train.parquet").replace("train", "test")}},
    {"tool_name": "Grep", "tool_input": {"pattern": "x", "path": "data/depmap"}},
])
def test_hook_denies(event):
    assert _run(event) == "deny"


@pytest.mark.parametrize("event", [
    {"tool_name": "Read", "tool_input": {"file_path": str(toy.PROJECT / "DESIGN.md")}},
    {"tool_name": "Bash", "tool_input": {"command": ".venv/bin/python scripts/make_splits.py --help"}},
])
def test_hook_allows(event):
    assert _run(event) == "allow"
