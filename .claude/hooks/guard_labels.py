#!/usr/bin/env python3
"""PreToolUse guard: block direct agent access to raw DepMap label files.

Labels must be read only through src/nextgem/labels.py (run inside project scripts),
never by an agent's own Read/Grep/Glob/Bash calls. See DESIGN.md (leakage controls) and
notes/design_review_r1/PI_decisions.md item 15.

Checks: Read/NotebookRead file paths, Grep/Glob path/glob/pattern, Bash command text.
Edit/Write are not checked (writing about the files is not reading them).
Bash commands that are git operations or the checksum script are allowed.
"""

import json
import re
import sys

# Raw label sources and the sealed test-label file.
PROTECTED = re.compile(
    r"CRISPRGeneEffect|CRISPRGeneDependency|ScreenGeneEffect|"
    r"CRISPRInferredCommonEssentials|AchillesCommonEssentialControls|"
    r"AchillesNonessentialControls|labels_test|labels_val|data/locked"
)
# Recursive searches over data/ reach the label files without naming them.
RECURSIVE_DATA_SEARCH = re.compile(r"\b(grep|rg|ag|ack|find)\b[^;&|\n]*\s(\S*/)?data(/depmap)?/?(\s|$)")
# Wildcards over the DepMap folder could reach the label files without naming them.
DEPMAP_WILDCARD = re.compile(r"depmap/[^\s'\"]*[*?\[]")
DEPMAP_DIR = re.compile(r"data/depmap/?$")
BASH_ALLOWED = re.compile(r"^\s*(git\s|\S*python\S*\s+scripts/checksum_data\.py)")


def deny(reason):
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": (
                f"Blocked by label guard: {reason}. Raw DepMap label files may only be read via "
                "src/nextgem/labels.py inside a project script (DESIGN.md leakage controls)."
            ),
        }
    }))
    sys.exit(0)


def main():
    event = json.load(sys.stdin)
    tool = event.get("tool_name", "")
    args = event.get("tool_input", {}) or {}

    if tool in ("Read", "NotebookRead"):
        path = str(args.get("file_path", "")) + str(args.get("notebook_path", ""))
        if PROTECTED.search(path):
            deny(f"direct read of {path}")
    elif tool in ("Grep", "Glob"):
        fields = " ".join(str(args.get(k, "")) for k in ("path", "glob", "pattern"))
        if PROTECTED.search(fields) or DEPMAP_WILDCARD.search(fields) or DEPMAP_DIR.search(str(args.get("path", ""))):
            deny(f"{tool} over raw label files ({fields.strip()})")
    elif tool == "Bash":
        # Check each chained segment, so an allowed prefix can't carry a read along.
        for segment in re.split(r"&&|\|\||[;|\n]", str(args.get("command", ""))):
            if BASH_ALLOWED.match(segment):
                if not segment.lstrip().startswith("git"):
                    continue
                # git: names may appear in commit messages, so only check unquoted text.
                segment = re.sub(r"\"[^\"]*\"|'[^']*'", "", segment)
            if PROTECTED.search(segment) or DEPMAP_WILDCARD.search(segment) or RECURSIVE_DATA_SEARCH.search(segment):
                deny("shell command references raw label files")


if __name__ == "__main__":
    # Fail closed: a guard that crashes must block, not let the call through.
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:
        deny(f"guard error ({type(exc).__name__}); refusing by default")
