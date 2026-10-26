"""Install the label tripwire into the project venv (DESIGN.md §11.1, PI decision r2-15).

Writes ``<venv>/lib/python3.11/site-packages/sitecustomize.py``, which loads
``src/nextgem/tripwire.py`` as the standalone module ``nextgem_tripwire`` at every interpreter launch.
If loading fails, the sitecustomize installs a minimal fail-closed audit hook that denies *every* open
of a protected label path (and says so on stderr).

  /data/tmpA/nextgem/.venv/bin/python scripts/install_tripwire.py [--check]
"""

from __future__ import annotations

import argparse
import sys
import sysconfig
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
TRIPWIRE = PROJECT / "src" / "nextgem" / "tripwire.py"
MARK = "# nextgem label tripwire"

TEMPLATE = '''{mark} (DESIGN.md 11.1); written by scripts/install_tripwire.py. Do not edit.
import sys as _sys


def _nextgem_tripwire():
    import importlib.util
    path = {path!r}
    spec = importlib.util.spec_from_file_location("nextgem_tripwire", path)
    mod = importlib.util.module_from_spec(spec)
    _sys.modules["nextgem_tripwire"] = mod
    spec.loader.exec_module(mod)
    mod.install()


def _nextgem_fail_closed(exc):
    import re
    pat = re.compile(r"CRISPRGeneEffect|CRISPRGeneDependency|ScreenGeneEffect|CRISPRInferredCommonEssentials|"
                     r"AchillesCommonEssentialControls|AchillesNonessentialControls|labels_(train|val|test)|"
                     r"label_derived_lists|data/locked")

    def hook(event, args):
        if event == "open" and args and isinstance(args[0], (str, bytes)):
            p = args[0].decode(errors="replace") if isinstance(args[0], bytes) else args[0]
            if pat.search(p):
                raise PermissionError("nextgem tripwire failed to load (" + repr(exc) + "); all label access denied")

    _sys.addaudithook(hook)
    _sys.stderr.write("nextgem: tripwire failed to load; fail-closed hook active: " + repr(exc) + "\\n")


try:
    _nextgem_tripwire()
except Exception as _e:  # fail closed
    _nextgem_fail_closed(_e)
'''


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="only report whether the installed file is current")
    a = ap.parse_args(argv)
    if not sys.prefix.startswith(str(PROJECT / ".venv")):
        print(f"run with the project venv interpreter ({PROJECT}/.venv/bin/python), not {sys.executable}", file=sys.stderr)
        return 1
    target = Path(sysconfig.get_paths()["purelib"]) / "sitecustomize.py"
    text = TEMPLATE.format(mark=MARK, path=str(TRIPWIRE))
    current = target.read_text() if target.exists() else None
    if a.check:
        print(f"{target}: {'current' if current == text else 'missing or stale'}")
        return 0 if current == text else 1
    if current is not None and MARK not in current:
        print(f"{target} exists and is not ours; refusing to overwrite", file=sys.stderr)
        return 1
    target.write_text(text)
    print(f"installed {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
