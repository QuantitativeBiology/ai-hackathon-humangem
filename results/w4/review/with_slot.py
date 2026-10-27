"""Review helper: hold one slot of the shared Gurobi pool (nextgem.gurobi_pool, n=4) while running a command.
Usage: .venv/bin/python results/w4/review/with_slot.py <cmd> [args...]"""
import subprocess
import sys
import time

from nextgem import gurobi_pool as gp

t0 = time.time()
slot = gp.hold_slot(timeout=3600, poll=10)
print(f"[with_slot] holding slot {slot.index} after {time.time() - t0:.0f}s", file=sys.stderr, flush=True)
try:
    rc = subprocess.run(sys.argv[1:]).returncode
finally:
    slot.release()
sys.exit(rc)
