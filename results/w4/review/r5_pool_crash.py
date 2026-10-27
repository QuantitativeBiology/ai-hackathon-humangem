"""W4 review, part 5: Gurobi slot pool semantics not covered by tests/test_gurobi_pool.py (no Gurobi session is
opened; a private lock directory under results/w4/review/tmp/ is used).
1. a holder killed with SIGKILL releases its slot (kernel drops the flock);
2. leftover lock files are not stale locks;
3. callers with different n share slot files 0..n-1 (the pool size is per caller, not global)."""
import json, os, signal, subprocess, sys, time
from pathlib import Path
from nextgem import gurobi_pool as gp
d = Path(__file__).parent / "tmp" / "poolcheck"; d.mkdir(parents=True, exist_ok=True)
res = {}
code = f"import time; from nextgem import gurobi_pool as gp; s = gp.hold_slot(1, {str(d)!r}, 'k', timeout=0); print('held', flush=True); time.sleep(60)"
p = subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.PIPE, text=True)
assert p.stdout.readline().strip() == "held"
res["busy_while_child_holds"] = gp.try_slot(1, d, "k") is None
os.kill(p.pid, signal.SIGKILL); p.wait()
s = gp.try_slot(1, d, "k")
res["free_after_sigkill"] = s is not None
res["lock_file_left_behind"] = (d / "k_0.lock").exists()
s.release()
a = [gp.hold_slot(4, d, "n", timeout=0) for _ in range(4)]
b = gp.try_slot(6, d, "n")
res["caller_with_n6_gets_5th_slot_while_4_held"] = b is not None and b.index == 4
for x in a + ([b] if b else []): x.release()
(Path(__file__).parent / "r5_pool.json").write_text(json.dumps(res, indent=1)); print(res)
