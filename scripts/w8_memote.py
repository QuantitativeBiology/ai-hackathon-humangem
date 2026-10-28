"""W8: MEMOTE on a sample of resource models (DESIGN W8: "MEMOTE (separate venv) and FROG are run on 10 models and
reported, not gated"). Label-free.

Runs ONLY under the separate venv /data/tmpA/nextgem/.venv_memote (memote 0.17.0 + its own cobra), never the main
.venv; it imports neither nextgem nor any project code, only memote/cobra and the standard library.

Expensive tests are skipped (``SKIP``): FVA- or per-metabolite LP sweeps over the whole Human-GEM network (tens of
thousands of LPs, or a MILP; several CPU-h per model) and the unscored matrix module (dense rank / nullspace; ~25
CPU-min per model in a timing run). Everything else (annotation, SBO,
basic, biomass, the remaining consistency tests: stoichiometric consistency, charge/mass balance, energy-
generating cycles, orphans/dead-ends/disconnected) runs. The score is MEMOTE's own weighted score over the tests
that ran; skipped tests count as 0 in MEMOTE's score, so the score is reported together with the skip list.

Usage:
  .venv_memote/bin/python scripts/w8_memote.py MODEL.xml.gz [...] --out results/w8/memote [--solver-timeout 300]
      [--summary memote_summary.json]
"""

from __future__ import annotations

import argparse
import gzip
import json
import resource as _res
import sys
import time
from pathlib import Path

SKIP = (
    "test_find_stoichiometrically_balanced_cycles",  # FVA over all reactions
    "test_blocked_reactions",  # FVA over all reactions
    "test_find_reactions_unbounded_flux_default_condition",  # FVA over all reactions
    "test_find_metabolites_not_produced_with_open_bounds",  # one LP per metabolite (~8.4k)
    "test_find_metabolites_not_consumed_with_open_bounds",  # one LP per metabolite (~8.4k)
    "test_unconserved_metabolites",  # one LP per metabolite
    "test_inconsistent_min_stoichiometry",  # MILP enumeration
    "test_matrix",  # rank / nullspace of the dense 8.4k x 13k matrix: ~25 CPU-min per model; unscored by MEMOTE
)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("models", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--solver-timeout", type=int, default=300)
    ap.add_argument("--summary", default="memote_summary.json", help="summary file name under --out")
    a = ap.parse_args(argv)
    if "nextgem" in sys.modules or ".venv_memote" not in sys.prefix:
        raise SystemExit("run with /data/tmpA/nextgem/.venv_memote/bin/python (separate venv)")
    import cobra
    import memote
    from memote.suite.api import snapshot_report, test_model, validate_model

    cobra.Configuration().solver = "glpk"
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for mp in a.models:
        mp = Path(mp)
        stem = mp.name.replace(".xml.gz", "")
        t0, c0 = time.time(), time.process_time()
        raw = gzip.decompress(mp.read_bytes())
        tmp = out / f".{stem}.xml"
        tmp.write_bytes(raw)
        try:
            model, sbml_ver, notifications = validate_model(str(tmp))
        finally:
            tmp.unlink()
        if model is None:
            rows.append({"model": mp.name, "status": "invalid", "errors": notifications.get("errors", [])[:5]})
            continue
        model.solver = "glpk"
        code, result = test_model(model, sbml_version=sbml_ver, results=True, skip=list(SKIP),
                                  solver_timeout=a.solver_timeout)
        rep = json.loads(snapshot_report(result, html=False))
        score = rep.get("score", {})
        sections = {s["section"]: round(s["score"], 4) for s in score.get("sections", [])}
        tests = {}
        for name, case in rep.get("tests", {}).items():
            r = case.get("result")
            if isinstance(r, dict):  # parametrised test
                vals = list(r.values())
                r = "passed" if all(v == "passed" for v in vals) else ",".join(sorted(set(vals)))
            tests[name] = r
        with gzip.open(out / f"{stem}.memote.json.gz", "wt") as fh:
            json.dump(rep, fh, default=str)
        row = {"model": mp.name, "status": "ok", "pytest_code": int(code),
               "total_score": round(score.get("total_score", float("nan")), 4), "sections": sections,
               "n_tests": len(tests), "n_passed": sum(v == "passed" for v in tests.values()),
               "n_failed": sum(v == "failed" for v in tests.values()),
               "n_skipped": sum(v == "skipped" for v in tests.values()),
               "failed": sorted(k for k, v in tests.items() if v == "failed"),
               "wall_s": round(time.time() - t0, 1), "cpu_s": round(time.process_time() - c0, 1)}
        rows.append(row)
        print(json.dumps(row), flush=True)
    summ = {"memote_version": memote.__version__, "cobra_version": cobra.__version__, "solver": "glpk",
            "solver_timeout_s": a.solver_timeout, "skipped_tests": list(SKIP), "models": rows,
            "peak_rss_mb": round(_res.getrusage(_res.RUSAGE_SELF).ru_maxrss / 1024, 1)}
    (out / a.summary).write_text(json.dumps(summ, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
