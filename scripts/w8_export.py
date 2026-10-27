"""W8 export stub / demo: SBML L3V1 FBC v2 export + validation gates (DESIGN.md section 9 W8). Label-free.

Demo mode (default) writes under results/w8/demo/ and tags every model "demo, not a release". It exports the two
generic parents on M0 (Ham's) and N per-line models reconstructed from the stored ``open_mask`` of Stage-A F0
builds (results/w4/stage_a/builds.tsv). The real resource (resource/v{X.Y.Z}/ + resource/manifest.tsv) is built
after T1 with ``--release X.Y.Z``; this stub refuses ``--release`` until the build driver exists.

GLPK only, one process. Usage: .venv/bin/python scripts/w8_export.py [--lines ACH-... ACH-...] [--out DIR]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"  # before any model is loaded

import pandas as pd  # noqa: E402

from nextgem import paths, resource as rs  # noqa: E402
from nextgem.knockouts import configure  # noqa: E402
from nextgem.model import GEMS, load_model  # noqa: E402
from nextgem.run import git_commit, provenance, sha256_file  # noqa: E402

DEMO_TAG = "demo, not a release"
DEFAULT_LINES = ["ACH-001283", "ACH-003383"]  # Bone / Bowel A6' lines of the W4 smoke set (descriptive, A7)
FAMILY, RULE = "F0", "R1"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lines", nargs="*", default=DEFAULT_LINES)
    ap.add_argument("--out", default=str(paths.results_dir() / "w8" / "demo"))
    ap.add_argument("--stage-a", default=str(paths.results_dir() / "w4" / "stage_a"))
    ap.add_argument("--release", default=None, help="X.Y.Z: real resource build (not implemented in the stub)")
    a = ap.parse_args(argv)
    if a.release:
        print("real resource build is not implemented (stub; builds start after T1)", file=sys.stderr)
        return 2
    out = Path(a.out)
    builds_path = Path(a.stage_a) / "builds.tsv"
    inputs = [builds_path, GEMS["human2"].sbml, GEMS["human1"].sbml, GEMS["human2"].essential_tasks,
              GEMS["human1"].essential_tasks]
    with provenance(paths.results_dir() / "w8" / "PROVENANCE.tsv", inputs=inputs, seeds=[], solver="glpk",
                    workers=1, allow_dirty=True, author="W8 resource engineer agent") as prov:
        t_all = time.perf_counter()
        commit = git_commit()
        shas = {Path(p).name: sha256_file(p) for p in inputs}
        rows, report = [], []
        builds = pd.read_csv(builds_path, sep="\t", usecols=["line", "family", "rule", "variant", "genomic", "failed",
                                                              "n_reopened", "open_mask", "wt_growth"])
        for parent_name in ("human2", "human1"):
            gem = GEMS[parent_name]
            t0 = time.perf_counter()
            parent = load_model(parent_name)
            configure(parent)  # M0 = Ham's medium (GR task inputs), biomass objective
            raw_errors = rs.libsbml_errors(gem.sbml)
            ver = gem.version.lstrip("v")
            variants = [("generic", "none", parent, 0, 0, None)]
            if parent_name == "human2":
                for line in a.lines:
                    r = builds[(builds.line == line) & (builds.family == FAMILY) & (builds.rule == RULE)
                               & (builds.variant == "real") & (builds.genomic == "none") & (~builds.failed)]
                    if len(r) != 1:
                        raise SystemExit(f"{line}: expected 1 {FAMILY}/{RULE} real/none build, found {len(r)}")
                    m, n_closed = rs.line_model_from_mask(parent, r.iloc[0].open_mask)
                    variants.append((line, f"{FAMILY}{RULE}", m, n_closed, int(r.iloc[0].n_reopened),
                                     float(r.iloc[0].wt_growth)))
            for model_id, pipeline, model, n_closed, n_rep, stage_wt in variants:
                t1 = time.perf_counter()
                meta = rs.Meta(model_id=model_id, parent=parent_name, parent_version=ver, pipeline_id=pipeline,
                               git_commit=commit, input_shas=shas, n_repaired_reactions=n_rep,
                               medium_source_tier="generic Ham's (GR task); line medium tiers not applied in demo",
                               release_status=DEMO_TAG, benchmarked=False,
                               flags="demo;mask-reconstructed;repaired_list_not_stored" if model_id != "generic" else "demo")
                path = rs.export_model(model, meta, out / "v0.0.0-demo")
                meta.build_time_s = time.perf_counter() - t1
                tasks = gem.essential_tasks
                val = rs.validate_model(model, path, gem.sbml, tasks, parent_errors=raw_errors)
                # re-export is not needed for build_time; manifest records export+validation time
                meta.build_time_s = time.perf_counter() - t1
                row = rs.manifest_row(path, meta, model, val, n_closed)
                rows.append(row)
                gates = {k: (v[0] if not isinstance(v[0], bool) else ("pass" if v[0] else "FAIL"))
                         for k, v in val.gates.items()}
                grow = val.gates.get("wt_growth_match", (None, (None, None)))[1]
                report.append((path.name, row["size_bytes"], val.status, gates, grow, stage_wt,
                               round(meta.build_time_s, 1)))
                print(path.name, row["size_bytes"], val.status, gates, "growth", grow, "stageA", stage_wt, flush=True)
            print(f"{parent_name} done in {time.perf_counter() - t0:.0f} s", flush=True)
            del parent
        man = out / "manifest.tsv"
        rs.write_manifest(rows, man)
        prov.output(man, *[Path(r["file"]) for r in rows])
        prov.note(release_status=DEMO_TAG, n_models=len(rows), all_pass=all(r["validation_status"] == "pass" for r in rows),
                  total_wall_s=round(time.perf_counter() - t_all, 1))
    print("total", round(time.perf_counter() - t_all, 1), "s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
