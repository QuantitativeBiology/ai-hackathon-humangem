"""W8: FROG reproducibility report (fbc_curation, cobrapy curator) for resource models (DESIGN W8: "MEMOTE (separate
venv) and FROG are run ... and reported, not gated"). Label-free.

Runs ONLY under the separate venv /data/tmpA/nextgem/.venv_memote (fbc_curation 0.3.2 + its own cobra/GLPK, single
process: fbc_curation sets cobra's ``processes = 1``); imports no project code.

FROG = objective value, FVA at 100% of the optimum on every reaction, every single-gene and every single-reaction
deletion. Only the cobrapy curator is run (``runfrog`` would also run the cameo curator: a second, independent
implementation of the same four analyses, i.e. twice the cost). Each part is timed separately.

Outputs per model under --out/<stem>/: FROG TSVs + JSON (``FrogReport.to_tsv`` / ``to_json``) and timing.json; a
frog_summary.json over all models.

Usage:
  .venv_memote/bin/python scripts/w8_frog.py MODEL.xml.gz [...] --out results/w8/frog [--parts objective,gene,...]
"""

from __future__ import annotations

import argparse
import gzip
import json
import resource as _res
import sys
import tempfile
import time
from pathlib import Path

PARTS = ("objective", "gene_deletions", "reaction_deletions", "fva")


def _cpu() -> float:
    r = _res.getrusage(_res.RUSAGE_SELF)
    c = _res.getrusage(_res.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime + c.ru_utime + c.ru_stime


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("models", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--parts", default=",".join(PARTS))
    a = ap.parse_args(argv)
    if "nextgem" in sys.modules or ".venv_memote" not in sys.prefix:
        raise SystemExit("run with /data/tmpA/nextgem/.venv_memote/bin/python (separate venv)")
    parts = [p for p in a.parts.split(",") if p]
    if set(parts) - set(PARTS):
        raise SystemExit(f"unknown parts {set(parts) - set(PARTS)}")
    import cobra
    import fbc_curation
    from fbc_curation.curator.cobrapy_curator import CuratorCobrapy
    from fbc_curation.frog import FrogReport

    cobra.Configuration().solver = "glpk"
    cobra.Configuration().processes = 1
    out = Path(a.out)
    rows = []
    for mp in a.models:
        mp = Path(mp)
        stem = mp.name.replace(".xml.gz", "")
        odir = out / stem
        odir.mkdir(parents=True, exist_ok=True)
        timing = {}
        with tempfile.TemporaryDirectory() as td:
            x = Path(td) / f"{stem}.xml"
            x.write_bytes(gzip.decompress(mp.read_bytes()))
            cur = CuratorCobrapy(model_path=x, frog_id="cobrapy", curators=[])
            res = {"metadata": cur.set_metadata()}
            calls = {"objective": ("objectives", cur.objectives), "gene_deletions": ("gene_deletions", cur.gene_deletions),
                     "reaction_deletions": ("reaction_deletions", cur.reaction_deletions), "fva": ("fva", cur.fva)}
            for p in parts:
                t0, c0 = time.time(), _cpu()
                key, fn = calls[p]
                res[key] = fn()
                timing[p] = {"wall_s": round(time.time() - t0, 1), "cpu_s": round(_cpu() - c0, 1)}
                print(stem, p, timing[p], flush=True)
        if set(parts) == set(PARTS):
            rep = FrogReport(**res)
            rep.to_tsv(odir)
            dfs = rep.to_dfs()
        else:
            dfs = {k: v.to_df() for k, v in res.items() if k != "metadata"}
            for k, df in dfs.items():
                df.to_csv(odir / f"{k}.tsv", sep="\t", index=False)
        summ = {"model": mp.name, "parts": parts, "timing": timing,
                "objective": float(dfs["objectives"].value.iloc[0]) if "objectives" in dfs else None,
                "n_gene_deletions": int(len(dfs["gene_deletions"])) if "gene_deletions" in dfs else None,
                "n_reaction_deletions": int(len(dfs["reaction_deletions"])) if "reaction_deletions" in dfs else None,
                "n_fva": int(len(dfs["fva"])) if "fva" in dfs else None,
                "status_counts": {k: df.status.astype(str).value_counts().to_dict() for k, df in dfs.items()
                                  if "status" in df}}
        (odir / "timing.json").write_text(json.dumps(summ, indent=1, default=str) + "\n")
        rows.append(summ)
    (out / "frog_summary.json").write_text(json.dumps(
        {"fbc_curation_version": fbc_curation.__version__, "cobra_version": cobra.__version__, "solver": "glpk",
         "curator": "cobrapy", "models": rows}, indent=1, default=str) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
