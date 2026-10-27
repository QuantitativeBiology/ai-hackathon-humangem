"""W4 re-check (B2): the history-dependence repro through the author's one-job-per-process Runner, plus the
watchdog. Label-free. Gurobi lane: 1 worker (holds 1 pool slot at a time); GLPK lane: 1 worker.
 A: F3 ACH-002024 R1 run after F3 ACH-000957 in the same lane;  B: F3 ACH-002024 alone.
 C: F0 ACH-000957 run after F0 ACH-002024;                       D: F0 ACH-000957 alone.
 E: watchdog: an F0 job with a 10 s CPU limit followed by a normal job in the same lane."""
import hashlib, json, sys, time
import cobra
cobra.Configuration().solver = "glpk"
import pandas as pd
from nextgem import context_inputs as ci, context_run as cr, paths
W4 = paths.results_dir() / "w4"


def main():
    global genes, rna, st
    m = cobra.io.read_sbml_model(str(paths.data_root() / "Human-GEM-main/model/Human-GEM.xml"))
    genes = ci.gene_table(m); del m
    lines = ["ACH-000957", "ACH-002024"]
    rna = ci.read_expression(genes, lines)
    st = ci.ThresholdStats("rna", 0, 0, pd.Series(dtype=float), 0, 0)  # R1 is a constant threshold

    def job(line, fam, lim=None):
        lane = "gurobi" if fam == "F3" else "glpk"
        meta = {"row": "review", "solver": lane}
        if lim: meta["cpu_limit_s"] = lim
        return cr.Job(build_id=f"{line}|{fam}", line=line, family=fam, rule="R1", modality="rna", variant="real",
                      scores=ci.signed_scores(line, "R1", "rna", genes, rna, st), genomic=("none",), tasks=False, meta=meta)

    def run(items, workers, timeout=1800):
        out = {}
        def done(key, status, res, payload, wall):
            if status == "ok":
                r = res["rows"][0]; ko = res["ko"]["none"]
                out[key] = {"status": status, "open_sha": hashlib.sha256(r["open_mask"].encode()).hexdigest()[:16],
                            "n_open": r["n_open"], "ko_sha": hashlib.sha256(pd.util.hash_pandas_object(ko[["gene", "ratio_r"]], index=False).values.tobytes()).hexdigest()[:16],
                            "pid": r["pid"], "wall": round(wall)}
            else:
                out[key] = {"status": status, "wall": round(wall), "detail": (res or "")[-300:] if isinstance(res, str) else None}
        cr.Runner(workers, timeout, str(W4 / "fastcc_human2_M0.json")).run_all(items, done)
        return out

    t0 = time.time(); res = {}
    res["A_B_C_D"] = run([("A_prev", "gurobi", job("ACH-000957", "F3"), None, None), ("A", "gurobi", job("ACH-002024", "F3"), None, None),
                          ("C_prev", "glpk", job("ACH-002024", "F0"), None, None), ("C", "glpk", job("ACH-000957", "F0"), None, None),
                          ("B", "gurobi", job("ACH-002024", "F3"), None, None), ("D", "glpk", job("ACH-000957", "F0"), None, None)],
                         {"gurobi": 1, "glpk": 1})
    r = res["A_B_C_D"]
    res["F3_identical_A_B"] = r["A"]["open_sha"] == r["B"]["open_sha"] and r["A"]["ko_sha"] == r["B"]["ko_sha"]
    res["F0_identical_C_D"] = r["C"]["open_sha"] == r["D"]["open_sha"] and r["C"]["ko_sha"] == r["D"]["ko_sha"]
    res["distinct_pids"] = len({v["pid"] for v in r.values()}) == len(r)
    res["E"] = run([("E_timeout", "glpk", job("ACH-000957", "F0", lim=10), None, None), ("E_next", "glpk", job("ACH-002024", "F0"), None, None)], {"glpk": 1})
    res["wall_s"] = round(time.time() - t0)
    (W4 / "review" / "r10_runner_repro.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
