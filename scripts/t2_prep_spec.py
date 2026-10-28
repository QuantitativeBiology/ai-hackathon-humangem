"""T2 prep: assemble results/m6/p2_spec.json (template results/m6/p2_spec.template.json) for P2 = W7-05
(RESEARCH_TREE D12, D14: B6 with a monotone increasing constraint on b2/b2m, non-GEM, label-fitted).

P2 is label-fitted, so the harness (scripts/evaluate_test.py) does NOT read P2 from a score file: it refits W7-05
once on training ∪ validation through the reviewed ``w7_idea`` plug-in (``refit = {w7_idea, W7-05, P2}``). The
spec's test score file therefore carries only the generic columns H1 and H2 (gene layout, line-constant), which the
harness checks equal to the generic scores of the W1 look.

Three outputs, all below results/m6/t2/ except the spec (tracked; no line-keyed or label-shaped content):

1. ``p2_scores_test.parquet``: entrez, H1 = e_H1, H2 = e_H2 from ``nextgem.w1.generic_table`` over the gene
   universe U (label-free).
2. ``skew.json``: the §6.5 check that fixes ``bootstrap_method`` for C4 / C5. Rule: BCa iff |skew| > 1 of the
   training-tier line ΔtAP. Decided on C4 (P2 − G*, first in the chain), C5 (P2 − B6) reported alongside and a
   disagreement flagged for the PI (as T1 decided on C1 and reported C1b). Source: the W7 training screen
   (nested 12-fold LOLO on the 686 W1 primary training lines, L_prob, scripts/w7_screen.py):
   results/w7/tmp/perline_r2b_W7-05.parquet (W7-05) and perline_r1_W7-01.parquet (G* = H1 as ``Gstar``; B6 =
   W7-01, the B6 ablation), seed SEEDS[0]. These are per-line tAP scalars on training; no label file is opened
   here and no validation or test quantity is read. Input SHA-256s are recorded in skew.json.
3. ``results/m6/p2_spec.json``, checked with ``evaluate_test.check_p2_spec`` and the T2 config rules
   (primary label L_prob, b6.hgb == nextgem.w7.HGB_M2) before it is written.

It is NOT committed or frozen here: the PI commits it, adds the FROZEN.tsv rows listed by
``scripts/evaluate_test.py --look T2 --print-frozen-candidates``, then appends the W7 stop row.

Usage: .venv/bin/python scripts/t2_prep_spec.py [--check]   (--check: rebuild in memory, compare, write nothing)
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
M6 = PROJECT / "results" / "m6"
T2 = M6 / "t2"  # tracked
SPEC_PATH = M6 / "p2_spec.json"
SCORES_PATH = T2 / "p2_scores_test.parquet"
SKEW_PATH = T2 / "skew.json"
MANIFEST = T2 / "MANIFEST.tsv"
AUTHOR = "T2-preparation engineer agent (Claude Opus)"

SPEC_ID = "P2_M7"
MILESTONE = "M7"
PRIMARY = "P2"
IDEA = "W7-05"
METHODS = ["H1", "H2", PRIMARY]
GRID = {"labels": ["L_prob", "L_ge05", "L_ge10"], "genesets": ["all", "noCE", "noDepMapCE"],
        "universes": ["U", "Ucap"], "strata": ["none", "expressed_pairs"]}  # the full W1 grid (§6.6)
W7_TMP = PROJECT / "results" / "w7" / "tmp"
SKEW_SOURCES = {"p2": (W7_TMP / "perline_r2b_W7-05.parquet", IDEA),
                "gstar": (W7_TMP / "perline_r1_W7-01.parquet", "Gstar"),
                "b6": (W7_TMP / "perline_r1_W7-01.parquet", "W7-01")}
SELECTION = ("results/w7/val/look_spec_W7_P2_select.json", "results/w7/val/val_look_W7_P2_select_decision.json")
PROVENANCE_FILES = ("src/nextgem/w7.py", "scripts/w1_b6.py", "scripts/t2_prep_spec.py",
                    "results/m6/t2/skew.json") + SELECTION


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(p: Path) -> str:
    return str(Path(p).resolve().relative_to(PROJECT))


@functools.lru_cache(maxsize=1)
def _evaluate_test():
    s = importlib.util.spec_from_file_location("evaluate_test", PROJECT / "scripts" / "evaluate_test.py")
    et = importlib.util.module_from_spec(s)
    s.loader.exec_module(et)
    return et


# ============================================================================= pure helpers (tested)
def generic_scores(gen: pd.DataFrame) -> pd.DataFrame:
    """The P2 spec's test score file: entrez + generic H1/H2 (``e_H1``/``e_H2`` of ``w1.generic_table``)."""
    out = pd.DataFrame({"entrez": gen.index.astype("int64"),
                        "H1": gen["e_H1"].to_numpy(dtype="float64"), "H2": gen["e_H2"].to_numpy(dtype="float64")})
    return out.sort_values("entrez", kind="mergesort").reset_index(drop=True)


def line_deltas(perline: dict[str, pd.DataFrame], a: str, b: str, seed: int) -> pd.Series:
    """Per-line ΔtAP (a − b) on the W7 screen's training lines. ``perline[name]`` holds the per-line table of
    method ``name`` (columns seed, ModelID, tap); the fitted method uses ``seed``; a fixed comparator that has
    only one seed (G*) uses it whatever its value. Lines must match exactly."""
    def tap(name):
        t = perline[name]
        if t.seed.nunique() > 1:
            t = t[t.seed == seed]
        if t.ModelID.duplicated().any():
            raise ValueError(f"{name}: duplicate lines at seed {seed}")
        return t.set_index("ModelID").tap.astype(float)

    x, y = tap(a), tap(b)
    if set(x.index) != set(y.index):
        raise ValueError(f"{a} and {b} cover different lines")
    return (x - y.reindex(x.index)).sort_index()


def skew_decision(deltas: dict[str, pd.Series], decide_on: str = "C4") -> dict:
    """§6.5: BCa iff |skew| of the training line ΔtAP > 1 (``metrics.skew_triggers_bca``), decided on
    ``decide_on``; every contrast reported, a disagreement flagged."""
    from scipy import stats

    from nextgem import metrics

    out = {}
    for k, d in deltas.items():
        v = pd.Series(d, dtype=float).dropna()
        bca = metrics.skew_triggers_bca(v)
        out[k] = {"skew": float(stats.skew(v.to_numpy())), "n_lines": int(v.size), "bca": bca,
                  "bootstrap_method": "bca" if bca else "percentile", "mean_delta": float(v.mean())}
    method = out[decide_on]["bootstrap_method"]
    return {"contrasts": out, "decided_on": decide_on, "bootstrap_method": method,
            "disagrees": sorted(k for k, v in out.items() if v["bootstrap_method"] != method)}


def build_spec(scores_entry: dict, bootstrap_method: str, provenance: dict) -> dict:
    from nextgem import w7

    notes = (
        "P2 (RESEARCH_TREE D14, look 5/40 W7_P2_select, W7 2/20) = W7-05: the pre-registered B6 (W1 feature set, "
        "HistGradientBoosting at w7.HGB_M2) with a monotone increasing constraint on b2 and b2m; non-GEM, label-"
        "fitted. It is refitted once on training ∪ validation inside the harness through the reviewed w7_idea "
        "plug-in (w7.fit / w7.predict, seed SEEDS[0]); the score file carries only the generic H1/H2 columns "
        "(gene layout). No failed build (P2 builds no model). D12 wording applies in REPORT.md: C4 (P2 − G*) "
        "compares gene-prior ML with generic FBA and says nothing about GEMs; C5 (P2 − B6, the pre-registered B6) "
        "measures only the monotone constraint; W7's GEM result is negative (no idea with a GEM feature was "
        "promoted); P2 is reported next to B2, B2m and B8±. C3 was not confirmed at T1 (D15), so C4/C5 are "
        f"estimates only. bootstrap_method from the §6.5 training skew check (results/m6/t2/skew.json): "
        f"{bootstrap_method}."
    )
    return {
        "spec_id": SPEC_ID,
        "kind": "p2",
        "milestone": MILESTONE,
        "pipeline": (f"{IDEA} (B6 with a monotone increasing constraint on b2 and b2m; nextgem.w7.IDEAS['{IDEA}']), "
                     "refitted on training ∪ validation (selected in the W7 validation look W7_P2_select, D14)."),
        "primary": PRIMARY,
        "methods": list(METHODS),
        "scores": {"test": dict(scores_entry)},
        "failed_builds": {"test": None},
        **{k: list(v) for k, v in GRID.items()},
        "b6": {"features": "w1", "hgb": dict(w7.HGB_M2)},
        "bootstrap_method": bootstrap_method,
        "fitted_on_labels": True,
        "refit": {"recipe": "w7_idea", "idea": IDEA, "method": PRIMARY},
        "provenance": dict(provenance),
        "notes": notes,
    }


def check(spec: dict) -> None:
    """Static spec checks plus the T2 config rules that need no label (check_p2_spec; hgb == w7.HGB_M2)."""
    from nextgem import w7

    et = _evaluate_test()
    et.check_p2_spec(spec)
    if spec["b6"]["hgb"] != dict(w7.HGB_M2):
        raise et.PreflightError("b6.hgb differs from w7.HGB_M2")
    idea = w7.IDEAS[spec["refit"]["idea"]]
    if idea.kind != "clf" or idea.extra or idea.drop or idea.mono_extra or not idea.mono:
        raise et.PreflightError(f"{idea.id} is not the monotone B6-feature classifier")


# ============================================================================= real-data steps
def make_scores() -> pd.DataFrame:
    from nextgem import labels, w1

    U = labels.gene_universe()
    return generic_scores(w1.generic_table(U.entrez))


def make_skew() -> dict:
    from nextgem.config import SEEDS

    tabs, ins = {}, {}
    for name, (path, method) in SKEW_SOURCES.items():
        df = pd.read_parquet(path, columns=["method", "seed", "ModelID", "tap"])
        tabs[name] = df[df.method == method]
        if tabs[name].empty:
            raise RuntimeError(f"{rel(path)} has no rows for {method}")
        ins[rel(path)] = sha256(path)
    d = {"C4": line_deltas(tabs, "p2", "gstar", SEEDS[0]), "C5": line_deltas(tabs, "p2", "b6", SEEDS[0])}
    dec = skew_decision(d)
    seeds = sorted(set(tabs["p2"].seed) & set(tabs["b6"].seed))
    sens = {f"seed_{s}": {k: round(float(v["skew"]), 6) for k, v in
                          skew_decision({"C4": line_deltas(tabs, "p2", "gstar", s),
                                         "C5": line_deltas(tabs, "p2", "b6", s)})["contrasts"].items()}
            for s in seeds}
    return {"rule": "DESIGN §6.5: BCa for C4/C5 iff |skew| of the training-tier line ΔtAP > 1 "
                    "(metrics.skew_triggers_bca); decided on C4 (P2 − G*), C5 (P2 − B6) reported",
            "label": "L_prob", "lines": "W1 primary training set, W7 screen nested 12-fold LOLO (held-out lines)",
            "p2": f"{IDEA} (seed {SEEDS[0]})", "gstar": "H1 (Gstar)", "b6": f"W7-01 = B6 (seed {SEEDS[0]})",
            "synthetic_labels": False, **dec, "skew_by_seed": sens, "inputs": ins}


def build_all() -> tuple[pd.DataFrame, dict]:
    return make_scores(), make_skew()


def write_parquet(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_parquet(tmp, index=False)
    os.replace(tmp, path)
    return path


def write_json(obj, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare with the files on disk")
    a = ap.parse_args(argv)
    from nextgem.config import SEEDS
    from nextgem.run import provenance

    scores, skew = build_all()
    if a.check:
        et = _evaluate_test()
        spec = json.loads(SPEC_PATH.read_text())
        et.check_p2_spec(spec)
        disk = pd.read_parquet(SCORES_PATH)
        ok = {"scores": disk.equals(scores), "skew": json.loads(SKEW_PATH.read_text()) == json.loads(json.dumps(skew)),
              "spec_pins": all(sha256(PROJECT / p) == s for p, s in spec["provenance"].items())
              and sha256(SCORES_PATH) == spec["scores"]["test"]["sha256"]}
        print(json.dumps(ok, indent=1))
        return 0 if all(ok.values()) else 1
    inputs = [PROJECT / p for p in ("src/nextgem/w7.py", "scripts/w1_b6.py", "scripts/t2_prep_spec.py",
                                    "scripts/evaluate_test.py", *SELECTION)] + sorted({p for p, _ in SKEW_SOURCES.values()})
    with provenance(MANIFEST, inputs=inputs, seeds=[SEEDS[0]], allow_dirty=True, author=AUTHOR) as prov:
        write_parquet(scores, SCORES_PATH)
        write_json(skew, SKEW_PATH)
        entry = {"path": rel(SCORES_PATH), "sha256": sha256(SCORES_PATH), "layout": "gene"}
        spec = build_spec(entry, skew["bootstrap_method"], {p: sha256(PROJECT / p) for p in PROVENANCE_FILES})
        check(spec)
        write_json(spec, SPEC_PATH)
        prov.output(SCORES_PATH, SKEW_PATH, SPEC_PATH)
    print(json.dumps({"spec": rel(SPEC_PATH), "sha256": sha256(SPEC_PATH), "scores_sha256": entry["sha256"],
                      "n_genes": len(scores), "bootstrap_method": skew["bootstrap_method"],
                      "skew": {k: round(v["skew"], 4) for k, v in skew["contrasts"].items()},
                      "skew_disagrees": skew["disagrees"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
