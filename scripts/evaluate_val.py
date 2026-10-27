"""Validation looks (DESIGN.md §4.2, §10, §11.1, §11.3, §11.4). The only entry point allowed to read validation labels.

Order of a look (everything that does not need validation labels happens BEFORE the ledger ``started`` row,
so that no label-free failure can strand a look id after ``started`` is committed; review M2_W1 A5):

1. Pre-flight (refuses to start, nothing written):
   - every file in ``required_files(cfg)`` (this script, labels.py, w1.py, metrics.py, run.py and every input of
     the look) has a ``file`` row in the **committed** ``HEAD:FROZEN.tsv`` and its current SHA-256 equals it;
   - those files are clean in git (tracked and unmodified; git-ignored inputs such as ``data/`` and
     ``results/**/tmp/`` are checked by SHA only);
   - ``scripts/check_frozen.py`` passes **with** the test-run check (no ``--no-test-runs``; §11.3).
2. ``nextgem.run.provenance`` (refuses a dirty src/scripts/tests; writes a MANIFEST row on success).
3. Label-free and training-only work: the label-free context (scores' inputs), fitting on ``labels.train()``,
   all validation-line score matrices. The config JSON (input/code SHA-256 and its own SHA) is written.
4. Ledger ``started`` (committed through ``nextgem.ledger``; look cap 40, <= 20 for W7; a look id is never
   repeated). A resumed crashed look (same config SHA) appends an ``output`` row whose detail is
   ``{"event": "resumed", ...}`` before any label is read (``ledger.py`` is frozen and has no ``resumed`` event;
   the marker row does not create a new look, so the cap counts the look once).
5. Only now: ``labels.validation()``, scoring, scalar-only outputs (guarded by :func:`assert_scalar_only`),
   ledger ``output`` (SHA-256 of the outputs) and ``completed``.

Look implemented: ``W1_M2_baseline`` (M2, W1). Generic Human1/Human2 and B0–B3, B8±/B8′ (B6 is excluded:
``include_b6 = False``; it is untuned and W7 owns its validation looks) on the validation expression lines:
tAP grid scalars (L_prob/L_ge05/L_ge10 x {all, minus CE_fit, minus DepMap CE} x {U, U∩} x {none, expressed
pairs}) with paired Δ̄_lin vs H1 and vs H2; §6.8 binary endpoints for generic H1/H2 and the paired H2 − H1 MCC;
and G-star: **Human1 iff Δ̄_lin(H2 − H1) < 0** (L_prob, all genes, U), otherwise Human2 (a tie goes to Human2,
the default parent, DESIGN §3.1). CE_fit, B2/B2m/B2c come from all training lines (§5.3).

Usage:
  .venv/bin/python scripts/evaluate_val.py --look W1_M2_baseline                       # the look (once)
  .venv/bin/python scripts/evaluate_val.py --look W1_M2_baseline --print-config        # reads no label
  .venv/bin/python scripts/evaluate_val.py --look W1_M2_baseline --print-frozen-candidates   # reads no label
  .venv/bin/python scripts/evaluate_val.py --look W1_M2_baseline --preflight-only      # reads no label
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
LABEL_COLS = ["gene_effect", "dep_prob", "L_prob", "L_ge05", "L_ge06", "L_ge10"]
W1_STRATA = ("none", "expressed_pairs")
W1_GENESETS = ("all", "noCE", "noDepMapCE")
W1_UNIVERSES = ("U", "Ucap")
INCLUDE_B6 = False  # review M2_W1 item 4: B6 is untuned/provisional; W7 owns its validation looks
CODE_FILES = ("scripts/evaluate_val.py", "src/nextgem/labels.py", "src/nextgem/w1.py", "src/nextgem/metrics.py",
              "src/nextgem/run.py")
KNOWN_ANSWER_TESTS = ("tests/test_metrics.py", "tests/test_metrics_binary.py", "tests/test_evaluate_val.py",
                      "tests/test_run.py")
BINARY_PROTOCOLS = (("Human1-paper", "P57", "L_ge06"), ("Human1-paper", "PGR", "L_ge06"), ("Human2-paper", "P57", "L_ge06"),
                    ("fixed-cutoff", "FBA<0.01", "L_prob"), ("fixed-cutoff", "FBA<0.01", "L_ge06"))
PAIRED_PROTOCOLS = BINARY_PROTOCOLS[:4]
GRID_KEY = ["section", "label", "geneset", "universe", "stratum", "method"]
BINARY_KEY = ["section", "protocol", "method", "call_rule", "label", "gene_set", "metric", "nan_policy"]
ALLOWED_NUMERIC = {
    "n_lines", "n_lines_estimable", "na", "n_lineages", "tap_lin", "tap_lin_lo", "tap_lin_hi", "e_random_lin",
    "lift_lin", "delta_null_lin", "tap_line", "e_random_line", "lift_line", "prevalence_lin", "n_genes_estimable",
    "m_delta_ap_gene", "pooled_tap", "pooled_e_random", "pooled_lift", "pooled_n", "pooled_n_pos", "auroc_lin",
    "ap_sklearn_lin", "lineage_equal_mean", "line_mean", "median", "n_lines_used", "n_nan", "ci_lo", "ci_hi",
    *[f"{k}_vs_{r}" for k in ("delta_lin", "delta_lo", "delta_hi", "delta_line") for r in ("H1", "H2")],
}
ID_RE = re.compile(r"^ACH-\d+$")


class PreflightError(RuntimeError):
    pass


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------------------- configuration
def input_files(include_b6: bool = INCLUDE_B6) -> dict[str, str]:
    """Inputs of the W1 look, as paths relative to the project."""
    d = {
        "human1_ko": "results/baseline/human1_ko.parquet", "human2_ko": "results/baseline/human2_ko.parquet",
        "p57_human1": "results/w1/p57_human1.parquet", "p57_human2": "results/w1/p57_human2.parquet",
        "line_features": "results/w1/line_features.tsv", "gene_features": "results/w1/gene_features.tsv",
        "expr": "results/w1/tmp/expr.parquet", "cn": "results/w1/tmp/cn.parquet",
        "gene_panels": "data/controls/gene_panels.tsv", "sesoi_c1b": "results/w1/sesoi_c1b.json",
    }
    if include_b6:
        d.update({"b6_full": "results/w1/tmp/b6_full.pkl", "gene_neighbours": "results/w1/tmp/gene_neighbours.json",
                  "mut_dam": "results/w1/tmp/mut_dam.parquet", "mut_hot": "results/w1/tmp/mut_hot.parquet"})
    return d


def required_files(cfg: dict) -> list[str]:
    """Files that must be frozen (committed HEAD:FROZEN.tsv) and clean before the look starts."""
    return list(CODE_FILES) + list(cfg["_inputs_rel"].values())


def w1_config(include_b6: bool = INCLUDE_B6) -> dict:
    """Everything that defines the W1 look; its SHA is the ledger config_sha. Label-free."""
    from nextgem import metrics, paths, w1
    from nextgem.config import SEEDS

    root = paths.state_root()
    rel = input_files(include_b6)
    return {
        "look_id": "W1_M2_baseline", "workstream": "W1", "milestone": "M2", "include_b6": include_b6,
        "labels": list(w1.GRID_LABELS), "binary_label": "L_ge06", "genesets": list(W1_GENESETS),
        "universes": list(W1_UNIVERSES), "strata": list(W1_STRATA), "B": metrics.B_OTHER, "seed": SEEDS[0],
        "random_seeds": SEEDS[0:5],
        "gstar_rule": "Human1 iff delta_lin(H2 - H1) < 0 on L_prob/all genes/U/validation expression lines; else Human2 "
                      "(a tie goes to Human2, the default parent)",
        "inputs": {k: _sha(root / v) for k, v in rel.items()},
        "code": {Path(c).name: _sha(paths.PROJECT_ROOT / c) for c in CODE_FILES},
        "_inputs_rel": rel,
    }


def config_sha(cfg: dict) -> str:
    c = {k: v for k, v in cfg.items() if not k.startswith("_")}
    return hashlib.sha256(json.dumps(c, sort_keys=True).encode()).hexdigest()


def frozen_candidates(cfg: dict) -> list[str]:
    """Candidate FROZEN.tsv rows (item, kind, path, sha256_or_value, frozen_on, amendment, note) for the PI:
    every required file and known-answer test not already frozen with its current SHA in HEAD:FROZEN.tsv."""
    from nextgem import frozen, paths

    head = subprocess.run(["git", "-C", str(paths.state_root()), "show", "HEAD:FROZEN.tsv"], capture_output=True, text=True)
    have = {r["path"]: r["sha256_or_value"] for r in frozen.parse(head.stdout) if r["kind"] == "file"} if head.returncode == 0 else {}
    day = time.strftime("%Y-%m-%d")
    rows = []
    for rel in required_files(cfg) + list(KNOWN_ANSWER_TESTS):
        p = paths.state_root() / rel
        if have.get(rel) == _sha(p):
            continue
        note = "W1 validation look: code" if rel in CODE_FILES else (
            "known-answer test (DESIGN 11.3)" if rel in KNOWN_ANSWER_TESTS else "W1 validation look: input")
        rows.append("\t".join([rel, "file", rel, _sha(p), day, "", note]))
    s = json.loads((paths.results_dir() / "w1" / "sesoi_c1b.json").read_text())
    rows.append("\t".join(["sesoi.C1b", "value", "", repr(s["SESOI_C1b"]), day, "",
                           "DESIGN 6.4, frozen from training at M2 (results/w1/sesoi_c1b.json)"]))
    return rows


# ----------------------------------------------------------------------------- pre-flight
def preflight(required: list[str], root: Path | None = None) -> list[str]:
    """Problems that block the look (empty list = OK). Checks committed HEAD:FROZEN.tsv and git cleanliness."""
    from nextgem import frozen, paths

    root = Path(root or paths.state_root())

    def git(*a):
        return subprocess.run(["git", "-C", str(root), *a], capture_output=True, text=True)

    head = git("show", "HEAD:FROZEN.tsv")
    rows = {r["path"]: r for r in frozen.parse(head.stdout) if r["kind"] == "file"} if head.returncode == 0 else {}
    problems = []
    for rel in required:
        p = root / rel
        if not p.exists():
            problems.append(f"{rel}: missing")
            continue
        r = rows.get(rel)
        if r is None:
            problems.append(f"{rel}: no file row in the committed HEAD:FROZEN.tsv")
        elif r["sha256_or_value"] != _sha(p):
            problems.append(f"{rel}: SHA-256 differs from HEAD:FROZEN.tsv")
        ignored = git("check-ignore", "-q", "--", rel).returncode == 0
        if not ignored:
            if git("ls-files", "--error-unmatch", "--", rel).returncode != 0:
                problems.append(f"{rel}: not tracked in git")
            elif git("status", "--porcelain", "--", rel).stdout.strip():
                problems.append(f"{rel}: uncommitted changes")
    return problems


def check_frozen_full() -> list[str]:
    """scripts/check_frozen.py WITH the test-run check (§11.3)."""
    r = subprocess.run([sys.executable, str(PROJECT / "scripts" / "check_frozen.py")], capture_output=True, text=True)
    return [] if r.returncode == 0 else [f"check_frozen.py failed: {(r.stdout + r.stderr).strip()[:2000]}"]


# ----------------------------------------------------------------------------- label-free + training-only preparation
def prepare_look(cfg: dict, ctx: dict, train_long: pd.DataFrame) -> dict:
    """Fit on training, then compute every validation-line score matrix (label-free). Runs BEFORE ``started``."""
    from nextgem import w1

    fitted = w1.fit_components(train_long, ctx["line_class"])
    lines, genes = ctx["val_expr_lines"], ctx["genes"]
    b6 = ctx["b6"](lines, fitted) if (cfg["include_b6"] and ctx.get("b6")) else None
    if cfg["include_b6"] and b6 is None:
        raise RuntimeError("include_b6 is set but no B6 scorer is available")
    S = w1.score_matrices(lines, genes, fitted, ctx["gen"], ctx["line_class"], ctx["omics"]["expr"], b6=b6,
                          random_seeds=cfg["random_seeds"])
    for n, df in S.items():
        if not np.isfinite(df.to_numpy()).all():
            raise RuntimeError(f"non-finite scores for {n}")
    return {"fitted": fitted, "S": S}


# ----------------------------------------------------------------------------- label-dependent scoring (after started)
def w1_scalars(val_long: pd.DataFrame, prep: dict, ctx: dict, cfg: dict) -> tuple[pd.DataFrame, dict]:
    """Scalar rows only, from precomputed score matrices ``prep['S']``."""
    from nextgem import metrics, w1

    fitted, S = prep["fitted"], prep["S"]
    lines, genes, lineage = ctx["val_expr_lines"], ctx["genes"], ctx["lineage"]
    Yall = {lab: w1.wide(val_long, lab, lines, genes) for lab in cfg["labels"] + [cfg["binary_label"]]}
    rows = []
    B, seed = cfg["B"], cfg["seed"]
    for lab in cfg["labels"]:
        names = w1.methods_for(lab, S.keys())
        for gs in cfg["genesets"]:
            for u in cfg["universes"]:
                for st in cfg["strata"]:
                    Y = Yall[lab].loc[:, genes if u == "U" else ctx["ucap"]]
                    if gs == "noCE":
                        Y = Y.loc[:, [g for g in Y.columns if g not in fitted.ce[lab]]]
                    elif gs == "noDepMapCE":
                        Y = Y.loc[:, [g for g in Y.columns if g not in ctx["depmap_ce"]]]
                    Y = w1.apply_stratum(Y, ctx["strata"][st], ctx["lf"], ctx["omics"]["expr"], ctx["omics"]["cn"], None)
                    lts, cell = {}, []
                    for n in names:
                        r, lt, _ = w1.method_summary(S[n].loc[Y.index, Y.columns], Y, lineage, seed=seed, B=B)
                        lts[n] = lt
                        r.update(label=lab, geneset=gs, universe=u, stratum=st, method=n)
                        cell.append(r)
                    for ref in ("H1", "H2"):
                        for r in cell:
                            if r["method"] != ref:
                                d = w1.delta_summary(lts[r["method"]], lts[ref], lineage, seed=seed, B=B)
                                r.update({f"{k}_vs_{ref}": v for k, v in d.items()})
                    rows += cell
    grid = pd.DataFrame(rows)
    brows, perline = [], {}
    for s in ("H1", "H2"):
        pg = w1.parent_gene_set(ctx["gen"], s)
        for proto, rule, lab in BINARY_PROTOCOLS:
            col = {"P57": f"p57_{s}", "PGR": f"pgr_{s}", "FBA<0.01": f"fba_call_{s}"}[rule]
            summ, pl = w1.binary_endpoints(ctx["gen"][col], Yall[lab], lineage, pg, seed=seed, B=B)
            perline[(proto, rule, lab, s)] = pl
            brows.append(summ.assign(protocol=proto, method=s, call_rule=rule, label=lab, gene_set="parent"))
    for proto, rule, lab in PAIRED_PROTOCOLS:
        d = (perline[(proto, rule, lab, "H2")]["mcc"] - perline[(proto, rule, lab, "H1")]["mcc"]).dropna()
        if len(d):
            bt = metrics.bootstrap_mean(d, lineage.reindex(d.index), B=B, seed=seed)
            brows.append(pd.DataFrame([{"protocol": proto, "method": "H2 - H1", "call_rule": rule, "label": lab,
                                        "gene_set": "parent", "metric": "mcc", "nan_policy": "exclude",
                                        "lineage_equal_mean": bt.estimate, "ci_lo": bt.ci_low, "ci_hi": bt.ci_high,
                                        "n_lines_used": len(d)}]))
    binary = pd.concat(brows, ignore_index=True)
    out = pd.concat([grid.assign(section="tAP_grid"), binary.assign(section="binary_6.8")], ignore_index=True)
    out = out[["section"] + [c for c in out.columns if c != "section"]]
    gstar = decide_gstar(grid)
    assert_scalar_only(out, cfg, ctx, n_methods={lab: len(w1.methods_for(lab, S.keys())) for lab in cfg["labels"]})
    return out, {**gstar, "rule": cfg["gstar_rule"]}


def decide_gstar(grid: pd.DataFrame) -> dict:
    """DESIGN §6.3: G-star = the better generic parent on validation Δ̄_lin. Decided on the SIGN of
    Δ̄_lin(H2 − H1) (L_prob, all genes, U, no stratum): Human1 iff it is < 0; a tie (0) goes to Human2."""
    prim = grid[(grid.label == "L_prob") & (grid.geneset == "all") & (grid.universe == "U") & (grid.stratum == "none")]
    prim = prim.set_index("method")
    d = float(prim.loc["H2", "delta_lin_vs_H1"])
    return {"G_star": "Human1" if d < 0 else "Human2", "delta_lin_H2_vs_H1": d,
            "delta_ci": [float(prim.loc["H2", "delta_lo_vs_H1"]), float(prim.loc["H2", "delta_hi_vs_H1"])],
            "tap_lin_H1": float(prim.loc["H1", "tap_lin"]), "tap_lin_H2": float(prim.loc["H2", "tap_lin"])}


def assert_scalar_only(out: pd.DataFrame, cfg: dict, ctx: dict, n_methods: dict) -> None:
    """PI decision 21: a validation look reveals scalars only. Raises RuntimeError unless
    (a) every row is one scalar family (unique key per section) and no section exceeds its expected row count;
    (b) no column name or string cell is a line ID (ACH-...), a gene (Entrez or symbol) or a lineage name;
    (c) every numeric column is a known aggregate (whitelist ``ALLOWED_NUMERIC``)."""
    if "section" not in out.columns or not set(out.section) <= {"tAP_grid", "binary_6.8"}:
        raise RuntimeError("unknown output section")
    n_grid = sum(n_methods[lab] for lab in cfg["labels"]) * len(cfg["genesets"]) * len(cfg["universes"]) * len(cfg["strata"])
    from nextgem import metrics

    n_bin = 2 * len(BINARY_PROTOCOLS) * len(metrics.BINARY_METRICS) * 2 + len(PAIRED_PROTOCOLS)
    for sec, key, cap in (("tAP_grid", GRID_KEY, n_grid), ("binary_6.8", BINARY_KEY, n_bin)):
        sub = out[out.section == sec]
        if len(sub) > cap:
            raise RuntimeError(f"{sec}: {len(sub)} rows > expected {cap}")
        if sub[key].astype(str).duplicated().any():
            raise RuntimeError(f"{sec}: more than one row per scalar family")
    lineages = {str(x) for x in ctx["lineage"].dropna().unique()}
    genes = {str(g) for g in ctx["genes"]} | {str(s) for s in ctx.get("gene_symbols", ())}
    ids = {str(i) for i in ctx["all_ids"]}

    def bad(x: str) -> str | None:
        if ID_RE.match(x) or x in ids:
            return "line ID"
        if x in lineages:
            return "lineage name"
        if x in genes:
            return "gene"
        return None

    for c in out.columns:
        why = bad(str(c))
        if why:
            raise RuntimeError(f"column {c!r} is a {why}")
        if out[c].dtype == object:
            for v in out[c].dropna().astype(str).unique():
                why = bad(v)
                if why:
                    raise RuntimeError(f"column {c!r} holds a {why} value")
        elif c not in ALLOWED_NUMERIC:
            raise RuntimeError(f"numeric column {c!r} is not a known scalar aggregate (per-unit columns are refused)")


# ----------------------------------------------------------------------------- the look
def write_config(cfg: dict, out_dir: Path) -> Path:
    path = out_dir / f"val_look_{cfg['look_id']}_config.json"
    out_dir.mkdir(parents=True, exist_ok=True)
    doc = {k: v for k, v in cfg.items() if not k.startswith("_")}
    path.write_text(json.dumps({**doc, "config_sha": config_sha(cfg)}, indent=1, sort_keys=True))
    return path


def run_look(cfg: dict, *, prep: dict, load_val, ctx: dict, out_dir: Path, author: str) -> dict:
    """Ledger part. Everything label-free (``prep``) and the config JSON exist before this is called."""
    from nextgem import ledger

    csha = config_sha(cfg)
    cpath = out_dir / f"val_look_{cfg['look_id']}_config.json"
    if not cpath.exists() or json.loads(cpath.read_text()).get("config_sha") != csha:
        raise RuntimeError("config JSON missing or stale: write_config() must run before the look")
    prior = ledger.looks(ledger.read("validation")).get(cfg["look_id"])
    resumed = prior is not None and prior["started"] is not None and not prior["completed"]
    detail = json.dumps({"milestone": cfg["milestone"], "scalars_only": True, "include_b6": cfg["include_b6"]}, sort_keys=True)
    with ledger.look("validation", cfg["look_id"], workstream=cfg["workstream"], config_sha=csha, detail=detail,
                     author=author) as lk:
        if resumed:
            n_prev = sum(1 for r in ledger.read("validation") if r["look_id"] == cfg["look_id"] and r["event"] == "output"
                         and '"resumed"' in r["detail"])
            marker = hashlib.sha256(f"{csha}|resumed|{n_prev + 1}".encode()).hexdigest()
            lk.output(marker, detail=json.dumps({"event": "resumed", "attempt": n_prev + 1}), author=author)
        val_long = load_val()  # first and only validation-label read
        out, gstar = w1_scalars(val_long, prep, ctx, cfg)
        path = out_dir / f"val_look_{cfg['look_id']}.tsv"
        out.to_csv(path, sep="\t", index=False, float_format="%.6g")
        gpath = out_dir / f"val_look_{cfg['look_id']}_gstar.json"
        gpath.write_text(json.dumps(gstar, indent=1))
        h = hashlib.sha256(path.read_bytes() + gpath.read_bytes()).hexdigest()
        lk.output(h, detail=json.dumps({"files": [path.name, gpath.name]}), author=author)
    return {"output": str(path), "gstar_file": str(gpath), "gstar": gstar, "config_sha": csha, "output_sha": h,
            "resumed": resumed}


def real_context(cfg: dict) -> dict:
    """Label-free context on the real data (runs before ``started``)."""
    from nextgem import labels, paths, w1

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    val = sorted(split.index[split.tier == "validation"])
    omics = w1.load_omics()
    gf = w1.gene_features()
    syms = set()
    for p in ("human1", "human2"):
        syms |= set(pd.read_csv(paths.baseline_dir() / f"{p}_gene_mapping.tsv", sep="\t", usecols=["symbol"]).symbol.dropna())
    ctx = {"genes": genes, "ucap": [g for g, a, b in zip(U.entrez, U.in_human1, U.in_human2) if a and b],
           "gen": w1.generic_table(U.entrez), "omics": omics, "lf": lf, "lineage": split.lineage,
           "line_class": lf.base_medium_class, "val_lines": val,
           "val_expr_lines": [l for l in val if lf.loc[l, "has_expression_matrix"]], "all_ids": list(split.index),
           "gene_symbols": syms,
           "depmap_ce": labels.label_derived_list("depmap_inferred_common_essentials", purpose="reporting"),
           "strata": w1.stratum_specs(w1.gene_panels(), gf), "b6": None}
    if cfg["include_b6"]:
        root = paths.state_root()
        nb = json.loads((root / cfg["_inputs_rel"]["gene_neighbours"]).read_text())
        model = pickle.loads((root / cfg["_inputs_rel"]["b6_full"]).read_bytes())
        aggs = w1.expression_aggregates(omics["expr"], genes, nb["neighbours_h2"], nb["paralogs"])

        def b6_scores(lines, fitted):
            b2, b2m = w1.heldout_priors(fitted, lines, ctx["line_class"], genes)
            X = w1.b6_matrix(lines, genes, b2, b2m, omics, aggs, gf, lf, model["class_codes"])
            mi = w1.B6_FEATURES.index("medium_class")
            X[:, mi] = np.where(X[:, mi] < 0, np.nan, X[:, mi])
            p = model["clf"].predict_proba(X)[:, 1].reshape(len(lines), len(genes))
            return pd.DataFrame(np.round(p, 8), index=lines, columns=genes)

        ctx["b6"] = b6_scores
    return ctx


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--look", required=True, choices=["W1_M2_baseline"])
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--print-config", action="store_true")
    g.add_argument("--print-frozen-candidates", action="store_true")
    g.add_argument("--preflight-only", action="store_true")
    a = ap.parse_args(argv)
    t0 = time.time()
    cfg = w1_config()
    if a.print_config:
        print(json.dumps({**{k: v for k, v in cfg.items() if not k.startswith("_")}, "config_sha": config_sha(cfg)}, indent=1))
        return 0
    if a.print_frozen_candidates:
        print("\n".join(frozen_candidates(cfg)))
        return 0
    problems = preflight(required_files(cfg)) + check_frozen_full()
    if problems:
        print("PRE-FLIGHT FAILED; the look was not started (nothing written to the ledger):", file=sys.stderr)
        print("\n".join(f"  - {p}" for p in problems), file=sys.stderr)
        return 2
    if a.preflight_only:
        print("pre-flight OK")
        return 0
    from nextgem import labels, paths
    from nextgem.run import provenance

    out_dir = paths.results_dir() / "w1"
    root = paths.state_root()
    with provenance(out_dir / "MANIFEST.tsv", inputs=[root / r for r in required_files(cfg)], seeds=[cfg["seed"]],
                    allow_dirty=False, author="scripts/evaluate_val.py") as prov:
        ctx = real_context(cfg)
        prep = prepare_look(cfg, ctx, labels.train(columns=LABEL_COLS))
        cpath = write_config(cfg, out_dir)
        res = run_look(cfg, prep=prep, load_val=labels.validation, ctx=ctx, out_dir=out_dir, author="scripts/evaluate_val.py")
        prov.output(cpath, res["output"], res["gstar_file"])
        prov.note(look_id=cfg["look_id"], config_sha=res["config_sha"], output_sha=res["output_sha"])
    print(json.dumps({**res, "wall_s": round(time.time() - t0, 1)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
