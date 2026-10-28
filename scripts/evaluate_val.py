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

Looks implemented:

- ``W1_M2_baseline`` (M2, W1). Generic Human1/Human2 and B0–B3, B8±/B8′ (B6 is excluded:
  ``include_b6 = False``; it is untuned and W7 owns its validation looks) on the validation expression lines:
  tAP grid scalars (L_prob/L_ge05/L_ge10 x {all, minus CE_fit, minus DepMap CE} x {U, U∩} x {none, expressed
  pairs}) with paired Δ̄_lin vs H1 and vs H2; §6.8 binary endpoints for generic H1/H2 and the paired H2 − H1 MCC;
  and G-star: **Human1 iff Δ̄_lin(H2 − H1) < 0** (L_prob, all genes, U), otherwise Human2 (a tie goes to Human2,
  the default parent, DESIGN §3.1). CE_fit, B2/B2m/B2c come from all training lines (§5.3).
- **shortlist looks** (A6: a new look = a reviewed, frozen look spec, no harness change). A spec
  ``results/<ws>/val/look_spec_<look_id>.json`` (``kind: "shortlist"``; see :func:`shortlist_config`) names a
  label-free scores parquet pinned by SHA-256 whose columns are the look's methods (line-constant ``entrez`` x
  method, or ``ModelID`` + ``entrez`` x method). They are scored next to generic H1/H2 with the W1 machinery
  (same tAP grid columns, paired Δ̄_lin vs H1 and vs H2, scalar guard, ledger). The decision (§8/§9) is written
  to ``val_look_<look_id>_decision.json``: a candidate is eligible iff its paired Δ̄_lin vs the reference (which
  must be G-star from the W1 look) has CI lower bound > 0 on the primary cell (first label, all genes, U, no
  stratum); among eligible candidates the paired-difference one-SE rule picks the first candidate (the spec
  lists candidates simplest first) whose Δ̄_lin(best − candidate) <= its bootstrap SE. The §9 winner's-curse
  term σ·√(2 ln m) of the best candidate is reported alongside (it does not change the gate). Spec rules
  (A11g): the workstream equals the spec's ``results/<ws>/val`` directory (look caps), labels start with L_prob,
  B and seed are the integers B_OTHER and SEEDS[0], the scores file lives in ``results/<ws>/val/``, and a
  look_id has exactly one spec.

Usage (``<look>`` = W1_M2_baseline or a spec's look_id):
  .venv/bin/python scripts/evaluate_val.py --look <look>                       # the look (once)
  .venv/bin/python scripts/evaluate_val.py --look <look> --print-config        # reads no label
  .venv/bin/python scripts/evaluate_val.py --look <look> --print-frozen-candidates   # reads no label
  .venv/bin/python scripts/evaluate_val.py --look <look> --preflight-only      # reads no label
  .venv/bin/python scripts/evaluate_val.py --look <look> --dry-run             # synthetic validation labels;
                                                       # reads training labels only, writes nothing (A6)
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
                      "tests/test_run.py", "tests/test_evaluate_val_shortlist.py")
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
W1_LOOK = "W1_M2_baseline"
W1_GSTAR = "results/w1/val_look_W1_M2_baseline_gstar.json"  # G-star decided by the W1 look
SPEC_GLOB = "results/*/val/look_spec_*.json"
GSTAR_SHORT = {"Human1": "H1", "Human2": "H2"}
# a method name or scores column that looks like a label, a line key or a raw DepMap label source is refused
LABEL_LIKE_RE = re.compile(r"^(L_\w+|gene_effect|dep_prob|lineage|tier|ModelID|entrez|label\w*)$|CRISPR|Chronos|"
                           r"GeneEffect|Dependency", re.I)
SPEC_KEYS = {"look_id", "kind", "workstream", "milestone", "scores", "methods", "candidates", "reference", "labels",
             "genesets", "universes", "strata", "B", "seed", "decision_rule", "provenance", "notes"}


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
    return list(CODE_FILES) + list(cfg["_inputs_rel"].values()) + list(cfg.get("provenance", {}))


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


# ----------------------------------------------------------------------------- shortlist looks (look specs)
def look_specs() -> dict[str, Path]:
    """Registered shortlist looks: look_id -> ``results/<ws>/val/look_spec_<look_id>.json`` under the state root.
    A look_id with a spec in more than one workstream directory is refused (PreflightError)."""
    from nextgem import paths

    specs = sorted(paths.state_root().glob(SPEC_GLOB))
    ids = [p.stem.removeprefix("look_spec_") for p in specs]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        raise PreflightError(f"look_id(s) {dup} have more than one look spec ({SPEC_GLOB})")
    return dict(zip(ids, specs))


def check_spec(spec: dict, look_id: str, ws_dir: str) -> None:
    """Static, label-free checks of a look spec found in ``results/<ws_dir>/val/``. Raises PreflightError."""
    from nextgem import metrics, w1
    from nextgem.config import SEEDS

    def need(ok, msg):
        if not ok:
            raise PreflightError(f"look spec {look_id}: {msg}")

    def under(p, *tops):  # a relative path, without '..', strictly below one of the given directories
        if not isinstance(p, str) or Path(p).is_absolute() or ".." in Path(p).parts:
            return False
        parts = Path(p).parts
        return any(len(parts) > len(t) and parts[:len(t)] == t for t in tops)

    need(SPEC_KEYS - {"provenance", "notes"} <= set(spec) <= SPEC_KEYS, f"keys must be {sorted(SPEC_KEYS)} "
         "(provenance and notes optional)")
    need(spec["look_id"] == look_id != W1_LOOK and spec["kind"] == "shortlist",
         "look_id must match the file name (and not be the W1 look); kind must be 'shortlist'")
    ws = spec["workstream"]
    need(isinstance(ws, str) and re.fullmatch(r"[WM][1-9]\d*", ws) is not None and ws.lower() == ws_dir
         and (ws[0] == "W" or spec.get("milestone") == ws),
         f"workstream must be W<n> or M<n> (a milestone-level look whose milestone equals it, A13), without leading "
         f"zeros, and match the spec's directory results/{ws_dir}/val (DESIGN §17 A11g; look caps)")
    methods, cands = spec["methods"], spec["candidates"]
    need(isinstance(methods, list) and len(set(methods)) == len(methods) and {"H1", "H2"} <= set(methods),
         "methods must be a list of unique names that includes H1 and H2")
    for m in methods:
        need(isinstance(m, str) and m.strip() == m and m and not ID_RE.match(m) and not m.isdigit()
             and not LABEL_LIKE_RE.search(m), f"method name {m!r} is refused (line-ID-, gene- or label-like)")
    need(isinstance(cands, list) and 1 <= len(cands) <= 3 and len(set(cands)) == len(cands)
         and set(cands) <= set(methods) - {"H1", "H2"},
         "candidates: 1-3 distinct non-generic methods, listed simplest first (DESIGN §9)")
    need(spec["reference"] in ("H1", "H2"), "reference must be H1 or H2 (G-star)")
    labs = spec["labels"]
    need(isinstance(labs, list) and labs and labs[0] == "L_prob" and set(labs) <= set(w1.GRID_LABELS),
         f"labels must start with the primary label L_prob (§8 gate cell, A11g) and be a subset of {w1.GRID_LABELS}")
    for k, allowed, prim in (("genesets", W1_GENESETS, "all"), ("universes", W1_UNIVERSES, "U"), ("strata", W1_STRATA, "none")):
        need(prim in spec[k] and set(spec[k]) <= set(allowed), f"{k} must include {prim!r} and be a subset of {allowed}")
    need(type(spec["B"]) is int and spec["B"] == metrics.B_OTHER and type(spec["seed"]) is int and spec["seed"] == SEEDS[0],
         f"B must be the integer {metrics.B_OTHER} (§6.5) and seed the integer SEEDS[0] = {SEEDS[0]} (§7)")
    sc = spec["scores"]
    sha_ok = lambda s: isinstance(s, str) and re.fullmatch(r"[0-9a-f]{64}", s) is not None  # noqa: E731
    need(isinstance(sc, dict) and set(sc) == {"path", "sha256", "layout"} and sc["layout"] in ("gene", "line_gene")
         and sha_ok(sc["sha256"]), "scores must be {path, sha256 (hex), layout: gene | line_gene}")
    need(under(sc["path"], ("results", ws_dir, "val")) and len(Path(sc["path"]).parts) == 4
         and sc["path"].endswith(".parquet"), f"scores.path must be a .parquet file in results/{ws_dir}/val/")
    prov = spec.get("provenance", {})
    need(isinstance(prov, dict) and all(under(p, ("results",), ("src",), ("scripts",)) and sha_ok(v) for p, v in prov.items()),
         "provenance must map paths under results/, src/ or scripts/ to SHA-256 pins")


def load_scores(path: Path, layout: str, methods) -> pd.DataFrame:
    """Read and check a look's scores parquet (label-free). Only key columns (``entrez``; plus ``ModelID`` for
    ``line_gene``) and one float score column per spec method (H1 and H2 included; they must equal the W1 look's
    generic scores) are allowed: no more, no fewer. Raises PreflightError."""
    def refuse(msg):
        raise PreflightError(f"scores {path.name}: {msg}")

    df = pd.read_parquet(path)
    keys = ["entrez"] if layout == "gene" else ["ModelID", "entrez"]
    if [k for k in keys if k not in df.columns]:
        refuse(f"key columns must be {keys}")
    extra = [str(c) for c in df.columns if c not in keys and c not in methods]
    if extra:
        refuse(f"non-score column(s) {extra}: only the spec's methods are allowed")
    absent = [m for m in methods if m not in df.columns]
    if absent:
        refuse(f"no scores column for method(s) {absent}")
    if not pd.api.types.is_integer_dtype(df["entrez"]):
        refuse("entrez must be an integer column")
    if df.duplicated(keys).any():
        refuse(f"duplicate {keys} rows")
    for c in (c for c in df.columns if c not in keys):
        if not pd.api.types.is_float_dtype(df[c]) or not np.isfinite(df[c].to_numpy()).all():
            refuse(f"column {c!r} must be finite floats")
        if layout == "line_gene" and set(np.unique(df[c])) <= {0.0, 1.0}:
            refuse(f"column {c!r} is a binary line x gene matrix (label-shaped); submit continuous scores")
    if layout == "line_gene" and not df.ModelID.map(lambda x: isinstance(x, str) and ID_RE.match(x) is not None).all():
        refuse("ModelID values must be line IDs")
    return df


def shortlist_config(look_id: str) -> dict:
    """Config of a shortlist look from its spec; its SHA is the ledger config_sha. Label-free. The scores file
    must carry the SHA-256 pinned in the spec and pass :func:`load_scores`, and the reference must be the G-star
    decided by the W1 look, or PreflightError is raised before anything else (no label is read)."""
    from nextgem import paths

    root = paths.state_root()
    spath = look_specs().get(look_id)
    if spath is None:
        raise PreflightError(f"no look spec for {look_id!r} ({SPEC_GLOB})")
    spec = json.loads(spath.read_text())
    check_spec(spec, look_id, spath.parent.parent.name)
    sc = spec["scores"]
    sp = root / sc["path"]
    if not sp.is_file() or _sha(sp) != sc["sha256"]:
        raise PreflightError(f"{sc['path']}: missing or SHA-256 differs from the pin in {spath.name}")
    load_scores(sp, sc["layout"], spec["methods"])
    try:
        gstar = json.loads((root / W1_GSTAR).read_text())["G_star"]
    except (OSError, ValueError, KeyError, TypeError) as e:
        raise PreflightError(f"G-star file {W1_GSTAR} missing or unreadable ({type(e).__name__})") from None
    if gstar not in GSTAR_SHORT:
        raise PreflightError(f"G-star file {W1_GSTAR}: G_star must be one of {sorted(GSTAR_SHORT)}")
    if spec["reference"] != GSTAR_SHORT[gstar]:
        raise PreflightError(f"reference {spec['reference']!r} is not G-star = {gstar} ({W1_GSTAR})")
    rel = {k: v for k, v in input_files(False).items() if k != "sesoi_c1b"}  # context of the W1 look (real_context)
    rel.update(look_spec=str(spath.relative_to(root)), scores=sc["path"], w1_gstar=W1_GSTAR)
    cfg = {k: spec[k] for k in ("look_id", "kind", "workstream", "milestone", "labels", "genesets", "universes",
                                "strata", "B", "seed", "methods", "candidates", "reference", "decision_rule")}
    cfg.update({"scores_layout": sc["layout"], "G_star": gstar, "provenance": spec.get("provenance", {}),
                "notes": spec.get("notes", ""), "inputs": {k: _sha(root / v) for k, v in rel.items()},
                "code": {Path(c).name: _sha(paths.PROJECT_ROOT / c) for c in CODE_FILES}, "_inputs_rel": rel})
    return cfg


def look_config(look_id: str) -> dict:
    return w1_config() if look_id == W1_LOOK else shortlist_config(look_id)


def look_out_dir(cfg: dict) -> Path:
    """W1: results/w1; a shortlist look writes next to its spec."""
    from nextgem import paths

    if cfg.get("kind") == "shortlist":
        return paths.state_root() / Path(cfg["_inputs_rel"]["look_spec"]).parent
    return paths.results_dir() / "w1"


def pinned_problems(cfg: dict) -> list[str]:
    """Provenance files of a look spec must still carry their pinned SHA-256."""
    from nextgem import paths

    root = paths.state_root()
    return [f"{p}: SHA-256 differs from the look spec's provenance pin" for p, s in cfg.get("provenance", {}).items()
            if (root / p).is_file() and _sha(root / p) != s]


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
        if not p.is_file():
            rows.append(f"# MISSING\t{rel}\t(no row possible until the file exists)")
            continue
        if have.get(rel) == _sha(p):
            continue
        ws = cfg["workstream"]
        note = f"{ws} validation look: code" if rel in CODE_FILES else (
            "known-answer test (DESIGN 11.3)" if rel in KNOWN_ANSWER_TESTS else f"{ws} validation look: input")
        rows.append("\t".join([rel, "file", rel, _sha(p), day, "", note]))
    if cfg.get("kind") == "shortlist":
        return rows
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


def prepare_shortlist(cfg: dict, ctx: dict, train_long: pd.DataFrame) -> dict:
    """Shortlist look: CE_fit from all training lines (§5.3), generic H1/H2 exactly as in the W1 look, and the
    spec's methods from the pinned scores file (re-checked here). Label-free; runs BEFORE ``started``, so every
    label-free error that the scalar guard would otherwise raise after ``started`` (which burns the look id) is
    raised here: method names that are line IDs, lineage names or genes, and non-integer B / seed."""
    from nextgem import paths, w1

    bad = unit_key_checker(ctx)
    for m in cfg["methods"]:
        if bad(m):
            raise PreflightError(f"method name {m!r} is a {bad(m)}; the scalar guard would refuse the output")
    if type(cfg["B"]) is not int or type(cfg["seed"]) is not int:
        raise PreflightError("B and seed must be integers")
    sp = paths.state_root() / cfg["_inputs_rel"]["scores"]
    if _sha(sp) != cfg["inputs"]["scores"]:
        raise RuntimeError("scores file changed since the config was built")
    df = load_scores(sp, cfg["scores_layout"], cfg["methods"])
    fitted = w1.fit_components(train_long, ctx["line_class"])
    lines, genes = list(ctx["val_expr_lines"]), list(ctx["genes"])
    gen = w1.score_matrices(lines, genes, fitted, ctx["gen"], ctx["line_class"], ctx["omics"]["expr"], random_seeds=())
    S = {s: gen[s] for s in ("H1", "H2")}
    if set(df.entrez) != set(genes):
        raise RuntimeError(f"scores cover {df.entrez.nunique()} genes; the look needs exactly the {len(genes)} genes of U")
    if cfg["scores_layout"] == "line_gene" and set(df.ModelID) != set(lines):
        raise RuntimeError("line_gene scores must cover exactly the validation expression lines")
    for m in cfg["methods"]:
        if cfg["scores_layout"] == "gene":
            v = df.set_index("entrez")[m].reindex(genes).to_numpy(dtype=float)
            M = pd.DataFrame(np.tile(v, (len(lines), 1)), index=lines, columns=genes)
        else:
            M = df.pivot(index="ModelID", columns="entrez", values=m).reindex(index=lines, columns=genes).astype(float)
        if m in S:
            if not np.allclose(M.to_numpy(), S[m].to_numpy(), rtol=0, atol=1e-9):
                raise RuntimeError(f"scores column {m} differs from the generic {m} of the W1 look")
            continue
        S[m] = M
    for n, M in S.items():
        if not np.isfinite(M.to_numpy()).all():
            raise RuntimeError(f"non-finite scores for {n}")
    return {"fitted": fitted, "S": S}


# ----------------------------------------------------------------------------- label-dependent scoring (after started)
def tap_grid(Yall: dict, S: dict, names_for, fitted, ctx: dict, cfg: dict) -> tuple[list[dict], dict]:
    """tAP grid rows (label x geneset x universe x stratum x method) with paired Δ vs H1 and vs H2, and the
    per-line tables of every cell, keyed (label, geneset, universe, stratum)."""
    from nextgem import w1

    genes, lineage = ctx["genes"], ctx["lineage"]
    rows, cells = [], {}
    B, seed = cfg["B"], cfg["seed"]
    for lab in cfg["labels"]:
        names = names_for(lab)
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
                    cells[(lab, gs, u, st)] = lts
    return rows, cells


def w1_scalars(val_long: pd.DataFrame, prep: dict, ctx: dict, cfg: dict) -> tuple[pd.DataFrame, dict]:
    """Scalar rows only, from precomputed score matrices ``prep['S']``."""
    from nextgem import metrics, w1

    fitted, S = prep["fitted"], prep["S"]
    lines, genes, lineage = ctx["val_expr_lines"], ctx["genes"], ctx["lineage"]
    Yall = {lab: w1.wide(val_long, lab, lines, genes) for lab in cfg["labels"] + [cfg["binary_label"]]}
    B, seed = cfg["B"], cfg["seed"]
    rows, _ = tap_grid(Yall, S, lambda lab: w1.methods_for(lab, S.keys()), fitted, ctx, cfg)
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


def shortlist_scalars(val_long: pd.DataFrame, prep: dict, ctx: dict, cfg: dict) -> tuple[pd.DataFrame, dict]:
    """Shortlist look: the tAP grid of the spec's methods (no binary section) and the §8/§9 decision."""
    from nextgem import w1

    lines, genes = ctx["val_expr_lines"], ctx["genes"]
    Yall = {lab: w1.wide(val_long, lab, lines, genes) for lab in cfg["labels"]}
    rows, cells = tap_grid(Yall, prep["S"], lambda lab: list(cfg["methods"]), prep["fitted"], ctx, cfg)
    grid = pd.DataFrame(rows)
    out = grid.assign(section="tAP_grid")
    out = out[["section"] + [c for c in out.columns if c != "section"]]
    decision = decide_shortlist(grid, cells[(cfg["labels"][0], "all", "U", "none")], ctx["lineage"], cfg)
    assert_scalar_only(out, cfg, ctx, n_methods={lab: len(cfg["methods"]) for lab in cfg["labels"]}, n_bin=0)
    return out, decision


def decide_shortlist(grid: pd.DataFrame, lts: dict, lineage: pd.Series, cfg: dict) -> dict:
    """DESIGN §8 gate and §9 paired-difference one-SE rule on the primary cell (first label, all genes, U, no
    stratum). Eligible: CI lower bound of the paired Δ̄_lin vs the reference (G-star) > 0. Among eligible
    candidates, best = the largest Δ̄_lin; selected = the first candidate in spec order (simplest first) with
    Δ̄_lin(best − candidate) <= SE, the SD of its lineage-stratified bootstrap replicates. The winner's-curse
    term of §9 is added for reporting (:func:`winners_curse`)."""
    from nextgem import metrics

    lab, ref = cfg["labels"][0], cfg["reference"]
    prim = grid[(grid.label == lab) & (grid.geneset == "all") & (grid.universe == "U") & (grid.stratum == "none")]
    prim = prim.set_index("method")
    cand = []
    for c in cfg["candidates"]:
        d, lo, hi = (float(prim.loc[c, f"delta_{k}_vs_{ref}"]) for k in ("lin", "lo", "hi"))
        cand.append({"method": c, "delta_lin": d, "delta_ci": [lo, hi], "eligible": bool(lo > 0)})
    elig = [r for r in cand if r["eligible"]]
    out = {"look_id": cfg["look_id"], "reference": ref, "G_star": cfg["G_star"],
           "primary_cell": {"label": lab, "geneset": "all", "universe": "U", "stratum": "none"},
           "candidates": cand, "gate_passed": bool(elig), "best_eligible": None, "selected": None}
    if elig:
        best = max(elig, key=lambda r: r["delta_lin"])  # a tie keeps the earlier (simpler) candidate
        for r in elig:
            d = metrics.paired_line_delta(lts[best["method"]], lts[r["method"]])
            bt = metrics.bootstrap_mean(d, lineage.reindex(d.index), B=cfg["B"], seed=cfg["seed"])
            se = float(np.std(bt.replicates, ddof=1))
            r.update(diff_vs_best=bt.estimate, se_diff=se, within_1se=bool(bt.estimate <= se))
        out["best_eligible"] = best["method"]
        out["selected"] = next(r["method"] for r in elig if r["within_1se"])
    top = max(cand, key=lambda r: -np.inf if np.isnan(r["delta_lin"]) else r["delta_lin"])
    if not elig:
        out["best_point_estimate"] = top["method"]
    out["winners_curse"] = winners_curse(top, lts, lineage, cfg)
    out["rule"] = cfg["decision_rule"]
    return out


def winners_curse(top: dict, lts: dict, lineage: pd.Series, cfg: dict) -> dict:
    """DESIGN §9, reporting only (the gate is unchanged): expected optimism σ·√(2 ln m) of the best candidate
    (largest Δ̄_lin vs the reference), m = the spec's number of candidates, σ = the SD of the lineage-stratified
    bootstrap replicates of that candidate's paired Δ̄_lin (same B and seed as the grid)."""
    from nextgem import metrics

    m = len(cfg["candidates"])
    d = metrics.paired_line_delta(lts[top["method"]], lts[cfg["reference"]])
    sigma = float(np.std(metrics.bootstrap_mean(d, lineage.reindex(d.index), B=cfg["B"], seed=cfg["seed"]).replicates, ddof=1)) \
        if len(d.dropna()) else float("nan")
    term = 0.0 if m == 1 else sigma * float(np.sqrt(2 * np.log(m)))  # one candidate: no selection, no optimism
    return {"method": top["method"], "m": m, "sigma": sigma, "term": term, "best_delta_lin": top["delta_lin"],
            "best_delta_lin_minus_term": top["delta_lin"] - term}


def unit_key_checker(ctx: dict):
    """bad(x) -> 'line ID' | 'lineage name' | 'gene' if the string x names a unit, else None (scalar guard)."""
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

    return bad


def assert_scalar_only(out: pd.DataFrame, cfg: dict, ctx: dict, n_methods: dict, n_bin: int | None = None) -> None:
    """PI decision 21: a validation look reveals scalars only. Raises RuntimeError unless
    (a) every row is one scalar family (unique key per section) and no section exceeds its expected row count;
    (b) no column name or string cell is a line ID (ACH-...), a gene (Entrez or symbol) or a lineage name;
    (c) every numeric column is a known aggregate (whitelist ``ALLOWED_NUMERIC``).
    ``n_bin`` = expected binary_6.8 rows (default: the W1 look's; 0 for shortlist looks, which have none)."""
    if "section" not in out.columns or not set(out.section) <= {"tAP_grid", "binary_6.8"}:
        raise RuntimeError("unknown output section")
    n_grid = sum(n_methods[lab] for lab in cfg["labels"]) * len(cfg["genesets"]) * len(cfg["universes"]) * len(cfg["strata"])
    from nextgem import metrics

    if n_bin is None:
        n_bin = 2 * len(BINARY_PROTOCOLS) * len(metrics.BINARY_METRICS) * 2 + len(PAIRED_PROTOCOLS)
    for sec, key, cap in (("tAP_grid", GRID_KEY, n_grid), ("binary_6.8", BINARY_KEY, n_bin)):
        sub = out[out.section == sec]
        if len(sub) > cap:
            raise RuntimeError(f"{sec}: {len(sub)} rows > expected {cap}")
        if cap == 0:
            continue
        if sub[key].astype(str).duplicated().any():
            raise RuntimeError(f"{sec}: more than one row per scalar family")
    bad = unit_key_checker(ctx)
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
    _, scalars, side = look_kind(cfg)
    if cfg.get("kind") == "shortlist":
        detail = json.dumps({"milestone": cfg["milestone"], "scalars_only": True, "kind": "shortlist"}, sort_keys=True)
    else:
        detail = json.dumps({"milestone": cfg["milestone"], "scalars_only": True, "include_b6": cfg["include_b6"]}, sort_keys=True)
    with ledger.look("validation", cfg["look_id"], workstream=cfg["workstream"], config_sha=csha, detail=detail,
                     author=author) as lk:
        if resumed:
            n_prev = sum(1 for r in ledger.read("validation") if r["look_id"] == cfg["look_id"] and r["event"] == "output"
                         and '"resumed"' in r["detail"])
            marker = hashlib.sha256(f"{csha}|resumed|{n_prev + 1}".encode()).hexdigest()
            lk.output(marker, detail=json.dumps({"event": "resumed", "attempt": n_prev + 1}), author=author)
        val_long = load_val()  # first and only validation-label read
        out, gstar = scalars(val_long, prep, ctx, cfg)
        path = out_dir / f"val_look_{cfg['look_id']}.tsv"
        out.to_csv(path, sep="\t", index=False, float_format="%.6g")
        gpath = out_dir / f"val_look_{cfg['look_id']}_{side}.json"
        gpath.write_text(json.dumps(gstar, indent=1))
        h = hashlib.sha256(path.read_bytes() + gpath.read_bytes()).hexdigest()
        lk.output(h, detail=json.dumps({"files": [path.name, gpath.name]}), author=author)
    return {"output": str(path), f"{side}_file": str(gpath), side: gstar, "config_sha": csha, "output_sha": h,
            "resumed": resumed}


def look_kind(cfg: dict):
    """(prepare, scalars, side-file name) of the look's kind: the W1 look (G-star) or a shortlist look (decision)."""
    if cfg.get("kind") == "shortlist":
        return prepare_shortlist, shortlist_scalars, "decision"
    return prepare_look, w1_scalars, "gstar"


def synthetic_validation(lines, genes, rates: dict, seed) -> pd.DataFrame:
    """Dry-run labels: independent Bernoulli(rate) per (line, gene) and label; no label file is read."""
    rng = np.random.default_rng(seed)
    df = pd.MultiIndex.from_product([list(lines), list(genes)], names=["ModelID", "entrez"]).to_frame(index=False)
    for lab, p in rates.items():
        df[lab] = (rng.random(len(df)) < p).astype(float)
    return df


def dry_run(cfg: dict, ctx: dict, train_long: pd.DataFrame) -> dict:
    """A6 label-free dry run: the look's preparation and scalars on SYNTHETIC validation labels (Bernoulli at
    each label's training positive rate, SEEDS[0]). No ledger row, no file; labels.validation() is never called."""
    prepare, scalars, side = look_kind(cfg)
    prep = prepare(cfg, ctx, train_long)
    labs = cfg["labels"] + ([cfg["binary_label"]] if "binary_label" in cfg else [])
    rates = {lab: prep["fitted"].rate[lab] for lab in labs}
    out, res = scalars(synthetic_validation(ctx["val_expr_lines"], ctx["genes"], rates, cfg["seed"]), prep, ctx, cfg)
    # the training positive rates are not reported: dry-run readers stay at synthetic-only exposure
    return {"look_id": cfg["look_id"], "config_sha": config_sha(cfg),
            "synthetic_labels": f"Bernoulli at each label's training positive rate (not shown), seed {cfg['seed']}",
            "n_val_expr_lines": len(ctx["val_expr_lines"]), "methods": sorted(prep["S"]),
            "n_rows": len(out), "scalar_guard": "pass", "validation_labels_read": False, "ledger_rows_written": 0,
            f"synthetic_{side}": res}


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
    if cfg.get("include_b6"):
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
    ap.add_argument("--look", required=True, help=f"{W1_LOOK} or the look_id of a look spec ({SPEC_GLOB})")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--print-config", action="store_true")
    g.add_argument("--print-frozen-candidates", action="store_true")
    g.add_argument("--preflight-only", action="store_true")
    g.add_argument("--dry-run", action="store_true", help="synthetic validation labels; no ledger row, no file")
    a = ap.parse_args(argv)
    t0 = time.time()
    try:
        looks = [W1_LOOK, *look_specs()]
        if a.look not in looks:
            ap.error(f"--look must be one of {looks}")
        cfg = look_config(a.look)
    except PreflightError as e:
        print(f"LOOK REFUSED (no label read, nothing written): {e}", file=sys.stderr)
        return 2
    if a.print_config:
        print(json.dumps({**{k: v for k, v in cfg.items() if not k.startswith("_")}, "config_sha": config_sha(cfg)}, indent=1))
        return 0
    if a.print_frozen_candidates:
        print("\n".join(frozen_candidates(cfg)))
        return 0
    if a.dry_run:
        from nextgem import labels

        try:
            res = dry_run(cfg, real_context(cfg), labels.train(columns=LABEL_COLS))
        except PreflightError as e:
            print(f"DRY RUN REFUSED (no validation label read, nothing written): {e}", file=sys.stderr)
            return 2
        print(json.dumps({**res, "wall_s": round(time.time() - t0, 1)}, indent=1))
        return 0
    problems = preflight(required_files(cfg)) + pinned_problems(cfg) + check_frozen_full()
    if problems:
        print("PRE-FLIGHT FAILED; the look was not started (nothing written to the ledger):", file=sys.stderr)
        print("\n".join(f"  - {p}" for p in problems), file=sys.stderr)
        return 2
    if a.preflight_only:
        print("pre-flight OK")
        return 0
    from nextgem import labels, paths
    from nextgem.run import provenance

    out_dir = look_out_dir(cfg)
    root = paths.state_root()
    prepare, _, side = look_kind(cfg)
    with provenance(out_dir / "MANIFEST.tsv", inputs=[root / r for r in required_files(cfg)], seeds=[cfg["seed"]],
                    allow_dirty=False, author="scripts/evaluate_val.py") as prov:
        ctx = real_context(cfg)
        prep = prepare(cfg, ctx, labels.train(columns=LABEL_COLS))
        cpath = write_config(cfg, out_dir)
        res = run_look(cfg, prep=prep, load_val=labels.validation, ctx=ctx, out_dir=out_dir, author="scripts/evaluate_val.py")
        prov.output(cpath, res["output"], res[f"{side}_file"])
        prov.note(look_id=cfg["look_id"], config_sha=res["config_sha"], output_sha=res["output_sha"])
    print(json.dumps({**res, "wall_s": round(time.time() - t0, 1)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
