"""Test looks (DESIGN.md §6.3-6.5, §6.8, §10.1-10.3, §11.1, §11.3, §11.4). The only entry point allowed to read
test labels (tripwire ``ENTRY_PERMS``: read val/test/sealed, write sealed).

It reuses the reviewed validation harness ``scripts/evaluate_val.py`` (imported as a module, unchanged): the
pre-flight (committed ``HEAD:FROZEN.tsv`` + clean git), ``check_frozen.py`` with the test-run check, the
score-file loader, ``tap_grid``, the unit-key checker of the scalar guard and the W1 label-free context. The
statistics are the frozen ``nextgem.metrics`` functions; this script only orchestrates them.

Driven by committed specs (no harness change per look):

- **P1 spec** ``results/m6/p1_spec.json`` (template ``results/m6/p1_spec.template.json``): P1's score sources per
  tier (validation for the rehearsal, test for T1) pinned by SHA-256, layout ``gene`` (line-constant) or
  ``line_gene`` (line-specific); the §6.8 binary calls (P57, PGR, FBA<0.01) of the P1 arms (P1 on H1 and on H2);
  failed builds; the attribution contrasts (i)-(v) of §10.1 ("lever not selected" allowed, A12c); optional
  per-task calls and failure-account gene annotations.
- **P2 spec** ``results/m6/p2_spec.json`` (T2): P2's test scores and the B6 refit settings.

Baselines are computed here, not read from files: generic H1/H2 (= B4/B5, G*) from the baseline KO parquets;
B0-B3, B8+/B8-/B8' from ``w1.score_matrices``; B2/B2m/B2c and CE_fit refitted on the fitting tiers (training for
the rehearsal; training ∪ validation for test looks, §4.2/§5.3); B6 refitted on training ∪ validation at T2.

Modes
-----
  --print-config [--look T1|T2|T3 | --rehearsal | --unseal]   config JSON; reads no label
  --print-frozen-candidates [...]                              FROZEN.tsv rows the PI must add; reads no label
  --preflight-only [...]                                       pre-flight + look-state checks; reads no label
  --dry-run [--look T1|T2] [--tier test|validation]            SYNTHETIC validation and test labels (training
                                                               labels are real); no ledger row, no file
  --rehearsal        (§10.2) stage 1 = the dry run on test-tier shapes; stage 2 = the T1 pipeline on validation
                     lines, a VALIDATION LOOK (look id REHEARSAL_<spec_id>, workstream M6, validation ledger)
  --look T1          reveals scalars and CIs (committed under results/m6/test/); seals the rest
  --look T2          after the W7 stop row: C4, C5, P2's scalar grid
  --look T3          the third look (§10.2): only with a committed authorisation; reported as post hoc
  --unseal           after a completed T2: ledger ``unseal`` row, then sealed outputs re-derived, verified
                     against their ledger SHA-256 and published under results/m6/test/unsealed/<look>/

  A look must be named: a bare ``evaluate_test.py`` (no --look / --rehearsal / --unseal) is refused.
  --workers N        processes for the C3 bootstrap and permutation replicates (1-8, default 8); the result is
                     bitwise identical to the serial ``metrics.c3_decision`` (same batches, same spawned seeds).

Order inside every look (as in evaluate_val.py): look-state checks (a completed look id, or a crashed one with another
config, is refused here, before any label read or file write) → pre-flight → provenance → label-free and fitting-tier
work (contexts, score matrices, refits) → config JSON → ledger ``started`` (committed) → the evaluation-tier label
read → revealed outputs (scalar guard) → sealed outputs (``data/locked/sealed/<look>/`` only, written through
``labels._write_protected`` under the tripwire's ``sealed`` scope) → ledger ``output`` rows (SHA-256 only) →
``completed``. stdout carries scalars, file names and SHA-256 values only. BLAS/OpenMP are pinned to one thread
before numpy is imported; the library versions and the thread count are part of the config, and so of the
config_sha that the unseal re-checks against the ledger before it reads any label.

Revealed at T1 (and at the rehearsal's validation stage): ``test_look_T1.tsv``, by section —

- ``confirmatory``: C1 (all genes) and C1b (U minus CE_fit): estimate (Δ̄_lin), ci_lo/ci_hi (lineage-stratified
  bootstrap, B = 10,000), p (sign flip within lineage, 10,000), passed, sesoi, line_mean (Δ̄_line), n_lines,
  n_lineages, status, claim. C3 and C3_expressed_pairs: estimate, ci_lo/ci_hi (line bootstrap of mΔAP_gene(P1) −
  mΔAP_gene(B3)), p (permutation), passed, m_delta_ap_gene_p1, null_threshold (97.5th percentile of 2,000
  within-lineage label permutations), bootstrap_passed, permutation_passed, n_genes_eligible, mean_genes_dropped,
  status, claim (the expressed-pairs row feeds the §6.3 C3 claim condition).
- ``transfer`` (C1): t_mean, t_low, t_high (t-interval, df = K − 1), exact_sign_flip_p (2^K), k_positive,
  k_required (ceil 0.7 K), n_lineages, generalises.
- ``estimates`` (B = 2,000, percentile): C2:B2, C2:B2m, C2:B2c; B8+, B8-, B8' (built on G*); within_parent:H1,
  within_parent:H2; vs_G*:<non-primary arm> — each estimate, ci_lo, ci_hi, line_mean, n_lines, n_lineages;
  headroom_closed (estimate/CI only if reportable; denominator, den_ci_lo, den_ci_hi, reportable);
  tap:<P1 arms and G*> (tap_lin, tap_lin_lo, tap_lin_hi, e_random_lin, lift_lin, n_lines, n_lineages).
- ``binary_6.8``: H1, H2, P1[H1], P1[H2] × {Human1-paper P57 and PGR, Human2-paper P57} × L_ge06 and fixed-cutoff
  FBA<0.01 × {L_prob, L_ge06} × {mcc, f1, precision, recall, specificity, accuracy} × nan_policy {exclude, zero}:
  lineage_equal_mean, line_mean, median, n_lines_used, n_nan, ci_lo, ci_hi; plus paired MCC rows "<arm> - <G*>" and
  "<arm> - <its parent>". A failed build has no model, so all of its calls are non-essential: TP + FP = 0 makes its
  MCC, precision and F1 undefined (NaN). Under nan_policy=exclude such a line is dropped and counted in n_nan; under
  nan_policy=zero it counts as 0.
- ``failed_builds``: per method × lineage (and 'all'): n_failed, n_lines, frac_failed, flag_failing (> 5%). This is
  the only section allowed to name lineages (it is label-free).
- ``counts``: n_lines (evaluation lines), n_lines_no_expression, n_lineages.

``test_look_T1_gates.json``: the gate table (C1 → C1b → C3 → C4 → C5); the C1/C1b/C3 results; the C3 claim
condition; transfer; the headline (the first T1 number, §10.2); and the pre-declared claim rule "adds information
beyond generic + expression" = the paired Δ̄_lin CI lower bound > 0 for all three of B8+, B8−, B8′ (the
all-point-estimates-positive version is reported as descriptive only).

Runtime: one C3 call on the real test shape is ≈ 88 CPU-min serially and T1 runs two (all pairs, expressed pairs);
with --workers 8 they run in parallel replicate batches. A look over 1 CPU-h runs from a pinned worktree (§11.4).
"""

from __future__ import annotations

import os

if __name__ == "__main__":  # one BLAS/OpenMP thread, set before numpy is imported (review F2)
    for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[_v] = "1"

import argparse  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import io  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import platform  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PROJECT = Path(__file__).resolve().parents[1]


def _load_evaluate_val():
    spec = importlib.util.spec_from_file_location("evaluate_val", PROJECT / "scripts" / "evaluate_val.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ev = _load_evaluate_val()
PreflightError = ev.PreflightError
_sha = ev._sha

CODE_FILES = ("scripts/evaluate_test.py", "scripts/evaluate_val.py", "src/nextgem/labels.py", "src/nextgem/w1.py",
              "src/nextgem/metrics.py", "src/nextgem/run.py", "src/nextgem/ledger.py")
B6_CODE = ("scripts/w1_b6.py", "scripts/w1_common.py")
KNOWN_ANSWER_TESTS = tuple(ev.KNOWN_ANSWER_TESTS) + ("tests/test_evaluate_test.py",)
P1_SPEC = "results/m6/p1_spec.json"
P2_SPEC = "results/m6/p2_spec.json"
SPEC_TOP = ("results", "m6")
REVEAL_DIR = "results/m6/test"
REHEARSAL_DIR = "results/m6/rehearsal"
UNSEALED_DIR = "results/m6/test/unsealed"
TEST_LOOKS = ("T1", "T2", "T3")
SEALED_LOOKS = ("T1", "T3")  # looks with a sealed bundle (T2 seals nothing, §10.2)
W7_STOP = "W7_stop"
REHEARSAL_PREFIX = "REHEARSAL_"
REHEARSAL_WS = "M6"
THIRD_LOOK_AUTH = "results/authorisations/test_third_look.md"
SKEW_REL = "results/m6/t1/skew.json"  # §6.5 training skew check of P1 − G*; pins bootstrap_method (review F9)
MAX_WORKERS = 8
C3_WORKERS = 1  # set by main(--workers); the C3 result does not depend on it (bitwise)
CALL = 0.99  # e >= 0.99 <=> ratio_r < 0.01 or a task failure: the §6.1 call (as in W3's per-task table)
# Reviewed refit recipes for a label-fitted P2 (review F10): a W7 idea re-fitted on training ∪ validation with
# nextgem.w7.fit/predict. Only ideas whose features are the W1 B6 set (no extra, dropped or kNN features) qualify.
REFIT_RECIPES = {"w7_idea": ("W7-01", "W7-05")}
REFIT_CODE = ("src/nextgem/w7.py", "scripts/w1_b6.py")  # must be pinned in the P2 spec's provenance
LEVER_NOT_SELECTED = "lever not selected"
ATTRIBUTION = {"i": "P1 vs G*", "ii": "medium alone", "iii": "curation alone", "iv": "constraints alone",
               "v": "curation given constraints"}
CALL_COLS = {"P57": "p57", "PGR": "pgr", "FBA<0.01": "fba"}
FAILURE_GENE_CLASSES = ("paralog_buffered", "medium_sensitive", "gene_less_bypass", "complex_subunit",
                        "housekeeping_adjacent", "biomass_artefact_lethal", "model_contradicted_loss",
                        "expected_fba_failure")
FAILURE_HARNESS_CLASSES = ("common_essential", "cn_amplified", "unexpressed", "missing_from_model")
N_WORST_CELLS = 20  # §10.3
GRID_GENESETS = ("all", "noCE", "noDepMapCE")
GRID_UNIVERSES = ("U", "Ucap")
GRID_STRATA = ("none", "expressed_pairs")
SPEC_ID_RE = re.compile(r"^[A-Za-z0-9_]{1,40}$")
SHA_RE = re.compile(r"^[0-9a-f]{64}$")
P1_KEYS = {"spec_id", "kind", "milestone", "primary", "methods", "arms", "scores", "calls", "failed_builds",
           "attribution", "bootstrap_method", "fitted_on_labels"}
P1_OPTIONAL = {"pipeline", "task_calls", "gene_annotations", "provenance", "notes"}
P2_KEYS = {"spec_id", "kind", "milestone", "primary", "methods", "scores", "failed_builds", "labels", "genesets",
           "universes", "strata", "b6", "bootstrap_method", "fitted_on_labels"}
P2_OPTIONAL = {"pipeline", "provenance", "notes", "refit"}
B6_HGB_KEYS = {"max_iter", "learning_rate", "max_leaf_nodes", "early_stopping"}

# ---- revealed-output guard (T1 / rehearsal / T2 confirmatory part)
REVEALED_SECTIONS = ("confirmatory", "transfer", "estimates", "binary_6.8", "failed_builds", "counts")
STR_COLS = {"section", "name", "label", "geneset", "stratum", "method", "method_a", "method_b", "protocol", "call_rule",
            "gene_set", "metric", "nan_policy", "status", "claim", "passed", "bootstrap_passed", "permutation_passed",
            "generalises", "reportable", "flag_failing", "lineage"}
NUM_COLS = {"estimate", "ci_lo", "ci_hi", "p", "sesoi", "line_mean", "n_lines", "n_lineages", "m_delta_ap_gene_p1",
            "null_threshold", "n_genes_eligible", "mean_genes_dropped", "t_mean", "t_low", "t_high",
            "exact_sign_flip_p", "k_positive", "k_required", "denominator", "den_ci_lo", "den_ci_hi",
            "lineage_equal_mean", "median", "n_lines_used", "n_nan", "n_failed", "frac_failed", "n_lines_no_expression",
            "tap_lin", "tap_lin_lo", "tap_lin_hi", "e_random_lin", "lift_lin"}


# ============================================================================= specs
def _under(p, top=SPEC_TOP) -> bool:
    if not isinstance(p, str) or Path(p).is_absolute() or ".." in Path(p).parts:
        return False
    parts = Path(p).parts
    return len(parts) > len(top) and parts[:len(top)] == top


def _need(ok, msg):
    if not ok:
        raise PreflightError(msg)


def _check_method_names(names, what):
    _need(isinstance(names, list) and len(set(names)) == len(names) and {"H1", "H2"} <= set(names),
          f"{what}: methods must be a list of unique names that includes H1 and H2")
    for m in names:
        _need(isinstance(m, str) and m and m.strip() == m and not ev.ID_RE.match(m) and not m.isdigit()
              and not ev.LABEL_LIKE_RE.search(m) and m != "G*", f"{what}: method name {m!r} is refused")


def _check_file(entry, what, layouts=None):
    """{path, sha256[, layout]} with path below results/m6/ (results/m6/tmp/ for git-ignored, line-keyed files)."""
    keys = {"path", "sha256"} | ({"layout"} if layouts else set())
    _need(isinstance(entry, dict) and set(entry) == keys, f"{what} must be {sorted(keys)}")
    _need(isinstance(entry["sha256"], str) and SHA_RE.match(entry["sha256"]) is not None, f"{what}.sha256 must be hex")
    _need(_under(entry["path"]) and entry["path"].endswith(".parquet"), f"{what}.path must be a .parquet file below results/m6/")
    if layouts:
        _need(entry["layout"] in layouts, f"{what}.layout must be one of {layouts}")


def _check_provenance(prov, what):
    _need(isinstance(prov, dict) and all(
        isinstance(p, str) and (_under(p, ("results",)) or _under(p, ("src",)) or _under(p, ("scripts",)))
        and isinstance(v, str) and SHA_RE.match(v) for p, v in prov.items()),
        f"{what}: provenance must map paths under results/, src/ or scripts/ to SHA-256 pins")


def check_p1_spec(spec: dict) -> None:
    """Static, label-free checks of the P1 spec. Raises PreflightError."""
    w = "P1 spec"
    _need(isinstance(spec, dict) and P1_KEYS <= set(spec) <= P1_KEYS | P1_OPTIONAL,
          f"{w}: keys must be {sorted(P1_KEYS)} (optional {sorted(P1_OPTIONAL)})")
    _need(spec["kind"] == "p1", f"{w}: kind must be 'p1'")
    _need(isinstance(spec["spec_id"], str) and SPEC_ID_RE.match(spec["spec_id"]) is not None,
          f"{w}: spec_id must match {SPEC_ID_RE.pattern}")
    _check_method_names(spec["methods"], w)
    arms = spec["arms"]
    _need(isinstance(arms, dict) and arms and set(arms) <= set(spec["methods"]) - {"H1", "H2"}
          and set(arms.values()) <= {"H1", "H2"} and len(set(arms.values())) == len(arms),
          f"{w}: arms maps each P1 arm (a non-generic method) to its parent H1 or H2, one arm per parent")
    _need(spec["primary"] in arms, f"{w}: primary must be one of the arms")
    for tier in ("validation", "test"):
        _need(isinstance(spec["scores"], dict) and set(spec["scores"]) == {"validation", "test"},
              f"{w}: scores must have exactly the tiers validation and test")
        _check_file(spec["scores"][tier], f"{w}: scores.{tier}", ("gene", "line_gene"))
        _need(isinstance(spec["calls"], dict) and set(spec["calls"]) == {"validation", "test"},
              f"{w}: calls must have exactly the tiers validation and test")
        _check_file(spec["calls"][tier], f"{w}: calls.{tier}", ("gene", "line_gene"))
        _need(isinstance(spec["failed_builds"], dict) and set(spec["failed_builds"]) == {"validation", "test"},
              f"{w}: failed_builds must have exactly the tiers validation and test (null = no failed build)")
        if spec["failed_builds"][tier] is not None:
            _check_file(spec["failed_builds"][tier], f"{w}: failed_builds.{tier}")
    att = spec["attribution"]
    _need(isinstance(att, dict) and set(att) == set(ATTRIBUTION), f"{w}: attribution keys must be {list(ATTRIBUTION)}")
    names = set(spec["methods"]) | {"G*"}
    for k, v in att.items():
        _need(v == LEVER_NOT_SELECTED or (isinstance(v, list) and len(v) == 2 and set(v) <= names and v[0] != v[1]),
              f"{w}: attribution.{k} must be [method_a, method_b] (methods or 'G*') or {LEVER_NOT_SELECTED!r}")
    _need(att["i"] == [spec["primary"], "G*"], f"{w}: attribution.i must be [primary, 'G*'] (§10.1 (i) P1 vs G*)")
    for k in ("task_calls", "gene_annotations"):
        if spec.get(k) is not None:
            _check_file(spec[k], f"{w}: {k}")
    _need(spec["bootstrap_method"] in ("percentile", "bca"), f"{w}: bootstrap_method must be percentile or bca (§6.5)")
    _need(spec["fitted_on_labels"] is False,
          f"{w}: fitted_on_labels must be false: a P1 with label-fitted components needs its refit on training ∪ "
          "validation inside this harness (§4.2), which is not implemented")
    _check_provenance(spec.get("provenance", {}), w)
    _need(SKEW_REL in spec.get("provenance", {}),
          f"{w}: provenance must pin {SKEW_REL} (the §6.5 skew check that decides bootstrap_method)")


def check_p2_spec(spec: dict) -> None:
    w = "P2 spec"
    _need(isinstance(spec, dict) and P2_KEYS <= set(spec) <= P2_KEYS | P2_OPTIONAL,
          f"{w}: keys must be {sorted(P2_KEYS)} (optional {sorted(P2_OPTIONAL)})")
    _need(spec["kind"] == "p2", f"{w}: kind must be 'p2'")
    _need(isinstance(spec["spec_id"], str) and SPEC_ID_RE.match(spec["spec_id"]) is not None, f"{w}: bad spec_id")
    _check_method_names(spec["methods"], w)
    _need(spec["primary"] in set(spec["methods"]) - {"H1", "H2"}, f"{w}: primary must be a non-generic method")
    _need(isinstance(spec["scores"], dict) and set(spec["scores"]) == {"test"}, f"{w}: scores must have the tier test only")
    _check_file(spec["scores"]["test"], f"{w}: scores.test", ("gene", "line_gene"))
    _need(isinstance(spec["failed_builds"], dict) and set(spec["failed_builds"]) == {"test"}, f"{w}: failed_builds.test")
    if spec["failed_builds"]["test"] is not None:
        _check_file(spec["failed_builds"]["test"], f"{w}: failed_builds.test")
    from nextgem import w1

    _need(isinstance(spec["labels"], list) and spec["labels"] and set(spec["labels"]) <= set(w1.GRID_LABELS),
          f"{w}: labels must be a subset of {w1.GRID_LABELS}")
    for k, allowed, prim in (("genesets", GRID_GENESETS, "all"), ("universes", GRID_UNIVERSES, "U"),
                             ("strata", GRID_STRATA, "none")):
        _need(isinstance(spec[k], list) and prim in spec[k] and set(spec[k]) <= set(allowed),
              f"{w}: {k} must include {prim!r} and be a subset of {allowed}")
    b6 = spec["b6"]
    _need(isinstance(b6, dict) and set(b6) == {"features", "hgb"} and b6["features"] == "w1"
          and isinstance(b6["hgb"], dict) and set(b6["hgb"]) == B6_HGB_KEYS,
          f"{w}: b6 must be {{features: 'w1', hgb: {sorted(B6_HGB_KEYS)}}} (the W1 feature set refitted on training ∪ validation)")
    _need(spec["bootstrap_method"] in ("percentile", "bca"), f"{w}: bootstrap_method must be percentile or bca")
    _check_provenance(spec.get("provenance", {}), w)
    refit = spec.get("refit")
    if spec["fitted_on_labels"] is False:
        _need(refit is None, f"{w}: refit is only allowed with fitted_on_labels: true")
        return
    _need(spec["fitted_on_labels"] is True, f"{w}: fitted_on_labels must be a boolean")
    _need(isinstance(refit, dict) and set(refit) == {"recipe", "idea", "method"}
          and refit["recipe"] in REFIT_RECIPES and refit["idea"] in REFIT_RECIPES[refit["recipe"]],
          f"{w}: a label-fitted P2 needs refit = {{recipe, idea, method}} with a reviewed recipe "
          f"({ {k: list(v) for k, v in REFIT_RECIPES.items()} })")
    _need(refit["method"] == spec["primary"], f"{w}: refit.method must be the primary P2 method")
    _need(set(REFIT_CODE) <= set(spec.get("provenance", {})),
          f"{w}: the refit recipe's code {list(REFIT_CODE)} must be pinned in provenance")


def load_spec(rel: str, kind: str) -> tuple[dict, Path]:
    from nextgem import paths

    p = paths.state_root() / rel
    if not p.is_file():
        raise PreflightError(f"{kind.upper()} spec {rel} does not exist")
    try:
        spec = json.loads(p.read_text())
    except ValueError:
        raise PreflightError(f"{rel} is not valid JSON") from None
    (check_p1_spec if kind == "p1" else check_p2_spec)(spec)
    return spec, p


# ============================================================================= label-free file loaders
def _read(rel: str) -> pd.DataFrame:
    from nextgem import paths

    return pd.read_parquet(paths.state_root() / rel)


def load_calls(rel: str, layout: str, arms) -> pd.DataFrame:
    """§6.8 binary calls of the P1 arms (label-free): keys + ``method`` + bool ``p57``, ``pgr``, ``fba``."""
    df = _read(rel)
    keys = ["entrez"] if layout == "gene" else ["ModelID", "entrez"]
    cols = keys + ["method", *CALL_COLS.values()]
    if sorted(df.columns) != sorted(cols):
        raise PreflightError(f"calls {Path(rel).name}: columns must be exactly {cols}")
    if not pd.api.types.is_integer_dtype(df["entrez"]):
        raise PreflightError(f"calls {Path(rel).name}: entrez must be integer")
    if set(df["method"]) != set(arms):
        raise PreflightError(f"calls {Path(rel).name}: method values must be exactly the arms {sorted(arms)}")
    for c in CALL_COLS.values():
        if not pd.api.types.is_bool_dtype(df[c]):
            raise PreflightError(f"calls {Path(rel).name}: column {c!r} must be boolean")
    if df.duplicated(keys + ["method"]).any():
        raise PreflightError(f"calls {Path(rel).name}: duplicate rows")
    if layout == "line_gene" and not df.ModelID.map(lambda x: isinstance(x, str) and ev.ID_RE.match(x) is not None).all():
        raise PreflightError(f"calls {Path(rel).name}: ModelID values must be line IDs")
    return df


def load_failed(rel: str | None, methods) -> pd.DataFrame:
    if rel is None:
        return pd.DataFrame({"ModelID": pd.Series(dtype=str), "method": pd.Series(dtype=str)})
    df = _read(rel)
    if sorted(df.columns) != ["ModelID", "method"]:
        raise PreflightError(f"failed builds {Path(rel).name}: columns must be exactly ModelID, method")
    bad = set(df["method"]) - (set(methods) - {"H1", "H2"})
    if bad:
        raise PreflightError(f"failed builds {Path(rel).name}: unknown or generic method(s) {sorted(bad)}")
    if not df.ModelID.map(lambda x: isinstance(x, str) and ev.ID_RE.match(x) is not None).all() or df.duplicated().any():
        raise PreflightError(f"failed builds {Path(rel).name}: ModelID must be line IDs, no duplicate rows")
    return df


def load_task_calls(rel: str | None) -> pd.DataFrame | None:
    if rel is None:
        return None
    df = _read(rel)
    if sorted(df.columns) != ["entrez", "fails", "task"] or not pd.api.types.is_bool_dtype(df["fails"]) \
            or not pd.api.types.is_integer_dtype(df["entrez"]) or df.duplicated(["entrez", "task"]).any():
        raise PreflightError(f"task calls {Path(rel).name}: columns entrez (int), task (str), fails (bool); no duplicates")
    return df


def load_annotations(rel: str | None) -> pd.DataFrame | None:
    if rel is None:
        return None
    df = _read(rel)
    extra = set(df.columns) - {"entrez", *FAILURE_GENE_CLASSES}
    if "entrez" not in df.columns or extra or not pd.api.types.is_integer_dtype(df["entrez"]) or df.entrez.duplicated().any():
        raise PreflightError(f"gene annotations {Path(rel).name}: entrez (int, unique) + boolean columns among "
                             f"{FAILURE_GENE_CLASSES}; unknown {sorted(extra)}")
    for c in df.columns.drop("entrez"):
        if not pd.api.types.is_bool_dtype(df[c]):
            raise PreflightError(f"gene annotations {Path(rel).name}: {c!r} must be boolean")
    return df.set_index("entrez")


# ============================================================================= configuration
def _frozen_value(item: str):
    from nextgem import frozen

    for r in frozen.read_rows():
        if r["item"] == item and r["kind"] == "value":
            return r["sha256_or_value"]
    raise PreflightError(f"FROZEN.tsv has no value row {item!r}")


def _gstar() -> str:
    from nextgem import paths

    try:
        g = json.loads((paths.state_root() / ev.W1_GSTAR).read_text())["G_star"]
    except (OSError, ValueError, KeyError, TypeError):
        raise PreflightError(f"G-star file {ev.W1_GSTAR} missing or unreadable") from None
    if g not in ev.GSTAR_SHORT:
        raise PreflightError(f"G-star file {ev.W1_GSTAR}: G_star must be one of {sorted(ev.GSTAR_SHORT)}")
    return g


def runtime_info() -> dict:
    """Library versions and the BLAS/OpenMP thread count: part of the config (and so of the config_sha), because
    the unseal re-derives the sealed bytes and must run under the same numerics (review F2)."""
    import scipy
    import sklearn

    try:
        from threadpoolctl import threadpool_info

        threads = max([int(p.get("num_threads", 1)) for p in threadpool_info()] or [1])
    except Exception:  # noqa: BLE001 - threadpoolctl missing: fall back to the environment
        threads = int(os.environ.get("OMP_NUM_THREADS", "0") or 0)
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "pandas": pd.__version__, "sklearn": sklearn.__version__, "blas_threads": threads}


def _common_config(spec: dict, spec_rel: str, look: str, tier: str) -> dict:
    from nextgem import metrics, paths, w1
    from nextgem.config import SEEDS

    root = paths.state_root()
    gstar = _gstar()
    primary_label = _frozen_value("label.primary")
    if primary_label not in w1.GRID_LABELS:
        raise PreflightError(f"label.primary {primary_label!r} is not a grid label")
    rel = dict(ev.input_files(False))
    rel.update(w1_gstar=ev.W1_GSTAR)
    return {
        "look": look, "tier": tier, "fit_tiers": ["train"] if tier == "validation" else ["train", "validation"],
        "spec_path": spec_rel, "spec_id": spec["spec_id"], "kind": spec["kind"], "milestone": spec["milestone"],
        "primary": spec["primary"], "methods": list(spec["methods"]), "score_methods": list(spec["methods"]),
        "bootstrap_method": spec["bootstrap_method"], "runtime": runtime_info(),
        "G_star": gstar, "G_star_short": ev.GSTAR_SHORT[gstar], "primary_label": primary_label,
        "sesoi": metrics.SESOI, "sesoi_c1b": float(_frozen_value("sesoi.C1b")),
        "B_conf": metrics.B_CONFIRMATORY, "n_flip": metrics.N_SIGN_FLIPS, "n_perm": metrics.N_C3_PERMUTATIONS,
        "B": metrics.B_OTHER, "seed": SEEDS[0], "random_seeds": SEEDS[0:5], "alpha": metrics.ALPHA,
        "labels": list(w1.GRID_LABELS), "binary_label": "L_ge06", "genesets": list(GRID_GENESETS),
        "universes": list(GRID_UNIVERSES), "strata": list(GRID_STRATA),
        "provenance": spec.get("provenance", {}), "notes": spec.get("notes", ""),
        "_inputs_rel": rel,
    }


def _finish_config(cfg: dict, code_files) -> dict:
    from nextgem import paths

    root = paths.state_root()
    missing = [v for v in cfg["_inputs_rel"].values() if not (root / v).is_file()]
    if missing:
        raise PreflightError(f"input file(s) missing: {missing}")
    cfg["inputs"] = {k: _sha(root / v) for k, v in cfg["_inputs_rel"].items()}
    cfg["code"] = {c: _sha(paths.PROJECT_ROOT / c) for c in code_files if (paths.PROJECT_ROOT / c).is_file()}
    cfg["_code_files"] = list(code_files)
    return cfg


def p1_config(look: str, spec_rel: str = P1_SPEC) -> dict:
    """Config of the rehearsal (look = REHEARSAL_<spec_id>, tier validation) or of T1/T3 (tier test).
    Label-free. The score/call files must carry their pinned SHA-256 and pass the loaders."""
    from nextgem import paths

    spec, spath = load_spec(spec_rel, "p1")
    tier = "validation" if look == "REHEARSAL" else "test"
    look_id = REHEARSAL_PREFIX + spec["spec_id"] if look == "REHEARSAL" else look
    cfg = _common_config(spec, spec_rel, look_id, tier)
    root = paths.state_root()
    rel = cfg["_inputs_rel"]
    rel["spec"] = spec_rel
    src = {}
    for t in ("validation", "test"):
        sc, ca, fb = spec["scores"][t], spec["calls"][t], spec["failed_builds"][t]
        for name, e in ((f"scores_{t}", sc), (f"calls_{t}", ca), (f"failed_{t}", fb)):
            if e is None:
                continue
            p = root / e["path"]
            if not p.is_file() or _sha(p) != e["sha256"]:
                raise PreflightError(f"{e['path']}: missing or SHA-256 differs from the pin in {spath.name}")
            rel[name] = e["path"]
        ev.load_scores(root / sc["path"], sc["layout"], spec["methods"])
        load_calls(ca["path"], ca["layout"], spec["arms"])
        load_failed(fb["path"] if fb else None, spec["methods"])
        src[t] = {"scores": sc["path"], "scores_layout": sc["layout"], "calls": ca["path"], "calls_layout": ca["layout"],
                  "failed_builds": fb["path"] if fb else None}
    for k in ("task_calls", "gene_annotations"):
        e = spec.get(k)
        if e is not None:
            p = root / e["path"]
            if not p.is_file() or _sha(p) != e["sha256"]:
                raise PreflightError(f"{e['path']}: missing or SHA-256 differs from the pin in {spath.name}")
            rel[k] = e["path"]
    load_task_calls(spec["task_calls"]["path"] if spec.get("task_calls") else None)
    load_annotations(spec["gene_annotations"]["path"] if spec.get("gene_annotations") else None)
    # §6.5: bootstrap_method must be the one decided by the pinned training skew check (review F9)
    sk = root / SKEW_REL
    if not sk.is_file() or _sha(sk) != spec["provenance"][SKEW_REL]:
        raise PreflightError(f"{SKEW_REL}: missing or SHA-256 differs from the provenance pin in {spath.name}")
    try:
        skew_method = json.loads(sk.read_text())["bootstrap_method"]
    except (ValueError, KeyError, TypeError):
        raise PreflightError(f"{SKEW_REL} has no bootstrap_method") from None
    if skew_method != spec["bootstrap_method"]:
        raise PreflightError(f"bootstrap_method {spec['bootstrap_method']!r} differs from {skew_method!r} in {SKEW_REL}")
    cfg.update(arms=dict(spec["arms"]), attribution=dict(spec["attribution"]), sources=src,
               task_calls=spec["task_calls"]["path"] if spec.get("task_calls") else None,
               gene_annotations=spec["gene_annotations"]["path"] if spec.get("gene_annotations") else None,
               post_hoc=look == "T3")
    if look == "T3" and (root / THIRD_LOOK_AUTH).is_file():  # its absence is reported by state_problems
        rel["third_look_authorisation"] = THIRD_LOOK_AUTH
    return _finish_config(cfg, CODE_FILES)


def p2_config(spec_rel: str = P2_SPEC, p1_rel: str = P1_SPEC) -> dict:
    from nextgem import paths

    from nextgem import w7

    spec, spath = load_spec(spec_rel, "p2")
    cfg = _common_config(spec, spec_rel, "T2", "test")
    root = paths.state_root()
    if cfg["primary_label"] != "L_prob":
        raise PreflightError("T2 needs the primary label L_prob (B6 and the W7 recipes are fitted on L_prob only)")
    if spec["b6"]["hgb"] != dict(w7.HGB_M2):
        raise PreflightError(f"b6.hgb must equal W7's B6 settings {dict(w7.HGB_M2)} (review F10)")
    refit = spec.get("refit")
    if refit is not None:
        cfg["refit"] = dict(refit)
        cfg["score_methods"] = [m for m in spec["methods"] if m != refit["method"]]
    rel = cfg["_inputs_rel"]
    rel["spec"] = spec_rel
    if (root / p1_rel).is_file():  # always present for a real T2 (T1 completed); absent only in a dry run
        rel["p1_spec"] = p1_rel
    sc, fb = spec["scores"]["test"], spec["failed_builds"]["test"]
    for name, e in (("scores_test", sc), ("failed_test", fb)):
        if e is None:
            continue
        p = root / e["path"]
        if not p.is_file() or _sha(p) != e["sha256"]:
            raise PreflightError(f"{e['path']}: missing or SHA-256 differs from the pin in {spath.name}")
        rel[name] = e["path"]
    ev.load_scores(root / sc["path"], sc["layout"], cfg["score_methods"])
    load_failed(fb["path"] if fb else None, spec["methods"])
    rel.update(mut_dam="results/w1/tmp/mut_dam.parquet", mut_hot="results/w1/tmp/mut_hot.parquet",
               gene_neighbours="results/w1/tmp/gene_neighbours.json")
    cfg.update(labels=list(spec["labels"]), genesets=list(spec["genesets"]), universes=list(spec["universes"]),
               strata=list(spec["strata"]), b6=dict(spec["b6"]), arms={}, post_hoc=False,
               sources={"test": {"scores": sc["path"], "scores_layout": sc["layout"], "calls": None, "calls_layout": None,
                                 "failed_builds": fb["path"] if fb else None}})
    code = CODE_FILES + B6_CODE + (tuple(c for c in REFIT_CODE if c not in B6_CODE) if refit else ())
    return _finish_config(cfg, code)


def look_config(mode: str, spec_rel: str | None = None) -> dict:
    """mode: REHEARSAL, T1, T2, T3 or UNSEAL (the T1 config: the unseal re-derives T1's sealed bundle)."""
    if mode == "T2":
        return p2_config(spec_rel or P2_SPEC)
    return p1_config("T1" if mode == "UNSEAL" else mode, spec_rel or P1_SPEC)


def config_sha(cfg: dict) -> str:
    return ev.config_sha(cfg)


def required_files(cfg: dict) -> list[str]:
    return list(cfg["_code_files"]) + list(cfg["_inputs_rel"].values()) + list(cfg.get("provenance", {}))


def pinned_problems(cfg: dict) -> list[str]:
    return ev.pinned_problems(cfg)


def frozen_candidates(cfg: dict) -> list[str]:
    from nextgem import frozen, paths

    head = subprocess.run(["git", "-C", str(paths.state_root()), "show", "HEAD:FROZEN.tsv"], capture_output=True, text=True)
    have = {r["path"]: r["sha256_or_value"] for r in frozen.parse(head.stdout) if r["kind"] == "file"} if head.returncode == 0 else {}
    day = time.strftime("%Y-%m-%d")
    rows = []
    for rel in dict.fromkeys(required_files(cfg) + list(KNOWN_ANSWER_TESTS)):
        p = paths.state_root() / rel
        if not p.is_file():
            rows.append(f"# MISSING\t{rel}\t(no row possible until the file exists)")
            continue
        if have.get(rel) == _sha(p):
            continue
        note = "test looks: code" if rel in cfg["_code_files"] else (
            "known-answer test (DESIGN 11.3)" if rel in KNOWN_ANSWER_TESTS else f"test looks ({cfg['look']}): input")
        rows.append("\t".join([rel, "file", rel, _sha(p), day, "", note]))
    return rows


def check_constants(cfg: dict) -> list[str]:
    """Real looks run with the frozen §6.3/§6.5 constants only (tests may shrink them on toy data)."""
    from nextgem import metrics
    from nextgem.config import SEEDS

    want = {"B_conf": metrics.B_CONFIRMATORY, "n_flip": metrics.N_SIGN_FLIPS, "n_perm": metrics.N_C3_PERMUTATIONS,
            "B": metrics.B_OTHER, "seed": SEEDS[0], "sesoi": metrics.SESOI, "alpha": metrics.ALPHA}
    return [f"{k} = {cfg.get(k)!r}, frozen value {v!r}" for k, v in want.items() if cfg.get(k) != v or type(cfg.get(k)) is not type(v)]


# ============================================================================= look-state checks (ledgers only)
def _looks(name: str) -> dict:
    from nextgem import ledger

    return ledger.looks(ledger.read(name))


def w7_stop_row(p2_sha: str | None) -> tuple[dict | None, str]:
    """The W7 stop: a ``non_look`` row with look_id ``W7_stop`` in the validation or test ledger whose detail
    JSON names the frozen P2 spec (``p2_spec_sha256``)."""
    from nextgem import ledger

    rows = [r for n in ("validation", "test") for r in ledger.read(n) if r["event"] == "non_look" and r["look_id"] == W7_STOP]
    if not rows:
        return None, "no W7 stop row (non_look, look_id W7_stop) in the validation or test ledger"
    r = rows[-1]
    if p2_sha is None or ledger.detail_json(r).get("p2_spec_sha256") != p2_sha:
        return None, "the W7 stop row does not name the current P2 spec (detail p2_spec_sha256)"
    return r, "ok"


def code_sha(cfg: dict) -> str:
    """SHA-256 of a config's harness-code map; pinned in every look's ``started`` row (review F5)."""
    return hashlib.sha256(json.dumps(cfg["code"], sort_keys=True).encode()).hexdigest()


def _started_detail(lk: dict | None) -> dict:
    try:
        return json.loads(lk["started"]["detail"] or "{}")
    except (TypeError, KeyError, ValueError):
        return {}


def _spec_at(commit: str, rel: str) -> dict | None:
    from nextgem import paths

    r = subprocess.run(["git", "-C", str(paths.state_root()), "show", f"{commit}:{rel}"], capture_output=True, text=True)
    try:
        return json.loads(r.stdout) if r.returncode == 0 else None
    except ValueError:
        return None


T3_AUTH_KEYS = ("regression_test", "fix_commit", "t1_spec_sha256", "t3_spec_sha256")
P1_IDENTITY = ("spec_id", "primary", "methods", "arms", "attribution")


def third_look_authorised(cfg: dict | None = None) -> list[str]:
    """§10.2 third look (review F4). The committed authorisation names: the label-free regression test (its file
    frozen with its current SHA-256 in FROZEN.tsv), the fix commit (an ancestor of HEAD), T1's P1 spec SHA-256
    (``t1_spec_sha256`` = the spec_sha of T1's started row) and this look's (``t3_spec_sha256`` = the current spec).
    If the spec changed since T1, its identity (spec_id, primary, methods, arms, attribution) must equal the spec
    committed at T1: the third look re-evaluates T1's P1 after a bug fix, never another pipeline."""
    from nextgem import frozen, ledger, paths

    root = paths.state_root()
    p = root / THIRD_LOOK_AUTH
    if not p.is_file() or not ledger.is_committed_clean(root, THIRD_LOOK_AUTH):
        return [f"the third look needs a committed {THIRD_LOOK_AUTH} (DESIGN §10.2)"]
    kv = dict(re.findall(rf"^({'|'.join(T3_AUTH_KEYS)}):\s*(\S+)\s*$", p.read_text(), re.M))
    if set(kv) != set(T3_AUTH_KEYS):
        return [f"{THIRD_LOOK_AUTH} must name {', '.join(T3_AUTH_KEYS)} (one 'key: value' line each)"]
    probs = []
    r = subprocess.run(["git", "-C", str(root), "merge-base", "--is-ancestor", kv["fix_commit"], "HEAD"], capture_output=True)
    if r.returncode != 0:
        probs.append("fix_commit is not an ancestor of HEAD")
    f, _, name = kv["regression_test"].partition("::")
    tp = paths.PROJECT_ROOT / f
    if not (f.startswith("tests/") and tp.is_file() and re.search(rf"^def {re.escape(name)}\(", tp.read_text(), re.M)):
        probs.append(f"regression test {kv['regression_test']} not found")
    elif {r["path"]: r["sha256_or_value"] for r in frozen.read_rows() if r["kind"] == "file"}.get(f) != _sha(tp):
        probs.append(f"the regression test file {f} is not frozen with its current SHA-256 in FROZEN.tsv")
    t1 = _looks("test").get("T1")
    if t1 is None or t1["started"] is None:
        return probs + ["T1 has no started row"]
    t1_sha = _started_detail(t1).get("spec_sha")
    spec_rel = cfg["spec_path"] if cfg else P1_SPEC
    cur = _sha(root / spec_rel) if (root / spec_rel).is_file() else None
    if kv["t1_spec_sha256"] != t1_sha:
        probs.append("t1_spec_sha256 is not the P1 spec SHA-256 recorded by the T1 look")
    if kv["t3_spec_sha256"] != cur:
        probs.append("t3_spec_sha256 is not the SHA-256 of the current P1 spec")
    if cur != t1_sha:
        old = _spec_at(t1["started"]["head_commit"], spec_rel)
        new = json.loads((root / spec_rel).read_text()) if cur else None
        if old is None or new is None or any(old.get(k) != new.get(k) for k in P1_IDENTITY):
            probs.append(f"the P1 spec differs from T1's in {'/'.join(P1_IDENTITY)}: the third look re-evaluates T1's P1 only")
    return probs


def state_problems(mode: str, cfg: dict | None = None) -> list[str]:
    """Look-order rules (§10.2, §11.3), from the ledgers and the config alone; no label is read and no file is
    written. ``cfg`` defaults to the mode's config (built label-free)."""
    test, val = _looks("test"), _looks("validation")
    done = lambda lk, k: k in lk and lk[k]["completed"]  # noqa: E731
    probs = []
    if cfg is None and mode != "UNSEAL":
        cfg = look_config(mode)
    if mode in ("REHEARSAL", "T1", "T2", "T3"):  # review F3: refuse before any label read or config write
        name, lks = ("validation", val) if mode == "REHEARSAL" else ("test", test)
        lk = lks.get(cfg["look"])
        if lk is not None and lk["completed"]:
            probs.append(f"look {cfg['look']} is already completed in the {name} ledger; looks are never repeated")
        elif lk is not None and lk["started"] is not None and lk["config_sha"] != config_sha(cfg):
            probs.append(f"look {cfg['look']} crashed earlier with another config; it resumes only with the same config_sha")
    if mode == "T1":
        rid = REHEARSAL_PREFIX + cfg["spec_id"]
        if not done(val, rid):
            probs.append(f"the rehearsal ({rid}) is not completed in the validation ledger (DESIGN §10.2: rehearsal "
                         "before T1)")
        else:
            d = _started_detail(val[rid])
            if d.get("spec_sha") != cfg["inputs"]["spec"]:
                probs.append("the completed rehearsal used a different P1 spec (spec_sha)")
            if d.get("code_sha") != code_sha(cfg):
                probs.append("the completed rehearsal ran another harness version (code_sha; review F5)")
    if mode in ("T2", "T3") and not done(test, "T1"):
        probs.append("T1 is not completed")
    if mode == "T2":
        _, why = w7_stop_row(cfg["inputs"]["spec"])
        if why != "ok":
            probs.append(f"T2 refused: {why} (DESIGN §10.2: T2 runs after the W7 stop)")
    if mode == "T3":
        if done(test, "T3"):
            probs.append("the third look was already used")
        probs += third_look_authorised(cfg)
    if mode == "UNSEAL":
        if not done(test, "T2"):
            probs.append("unseal refused: T2 is not completed (DESIGN §10.2: unseal after T2)")
    return probs


def unseal_problems(cfg_for) -> list[str]:
    """Before the unseal reads any label (review F2): for every completed sealed look, the current config must
    reproduce the config_sha of its ``started`` row (same code, inputs, constants, library versions, threads)."""
    lks = _looks("test")
    probs = []
    for look in SEALED_LOOKS:
        if look in lks and lks[look]["completed"]:
            try:
                c = cfg_for(look)
            except PreflightError as e:
                probs.append(f"{look}: config cannot be rebuilt ({e})")
                continue
            if config_sha(c) != lks[look]["config_sha"]:
                probs.append(f"{look}: the current config_sha differs from the look's started row (code, inputs, "
                             "constants, library versions or threads changed); the sealed bytes would not reproduce")
    return probs


# ============================================================================= contexts and preparation
def real_context(cfg: dict, tier: str) -> dict:
    """Label-free context on the real data for the evaluation tier (runs before ``started``)."""
    from nextgem import labels, paths, w1

    ctx = ev.real_context({})  # the reviewed W1 context (genes, U∩, generic table, omics, strata, DepMap CE, ...)
    split = labels.load_split().set_index("ModelID")
    lf = ctx["lf"]
    tl = sorted(split.index[split.tier == tier])
    ctx.update(tier=tier, tier_lines=tl, eval_lines=[l for l in tl if l in lf.index and bool(lf.loc[l, "has_expression_matrix"])],
               gene_features=w1.gene_features())
    ctx["tier_of"] = split.tier
    if cfg.get("kind") == "p2":
        nb = json.loads((paths.state_root() / "results/w1/tmp/gene_neighbours.json").read_text())
        ctx["b6_fit"] = lambda fit_long, lines, genes, b6, seed: b6_fit_predict(fit_long, lines, genes, ctx, nb, b6, seed)
        if cfg.get("refit"):
            ctx["refit_inputs"] = lambda: w7_inputs(ctx, nb)
    return ctx


def w7_inputs(ctx: dict, nb: dict):
    """Label-free ``nextgem.w7.Inputs`` for the reviewed refit recipes: the W1 B6 feature inputs only (as
    scripts/w7_common.build_inputs, without GEM, kNN or deviation features). ``tier_of`` is None: the train-only
    assertion of w7.fit is replaced by :func:`refit_fit_predict`'s training ∪ validation guard."""
    from nextgem import w1, w7

    genes = list(ctx["genes"])
    omics = ctx["omics"]
    aggs = w1.expression_aggregates(omics["expr"], genes, nb["neighbours_h2"], nb["paralogs"])
    lf = ctx["lf"]
    class_codes = {c: i for i, c in enumerate(sorted(lf.base_medium_class.dropna().unique()))}
    return w7.Inputs(genes=genes, omics=omics, aggs=aggs, gfeat=ctx["gene_features"], lf=lf, class_codes=class_codes,
                     tier_of=None, lineage_of=ctx["lineage"])


class RefitLeakageError(AssertionError):
    pass


def refit_fit_predict(refit: dict, inp, fit_long: pd.DataFrame, lines, tier_of: pd.Series, seed) -> pd.DataFrame:
    """The reviewed refit plug-in for a label-fitted P2 (review F10; §4.2 refit once on training ∪ validation).

    ``refit`` = {recipe: 'w7_idea', idea: one of REFIT_RECIPES['w7_idea'], method}. The idea is looked up in
    ``nextgem.w7.IDEAS`` (code pinned by SHA in the P2 spec) and must use the W1 B6 features only; it is fitted with
    ``w7.fit`` on ``fit_long`` and scored with ``w7.predict`` — the exact W7 code path. Guards: every fitting row is a
    training or validation line of the split (a test line raises RefitLeakageError), and no predicted line is a
    fitting line. There is no other code path."""
    from nextgem import w7

    if refit.get("recipe") not in REFIT_RECIPES or refit.get("idea") not in REFIT_RECIPES[refit["recipe"]]:
        raise PreflightError(f"refit recipe {refit!r} is not a reviewed recipe {REFIT_RECIPES}")
    idea = w7.IDEAS[refit["idea"]]
    if idea.kind != "clf" or idea.extra or idea.drop or idea.mono_extra:
        raise PreflightError(f"{idea.id} is not a B6-feature classifier; it needs its own reviewed recipe")
    fit_ids = pd.Index(pd.unique(fit_long.ModelID))
    t = tier_of.reindex(fit_ids)
    bad = fit_ids[~t.isin(["train", "validation"]).to_numpy()]
    if len(bad):
        raise RefitLeakageError(f"{len(bad)} fitting line(s) outside training ∪ validation (e.g. {list(bad[:2])})")
    overlap = set(lines) & set(fit_ids)
    if overlap:
        raise RefitLeakageError(f"{len(overlap)} predicted line(s) are fitting lines")
    f = w7.fit(inp, fit_long, idea, seed)
    return w7.predict(inp, f, list(lines))


def b6_fit_predict(fit_long, lines, genes, ctx, nb, b6: dict, seed) -> pd.DataFrame:
    """B6 (§8) refitted once on the fitting tiers (training ∪ validation at T2, §4.2) with the frozen
    hyperparameters of the P2 spec: the W1 feature set and stacking of scripts/w1_b6.py, L_prob only."""
    from sklearn.ensemble import HistGradientBoostingClassifier

    from nextgem import w1

    sys.path.insert(0, str(PROJECT / "scripts"))
    import w1_b6  # noqa: E402

    line_class, omics, gf, lf = ctx["line_class"], ctx["omics"], ctx["gene_features"], ctx["lf"]
    class_codes = {c: i for i, c in enumerate(sorted(line_class.dropna().unique()))}
    med = w1.B6_FEATURES.index("medium_class")
    aggs = w1.expression_aggregates(omics["expr"], genes, nb["neighbours_h2"], nb["paralogs"])
    fitted = w1.fit_components(fit_long, line_class, labels_=("L_prob",))
    b2s, b2ms = w1_b6.stacked_priors(fit_long, line_class, genes)
    fl = list(b2s.index)
    X = w1.b6_matrix(fl, genes, b2s, b2ms, omics, aggs, gf, lf, class_codes)
    Yf = w1.wide(fit_long, "L_prob", fl, genes).to_numpy().reshape(-1)
    ok = ~np.isnan(Yf)
    Xo = X[ok]
    Xo[:, med] = np.where(np.isnan(Xo[:, med]) | (Xo[:, med] < 0), np.nan, Xo[:, med])
    clf = HistGradientBoostingClassifier(**b6["hgb"], random_state=seed, categorical_features=[med])
    clf.fit(Xo, Yf[ok].astype(int))
    b2, b2m = w1.heldout_priors(fitted, list(lines), line_class, genes)
    Xp = w1.b6_matrix(list(lines), genes, b2, b2m, omics, aggs, gf, lf, class_codes)
    Xp[:, med] = np.where(Xp[:, med] < 0, np.nan, Xp[:, med])
    p = clf.predict_proba(Xp)[:, 1].reshape(len(lines), len(genes))
    return pd.DataFrame(np.round(p, 8), index=list(lines), columns=list(genes))


def _tile(v: pd.Series, lines, genes) -> pd.DataFrame:
    return pd.DataFrame(np.tile(v.reindex(list(genes)).to_numpy(), (len(lines), 1)), index=list(lines), columns=list(genes))


def _spec_matrices(cfg, ctx, tier, S: dict, gen_S: dict) -> None:
    """Spec methods from the pinned scores file into S (generic H1/H2 must equal the W1 look's)."""
    from nextgem import paths

    src = cfg["sources"][tier]
    lines, genes = list(ctx["eval_lines"]), list(ctx["genes"])
    sp = paths.state_root() / src["scores"]
    if _sha(sp) != cfg["inputs"][f"scores_{tier}"]:
        raise RuntimeError("scores file changed since the config was built")
    df = ev.load_scores(sp, src["scores_layout"], cfg["score_methods"])
    if set(df.entrez) != set(genes):
        raise RuntimeError(f"scores cover {df.entrez.nunique()} genes; the look needs exactly the {len(genes)} genes of U")
    if src["scores_layout"] == "line_gene" and set(df.ModelID) != set(lines):
        raise RuntimeError(f"line_gene scores must cover exactly the {tier} expression lines")
    for m in cfg["score_methods"]:
        if src["scores_layout"] == "gene":
            M = _tile(df.set_index("entrez")[m].astype(float), lines, genes)
        else:
            M = df.pivot(index="ModelID", columns="entrez", values=m).reindex(index=lines, columns=genes).astype(float)
        if m in ("H1", "H2"):
            if not np.allclose(M.to_numpy(), gen_S[m].to_numpy(), rtol=0, atol=1e-9):
                raise RuntimeError(f"scores column {m} differs from the generic {m} of the W1 look")
            continue
        S[m] = M


def _call_matrices(cfg, ctx, tier) -> dict:
    """(method, rule) -> lines x genes bool, for generic H1/H2 (line-constant) and the P1 arms (spec file)."""
    from nextgem import paths

    lines, genes, gen = list(ctx["eval_lines"]), list(ctx["genes"]), ctx["gen"]
    out = {}
    for s in ("H1", "H2"):
        for rule, col in (("P57", f"p57_{s}"), ("PGR", f"pgr_{s}"), ("FBA<0.01", f"fba_call_{s}")):
            out[(s, rule)] = _tile(gen[col].astype(bool), lines, genes).astype(bool)
    src = cfg["sources"][tier]
    if not cfg.get("arms"):
        return out
    if _sha(paths.state_root() / src["calls"]) != cfg["inputs"][f"calls_{tier}"]:
        raise RuntimeError("calls file changed since the config was built")
    df = load_calls(src["calls"], src["calls_layout"], cfg["arms"])
    for arm in cfg["arms"]:
        sub = df[df.method == arm]
        if set(sub.entrez) != set(genes):
            raise RuntimeError(f"calls of {arm} must cover exactly the genes of U (absent genes: False)")
        for rule, col in CALL_COLS.items():
            if src["calls_layout"] == "gene":
                M = _tile(sub.set_index("entrez")[col].astype(bool), lines, genes)
            else:
                if set(sub.ModelID) != set(lines):
                    raise RuntimeError(f"line_gene calls of {arm} must cover exactly the {tier} expression lines")
                M = sub.pivot(index="ModelID", columns="entrez", values=col).reindex(index=lines, columns=genes)
            out[(arm, rule)] = M.astype(bool)
    return out


def fit_sha(fitted) -> str:
    h = hashlib.sha256()
    for lab in sorted(fitted.prior):
        h.update(lab.encode())
        h.update(np.round(fitted.prior[lab].sort_index().to_numpy(dtype=float), 12).tobytes())
        h.update(json.dumps(sorted(int(g) for g in fitted.ce[lab])).encode())
    h.update(np.round(fitted.b2c.sort_index().to_numpy(dtype=float), 12).tobytes())
    h.update(json.dumps(sorted(map(str, fitted.fit_lines))).encode())
    return h.hexdigest()


def prepare(cfg: dict, ctx: dict, fit_long: pd.DataFrame, tier: str) -> dict:
    """Everything before ``started``: refit of the fitted baselines on the fitting tiers, every score matrix of
    the evaluation lines, the binary calls and the failed-build rule (§11.4). Reads no evaluation-tier label."""
    from nextgem import metrics, paths, w1

    bad = ev.unit_key_checker(ctx)
    for m in cfg["methods"]:
        if bad(m):
            raise PreflightError(f"method name {m!r} is a {bad(m)}; the scalar guard would refuse the output")
    lines, genes = list(ctx["eval_lines"]), list(ctx["genes"])
    if not lines:
        raise RuntimeError(f"no {tier} expression lines")
    # labels.validation() frames carry no lineage column: attach it from the (label-free) split for every fitting row
    fit_long = fit_long.assign(lineage=fit_long.ModelID.map(ctx["lineage"]))
    if fit_long.lineage.isna().any():
        raise RuntimeError("fitting rows without a lineage in the split")
    fitted = w1.fit_components(fit_long, ctx["line_class"])
    if not fitted.ce[cfg["primary_label"]]:  # review F11: C1b would silently equal C1
        raise RuntimeError("CE_fit is empty on the fitting tiers: C1b would equal C1 (check the fitting labels)")
    gen_S = w1.score_matrices(lines, genes, fitted, ctx["gen"], ctx["line_class"], ctx["omics"]["expr"],
                              random_seeds=cfg["random_seeds"] if cfg["kind"] == "p1" else ())
    clash = (set(gen_S) - {"H1", "H2"}) & set(cfg["methods"])
    if clash:
        raise PreflightError(f"spec method(s) {sorted(clash)} collide with baseline names")
    S = dict(gen_S) if cfg["kind"] == "p1" else {"H1": gen_S["H1"], "H2": gen_S["H2"], "B3": gen_S["B3"]}
    _spec_matrices(cfg, ctx, tier, S, gen_S)
    if cfg["kind"] == "p2":
        S["B6[L_prob]"] = ctx["b6_fit"](fit_long, lines, genes, cfg["b6"], cfg["seed"]).reindex(index=lines, columns=genes)
        if cfg.get("refit"):
            S[cfg["refit"]["method"]] = refit_fit_predict(cfg["refit"], ctx["refit_inputs"](), fit_long, lines,
                                                          ctx["tier_of"], cfg["seed"]).reindex(index=lines, columns=genes)
    for n, M in S.items():
        if not np.isfinite(M.to_numpy(dtype=float)).all():
            raise RuntimeError(f"non-finite scores for {n}")
    src = cfg["sources"][tier]
    fb = src.get("failed_builds")
    if fb and _sha(paths.state_root() / fb) != cfg["inputs"][f"failed_{tier}"]:
        raise RuntimeError("failed-builds file changed since the config was built")
    fdf = load_failed(fb, cfg["methods"])
    failed = {m: sorted(set(fdf.loc[fdf.method == m, "ModelID"]) & set(lines)) for m in cfg["methods"] if m not in ("H1", "H2")}
    S_line, S_gene = dict(S), dict(S)
    for m, fl in failed.items():
        if fl:
            S_line[m] = metrics.no_information_scores(S[m], fl, "line")
            S_gene[m] = metrics.no_information_scores(S[m], fl, "gene")
    calls = _call_matrices(cfg, ctx, tier) if cfg["kind"] == "p1" else {}
    for (m, _), C in calls.items():  # a failed build has no model: every gene is called non-essential (ours)
        if failed.get(m):
            C.loc[failed[m], :] = False
    n_failed_all = {m: sorted(set(fdf.loc[fdf.method == m, "ModelID"]) & set(ctx["tier_lines"])) for m in failed}
    return {"fitted": fitted, "S": S_line, "S_gene": S_gene, "S_raw": S, "calls": calls, "failed": failed,
            "failed_tier": n_failed_all, "fit_sha": fit_sha(fitted), "tier": tier}


# ============================================================================= statistics blocks
def _lt(S, m, Y):
    from nextgem import metrics

    return metrics.per_line(S[m].loc[Y.index, Y.columns], Y)


def contrast_est(lt_a, lt_b, lineage, *, B, seed, method="percentile") -> dict:
    """Paired Δ̄_lin(a − b) with the lineage-stratified line bootstrap CI (§6.5); estimates only."""
    from nextgem import metrics

    d = metrics.paired_line_delta(lt_a, lt_b).dropna()
    if d.empty:
        return {"estimate": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "line_mean": np.nan, "n_lines": 0, "n_lineages": 0}
    bt = metrics.bootstrap_mean(d, lineage, B=B, seed=seed, method=method)
    return {"estimate": bt.estimate, "ci_lo": bt.ci_low, "ci_hi": bt.ci_high, "line_mean": metrics.line_weighted_mean(d),
            "n_lines": bt.n_lines, "n_lineages": bt.n_lineages}


def confirmatory(delta, lineage, contrast: str, sesoi: float, cfg: dict) -> dict:
    """C1/C1b/C4/C5 (§6.3): sign-flip permutation within lineage (n_flip flips), p = (1 + #)/(1 + n_flip), and the
    lineage-stratified bootstrap CI (B_conf), both from the spawned seed; plus the §6.4 transfer statistics."""
    from nextgem import metrics

    return metrics.confirmatory_contrast(delta, lineage, contrast=contrast, sesoi=sesoi, seed=cfg["seed"],
                                         B=cfg["B_conf"], n_flip=cfg["n_flip"], method=cfg["bootstrap_method"],
                                         alpha=cfg["alpha"])


def c3(S_p1, S_b3, Y, lineage, cfg: dict, workers: int | None = None) -> dict:
    """C3 (§6.3): bootstrap CI of mΔAP_gene(P1) − mΔAP_gene(B3) and the 97.5th percentile of the within-lineage
    label-permutation null (n_perm permutations). With ``workers`` == 1 this is ``metrics.c3_decision`` itself;
    with more, :func:`c3_parallel` (bitwise identical; tested)."""
    from nextgem import metrics

    a, b = S_p1.loc[Y.index, Y.columns], S_b3.loc[Y.index, Y.columns]
    w = C3_WORKERS if workers is None else workers
    if w <= 1:
        return metrics.c3_decision(a, b, Y, lineage, seed=cfg["seed"], B=cfg["B_conf"], n_perm=cfg["n_perm"])
    return c3_parallel(a, b, Y, lineage, seed=cfg["seed"], B=cfg["B_conf"], n_perm=cfg["n_perm"], workers=w)


_C3G: dict = {}  # read-only state shared with forked workers (set before the pool is created)


def _c3_boot_batch(s0: int):
    """One bootstrap batch of metrics.gene_bootstrap, rows s0 .. s0 + bs of the shared count matrix (same arrays,
    same shapes and the same operations as the serial loop, so the replicate values are bitwise equal)."""
    from nextgem import metrics

    g = _C3G
    cb = g["C"][s0: s0 + g["bs"]]
    vals, keep = [], None
    for st, y in zip(g["structs"], g["ys"]):
        t, n, p = st.tap(cb, y, g["h"])
        keep = (p > 0) & (p < n)
        e = np.where(keep, metrics.expected_ap_random(np.where(keep, n, 1), np.where(keep, p, 1)), np.nan)
        d = np.where(keep, t - e, np.nan)
        with np.errstate(invalid="ignore"):
            vals.append(np.nanmean(d, axis=1))
    return s0, vals[0] - (vals[1] if len(vals) > 1 else 0.0), (~keep).sum(1)


def _c3_perm_chunk(idx):
    """Permutations ``idx`` of metrics.c3_permutation_null, each from its own spawned child seed."""
    from nextgem import metrics

    g = _C3G
    out = []
    for i in idx:
        Yp = metrics.permute_within_lineage(g["Yn"], g["lin"], np.random.default_rng(g["children"][i]))
        t, _, _ = g["st_null"].tap(g["ones"], g["st_null"].sorted_labels(Yp), g["h_null"])
        out.append((i, np.mean(t[0] - g["e"])))
    return out


def c3_parallel(scores_p1, scores_b3, labels, lineage, *, seed, B: int, n_perm: int, workers: int,
                conf: float = 0.95) -> dict:
    """``metrics.c3_decision`` with the bootstrap batches and the permutations spread over ``workers`` forked
    processes (≤ 8). Nothing in metrics.py changes: this re-runs its steps with its own private helpers —
    the bootstrap count matrix is drawn once from the spawned bootstrap stream and cut into the serial batches
    (``_batch_size``), and permutation i uses the i-th child of the spawned permutation stream — so every replicate,
    the null, the CI, the threshold and p are bitwise identical to the serial function
    (tests/test_evaluate_test.py::test_c3_parallel_bitwise_equals_serial)."""
    import multiprocessing

    from nextgem import metrics as M

    workers = max(1, min(int(workers), MAX_WORKERS))
    M.assert_same_pairs(labels, scores_p1, scores_b3)
    for s in (scores_p1, scores_b3):
        M._align(s, labels)
    ss_boot, ss_perm = M._seed_sequence(seed).spawn(2)
    # --- gene_bootstrap (P1 − B3)
    ga = M.per_gene(scores_p1, labels)
    elig = ga.index[ga["estimable"]]
    est = M.m_delta_ap_gene(ga) - M.m_delta_ap_gene(M.per_gene(scores_b3, labels))
    lin = M._lineage_of(labels.index, lineage)
    Yb = labels.loc[:, elig].to_numpy(dtype=float)
    Mb = ~np.isnan(Yb)
    structs = [M._SortedColumns(s.loc[labels.index, elig].to_numpy(float), Mb) for s in (scores_p1, scores_b3)]
    bs = M._batch_size(*Yb.shape)
    # --- c3_permutation_null (P1)
    g = M.per_gene(scores_p1, labels)
    elig_n = g.index[g["estimable"]]
    obs = M.m_delta_ap_gene(g)
    Yn = labels.loc[:, elig_n].to_numpy(dtype=float)
    _C3G.clear()
    _C3G.update(C=M.stratified_bootstrap_counts(lin, B, ss_boot), bs=bs, structs=structs,
                ys=[st.sorted_labels(Yb) for st in structs], h=M._harmonic_table(len(lin) + 1),
                Yn=Yn, lin=M._lineage_of(labels.index, lineage).to_numpy(), children=M._seed_sequence(ss_perm).spawn(n_perm),
                st_null=M._SortedColumns(scores_p1.loc[labels.index, elig_n].to_numpy(float), ~np.isnan(Yn)),
                h_null=M._harmonic_table(len(lin) + 1), ones=np.ones((1, len(lin))),
                e=g.loc[elig_n, "e_random"].to_numpy())
    if __name__ not in sys.modules:  # loaded by file path (tests): make the worker functions picklable by reference
        import types

        mod = types.ModuleType(__name__)
        mod.__dict__.update({"_c3_boot_batch": _c3_boot_batch, "_c3_perm_chunk": _c3_perm_chunk})
        sys.modules[__name__] = mod
    starts = list(range(0, B, bs))
    chunks = [list(range(i, min(i + 25, n_perm))) for i in range(0, n_perm, 25)]
    reps, dropped, null = np.empty(B), np.empty(B), np.empty(n_perm)
    try:
        with multiprocessing.get_context("fork").Pool(workers) as pool:
            for s0, r, dr in pool.imap_unordered(_c3_boot_batch, starts):
                reps[s0: s0 + bs] = r
                dropped[s0: s0 + bs] = dr
            for part in pool.imap_unordered(_c3_perm_chunk, chunks):
                for i, v in part:
                    null[i] = v
    finally:
        _C3G.clear()
    lo, hi = np.quantile(reps, [(1 - conf) / 2, 1 - (1 - conf) / 2])
    thr = float(np.percentile(null, M.C3_NULL_PERCENTILE))
    tol = 1e-12
    perm_passed = bool(obs > thr + tol)
    boot_ok = bool(float(lo) > 0)
    return {"contrast": "C3", "estimate": float(est), "ci_low": float(lo), "ci_high": float(hi),
            "bootstrap_passed": boot_ok, "m_delta_ap_gene_p1": obs, "null_threshold": thr,
            "permutation_p": float((1 + np.sum(null >= obs - tol)) / (1 + n_perm)), "permutation_passed": perm_passed,
            "passed": bool(boot_ok and perm_passed), "n_genes_eligible": int(len(elig)),
            "mean_genes_dropped": float(dropped.mean())}


def transfer(delta, lineage) -> dict:
    """§6.4 transfer, scalars only (the per-lineage means are sealed): lineage t-interval (df = K − 1), exact
    sign-flip p over 2^K, k-of-K with k_required = ceil(0.7 K)."""
    from nextgem import metrics

    t = metrics.transfer_summary(delta, lineage)
    return {"t_mean": t["mean"], "t_low": t["t_low"], "t_high": t["t_high"], "exact_sign_flip_p": t["exact_sign_flip_p"],
            "k_positive": t["k_positive"], "k_required": t["k_required"], "n_lineages": t["K"],
            "generalises": bool(t["generalises"])}


def gate_table(results: dict) -> list[dict]:
    """§6.3 fixed sequence and §6.4 wording over the contrasts run so far (a prefix of C1 → C1b → C3 → C4 → C5).
    Contrasts not yet run are 'pending' while every earlier one is confirmed, else 'estimate only'."""
    from nextgem import metrics

    run = [c for c in metrics.GATE_ORDER if results.get(c) is not None]
    if run != list(metrics.GATE_ORDER[:len(run)]):
        raise ValueError(f"contrasts {run} are not a prefix of {metrics.GATE_ORDER}")
    rows = metrics.gatekept_claims({c: results[c] for c in run}).to_dict("records") if run else []
    alive = all(r["status"] == "confirmed" for r in rows)
    for c in metrics.GATE_ORDER[len(run):]:
        rows.append({"contrast": c, "passed": None, "confirmatory": False, "status": "pending" if alive else "estimate only",
                     "estimate": np.nan, "ci_low": np.nan, "sesoi": None,
                     "claim": "pending (T2)" if alive else "estimate only (earlier gate failed)"})
    for r in rows:
        for k in ("estimate", "ci_low", "sesoi"):
            if r.get(k) is not None and not isinstance(r[k], str):
                r[k] = None if (isinstance(r[k], float) and math.isnan(r[k])) else float(r[k])
        if r.get("passed") is not None:
            r["passed"] = bool(r["passed"])
        r["confirmatory"] = bool(r["confirmatory"])
    return rows


def _binary(C: pd.DataFrame, T: pd.DataFrame, lineage, *, B, seed):
    """Per-line confusion and §6.8 summaries (w1.binary_endpoints with a line x gene call matrix)."""
    from nextgem import metrics

    pl = metrics.binary_per_line(C.loc[T.index, T.columns], T)
    lin = lineage.reindex(pl.index)
    summ = metrics.summarise_binary(pl, lin)
    cis = {}
    for m in metrics.BINARY_METRICS:
        v = pl[m].dropna()
        if len(v):
            bt = metrics.bootstrap_mean(v, lin, B=B, seed=seed)
            cis[m] = (bt.ci_low, bt.ci_high)
    summ["ci_lo"] = [cis.get(m, (np.nan,) * 2)[0] if p == "exclude" else np.nan for m, p in zip(summ.metric, summ.nan_policy)]
    summ["ci_hi"] = [cis.get(m, (np.nan,) * 2)[1] if p == "exclude" else np.nan for m, p in zip(summ.metric, summ.nan_policy)]
    return summ, pl


def _parent(cfg, m):
    return m if m in ("H1", "H2") else cfg["arms"][m]


def binary_section(cfg, ctx, prep, Yall) -> pd.DataFrame:
    """§6.8 endpoints for generic H1/H2 (G* among them) and the P1 arms under both paper protocols and the
    fixed cutoff; paired MCC P1-arm − G* and P1-arm − its generic parent (estimates)."""
    from nextgem import metrics, w1

    lineage, rows, pls = ctx["lineage"], [], {}
    methods = ["H1", "H2", *cfg["arms"]]
    for proto, rule, lab in ev.BINARY_PROTOCOLS:
        for m in methods:
            g = w1.parent_gene_set(ctx["gen"], _parent(cfg, m))
            summ, pl = _binary(prep["calls"][(m, rule)].loc[:, g], Yall[lab].loc[:, g], lineage, B=cfg["B"], seed=cfg["seed"])
            pls[(proto, rule, lab, m)] = pl
            rows.append(summ.assign(protocol=proto, method=m, call_rule=rule, label=lab, gene_set="parent"))
    G = cfg["G_star_short"]
    for proto, rule, lab in ev.PAIRED_PROTOCOLS:
        for arm in cfg["arms"]:
            for ref in dict.fromkeys([G, _parent(cfg, arm)]):
                d = (pls[(proto, rule, lab, arm)]["mcc"] - pls[(proto, rule, lab, ref)]["mcc"]).dropna()
                if not len(d):
                    continue
                bt = metrics.bootstrap_mean(d, lineage.reindex(d.index), B=cfg["B"], seed=cfg["seed"])
                rows.append(pd.DataFrame([{"protocol": proto, "method": f"{arm} - {ref}", "call_rule": rule, "label": lab,
                                           "gene_set": "parent", "metric": "mcc", "nan_policy": "exclude",
                                           "lineage_equal_mean": bt.estimate, "ci_lo": bt.ci_low, "ci_hi": bt.ci_high,
                                           "n_lines_used": len(d)}]))
    return pd.concat(rows, ignore_index=True).assign(section="binary_6.8")


def _cells(cfg, ctx, prep, Ylong):
    from nextgem import w1

    lab, lines, genes = cfg["primary_label"], ctx["eval_lines"], ctx["genes"]
    Y = w1.wide(Ylong, lab, lines, genes)
    Ynoce = Y.loc[:, [g for g in Y.columns if g not in prep["fitted"].ce[lab]]]
    Yexp = w1.apply_stratum(Y, ("pairs", "expressed"), ctx["lf"], ctx["omics"]["expr"], ctx["omics"]["cn"], None)
    return Y, Ynoce, Yexp


def _row(section, name, **kw):
    return {"section": section, "name": name, **kw}


def _boolstr(x):
    return "" if x is None else str(bool(x))


def t1_reveal(cfg: dict, ctx: dict, prep: dict, Ylong: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """The revealed part of T1 (and of the rehearsal on validation): C1, C1b, C3 (+ expressed-pairs C3 for the
    claim condition), transfer, C2, B8±/B8′, within-parent, headroom, tAP of P1/G*, §6.8 binary endpoints, failed
    builds per method and lineage, counts. Scalars only (lineage names only in the failed_builds section)."""
    from nextgem import metrics, w1

    lab, P, G, lineage = cfg["primary_label"], cfg["primary"], cfg["G_star_short"], ctx["lineage"]
    S, Sg = prep["S"], prep["S_gene"]
    Y, Ynoce, Yexp = _cells(cfg, ctx, prep, Ylong)
    res, rows = {}, []
    for name, YY, sesoi in (("C1", Y, cfg["sesoi"]), ("C1b", Ynoce, cfg["sesoi_c1b"])):
        d = metrics.paired_line_delta(_lt(S, P, YY), _lt(S, G, YY))
        res[name] = confirmatory(d, lineage, name, sesoi, cfg)
        if name == "C1":
            tr = transfer(d, lineage)
    res["C3"] = c3(Sg[P], S["B3"], Y, lineage, cfg)
    c3e = c3(Sg[P], S["B3"], Yexp, lineage, cfg)
    gates = gate_table({k: res[k] for k in ("C1", "C1b", "C3")})
    status = {r["contrast"]: r for r in gates}
    for name in ("C1", "C1b"):
        r = res[name]
        rows.append(_row("confirmatory", name, label=lab, geneset="all" if name == "C1" else "noCE", stratum="none",
                         method_a=P, method_b=G, estimate=r["estimate"], ci_lo=r["ci_low"], ci_hi=r["ci_high"],
                         p=r["p_sign_flip"], passed=_boolstr(r["passed"]), sesoi=r["sesoi"], line_mean=r["line_weighted"],
                         n_lines=r["n_lines"], n_lineages=r["n_lineages"], status=status[name]["status"],
                         claim=status[name]["claim"]))
    for name, r, st in (("C3", res["C3"], status["C3"]["status"]), ("C3_expressed_pairs", c3e, "estimate (C3 claim condition)")):
        rows.append(_row("confirmatory", name, label=lab, geneset="all", stratum="none" if name == "C3" else "expressed_pairs",
                         method_a=P, method_b="B3", estimate=r["estimate"], ci_lo=r["ci_low"], ci_hi=r["ci_high"],
                         p=r["permutation_p"], passed=_boolstr(r["passed"]), m_delta_ap_gene_p1=r["m_delta_ap_gene_p1"],
                         null_threshold=r["null_threshold"], bootstrap_passed=_boolstr(r["bootstrap_passed"]),
                         permutation_passed=_boolstr(r["permutation_passed"]), n_genes_eligible=r["n_genes_eligible"],
                         mean_genes_dropped=r["mean_genes_dropped"], status=st,
                         claim=status["C3"]["claim"] if name == "C3" else ""))
    rows.append(_row("transfer", "C1", label=lab, geneset="all", stratum="none", method_a=P, method_b=G,
                     **{k: (_boolstr(v) if k == "generalises" else v) for k, v in tr.items()}))
    # estimates (§6.3 list; B = B_OTHER, percentile)
    est = {}
    lts = {m: _lt(S, m, Y) for m in dict.fromkeys([P, G, *cfg["arms"], "H1", "H2", f"B2[{lab}]", f"B2m[{lab}]", "B2c",
                                                    f"B8+_{G}", f"B8-_{G}", f"B8'_{G}"])}
    pairs = [(f"C2:B2", P, f"B2[{lab}]"), ("C2:B2m", P, f"B2m[{lab}]"), ("C2:B2c", P, "B2c"),
             ("B8+", P, f"B8+_{G}"), ("B8-", P, f"B8-_{G}"), ("B8'", P, f"B8'_{G}")]
    pairs += [(f"within_parent:{par}", arm, par) for arm, par in cfg["arms"].items()]
    pairs += [(f"vs_G*:{arm}", arm, G) for arm in cfg["arms"] if arm != P]
    for name, a, b in pairs:
        est[name] = contrast_est(lts[a], lts[b], lineage, B=cfg["B"], seed=cfg["seed"])
        rows.append(_row("estimates", name, label=lab, geneset="all", stratum="none", method_a=a, method_b=b, **est[name]))
    hc = metrics.headroom_closed(metrics.paired_line_delta(lts[P], lts[G]), metrics.paired_line_delta(lts[f"B2[{lab}]"], lts[G]),
                                 lineage, B=cfg["B"], seed=cfg["seed"])
    rows.append(_row("estimates", "headroom_closed", label=lab, geneset="all", stratum="none", method_a=P,
                     method_b=f"B2[{lab}]", estimate=hc["estimate"] if hc["reportable"] else np.nan,
                     ci_lo=hc["ci_low"] if hc["reportable"] else np.nan, ci_hi=hc["ci_high"] if hc["reportable"] else np.nan,
                     denominator=hc["denominator"], den_ci_lo=hc["denominator_ci"][0], den_ci_hi=hc["denominator_ci"][1],
                     reportable=_boolstr(hc["reportable"])))
    for m in dict.fromkeys([P, G, *cfg["arms"]]):
        lt = lts[m]
        e = lt[lt.estimable]
        bt = metrics.bootstrap_mean(e.tap, lineage, B=cfg["B"], seed=cfg["seed"])
        el = metrics.lineage_equal_mean(e.e_random, lineage)
        rows.append(_row("estimates", f"tap:{m}", label=lab, geneset="all", stratum="none", method=m, tap_lin=bt.estimate,
                         tap_lin_lo=bt.ci_low, tap_lin_hi=bt.ci_high, e_random_lin=el, lift_lin=bt.estimate / el,
                         n_lines=bt.n_lines, n_lineages=bt.n_lineages))
    Yall = {lab_: w1.wide(Ylong, lab_, ctx["eval_lines"], ctx["genes"]) for lab_ in {lab, "L_ge06", "L_prob"}}
    binary = binary_section(cfg, ctx, prep, Yall)
    fb = failed_rows(cfg, ctx, prep)
    rows.append(_row("counts", "lines", n_lines=len(ctx["eval_lines"]),
                     n_lines_no_expression=len(set(ctx["tier_lines"]) - set(ctx["eval_lines"])),
                     n_lineages=int(lineage.reindex(ctx["eval_lines"]).nunique())))
    out = pd.concat([pd.DataFrame(rows), binary, fb], ignore_index=True)
    out = out[["section"] + [c for c in out.columns if c != "section"]]
    for c in out.columns:
        if c in STR_COLS:
            out[c] = out[c].astype(object)
    b8 = [est[k] for k in ("B8+", "B8-", "B8'")]
    doc = {"look": cfg["look"], "tier": cfg["tier"], "primary": P, "G_star": cfg["G_star"], "label": lab,
           "gates": gates,
           "contrasts": {k: {"contrast": k, "estimate": res[k]["estimate"], "ci_low": res[k]["ci_low"],
                             "ci_high": res[k]["ci_high"], "passed": bool(res[k]["passed"]),
                             "p": res[k].get("p_sign_flip", res[k].get("permutation_p")),
                             "sesoi": None if k == "C3" else res[k]["sesoi"]} for k in ("C1", "C1b", "C3")},
           "c3_claim_condition": metrics.c3_claim_condition(res["C3"], c3e),
           "transfer_C1": tr, "generalises_to_new_tissues": tr["generalises"],
           "adds_beyond_generic_expression": {
               "rule": "pre-declared claim rule: the paired Δ̄_lin(P1 − B8x) CI lower bound > 0 for all three of B8+, "
                       "B8−, B8′ (§6.3, §8); the all-point-estimates-positive version is descriptive only",
               "claim": bool(all(x["ci_lo"] > 0 for x in b8)),
               "ci_lower_all_gt_0": bool(all(x["ci_lo"] > 0 for x in b8)),
               "estimate_all_gt_0_descriptive": bool(all(x["estimate"] > 0 for x in b8))},
           "headline": {"contrast": "C1", "estimate": res["C1"]["estimate"], "ci": [res["C1"]["ci_low"], res["C1"]["ci_high"]],
                        "p_sign_flip": res["C1"]["p_sign_flip"], "claim": status["C1"]["claim"],
                        "conditional_on": f"the {tr['n_lineages']} {cfg['tier']} lineages"},
           "with_C1b": {"estimate": res["C1b"]["estimate"], "ci": [res["C1b"]["ci_low"], res["C1b"]["ci_high"]],
                        "claim": status["C1b"]["claim"]},
           "post_hoc": bool(cfg.get("post_hoc"))}
    assert_revealed_scalar_only(out, ctx)
    return out, _jsonable(doc)


def failed_rows(cfg, ctx, prep) -> pd.DataFrame:
    """§11.4: failed-build counts per method and lineage (revealed at T1; label-free), > 5% flagged as failing."""
    lin = ctx["lineage"]
    rows = []
    tl = pd.Index(ctx["tier_lines"])
    for m, fl in prep["failed_tier"].items():
        n_tot = len(tl)
        frac = len(fl) / n_tot if n_tot else np.nan
        rows.append({"section": "failed_builds", "name": "failed_builds", "method": m, "lineage": "all",
                     "n_failed": len(fl), "n_lines": n_tot, "frac_failed": frac, "flag_failing": _boolstr(frac > 0.05)})
        for k, idx in pd.Series(tl, index=tl).groupby(lin.reindex(tl).to_numpy()).groups.items():
            nk = len(idx)
            f = len(set(fl) & set(idx))
            rows.append({"section": "failed_builds", "name": "failed_builds", "method": m, "lineage": str(k),
                         "n_failed": f, "n_lines": nk, "frac_failed": f / nk if nk else np.nan, "flag_failing": ""})
    return pd.DataFrame(rows, columns=["section", "name", "method", "lineage", "n_failed", "n_lines", "frac_failed",
                                       "flag_failing"])


def assert_revealed_scalar_only(out: pd.DataFrame, ctx: dict) -> None:
    """Revealed outputs carry scalars only (PI decision 21, §10.2): known sections; every column is a known
    string or numeric aggregate column; no cell or column names a line ID or a gene; lineage names only in the
    ``lineage`` column of the label-free failed_builds section; one row per scalar family."""
    if "section" not in out.columns or not set(out.section) <= set(REVEALED_SECTIONS):
        raise RuntimeError("unknown output section")
    bad = ev.unit_key_checker(ctx)
    lineages = {str(x) for x in ctx["lineage"].dropna().unique()}
    for c in out.columns:
        if bad(str(c)):
            raise RuntimeError(f"column {c!r} is a {bad(str(c))}")
        if c in STR_COLS:
            vals = out[c].dropna().astype(str)
            if c == "lineage":
                fb = out.section == "failed_builds"
                if out.loc[~fb, c].notna().any():
                    raise RuntimeError("lineage names are allowed only in the failed_builds section")
                if not set(vals) <= lineages | {"all"}:
                    raise RuntimeError("failed_builds.lineage holds a non-lineage value")
                continue
            for v in vals.unique():
                if bad(v):
                    raise RuntimeError(f"column {c!r} holds a {bad(v)} value")
        elif c in NUM_COLS:
            if not (pd.api.types.is_numeric_dtype(out[c]) or out[c].isna().all()):
                raise RuntimeError(f"column {c!r} must be numeric")
        else:
            raise RuntimeError(f"column {c!r} is not a known scalar column (per-unit columns are refused)")
    key = [c for c in out.columns if c in STR_COLS and c not in ("status", "claim", "passed", "bootstrap_passed",
                                                                    "permutation_passed", "generalises", "reportable",
                                                                    "flag_failing")]
    if out[key].astype(str).duplicated().any():
        raise RuntimeError("more than one row per scalar family")


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    return o


def assert_stdout_scalar_only(text: str, ctx: dict) -> None:
    """stdout of a look: no line ID, gene or lineage anywhere (keys or values)."""
    bad = ev.unit_key_checker(ctx)
    if re.search(r"ACH-\d+", text):
        raise RuntimeError("stdout names a line ID")

    def walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                walk(k)
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
        elif isinstance(o, str) and bad(o):
            raise RuntimeError(f"stdout holds a {bad(o)} value")

    walk(json.loads(text))


# ============================================================================= sealed bundle (T1, T3; rehearsal)
def t1_sealed(cfg: dict, ctx: dict, prep: dict, Ylong: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """§10.2 sealed at T1: per-lineage and per-gene tables, the full W1 benchmark, generic models on every tier line,
    §5.4 strata (incl. library and organoid), attribution contrasts (i)-(v), the Human1-paper threshold grid, the
    per-task table, the failure account (§10.3) and its 20 worst lineage x gene cells, failed-build sensitivities."""
    from nextgem import metrics, w1

    lab, P, G, lineage = cfg["primary_label"], cfg["primary"], cfg["G_star_short"], ctx["lineage"]
    S, Sg, fitted = prep["S"], prep["S_gene"], prep["fitted"]
    Y, _, _ = _cells(cfg, ctx, prep, Ylong)
    B, seed = cfg["B"], cfg["seed"]
    out = {}
    focus = list(dict.fromkeys([P, *cfg["arms"], "H1", "H2", f"B2[{lab}]", f"B2m[{lab}]", "B2c", "B3",
                                f"B8+_{G}", f"B8-_{G}", f"B8'_{G}"]))
    lts = {m: _lt(S, m, Y) for m in focus}
    out["per_lineage"] = pd.concat([metrics.per_lineage_summary(lts[m], lineage, reference=lts[G]).reset_index().assign(method=m)
                                    for m in focus], ignore_index=True)
    out["per_gene"] = pd.concat([metrics.per_gene(Sg[m].loc[Y.index, Y.columns], Y).reset_index(names="entrez").assign(method=m)
                                 for m in focus], ignore_index=True)
    # the full W1 benchmark grid (every W1 method + every spec method; deltas vs H1 and H2)
    Yall = {l_: w1.wide(Ylong, l_, ctx["eval_lines"], ctx["genes"]) for l_ in cfg["labels"]}
    spec_m = [m for m in cfg["methods"] if m not in ("H1", "H2")]
    w1_names = [n for n in S if n not in spec_m]
    gcfg = {k: cfg[k] for k in ("labels", "genesets", "universes", "strata", "B", "seed")}
    rows, _ = ev.tap_grid(Yall, S, lambda l_: w1.methods_for(l_, w1_names) + spec_m, fitted, ctx, gcfg)
    out["w1_benchmark"] = pd.DataFrame(rows)
    # generic models and line-constant baselines on every tier line (incl. lines without expression, §4.3)
    tl = list(ctx["tier_lines"])
    Sall = w1.score_matrices(tl, ctx["genes"], fitted, ctx["gen"], ctx["line_class"], ctx["omics"]["expr"], random_seeds=())
    Yt = w1.wide(Ylong, lab, tl, ctx["genes"])
    gr = []
    for m in ("H1", "H2", "B0", f"B2[{lab}]", f"B2m[{lab}]", "B2c"):
        r, _, _ = w1.method_summary(Sall[m].loc[Yt.index, Yt.columns], Yt, lineage, seed=seed, B=B)
        gr.append({**r, "method": m, "label": lab, "lines": "all tier lines"})
    out["generic_all_lines"] = pd.DataFrame(gr)
    # §5.4 strata of the C1 contrast (estimates)
    sr = []
    for st, spec in ctx["strata"].items():
        if st in ("none", "selective_genes"):
            continue
        Ys = w1.apply_stratum(Y, spec, ctx["lf"], ctx["omics"]["expr"], ctx["omics"]["cn"], None)
        if Ys.empty or Ys.notna().sum().sum() == 0:
            sr.append({"stratum": st, "estimate": np.nan, "n_lines": 0})
            continue
        sr.append({"stratum": st, **contrast_est(_lt(S, P, Ys), _lt(S, G, Ys), lineage, B=B, seed=seed)})
    out["strata"] = pd.DataFrame(sr).assign(method_a=P, method_b=G, label=lab)
    # attribution (§10.1 (i)-(v)); 'lever not selected' is reported as such (A12c)
    ar = []
    for k, v in cfg["attribution"].items():
        if v == LEVER_NOT_SELECTED:
            ar.append({"contrast": k, "name": ATTRIBUTION[k], "status": LEVER_NOT_SELECTED})
            continue
        a, b = (G if x == "G*" else x for x in v)
        ar.append({"contrast": k, "name": ATTRIBUTION[k], "status": "estimated", "method_a": a, "method_b": b,
                   **contrast_est(_lt(S, a, Y), _lt(S, b, Y), lineage, B=B, seed=seed)})
    out["attribution"] = pd.DataFrame(ar)
    # Human1-paper robustness grid (21 gene-effect thresholds), P57 calls
    ge = w1.wide(Ylong, "gene_effect", ctx["eval_lines"], ctx["genes"])
    tg = []
    for m in ["H1", "H2", *cfg["arms"]]:
        g = w1.parent_gene_set(ctx["gen"], _parent(cfg, m))
        grid = metrics.binary_threshold_grid(prep["calls"][(m, "P57")].loc[ge.index, g], ge.loc[:, g])
        for thr, sub in grid.groupby("threshold"):
            s = metrics.summarise_binary(sub.drop(columns="threshold"), lineage.reindex(sub.index))
            tg.append(s.assign(threshold=thr, method=m, call_rule="P57"))
    out["binary_threshold_grid"] = pd.concat(tg, ignore_index=True)
    out["per_task"] = per_task_table(cfg, ctx, prep, Y)
    fa, worst = failure_account(cfg, ctx, prep, Y)
    out["failure_account"], out["worst_cells"] = fa, worst
    out["sensitivities"] = sensitivities(cfg, ctx, prep, Ylong, Y)
    return out


def per_task_table(cfg, ctx, prep, Y) -> pd.DataFrame:
    """The W3 M4 per-task table (scripts/w3_score_train.py ``per_task_table``) on the evaluation tier, so training
    and test are comparable (review F6). Per task of the primary arm's model, with base = the primary arm's score
    (W3's base was T0) and Y = the primary label:
    n_failing_in_U (KO-failing genes in U); n_failing_not_base_lethal (W3: failing genes with base e < 0.99; the base
    can be line-specific here, so the per-line count's median); precision_lin = lineage-equal mean over lines of the
    fraction of called genes that are positive; precision_new_calls_lin = the same over called pairs with base
    e < 0.99 (for a line-constant base identical to W3's gene list); tap_lin_base_plus_task = lineage-equal tAP_line
    of the base with the task-failing genes set to 1; delta_lin_vs_base / delta_lo / delta_hi = paired Δ̄_lin vs the
    base (w1.delta_summary, B = 2,000, SEEDS[0])."""
    from nextgem import metrics, w1

    tc = load_task_calls(cfg["task_calls"])
    if tc is None:
        return pd.DataFrame([{"task": "", "status": "not provided (spec task_calls is null)"}])
    lineage = ctx["lineage"]
    P = cfg["primary"]
    base = prep["S"][P].loc[Y.index, Y.columns]
    lt0 = metrics.per_line(base, Y)
    rows = []
    for task, sub in tc.groupby("task"):
        fails = set(sub.loc[sub.fails, "entrez"].astype(int))
        called = [g for g in Y.columns if g in fails]
        S = base.copy()
        if called:
            S.loc[:, called] = 1.0
        lt = metrics.per_line(S, Y)
        d = w1.delta_summary(lt, lt0, lineage, seed=cfg["seed"], B=cfg["B"])
        if called:
            Yc = Y.loc[:, called]
            newm = base.loc[:, called] < CALL
            prec = Yc.mean(axis=1)
            prec_new = Yc.where(newm).mean(axis=1)
            n_new = float(newm.sum(axis=1).median())
        else:
            prec = prec_new = pd.Series(np.nan, index=Y.index)
            n_new = 0.0
        rows.append({"parent": cfg["arms"][P], "task": task, "n_failing_in_U": len(called),
                     "n_failing_not_base_lethal": n_new,
                     "precision_lin": metrics.lineage_equal_mean(prec.dropna(), lineage) if prec.notna().any() else np.nan,
                     "precision_new_calls_lin": (metrics.lineage_equal_mean(prec_new.dropna(), lineage)
                                                 if prec_new.notna().any() else np.nan),
                     "tap_lin_base_plus_task": metrics.lineage_equal_mean(lt.tap[lt.estimable], lineage),
                     "delta_lin_vs_base": d["delta_lin"], "delta_lo": d["delta_lo"], "delta_hi": d["delta_hi"]})
    return pd.DataFrame(rows)


def failure_account(cfg, ctx, prep, Y) -> tuple[pd.DataFrame, pd.DataFrame]:
    """§10.3: P1's errors (fixed-cutoff FBA<0.01 calls of the primary arm vs the primary label) by class; the
    20 worst lineage x gene cells (most misclassified lines); the upper bound = the fraction of positives in U
    whose gene is addressable in P1's parent (the KO closes >= 1 reaction; proxy for 'has a metabolic mechanism')."""
    P = cfg["primary"]
    parent = cfg["arms"][P]
    C = prep["calls"][(P, "FBA<0.01")].loc[Y.index, Y.columns]
    m = Y.notna()
    fn = (~C) & (Y == 1) & m
    fp = C & (Y == 0) & m
    expr = ctx["omics"]["expr"].reindex(index=Y.index, columns=Y.columns)
    cn = ctx["omics"]["cn"].reindex(index=Y.index, columns=Y.columns)
    genes = Y.columns
    gmask = lambda s: pd.DataFrame(np.broadcast_to(genes.isin(list(s)), Y.shape), index=Y.index, columns=genes)  # noqa: E731
    classes = {
        "common_essential": gmask(prep["fitted"].ce[cfg["primary_label"]]),
        "cn_amplified": cn >= 4.0,
        "unexpressed": expr < 1.0,
        "missing_from_model": gmask(set(ctx["gen"].index[~ctx["gen"][f"in_{parent}"]])),
    }
    ann = load_annotations(cfg["gene_annotations"])
    for k in FAILURE_GENE_CLASSES:
        if ann is not None and k in ann.columns:
            classes[k] = gmask(set(ann.index[ann[k]]))
    nfn, nfp = int(fn.sum().sum()), int(fp.sum().sum())
    rows = [{"class": "all_errors", "n_fn": nfn, "n_fp": nfp, "frac_of_fn": 1.0, "frac_of_fp": 1.0, "annotated": True}]
    for k in FAILURE_HARNESS_CLASSES + FAILURE_GENE_CLASSES:
        if k not in classes:
            rows.append({"class": k, "annotated": False})
            continue
        a, b = int((fn & classes[k]).sum().sum()), int((fp & classes[k]).sum().sum())
        rows.append({"class": k, "n_fn": a, "n_fp": b, "frac_of_fn": a / nfn if nfn else np.nan,
                     "frac_of_fp": b / nfp if nfp else np.nan, "annotated": True})
    addr = ctx["gene_features"][f"addressable_{'human1' if parent == 'H1' else 'human2'}"]
    pos = (Y == 1) & m
    n_pos = int(pos.sum().sum())
    n_addr = int((pos & gmask(set(addr.index[addr.astype(bool)]))).sum().sum())
    rows.append({"class": "upper_bound_addressable_positives", "n_positive_pairs": n_pos,
                 "n_positive_pairs_addressable": n_addr, "upper_bound": n_addr / n_pos if n_pos else np.nan,
                 "annotated": True})
    fa = pd.DataFrame(rows)
    lin = ctx["lineage"].reindex(Y.index)
    err = (fn | fp).astype(int)
    cells = pd.DataFrame({"n_lines": m.groupby(lin).sum().stack(future_stack=True), "n_errors": err.groupby(lin).sum().stack(future_stack=True),
                          "n_fn": fn.astype(int).groupby(lin).sum().stack(future_stack=True), "n_fp": fp.astype(int).groupby(lin).sum().stack(future_stack=True)})
    cells.index.names = ["lineage", "entrez"]
    cells = cells.reset_index()
    cells["frac_errors"] = cells.n_errors / cells.n_lines.replace(0, np.nan)
    worst = cells.sort_values(["n_errors", "frac_errors", "lineage", "entrez"], ascending=[False, False, True, True]).head(N_WORST_CELLS)
    return fa, worst.reset_index(drop=True)


def sensitivities(cfg, ctx, prep, Ylong, Y) -> pd.DataFrame:
    """C1-type estimates: expressed pairs, the other labels, U∩, and the two failed-build sensitivities of §11.4
    ((a) lines failed by any method excluded for all; (b) generic fallback for failed lines)."""
    from nextgem import metrics, w1

    lab, P, G, lineage, S = cfg["primary_label"], cfg["primary"], cfg["G_star_short"], ctx["lineage"], prep["S"]
    B, seed = cfg["B"], cfg["seed"]
    rows = []
    _, _, Yexp = _cells(cfg, ctx, prep, Ylong)
    rows.append({"name": "C1_expressed_pairs", **contrast_est(_lt(S, P, Yexp), _lt(S, G, Yexp), lineage, B=B, seed=seed)})
    for l_ in cfg["labels"]:
        if l_ != lab:
            Yl = w1.wide(Ylong, l_, ctx["eval_lines"], ctx["genes"])
            rows.append({"name": f"C1_{l_}", **contrast_est(_lt(S, P, Yl), _lt(S, G, Yl), lineage, B=B, seed=seed)})
    Yu = Y.loc[:, ctx["ucap"]]
    rows.append({"name": "C1_Ucap", **contrast_est(_lt(S, P, Yu), _lt(S, G, Yu), lineage, B=B, seed=seed)})
    union = sorted(set().union(*prep["failed"].values())) if prep["failed"] else []
    if union:
        Ya = Y.drop(index=union)
        rows.append({"name": "C1_failed_excluded_all", **contrast_est(_lt(S, P, Ya), _lt(S, G, Ya), lineage, B=B, seed=seed)})
        fb = dict(S)
        fl = prep["failed"].get(P, [])
        if fl:
            M = prep["S_raw"][P].copy()
            M.loc[fl, :] = prep["S_raw"][G].loc[fl, :].to_numpy()
            fb[P] = M
        rows.append({"name": "C1_generic_fallback", **contrast_est(_lt(fb, P, Y), _lt(S, G, Y), lineage, B=B, seed=seed)})
    else:
        rows.append({"name": "failed_build_sensitivities", "estimate": np.nan, "n_lines": 0, "note": "no failed build"})
    return pd.DataFrame(rows).assign(method_a=P, method_b=G)


def tsv_bytes(df: pd.DataFrame) -> bytes:
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False, float_format="%.6g")
    return buf.getvalue().encode()


def bundle_bytes(bundle: dict[str, pd.DataFrame]) -> dict[str, bytes]:
    return {f"{k}.tsv": tsv_bytes(v) for k, v in sorted(bundle.items())}


def manifest_sha(hashes: dict[str, str]) -> str:
    return hashlib.sha256("".join(f"{k}\t{v}\n" for k, v in sorted(hashes.items())).encode()).hexdigest()


def write_sealed(look: str, name: str, data: bytes) -> str:
    """Write one sealed file to data/locked/sealed/<look>/<name> through labels.py's protected writer (tripwire
    scope 'sealed'; only evaluate_test.py may write it). Returns its SHA-256 (from the bytes; nothing is read)."""
    from nextgem import labels, paths

    if "/" in name or name.startswith("."):
        raise ValueError(f"bad sealed file name {name!r}")
    labels._write_protected(paths.sealed_dir(look) / name, data, "sealed")
    return hashlib.sha256(data).hexdigest()


# ============================================================================= T2
def t2_reveal(cfg: dict, ctx: dict, prep: dict, Ylong: pd.DataFrame, t1_doc: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """C4 = Δ̄_lin(P2 − G*), C5 = Δ̄_lin(P2 − B6) (§6.3), the full gatekeeping chain with T1's C1/C1b/C3, and
    P2's scalar tAP grid. Returns (confirmatory rows, grid rows, gates doc)."""
    from nextgem import metrics, w1

    lab, P, G, lineage, S = cfg["primary_label"], cfg["primary"], cfg["G_star_short"], ctx["lineage"], prep["S"]
    Y = w1.wide(Ylong, lab, ctx["eval_lines"], ctx["genes"])
    res, rows = {}, []
    for name, ref in (("C4", G), ("C5", "B6[L_prob]")):
        d = metrics.paired_line_delta(_lt(S, P, Y), _lt(S, ref, Y))
        res[name] = confirmatory(d, lineage, name, cfg["sesoi"], cfg)
        res[name + "_transfer"] = transfer(d, lineage)
    prior = {k: t1_doc["contrasts"][k] for k in ("C1", "C1b", "C3")}
    gates = gate_table({**prior, "C4": res["C4"], "C5": res["C5"]})
    status = {r["contrast"]: r for r in gates}
    for name, ref in (("C4", G), ("C5", "B6[L_prob]")):
        r = res[name]
        rows.append(_row("confirmatory", name, label=lab, geneset="all", stratum="none", method_a=P, method_b=ref,
                         estimate=r["estimate"], ci_lo=r["ci_low"], ci_hi=r["ci_high"], p=r["p_sign_flip"],
                         passed=_boolstr(r["passed"]), sesoi=r["sesoi"], line_mean=r["line_weighted"], n_lines=r["n_lines"],
                         n_lineages=r["n_lineages"], status=status[name]["status"], claim=status[name]["claim"]))
        rows.append(_row("transfer", name, label=lab, geneset="all", stratum="none", method_a=P, method_b=ref,
                         **{k: (_boolstr(v) if k == "generalises" else v) for k, v in res[name + "_transfer"].items()}))
    conf = pd.DataFrame(rows)
    for c in conf.columns:
        if c in STR_COLS:
            conf[c] = conf[c].astype(object)
    assert_revealed_scalar_only(conf, ctx)
    Yall = {l_: w1.wide(Ylong, l_, ctx["eval_lines"], ctx["genes"]) for l_ in cfg["labels"]}
    names = ["H1", "H2", *[m for m in cfg["methods"] if m not in ("H1", "H2")]]
    gcfg = {k: cfg[k] for k in ("labels", "genesets", "universes", "strata", "B", "seed")}
    grows, _ = ev.tap_grid(Yall, S, lambda l_: names + (["B6[L_prob]"] if l_ == "L_prob" else []), prep["fitted"], ctx, gcfg)
    grid = pd.DataFrame(grows).assign(section="tAP_grid")
    grid = grid[["section"] + [c for c in grid.columns if c != "section"]]
    nm = {l_: len(names) + (l_ == "L_prob") for l_ in cfg["labels"]}
    ev.assert_scalar_only(grid, gcfg, ctx, n_methods=nm, n_bin=0)
    doc = {"look": "T2", "primary": P, "G_star": cfg["G_star"], "label": lab, "gates": gates,
           "contrasts": {k: {"contrast": k, "estimate": res[k]["estimate"], "ci_low": res[k]["ci_low"],
                             "ci_high": res[k]["ci_high"], "passed": bool(res[k]["passed"]), "p": res[k]["p_sign_flip"],
                             "sesoi": res[k]["sesoi"]} for k in ("C4", "C5")},
           "transfer": {"C4": res["C4_transfer"], "C5": res["C5_transfer"]},
           "headline": t1_doc["headline"], "headline_rule": "the first T1 number stays the headline (§10.2)"}
    return conf, grid, _jsonable(doc)


def t1_record(look: str = "T1") -> tuple[dict, dict]:
    """The committed revealed T1 (or T3) outputs, verified against the test ledger's revealed output SHA-256."""
    from nextgem import ledger, paths

    root = paths.state_root()
    rows = [r for r in ledger.read("test") if r["look_id"] == look and r["event"] == "output"]
    rev = [r for r in rows if ledger.detail_json(r).get("kind") == "revealed"]
    if not rev:
        raise PreflightError(f"no revealed output row for {look} in the test ledger")
    files = ledger.detail_json(rev[-1])["files"]
    data = b"".join((root / REVEAL_DIR / f).read_bytes() for f in files)
    if hashlib.sha256(data).hexdigest() != rev[-1]["output_sha"]:
        raise PreflightError(f"the revealed {look} files differ from their ledger SHA-256")
    for f in files:
        if not ledger.is_committed_clean(root, f"{REVEAL_DIR}/{f}"):
            raise PreflightError(f"{REVEAL_DIR}/{f} is not committed")
    gates = json.loads((root / REVEAL_DIR / f"test_look_{look}_gates.json").read_text())
    sealed = [r for r in rows if ledger.detail_json(r).get("kind") == "sealed"]
    return gates, (ledger.detail_json(sealed[-1]) if sealed else {})


# ============================================================================= running looks
def write_config(cfg: dict, out_dir: Path, prefix: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{prefix}_{cfg['look']}_config.json"
    doc = {k: v for k, v in cfg.items() if not k.startswith("_")}
    p.write_text(json.dumps({**doc, "config_sha": config_sha(cfg)}, indent=1, sort_keys=True))
    return p


def _resume_marker(name: str, look_id: str, csha: str, lk, author: str) -> None:
    from nextgem import ledger

    n_prev = sum(1 for r in ledger.read(name) if r["look_id"] == look_id and r["event"] == "output" and '"resumed"' in r["detail"])
    lk.output(hashlib.sha256(f"{csha}|resumed|{n_prev + 1}".encode()).hexdigest(),
              detail=json.dumps({"event": "resumed", "attempt": n_prev + 1}), author=author)


def run_p1_look(cfg: dict, *, ctx: dict, prep: dict, load_labels, out_dir: Path, ledger_name: str, author: str,
                workstream: str = "") -> dict:
    """The ledger part of the rehearsal's validation stage (validation ledger) and of T1/T3 (test ledger).
    ``prep`` and the config JSON exist before this is called; the evaluation-tier labels are read only after
    the committed ``started`` row."""
    from nextgem import ledger

    csha = config_sha(cfg)
    prefix = "val_look" if ledger_name == "validation" else "test_look"
    cpath = out_dir / f"{prefix}_{cfg['look']}_config.json"
    if not cpath.exists() or json.loads(cpath.read_text()).get("config_sha") != csha:
        raise RuntimeError("config JSON missing or stale: write_config() must run before the look")
    prior = ledger.looks(ledger.read(ledger_name)).get(cfg["look"])
    resumed = prior is not None and prior["started"] is not None and not prior["completed"]
    detail = json.dumps({"milestone": cfg["milestone"], "kind": "test_harness", "tier": cfg["tier"],
                         "fit_tiers": cfg["fit_tiers"], "fit_sha": prep["fit_sha"], "spec_sha": cfg["inputs"]["spec"],
                         "scalars_only": True, "sealed_dir": f"data/locked/sealed/{cfg['look']}",
                         "post_hoc": bool(cfg.get("post_hoc")), "code_sha": code_sha(cfg)}, sort_keys=True)
    with ledger.look(ledger_name, cfg["look"], workstream=workstream or REHEARSAL_WS,
                     config_sha=csha, detail=detail, author=author) as lk:
        if resumed:
            _resume_marker(ledger_name, cfg["look"], csha, lk, author)
        Ylong = load_labels()  # first and only evaluation-tier label read
        out, doc = t1_reveal(cfg, ctx, prep, Ylong)
        tsv = out_dir / f"{prefix}_{cfg['look']}.tsv"
        gj = out_dir / f"{prefix}_{cfg['look']}_gates.json"
        tsv.write_bytes(tsv_bytes(out))
        gj.write_text(json.dumps(doc, indent=1, sort_keys=True))
        rsha = hashlib.sha256(tsv.read_bytes() + gj.read_bytes()).hexdigest()
        lk.output(rsha, detail=json.dumps({"kind": "revealed", "files": [tsv.name, gj.name],
                                           "headline": doc["headline"], "post_hoc": doc["post_hoc"]}), author=author)
        hashes = {name: write_sealed(cfg["look"], name, data) for name, data in bundle_bytes(t1_sealed(cfg, ctx, prep, Ylong)).items()}
        ssha = manifest_sha(hashes)
        lk.output(ssha, detail=json.dumps({"kind": "sealed", "dir": f"data/locked/sealed/{cfg['look']}", "files": hashes},
                                          sort_keys=True), author=author)
    return {"look": cfg["look"], "config_sha": csha, "resumed": resumed,
            "revealed": {"files": [str(tsv), str(gj)], "sha256": rsha},
            "sealed": {"dir": f"data/locked/sealed/{cfg['look']}", "manifest_sha256": ssha, "files": hashes},
            "gates": [{k: r[k] for k in ("contrast", "status", "claim", "estimate", "ci_low")} for r in doc["gates"]],
            "headline": doc["headline"], "post_hoc": doc["post_hoc"]}


def run_t2(cfg: dict, *, ctx: dict, prep: dict, load_labels, out_dir: Path, author: str) -> dict:
    from nextgem import ledger

    csha = config_sha(cfg)
    cpath = out_dir / "test_look_T2_config.json"
    if not cpath.exists() or json.loads(cpath.read_text()).get("config_sha") != csha:
        raise RuntimeError("config JSON missing or stale: write_config() must run before the look")
    t1_doc, _ = t1_record("T1")
    prior = ledger.looks(ledger.read("test")).get("T2")
    resumed = prior is not None and prior["started"] is not None and not prior["completed"]
    detail = json.dumps({"milestone": cfg["milestone"], "kind": "test_harness", "tier": "test", "fit_tiers": cfg["fit_tiers"],
                         "fit_sha": prep["fit_sha"], "spec_sha": cfg["inputs"]["spec"], "scalars_only": True,
                         "code_sha": code_sha(cfg), "refit": cfg.get("refit")}, sort_keys=True)
    with ledger.look("test", "T2", workstream="M7", config_sha=csha, detail=detail, author=author) as lk:
        if resumed:
            _resume_marker("test", "T2", csha, lk, author)
        Ylong = load_labels()
        conf, grid, doc = t2_reveal(cfg, ctx, prep, Ylong, t1_doc)
        tsv = out_dir / "test_look_T2.tsv"
        gt = out_dir / "test_look_T2_grid.tsv"
        gj = out_dir / "test_look_T2_gates.json"
        tsv.write_bytes(tsv_bytes(conf))
        gt.write_bytes(tsv_bytes(grid))
        gj.write_text(json.dumps(doc, indent=1, sort_keys=True))
        rsha = hashlib.sha256(tsv.read_bytes() + gt.read_bytes() + gj.read_bytes()).hexdigest()
        lk.output(rsha, detail=json.dumps({"kind": "revealed", "files": [tsv.name, gt.name, gj.name]}), author=author)
    return {"look": "T2", "config_sha": csha, "resumed": resumed, "revealed": {"files": [str(tsv), str(gt), str(gj)], "sha256": rsha},
            "gates": [{k: r[k] for k in ("contrast", "status", "claim", "estimate", "ci_low")} for r in doc["gates"]],
            "headline": doc["headline"]}


def run_unseal(*, cfg_for, ctx_for, fit_long_for, load_test, author: str) -> dict:
    """§10.2 unseal after T2: (1) the ledger ``unseal`` row (committed; the ledger itself refuses it without a
    completed T2); (2) every sealed bundle (T1, T3) is re-derived with the same harness and test labels and each
    file is published under results/m6/test/unsealed/<look>/ only if its SHA-256 equals the one in the ledger.
    The sealed files themselves are never read (the tripwire has no sealed reader in labels.py). Before the unseal
    row and before any label read, each sealed look's config must reproduce its started config_sha
    (:func:`unseal_problems`; review F2), otherwise PreflightError."""
    from nextgem import ledger, paths

    probs = unseal_problems(cfg_for)
    if probs:
        raise PreflightError("; ".join(probs))
    rows = ledger.read("test")
    if not any(r["event"] == "unseal" for r in rows):
        ledger.append("test", "unseal", look_id="unseal", detail=json.dumps({"looks": list(SEALED_LOOKS)}), author=author)
    looks = ledger.looks(ledger.read("test"))
    report = {}
    Ylong = None
    for look in SEALED_LOOKS:
        if look not in looks or not looks[look]["completed"]:
            continue
        _, sealed = t1_record(look)
        want = sealed.get("files", {})
        cfg = cfg_for(look)
        ctx = ctx_for(cfg)
        prep = prepare(cfg, ctx, fit_long_for(), "test")
        Ylong = load_test() if Ylong is None else Ylong
        got = bundle_bytes(t1_sealed(cfg, ctx, prep, Ylong))
        d = paths.state_root() / UNSEALED_DIR / look
        d.mkdir(parents=True, exist_ok=True)
        ok, bad = [], []
        for name, data in got.items():
            if want.get(name) == hashlib.sha256(data).hexdigest():
                (d / name).write_bytes(data)
                ok.append(name)
            else:
                bad.append(name)
        bad += sorted(set(want) - set(got))
        report[look] = {"published": ok, "sha_mismatch_not_published": bad, "config_sha_now": config_sha(cfg)}
    return report


# ============================================================================= dry run / synthetic labels
def synthetic_labels(lines, genes, priors: dict, seed) -> pd.DataFrame:
    """Synthetic labels of the evaluation shape (lines x U). ``priors``: label -> per-gene positive probability
    (a Series over Entrez: the training prior B2, review F11) or a single rate; each pair is an independent
    Bernoulli draw. gene_effect ~ N(−0.2, 0.5) and dep_prob ~ U(0, 1) for B2c and the threshold grid. No label file
    is read."""
    rng = np.random.default_rng(seed)
    lines, genes = list(lines), list(genes)
    df = pd.MultiIndex.from_product([lines, genes], names=["ModelID", "entrez"]).to_frame(index=False)
    df["gene_effect"] = rng.normal(-0.2, 0.5, len(df))
    df["dep_prob"] = rng.random(len(df))
    for lab, p in priors.items():
        pv = np.tile(p.reindex(genes).to_numpy(dtype=float), len(lines)) if isinstance(p, pd.Series) else float(p)
        df[lab] = (rng.random(len(df)) < pv).astype(float)
    return df


def synthetic_priors(train_long: pd.DataFrame, line_class: pd.Series, genes) -> dict:
    """Per-gene training priors (B2) of every label, a gene without training labels at the label's rate. Drawing
    synthetic validation labels at these priors keeps CE_fit (genes positive in >= 90% of the fitting lines)
    non-empty in the dry run, as with real validation labels (review F11)."""
    from nextgem import w1

    f = w1.fit_components(train_long, line_class)
    return {lab: f.prior[lab].reindex(list(genes)).fillna(f.rate[lab]) for lab in list(w1.GRID_LABELS) + ["L_ge06"]}


def dry_run(cfg: dict, ctx: dict, train_long: pd.DataFrame, tier: str) -> dict:
    """Synthetic run of a look on the evaluation-tier shape: training labels are real; the validation labels of
    the refit (test tier) and the evaluation labels are synthetic draws at the per-gene training priors. No ledger
    row, no file; returns scalars only (the priors themselves are never printed)."""
    from nextgem.config import SEEDS

    priors = synthetic_priors(train_long, ctx["line_class"], ctx["genes"])
    fit_long = train_long
    if tier == "test":
        val_lines = [l for l in ctx["lineage"].index if l not in set(train_long.ModelID) and l not in set(ctx["tier_lines"])]
        fit_long = pd.concat([train_long, synthetic_labels(val_lines, ctx["genes"], priors, SEEDS[1])], ignore_index=True)
    t0 = time.time()
    prep = prepare(cfg, ctx, fit_long, tier)  # raises if CE_fit is empty
    lab = cfg["primary_label"]
    ce = prep["fitted"].ce[lab]
    Ylong = synthetic_labels(ctx["eval_lines"], ctx["genes"], priors, cfg["seed"])
    extra = {"ce_fit_nonempty": bool(ce), "c1b_gene_set_smaller": bool(len(set(ctx["genes"]) - set(ce)) < len(ctx["genes"]))}
    if cfg["kind"] == "p1":
        out, doc = t1_reveal(cfg, ctx, prep, Ylong)
        hashes = {k: hashlib.sha256(v).hexdigest() for k, v in bundle_bytes(t1_sealed(cfg, ctx, prep, Ylong)).items()}
        n_rows = len(out)
        extra["c1b_differs_from_c1"] = bool(doc["contrasts"]["C1"]["estimate"] != doc["contrasts"]["C1b"]["estimate"])
    else:
        t1_doc = {"contrasts": {k: {"contrast": k, "passed": True, "estimate": 0.02, "ci_low": 0.01, "sesoi": s}
                                for k, s in (("C1", cfg["sesoi"]), ("C1b", cfg["sesoi_c1b"]), ("C3", None))},
                  "headline": {"contrast": "C1", "note": "synthetic T1 record (dry run)"}}
        conf, grid, doc = t2_reveal(cfg, ctx, prep, Ylong, t1_doc)
        hashes, n_rows = {}, len(conf) + len(grid)
    return {"look": cfg["look"], "tier": tier, "config_sha": config_sha(cfg), **extra,
            "synthetic_labels": f"Bernoulli per gene at the training prior of each label (not shown), seed {cfg['seed']}",
            "n_eval_lines": len(ctx["eval_lines"]), "methods": sorted(prep["S"]), "n_revealed_rows": n_rows,
            "scalar_guard": "pass", "validation_labels_read": False, "test_labels_read": False, "ledger_rows_written": 0,
            "files_written": 0, "sealed_files_in_memory": sorted(hashes), "synthetic_gates": doc["gates"],
            "compute_s": round(time.time() - t0, 1)}


# ============================================================================= CLI
def _print(obj, ctx=None) -> None:
    text = json.dumps(_jsonable(obj), indent=1, sort_keys=True, default=str)
    if ctx is not None:
        assert_stdout_scalar_only(text, ctx)
    print(text)


def _prov_inputs(cfg: dict) -> list[Path]:
    """Provenance inputs: code from the checkout, everything else from the state root (existing files only)."""
    from nextgem import paths

    root = paths.state_root()
    code = [paths.PROJECT_ROOT / c for c in cfg["_code_files"]]
    other = [root / r for r in list(cfg["_inputs_rel"].values()) + list(cfg.get("provenance", {}))]
    return [p for p in dict.fromkeys(code + other) if p.is_file()]


def preflight_problems(cfg: dict) -> list[str]:
    return ev.preflight(required_files(cfg)) + pinned_problems(cfg) + check_constants(cfg) + ev.check_frozen_full()


def main(argv=None) -> int:
    global C3_WORKERS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    act = ap.add_mutually_exclusive_group(required=True)  # review F9: no default real mode
    act.add_argument("--look", choices=TEST_LOOKS)
    act.add_argument("--rehearsal", action="store_true")
    act.add_argument("--unseal", action="store_true")
    mod = ap.add_mutually_exclusive_group()
    mod.add_argument("--print-config", action="store_true")
    mod.add_argument("--print-frozen-candidates", action="store_true")
    mod.add_argument("--preflight-only", action="store_true")
    mod.add_argument("--dry-run", action="store_true", help="synthetic validation/test labels; no ledger row, no file")
    ap.add_argument("--spec", help="spec path relative to the state root (dry modes only; default the canonical spec)")
    ap.add_argument("--tier", choices=("test", "validation"), default="test", help="--dry-run shape (default test)")
    ap.add_argument("--workers", type=int, default=MAX_WORKERS,
                    help=f"processes for the C3 replicates, 1-{MAX_WORKERS} (results are bitwise identical)")
    a = ap.parse_args(argv)
    if not 1 <= a.workers <= MAX_WORKERS:
        ap.error(f"--workers must be in 1..{MAX_WORKERS}")
    C3_WORKERS = a.workers
    try:  # one BLAS thread also when numpy was imported before this module (e.g. under a test runner)
        from threadpoolctl import threadpool_limits

        _limits = threadpool_limits(1)  # noqa: F841 - kept alive for the whole run
    except Exception:  # noqa: BLE001 - threadpoolctl missing: the environment variables set at import apply
        pass
    t0 = time.time()
    mode = "REHEARSAL" if a.rehearsal else "UNSEAL" if a.unseal else a.look
    dry = a.print_config or a.print_frozen_candidates or a.dry_run
    if a.spec and not dry:
        print("REFUSED: --spec is allowed only with --print-config, --print-frozen-candidates or --dry-run", file=sys.stderr)
        return 2
    if a.dry_run and mode not in ("T1", "T2"):
        print("REFUSED: --dry-run runs the T1 or T2 pipeline (--look T1|T2)", file=sys.stderr)
        return 2
    try:
        cfg = look_config(mode, a.spec)
        if a.dry_run and a.tier == "validation" and mode == "T1":
            cfg = look_config("REHEARSAL", a.spec)
    except PreflightError as e:
        print(f"LOOK REFUSED (no label read, nothing written): {e}", file=sys.stderr)
        return 2
    if a.print_config:
        _print({**{k: v for k, v in cfg.items() if not k.startswith("_")}, "config_sha": config_sha(cfg)})
        return 0
    if a.print_frozen_candidates:
        print("\n".join(frozen_candidates(cfg)))
        return 0
    from nextgem import labels

    if a.dry_run:
        try:
            ctx = real_context(cfg, cfg["tier"])
            res = dry_run(cfg, ctx, labels.train(columns=ev.LABEL_COLS), cfg["tier"])
        except PreflightError as e:
            print(f"DRY RUN REFUSED (no validation or test label read, nothing written): {e}", file=sys.stderr)
            return 2
        _print({**res, "wall_s": round(time.time() - t0, 1)}, ctx)
        return 0
    problems = state_problems(mode, cfg)
    if problems:
        print("REFUSED (look order, DESIGN §10.2; nothing read or written):", file=sys.stderr)
        print("\n".join(f"  - {p}" for p in problems), file=sys.stderr)
        return 2
    problems = preflight_problems(cfg)
    if problems:
        print("PRE-FLIGHT FAILED; nothing was started (no ledger row, no label read):", file=sys.stderr)
        print("\n".join(f"  - {p}" for p in problems), file=sys.stderr)
        return 2
    if a.preflight_only:
        print("pre-flight OK")
        return 0
    from nextgem import paths
    from nextgem.run import provenance

    root = paths.state_root()
    author = "scripts/evaluate_test.py"
    if mode == "REHEARSAL":
        # stage 1: synthetic labels of the test-tier shape (nothing written)
        tcfg = look_config("T1")
        tctx = real_context(tcfg, "test")
        stage1 = dry_run(tcfg, tctx, labels.train(columns=ev.LABEL_COLS), "test")
        out_dir = root / REHEARSAL_DIR
        with provenance(out_dir / "MANIFEST.tsv", inputs=_prov_inputs(cfg), seeds=[cfg["seed"]],
                        allow_dirty=False, author=author) as prov:
            ctx = real_context(cfg, "validation")
            prep = prepare(cfg, ctx, labels.train(columns=ev.LABEL_COLS), "validation")
            cp = write_config(cfg, out_dir, "val_look")
            res = run_p1_look(cfg, ctx=ctx, prep=prep, load_labels=labels.validation, out_dir=out_dir,
                              ledger_name="validation", author=author, workstream=REHEARSAL_WS)
            prov.output(cp, *res["revealed"]["files"])
            prov.note(look_id=cfg["look"], config_sha=res["config_sha"], sealed_manifest_sha=res["sealed"]["manifest_sha256"])
        _print({"stage1_synthetic": {k: stage1[k] for k in ("n_eval_lines", "scalar_guard", "n_revealed_rows",
                                                            "sealed_files_in_memory", "compute_s")},
                "stage2_validation_look": res, "wall_s": round(time.time() - t0, 1)}, ctx)
        return 0
    if mode == "UNSEAL":
        def fit_long_for():
            return pd.concat([labels.train(columns=ev.LABEL_COLS), labels.validation()], ignore_index=True)

        problems = unseal_problems(look_config)  # before the unseal row and before any label read (review F2)
        if problems:
            print("UNSEAL REFUSED (no ledger row, no label read):", file=sys.stderr)
            print("\n".join(f"  - {p}" for p in problems), file=sys.stderr)
            return 2
        with provenance(root / UNSEALED_DIR / "MANIFEST.tsv", inputs=_prov_inputs(cfg), seeds=[cfg["seed"]],
                        allow_dirty=False, author=author) as prov:
            rep = run_unseal(cfg_for=lambda look: look_config(look), ctx_for=lambda c: real_context(c, "test"),
                             fit_long_for=fit_long_for, load_test=labels.test, author=author)
            prov.output(*[root / UNSEALED_DIR / lk / f for lk, v in rep.items() for f in v["published"]])
            prov.note(unseal={lk: {"published": len(v["published"]), "mismatch": v["sha_mismatch_not_published"]}
                              for lk, v in rep.items()})
        _print({"unseal": rep, "wall_s": round(time.time() - t0, 1)})
        return 0 if all(not v["sha_mismatch_not_published"] for v in rep.values()) else 2
    out_dir = root / REVEAL_DIR
    with provenance(out_dir / "MANIFEST.tsv", inputs=_prov_inputs(cfg), seeds=[cfg["seed"]],
                    allow_dirty=False, author=author) as prov:
        ctx = real_context(cfg, "test")
        fit_long = pd.concat([labels.train(columns=ev.LABEL_COLS), labels.validation()], ignore_index=True)
        prep = prepare(cfg, ctx, fit_long, "test")
        cp = write_config(cfg, out_dir, "test_look")
        if mode == "T2":
            res = run_t2(cfg, ctx=ctx, prep=prep, load_labels=labels.test, out_dir=out_dir, author=author)
        else:
            res = run_p1_look(cfg, ctx=ctx, prep=prep, load_labels=labels.test, out_dir=out_dir, ledger_name="test",
                              author=author)
        prov.output(cp, *res["revealed"]["files"])  # sealed files are never hashed by reading them
        prov.note(look_id=res["look"], config_sha=res["config_sha"], revealed_sha=res["revealed"]["sha256"],
                  sealed_manifest_sha=res.get("sealed", {}).get("manifest_sha256", ""))
    _print({**res, "wall_s": round(time.time() - t0, 1)}, ctx)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
