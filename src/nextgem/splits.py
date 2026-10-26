"""Split draw function (DESIGN.md §4.1, steps 2-5). Label-free.

Moved here from ``scripts/design_numbers.py`` (§4.1 step 1); both ``design_numbers.py`` and
``scripts/make_splits.py`` import it. The behaviour is pinned three ways, checked by make_splits.py and
tests/test_splits.py:

1. :func:`draw` reproduces the reference draws for seeds 1-20 and the acceptance statistics for seeds
   1-1000 recorded in ``results/design/design_numbers.json`` at the freeze commit;
2. it equals, seed for seed, the ``draw``/``checks`` functions extracted verbatim from
   ``scripts/design_numbers.py`` at the freeze commit (:func:`frozen_reference_functions`);
3. ``FROZEN.tsv`` holds this file's SHA-256.
"""

from __future__ import annotations

import ast
import subprocess
from pathlib import Path

import numpy as np

from nextgem import config

MIN_LINES_ELIGIBLE = 15
STRATA = {"S1": (0, 8), "S2": (8, 16), "S3": (16, 23)}  # rank slices of eligible lineages
ALLOC_TEST = {"S1": 2, "S2": 2, "S3": 3}  # PI decision 3: 7 test lineages
ALLOC_VAL = {"S1": 1, "S2": 1, "S3": 2}  # 4 validation lineages
HEME = {"Lymphoid", "Myeloid"}
CHECK_NAMES = ("a_min10_expr_each", "b_test_ge250_expr_lines", "c_test_ge40_proteomic_lines", "d_not_both_heme_heldout")


class SeedRefused(ValueError):
    """The seed (or one of its attempt seeds) is not allowed (DESIGN §4.1, V-2)."""


def eligible_and_strata(lineage_counts: dict[str, int]) -> tuple[dict[str, int], dict[str, int], dict[str, list[str]]]:
    """(eligible, rare, strata): eligible = lineages with >= 15 screened lines, ranked by size (desc),
    ties alphabetical; S1 = ranks 1-8, S2 = 9-16, S3 = 17-23."""
    eligible = {k: v for k, v in lineage_counts.items() if v >= MIN_LINES_ELIGIBLE}
    rare = {k: v for k, v in lineage_counts.items() if v < MIN_LINES_ELIGIBLE}
    ranked = sorted(eligible.items(), key=lambda kv: (-kv[1], kv[0]))
    if len(ranked) != STRATA["S3"][1]:
        raise ValueError(f"expected 23 eligible lineages, got {len(ranked)}")
    strata = {s: [lin for lin, _ in ranked[a:b]] for s, (a, b) in STRATA.items()}
    return eligible, rare, strata


def draw(seed: int, strata: dict[str, list[str]]) -> tuple[list[str], list[str]]:
    """One draw: ``default_rng(seed)``; permute S1, S2, S3 in that order; test = first 2/2/3, validation =
    next 1/1/2. Identical to ``draw`` in design_numbers.py at the freeze commit."""
    rng = np.random.default_rng(seed)
    test, val = [], []
    for s in ("S1", "S2", "S3"):
        perm = list(rng.permutation(strata[s]))
        test += perm[:ALLOC_TEST[s]]
        val += perm[ALLOC_TEST[s]:ALLOC_TEST[s] + ALLOC_VAL[s]]
    return test, val


def checks(test: list[str], val: list[str], expr_by_lin: dict[str, int], prot_by_lin: dict[str, int]) -> dict[str, bool]:
    """Label-blind acceptance checks (a)-(d) of §4.1 step 5."""
    ntest_expr = sum(expr_by_lin[L] for L in test)
    return {
        "a_min10_expr_each": all(expr_by_lin[L] >= 10 for L in test + val),
        "b_test_ge250_expr_lines": ntest_expr >= 250,
        "c_test_ge40_proteomic_lines": sum(prot_by_lin[L] for L in test) >= 40,
        "d_not_both_heme_heldout": not (HEME <= set(test + val)),
    }


def retired_overlap(seed: int) -> bool:
    """True if any attempt seed S + a (a = 0..99) falls in the retired range 20261026..20261125."""
    lo, hi = config.RETIRED_SEED_RANGE
    return seed + config.MAX_SPLIT_ATTEMPTS - 1 >= lo and seed <= hi


def check_seed(seed: int) -> None:
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, bool) or seed <= 0:
        raise SeedRefused(f"the split seed must be a positive integer (got {seed!r})")
    if retired_overlap(int(seed)):
        lo, hi = config.RETIRED_SEED_RANGE
        raise SeedRefused(
            f"seed {seed}: its attempt range {seed}..{seed + config.MAX_SPLIT_ATTEMPTS - 1} overlaps the retired "
            f"range {lo}..{hi} (DESIGN §4.1; seeds {lo - config.MAX_SPLIT_ATTEMPTS + 1}..{hi} are refused)"
        )


def run_attempts(seed: int, strata, expr_by_lin, prot_by_lin):
    """The draw loop: attempts 0..99 with ``default_rng(seed + attempt)`` until all checks pass.

    Returns (test, val, attempt_log rows). Raises if no attempt passes.
    """
    check_seed(seed)
    lo, hi = config.RETIRED_SEED_RANGE
    log = []
    for attempt in range(config.MAX_SPLIT_ATTEMPTS):
        s = int(seed) + attempt
        if lo <= s <= hi:  # repeated per attempt inside the loop (V-2)
            raise SeedRefused(f"attempt seed {s} is in the retired range")
        test, val = draw(s, strata)
        c = checks(test, val, expr_by_lin, prot_by_lin)
        ok = all(c.values())
        log.append({"attempt": attempt, "rng_seed": s, "accepted": ok, **c,
                    "failed_checks": ",".join(k for k, v in c.items() if not v),
                    "test_lineages": ";".join(test), "validation_lineages": ";".join(val)})
        if ok:
            return test, val, log
    raise RuntimeError(f"no acceptable draw in {config.MAX_SPLIT_ATTEMPTS} attempts")


def acceptance_stats(seeds, strata, expr_by_lin, prot_by_lin) -> dict:
    """The simulation summary of design_numbers.py (``split_simulation``), for the cross-check.
    Aggregates only; individual draws are not returned."""
    rows = []
    for sd in seeds:
        t, v = draw(sd, strata)
        c = checks(t, v, expr_by_lin, prot_by_lin)
        ntest = sum(expr_by_lin[L] for L in t)
        rows.append({"ok": all(c.values()), **c, "test_expr": ntest, "val_expr": sum(expr_by_lin[L] for L in v),
                     "largest_share": max(expr_by_lin[L] for L in t) / ntest, "test_prot": sum(prot_by_lin[L] for L in t)})
    ok = np.array([r["ok"] for r in rows])
    acc = [r for r in rows if r["ok"]]

    def q(key, p):
        return float(np.quantile([r[key] for r in acc], p))

    return {
        "acceptance_rate": round(float(ok.mean()), 3),
        "constraint_fail_rate": {k: round(float(np.mean([not r[k] for r in rows])), 3) for k in CHECK_NAMES},
        "accepted_test_expr_lines_median_p05_p95": [q("test_expr", .5), q("test_expr", .05), q("test_expr", .95)],
        "accepted_val_expr_lines_median_p05_p95": [q("val_expr", .5), q("val_expr", .05), q("val_expr", .95)],
        "accepted_largest_lineage_share_median_p95": [round(q("largest_share", .5), 3), round(q("largest_share", .95), 3)],
        "accepted_test_proteomic_lines_median_p05": [q("test_prot", .5), q("test_prot", .05)],
    }


def frozen_reference_functions(repo: Path, commit: str, strata, expr_by_lin, prot_by_lin):
    """``draw(seed)`` and ``checks(test, val)`` extracted verbatim (AST) from
    ``scripts/design_numbers.py`` at ``commit``, bound to the given inputs. Running the whole script is
    not needed: only the two function definitions are executed."""
    src = subprocess.run(["git", "-C", str(repo), "show", f"{commit}:scripts/design_numbers.py"],
                         check=True, capture_output=True, text=True).stdout
    tree = ast.parse(src)
    ns: dict = {"np": np, "strata": strata, "ALLOC_TEST": dict(ALLOC_TEST), "ALLOC_VAL": dict(ALLOC_VAL),
                "expr_by_lin": expr_by_lin, "prot_by_lin": prot_by_lin}
    # constants the functions may read, taken from the frozen script itself
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in ("ALLOC_TEST", "ALLOC_VAL") for t in node.targets):
            exec(compile(ast.Module(body=[node], type_ignores=[]), "<frozen design_numbers.py>", "exec"), ns)
    found = {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in ("draw", "checks"):
            exec(compile(ast.Module(body=[node], type_ignores=[]), "<frozen design_numbers.py>", "exec"), ns)
            found[node.name] = ns[node.name]
    if set(found) != {"draw", "checks"}:
        raise RuntimeError(f"draw/checks not found in scripts/design_numbers.py at {commit}")
    return found["draw"], found["checks"]
