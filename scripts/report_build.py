"""REPORT.md builder STUB (DESIGN.md section 11.3 "Reports and citations", section 6.4 wording rules).

Structure only. It will assemble every number in REPORT.md from committed files under results/ (never from
label files; labels are reached only through labels.py inside the pre-registered pipeline) and apply the
section 6.4 claim-language rules and the flags "gate not passed", "method-dominated", "label-calibrated" and
"exploratory". The functions below are pure and unit-tested (tests/test_w8_resource.py); ``collect_numbers`` and
``render`` are placeholders to be filled once T1/T2 outputs exist.

Usage (stub): .venv/bin/python scripts/report_build.py --check   # lists which sources exist
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

SESOI = 0.01  # C1, C4, C5 (absolute, in delta-bar_lin); C3 has none
FLAGS = ("gate not passed", "method-dominated", "label-calibrated", "exploratory")
MIN_LINEAGES_7 = math.ceil(0.7 * 7)  # 5


def claim_tier(contrast: str, *, passes_test: bool, point: float, ci_lower: float, sesoi: float | None = SESOI) -> str:
    """Section 6.4 claim language for one contrast (C1, C1b, C4, C5 each with its own SESOI; C3 has none).

    - "detectable improvement": the one-sided alpha = 0.025 test passes;
    - "improves": passes and point >= SESOI;
    - "meaningfully improves": passes and the CI lower bound >= SESOI.
    A contrast without SESOI (C3) only ever gets "detectable improvement".
    """
    if not passes_test:
        return "no detectable improvement"
    if contrast.upper() == "C3" or sesoi is None:
        return "detectable improvement"
    if ci_lower >= sesoi:
        return "meaningfully improves"
    if point >= sesoi:
        return "improves"
    return "detectable improvement"


def sesoi_c1b(mean_ap_random_all: float, mean_ap_random_minus_ce: float, base: float = SESOI) -> float:
    """SESOI_C1b = 0.01 x [mean E[AP_random] on U minus CE_fit] / [mean E[AP_random] on U] (section 6.4)."""
    if mean_ap_random_all <= 0:
        raise ValueError("mean E[AP_random] on U must be positive")
    return base * mean_ap_random_minus_ce / mean_ap_random_all


def transfer_wording(lineage_t_lower: float, n_positive_lineages: int, n_lineages: int = 7) -> str:
    """'Generalises to new tissues' only if the lineage t-interval lower bound > 0 AND the descriptive
    criterion (>= ceil(0.7 n) positive lineages, 5 of 7) holds; the claim is conditional on the test lineages."""
    need = math.ceil(0.7 * n_lineages)
    if lineage_t_lower > 0 and n_positive_lineages >= need:
        return "generalises to new tissues (conditional on the 7 test lineages)"
    return "no claim of transfer to unseen tissue"


def headroom_closed(num: float, den: float, den_ci: tuple[float, float]) -> float | None:
    """Delta(P1-G*)/Delta(B2-G*), reported only if the denominator's CI excludes 0; else None (report raw contrasts)."""
    lo, hi = den_ci
    if lo <= 0 <= hi:
        return None
    return num / den


def apply_flags(sentence: str, flags: tuple[str, ...] | list[str]) -> str:
    """Append the mandatory flag tags to a result sentence; unknown flags raise."""
    bad = [f for f in flags if f not in FLAGS]
    if bad:
        raise ValueError(f"unknown flag(s): {bad}")
    if not flags:
        return sentence
    order = [f for f in FLAGS if f in flags]
    return sentence.rstrip() + " [" + "; ".join(order) + "]"


def derive_flags(*, gate_passed: bool, method_dominated: bool = False, label_calibrated: bool = False,
                 exploratory: bool = False) -> list[str]:
    out = []
    if not gate_passed:
        out.append("gate not passed")
    if method_dominated:
        out.append("method-dominated")
    if label_calibrated:
        out.append("label-calibrated")
    if exploratory:
        out.append("exploratory")
    return out


def render_claim(contrast: str, label: str, point: float, ci: tuple[float, float], *, passes_test: bool,
                 sesoi: float | None = SESOI, flags: list[str] | tuple[str, ...] = ()) -> str:
    """One result sentence with wording rule and flags applied."""
    tier = claim_tier(contrast, passes_test=passes_test, point=point, ci_lower=ci[0], sesoi=sesoi)
    s = f"{contrast} ({label}): {tier}; delta = {point:.4f} [{ci[0]:.4f}, {ci[1]:.4f}]"
    return apply_flags(s, flags)


# --- placeholders (to be filled after T1/T2) -------------------------------------------------------------
NUMBER_SOURCES = {
    # key: (path relative to results/, how it will be read). Training/validation aggregates already public in README.md
    "w1_generic_auprc": ("w1/", "TODO: W1 summary tables"),
    "w1_binary_mcc": ("w1/", "TODO: section 6.8 binary endpoints"),
    "t1_c1_c1b_c3": ("t1/", "TODO: after T1 (sealed until then)"),
    "t2_c4_c5": ("t2/", "TODO: after T2"),
    "resource_manifest": ("../resource/manifest.tsv", "TODO: after W8 build"),
}


def collect_numbers(results_dir: Path) -> dict:
    """Stub: report which number sources exist; never reads label files."""
    return {k: (Path(results_dir) / rel).exists() for k, (rel, _) in NUMBER_SOURCES.items()}


def render(template: str, numbers: dict) -> str:
    raise NotImplementedError("REPORT.md assembly is filled after T1/T2 (see REPORT.md TODO lines)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.parse_args()
    root = Path(__file__).resolve().parents[1] / "results"
    for k, ok in collect_numbers(root).items():
        print(f"{k}\t{'present' if ok else 'missing'}\t{NUMBER_SOURCES[k][1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
