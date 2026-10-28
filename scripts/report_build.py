"""REPORT.md builder and number checker (DESIGN.md section 11.3 "Reports and citations", section 6.4 wording rules).

Every number that REPORT.md quotes from a committed table is written behind an invisible marker::

    <!--n:m6.val.p1-->+0.0517 [0.0499, 0.0538] (results/m6/val/val_look_M6_factorial_decision.json)

``NUMBERS`` maps each marker key to its source file(s) and to a function that regenerates the exact string from that
file. ``--check`` regenerates every marked value and fails (exit 1) when REPORT.md and the file disagree, when a key is
unknown, or when the source file is not cited in the same paragraph, list item or table. ``--update`` rewrites the
marked values in place (a value written as ``?`` is a slot still to fill). Sources are read from the committed tree
(``git show HEAD:<path>``) by default, so uncommitted edits cannot leak into the report; ``--worktree`` reads the
working tree instead.

The script never reads label files: every source is a committed training or validation aggregate, a label-free table or
a provenance manifest. Test-tier numbers do not exist before T1; REPORT.md holds ``[T1: ...]`` / ``[T2: ...]``
placeholders for them, which ``--check`` counts and ``--require-filled`` turns into failures.

The section 6.4 helpers (``claim_tier``, ``sesoi_c1b``, ``transfer_wording``, ``headroom_closed``, flags) are unchanged
and unit-tested in tests/test_w8_resource.py; the marker machinery is tested in tests/test_report_build.py.

Usage:
    .venv/bin/python scripts/report_build.py --check            # verify REPORT.md against HEAD
    .venv/bin/python scripts/report_build.py --update           # rewrite marked values from HEAD
    .venv/bin/python scripts/report_build.py --sources          # list the T1/T2 sources and whether they exist
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import pandas as pd

SESOI = 0.01  # C1, C4, C5 (absolute, in delta-bar_lin); C3 has none
FLAGS = ("gate not passed", "method-dominated", "label-calibrated", "exploratory")
MIN_LINEAGES_7 = math.ceil(0.7 * 7)  # 5
ROOT = Path(__file__).resolve().parents[1]
MINUS = "−"


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


# --- T1/T2 sources (sealed until the looks run) ------------------------------------------------------------
NUMBER_SOURCES = {
    # key: (path relative to results/, how it will be read)
    "w1_generic_auprc": ("w1/", "W1 summary tables (now read through NUMBERS below)"),
    "w1_binary_mcc": ("w1/", "section 6.8 binary endpoints (now read through NUMBERS below)"),
    "t1_c1_c1b_c3": ("t1/", "TODO: after T1 (sealed until then)"),
    "t2_c4_c5": ("t2/", "TODO: after T2"),
    "resource_manifest": ("../resource/manifest.tsv", "read through NUMBERS below"),
}


def collect_numbers(results_dir: Path) -> dict:
    """Report which T1/T2 number sources exist; never reads label files."""
    return {k: (Path(results_dir) / rel).exists() for k, (rel, _) in NUMBER_SOURCES.items()}


def render(template: str, numbers: dict) -> str:
    raise NotImplementedError("REPORT.md test-tier assembly is filled after T1/T2 (see the [T1: ...] placeholders)")


# =========================================================================================================
# Source access (committed tree by default)
# =========================================================================================================
class Source:
    """Read committed files: ``git show <ref>:<path>`` when ``ref`` is set, else the working tree under ``root``."""

    def __init__(self, root: Path = ROOT, ref: str | None = "HEAD"):
        self.root = Path(root)
        self.ref = ref
        self._cache: dict = {}

    def text(self, rel: str) -> str:
        key = ("text", rel)
        if key not in self._cache:
            if self.ref:
                r = subprocess.run(["git", "show", f"{self.ref}:{rel}"], cwd=self.root, capture_output=True,
                                   text=True)
                if r.returncode != 0:
                    raise FileNotFoundError(f"{rel} is not in {self.ref}: {r.stderr.strip()}")
                self._cache[key] = r.stdout
            else:
                self._cache[key] = (self.root / rel).read_text()
        return self._cache[key]

    def json(self, rel: str):
        key = ("json", rel)
        if key not in self._cache:
            self._cache[key] = json.loads(self.text(rel).replace("NaN", "null"))
        return self._cache[key]

    def tsv(self, rel: str) -> pd.DataFrame:
        key = ("tsv", rel)
        if key not in self._cache:
            lines = [ln for ln in self.text(rel).splitlines() if not ln.startswith("#")]
            self._cache[key] = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t")
        return self._cache[key]

    def files(self, prefix: str) -> list[str]:
        """Committed (or, in worktree mode, existing) files under ``prefix``."""
        if self.ref:
            r = subprocess.run(["git", "ls-tree", "-r", "--name-only", self.ref, "--", prefix], cwd=self.root,
                               capture_output=True, text=True, check=True)
            return sorted(r.stdout.split())
        base = self.root / prefix
        return sorted(str(p.relative_to(self.root)) for p in base.rglob("*") if p.is_file()) if base.exists() else []

    def git(self, *args: str) -> str:
        return subprocess.run(["git", *args], cwd=self.root, capture_output=True, text=True, check=True).stdout


# --- getters -----------------------------------------------------------------------------------------------
def J(s: Source, path: str, *keys):
    v = s.json(path)
    for k in keys:
        v = v[k]
    return v


def T(s: Source, path: str, col: str, **where):
    """Exactly one row of a TSV matching ``where`` (keys may use ``__`` for characters not allowed in names)."""
    d = s.tsv(path)
    m = pd.Series(True, index=d.index)
    for k, v in where.items():
        m &= d[k] == v
    rows = d.loc[m, col]
    if len(rows) != 1:
        raise LookupError(f"{path}: {len(rows)} rows match {where} for column {col}")
    return rows.iloc[0]


def TW(s: Source, path: str, col: str, where: dict):
    return T(s, path, col, **where)


# --- formatters ----------------------------------------------------------------------------------------------
def fx(x: float, d: int = 3) -> str:
    """Fixed decimals with a typographic minus."""
    out = f"{float(x):.{d}f}"
    return out.replace("-", MINUS)


def sx(x: float, d: int = 4) -> str:
    """Signed fixed decimals (+0.0517, −0.0071)."""
    out = f"{float(x):+.{d}f}"
    return out.replace("-", MINUS)


def ci(lo: float, hi: float, d: int = 4) -> str:
    return f"[{fx(lo, d)}, {fx(hi, d)}]"


def sci(x: float, lo: float, hi: float, d: int = 4) -> str:
    return f"{sx(x, d)} {ci(lo, hi, d)}"


def ints(n) -> str:
    return f"{int(round(float(n))):,}"


def pct(x: float, d: int = 0) -> str:
    return f"{100 * float(x):.{d}f}%"


# =========================================================================================================
# Number registry: key -> (source paths, function(Source) -> exact string in REPORT.md)
# =========================================================================================================
SPLIT_META = "results/splits/split_meta.json"
SPLIT_SUM = "results/splits/split_summary.tsv"
DN = "results/design/design_numbers.json"
VLEDGER = "results/validation_ledger.tsv"
W1_SESOI = "results/w1/sesoi_c1b.json"
W1_GSTAR = "results/w1/val_look_W1_M2_baseline_gstar.json"
W1_VAL = "results/w1/val_look_W1_M2_baseline.tsv"
W1_GRID = "results/w1/train_grid.tsv"
W1_BIN = "results/w1/train_binary.tsv"
W1_ADDR = "results/w1/train_addressable.tsv"
W2_SL = "results/w2/train/shortlist.json"
W2_SHUF = "results/w2/train/shuffled_control.tsv"
W2_ATTR = "results/w2/train/attribution.json"
W2_STEP = "results/w2/train/step_scores.tsv"
W3_SL = "results/w3/train/shortlist.json"
W3_MODES = "results/w3/train/mode_scores.tsv"
W3_VAL = "results/w3/val/val_look_W3_M4_shortlist_decision.json"
W3_VALT = "results/w3/val/val_look_W3_M4_shortlist.tsv"
W4_SEL = "results/w4/stage_a_analysis/train_selection.json"
W4_H4 = "results/w4/stage_a_analysis/h4_tests.json"
W4_CORDA = "results/w4/stage_a_analysis/train_corda_a6prime.tsv"
W4_MAN = "results/w4/stage_a/MANIFEST.tsv"
W5_SUM = "results/w5/train_summary.json"
M5_DEC = "results/m5/train/train_decision.json"
M6_SUM = "results/m6/train/summary.json"
M6_CELLS = "results/m6/train/cells.tsv"
M6_INT = "results/m6/train/interaction.tsv"
M6_BIN = "results/m6/train/binary.tsv"
M6_BINP = "results/m6/train/binary_paired.tsv"
M6_DEC = "results/m6/val/val_look_M6_factorial_decision.json"
M6_VAL = "results/m6/val/val_look_M6_factorial.tsv"
W6_DEC = "results/w6/provisional/analysis/decisions.tsv"
W6_SEC = "results/w6/provisional/analysis/secondary.tsv"
W7_RES = "autoresearch/results.tsv"
W7_IDEAS = "autoresearch/ideas.tsv"
W7_GCV = "results/w7/screen_genecv.tsv"
W7_DEC = "results/w7/val/val_look_W7_P2_select_decision.json"
W7_VAL = "results/w7/val/val_look_W7_P2_select.tsv"
W8_SUM = "results/w8/p1/summary.json"
W8_MEMOTE = "results/w8/memote/memote_summary.json"
RES_MODELS = "resource/models.tsv"
T1_TSV = "results/m6/test/test_look_T1.tsv"
T1_GATES = "results/m6/test/test_look_T1_gates.json"
T2_GATES = "results/m6/test/test_look_T2_gates.json"
T2_GRID = "results/m6/test/test_look_T2_grid.tsv"
UN = "results/m6/test/unsealed/T1/"
UN_ATTR = UN + "attribution.tsv"
UN_LIN = UN + "per_lineage.tsv"
UN_GENE = UN + "per_gene.tsv"
UN_TASK = UN + "per_task.tsv"
UN_FAIL = UN + "failure_account.tsv"
UN_WORST = UN + "worst_cells.tsv"
UN_STRATA = UN + "strata.tsv"
UN_SENS = UN + "sensitivities.tsv"
UN_W1 = UN + "w1_benchmark.tsv"
UN_ALL = UN + "generic_all_lines.tsv"
H1_MAP = "results/baseline/human1_gene_mapping.tsv"
H2_MAP = "results/baseline/human2_gene_mapping.tsv"
W8_UNS = "results/w8/unscreened/build_summary.json"
W8_MEMOTE6 = "results/w8/memote/memote_summary_release6.json"
W8_FROG = "results/w8/frog/frog_summary.json"
W8_NOISE = "results/w8/ko_noise/human2_20cc933a4b852b6a_summary.json"
BUDGET = "results/budget.tsv"

PRIMARY = dict(label="L_prob", geneset="all", universe="U", stratum="none")


def _w1_train(s, method, geneset="all", col="tap_lin", stratum="none"):
    return TW(s, W1_GRID, col, dict(label="L_prob", geneset=geneset, universe="U", stratum=stratum,
                                    lineset="expr", method=method))


def _w1_val(s, method, geneset="all", col="tap_lin"):
    return TW(s, W1_VAL, col, dict(section="tAP_grid", label="L_prob", geneset=geneset, universe="U",
                                   stratum="none", method=method))


def _w1_bin(s, method, call, label="L_ge06", protocol="Human1-paper"):
    return TW(s, W1_BIN, "lineage_equal_mean", dict(protocol=protocol, method=method, call_rule=call, label=label,
                                                    gene_set="parent", lineset="expr", metric="mcc",
                                                    nan_policy="exclude"))


def _w1_valbin(s, method, call, label="L_ge06", protocol="Human1-paper"):
    return TW(s, W1_VAL, "lineage_equal_mean", dict(section="binary_6.8", protocol=protocol, method=method,
                                                    call_rule=call, label=label, gene_set="parent", metric="mcc",
                                                    nan_policy="exclude"))


def _val_tap(s, path, method, geneset="all", stratum="none", col="tap_lin"):
    return TW(s, path, col, dict(section="tAP_grid", label="L_prob", geneset=geneset, universe="U",
                                 stratum=stratum, method=method))


def _val_delta(s, path, method, geneset="all", stratum="none", d=4):
    w = dict(section="tAP_grid", label="L_prob", geneset=geneset, universe="U", stratum=stratum, method=method)
    return sci(TW(s, path, "delta_lin_vs_H1", w), TW(s, path, "delta_lo_vs_H1", w), TW(s, path, "delta_hi_vs_H1", w), d)


def _cand(s, path, method):
    for c in J(s, path, "candidates"):
        if c["method"] == method:
            return c
    raise LookupError(f"{path}: no candidate {method}")


def _split_list(s, tier, rare=None):
    d = s.tsv(SPLIT_SUM)
    d = d[d.tier == tier]
    if rare is not None:
        d = d[(d.stratum == "rare") == rare]
    d = d.sort_values(["n_lines", "lineage"], ascending=[False, True])
    return ", ".join(f"{r.lineage} ({ints(r.n_lines)})" for r in d.itertuples())


def _split_sum(s, tier, col, rare=None):
    d = s.tsv(SPLIT_SUM)
    d = d[d.tier == tier]
    if rare is not None:
        d = d[(d.stratum == "rare") == rare]
    return ints(d[col].sum())


def _w2_step(s, step, col):
    return TW(s, W2_STEP, col, dict(stratum="attribution", parent="human1", step=step))


def _w3_mode(s, parent, mode, col, medium="M0"):
    return TW(s, W3_MODES, col, dict(parent=parent, mode=mode, medium=medium))


def _w4_fam(s, fam, ref="vs_H1"):
    v = J(s, W4_SEL, "families", fam, ref)
    return v


def _w5_neg(s):
    vals = sorted(r["delta_lin"] for r in J(s, W5_SUM, "primary") if r["method"].startswith("H1_neg") and r["ref"] == "H1")
    return vals[0], vals[-1]


def _w5(s, method, ref):
    for r in J(s, W5_SUM, "primary"):
        if r["method"] == method and r["ref"] == ref:
            return r
    raise LookupError(f"{W5_SUM}: {method} vs {ref}")


def _w6(s, arm, hyp, analysis, col, contrast=None):
    d = s.tsv(W6_DEC)
    m = (d.arm == arm) & (d.hypothesis == hyp) & (d.analysis == analysis)
    if contrast:
        m &= d.contrast == contrast
    rows = d.loc[m, col]
    if len(rows) != 1:
        raise LookupError(f"{W6_DEC}: {len(rows)} rows for {arm}/{hyp}/{analysis}/{contrast}")
    return rows.iloc[0]


def _w7(s, idea, col):
    return TW(s, W7_RES, col, {"id": idea})


def _w7_num(text: str) -> float:
    return float(str(text).split()[0].replace(MINUS, "-"))


def _m6_bin(s, cell, protocol, call, label="L_ge06"):
    return TW(s, M6_BIN, "lineage_equal_mean", dict(metric="mcc", nan_policy="exclude", cell=cell,
                                                    protocol=protocol, call_rule=call, label=label))


def _m6_binp(s, cell, protocol, call, label="L_ge06"):
    w = dict(cell=cell, reference="H1", protocol=protocol, call_rule=call, label=label, metric="mcc")
    return sci(TW(s, M6_BINP, "delta_lin", w), TW(s, M6_BINP, "ci_lo", w), TW(s, M6_BINP, "ci_hi", w))


def _ledger_rows(s):
    return list(csv.DictReader(io.StringIO(s.text(VLEDGER)), delimiter="\t"))


def _manifest_rows(s: Source, rel: str) -> list[dict]:
    """Tolerant manifest reader: some early rows carry extra tab-separated output fields."""
    rows = list(csv.reader(io.StringIO(s.text(rel)), delimiter="\t"))
    head = rows[0]
    return [dict(zip(head, r)) for r in rows[1:] if r]


RUNTIME_GROUPS = ("w1", "w2", "w3", "w4", "w5", "w6", "w7", "w8", "m5", "m6")


def runtime_table(s: Source) -> dict[str, dict]:
    """Sum of recorded cpu_s and wall_s over committed results/**/{MANIFEST,PROVENANCE}.tsv, by top directory."""
    out: dict[str, dict] = {}
    for rel in s.files("results"):
        name = rel.rsplit("/", 1)[-1]
        if name not in ("MANIFEST.tsv", "PROVENANCE.tsv"):
            continue
        group = rel.split("/")[1]
        g = out.setdefault(group, {"cpu_s": 0.0, "wall_s": 0.0, "runs": 0, "files": []})
        for r in _manifest_rows(s, rel):
            try:
                g["cpu_s"] += float(r.get("cpu_s") or 0)
                g["wall_s"] += float(r.get("wall_s") or 0)
            except ValueError:
                continue
            g["runs"] += 1
        g["files"].append(rel)
    return out


def _rt(s, group, col="cpu_s"):
    t = runtime_table(s)
    if group == "total":
        return sum(g[col] for g in t.values()) / 3600
    if group == "other":
        return sum(g[col] for k, g in t.items() if k not in RUNTIME_GROUPS) / 3600
    return t[group][col] / 3600


def _rt_runs(s, group):
    t = runtime_table(s)
    if group == "total":
        return sum(g["runs"] for g in t.values())
    return t[group]["runs"]


def _elapsed_h(s):
    """Hours from the first commit to the last commit that touched results/ (stable while results are unchanged)."""
    last = s.git("log", "-1", "--format=%ct", s.ref or "HEAD", "--", "results").split()
    first = s.git("log", "--reverse", "--format=%ct", s.ref or "HEAD").split()
    return (int(last[0]) - int(first[0])) / 3600


def _pfmt(p: float) -> str:
    """Permutation p-values at the resolution floor 1/(1 + 10,000) are written as an inequality."""
    return "< 0.001" if p < 0.001 else fx(p, 3)


def _sci_e(x: float) -> str:
    return f"{x:.1e}".replace("-", MINUS)


def _gene_dn(s, method):
    d = s.tsv(UN_GENE)
    return d[(d.method == method) & (d.estimable.astype(bool))]["delta_null"]


def _frog_count(s, part):
    counts = J(s, W8_FROG, "models", 0, "status_counts")
    key = part if part in counts else part + "s"  # the committed summary uses singular keys; a later fix plural
    return counts[key]["optimal"]


def _worst_symbols(s):
    sym = {}
    for path in (H2_MAP, H1_MAP):  # Human1 symbols win where both parents map the gene
        m = s.tsv(path).dropna(subset=["entrez"])
        sym.update({str(int(float(e))): g for e, g in zip(m["entrez"], m["symbol"])})
    return ", ".join(sym.get(str(e), str(e)) for e in s.tsv(UN_WORST)["entrez"])


def _memote_scores(s, parent):
    vals = {round(m["total_score"], 4) for f in (W8_MEMOTE, W8_MEMOTE6) for m in J(s, f, "models") if parent in m["model"]}
    if len(vals) != 1:
        raise ValueError(f"MEMOTE scores differ across {parent} models: {sorted(vals)}")
    return vals.pop()


GIT_LOG = "git log"
RUNTIME_SRC = "results/*/**/MANIFEST.tsv and PROVENANCE.tsv"


@dataclass(frozen=True)
class Number:
    sources: tuple[str, ...]
    fn: Callable[[Source], str]


def _n(sources, fn) -> Number:
    return Number(tuple([sources] if isinstance(sources, str) else sources), fn)


def _registry() -> dict[str, Number]:
    R: dict[str, Number] = {}

    # --- cohort, labels, design constants from data -------------------------------------------------------
    R["split.n.total"] = _n(SPLIT_META, lambda s: ints(sum(J(s, SPLIT_META, "lines_per_tier").values())))
    for tier, k in (("test", "test"), ("validation", "val"), ("train", "train")):
        R[f"split.n.{k}"] = _n(SPLIT_META, lambda s, t=tier: ints(J(s, SPLIT_META, "lines_per_tier", t)))
        R[f"split.list.{k}"] = _n(SPLIT_SUM, lambda s, t=tier: _split_list(s, t))
        R[f"split.expr.{k}"] = _n(SPLIT_SUM, lambda s, t=tier: _split_sum(s, t, "n_expression"))
        R[f"split.prot.{k}"] = _n(SPLIT_SUM, lambda s, t=tier: _split_sum(s, t, "n_proteomics_any"))
    R["split.list.train"] = _n(SPLIT_SUM, lambda s: _split_list(s, "train", rare=False))
    R["split.rare.n"] = _n(SPLIT_SUM, lambda s: _split_sum(s, "train", "n_lines", rare=True))
    R["split.seed"] = _n(SPLIT_META, lambda s: str(J(s, SPLIT_META, "seed")))
    R["genes.U"] = _n(DN, lambda s: ints(J(s, DN, "genes", "U_union")))
    R["genes.Ucap"] = _n(DN, lambda s: ints(J(s, DN, "genes", "U_intersection")))
    R["labels.kappa"] = _n(VLEDGER, lambda s: fx(json.loads(_ledger_rows(s)[0]["detail"])["kappa"], 3))
    R["looks.used"] = _n(VLEDGER, lambda s: str(sum(r["event"] == "completed" for r in _ledger_rows(s))))
    R["sesoi.c1b"] = _n(W1_SESOI, lambda s: fx(J(s, W1_SESOI, "SESOI_C1b"), 4))
    R["cefit.n"] = _n(W1_SESOI, lambda s: ints(J(s, W1_SESOI, "CE_fit_size")))

    # --- W1 generic models and baselines --------------------------------------------------------------------
    for k, m in (("h1", "H1"), ("h2", "H2"), ("b0", "B0"), ("b2", "B2[L_prob]"), ("b2m", "B2m[L_prob]"),
                 ("b3", "B3"), ("b6", "B6[L_prob]"), ("b8p", "B8+_H1"), ("b8m", "B8-_H1"), ("b8q", "B8'_H1"),
                 ("h1task", "H1_task")):
        R[f"w1.tr.{k}"] = _n(W1_GRID, lambda s, m=m: fx(_w1_train(s, m)))
        R[f"w1.tr.{k}.noce"] = _n(W1_GRID, lambda s, m=m: fx(_w1_train(s, m, "noCE")))
    for k, m in (("h1", "H1"), ("h2", "H2"), ("b0", "B0"), ("b2", "B2[L_prob]"), ("b2m", "B2m[L_prob]"),
                 ("b3", "B3"), ("b8p", "B8+_H1"), ("b8m", "B8-_H1"), ("b8q", "B8'_H1"), ("h1task", "H1_task")):
        R[f"w1.val.{k}"] = _n(W1_VAL, lambda s, m=m: fx(_w1_val(s, m)))
        R[f"w1.val.{k}.noce"] = _n(W1_VAL, lambda s, m=m: fx(_w1_val(s, m, "noCE")))
    R["w1.tr.h2h1"] = _n(W1_GRID, lambda s: sci(_w1_train(s, "H2", col="delta_lin_vs_H1"),
                                                _w1_train(s, "H2", col="delta_lo_vs_H1"),
                                                _w1_train(s, "H2", col="delta_hi_vs_H1")))
    R["w1.val.h2h1"] = _n(W1_GSTAR, lambda s: sci(J(s, W1_GSTAR, "delta_lin_H2_vs_H1"),
                                                  *J(s, W1_GSTAR, "delta_ci")))
    R["w1.val.gstar"] = _n(W1_GSTAR, lambda s: J(s, W1_GSTAR, "G_star"))
    for k, m, c in (("h1.p57", "H1", "P57"), ("h2.p57", "H2", "P57"), ("h1.pgr", "H1", "PGR"),
                    ("h2.pgr", "H2", "PGR")):
        R[f"w1.tr.mcc.{k}"] = _n(W1_BIN, lambda s, m=m, c=c: fx(_w1_bin(s, m, c)))
        R[f"w1.val.mcc.{k}"] = _n(W1_VAL, lambda s, m=m, c=c: fx(_w1_valbin(s, m, c)))
    R["w1.addr.h2.out"] = _n(W1_ADDR, lambda s: pct(T(s, W1_ADDR, "frac_outside", label="L_prob", parent="human2")))
    R["w1.addr.h1.out"] = _n(W1_ADDR, lambda s: pct(T(s, W1_ADDR, "frac_outside", label="L_prob", parent="human1")))
    R["w1.addr.h2.n"] = _n(W1_ADDR, lambda s: ints(T(s, W1_ADDR, "n_addressable_U", label="L_prob",
                                                       parent="human2")))

    # --- W2 medium ----------------------------------------------------------------------------------------
    R["w2.m2"] = _n(W2_SL, lambda s: sci(*(J(s, W2_SL, "candidate_row", k)
                                           for k in ("delta_lin_vs_M0", "delta_lo_vs_M0", "delta_hi_vs_M0"))))
    R["w2.m2.folds"] = _n(W2_SL, lambda s: f"{int(J(s, W2_SL, 'candidate_row', 'folds_positive'))}/12")
    R["w2.m2.noce"] = _n(W2_SL, lambda s: sx(J(s, W2_SL, "candidate_row", "delta_noCE_vs_M0")))
    R["w2.nattr"] = _n(W2_SL, lambda s: ints(J(s, W2_SL, "candidate_row", "n_lines")))
    R["w2.shuf"] = _n(W2_SHUF, lambda s: sci(*(T(s, W2_SHUF, c, parent="human1", step="M2(f=1)", scope="pooled")
                                               for c in ("doc_minus_shuffled", "lo", "hi"))))
    R["w2.shuf.share"] = _n(W2_SHUF, lambda s: pct(
        T(s, W2_SHUF, "shuffled_minus_m0_mean", parent="human1", step="M2(f=1)", scope="pooled")
        / T(s, W2_SHUF, "doc_minus_m0", parent="human1", step="M2(f=1)", scope="pooled")))
    R["w2.mopen"] = _n(W2_STEP, lambda s: sx(_w2_step(s, "M-open", "delta_lin_vs_M0")))
    R["w2.mopen.noce"] = _n(W2_STEP, lambda s: sx(_w2_step(s, "M-open", "delta_noCE_vs_M0")))
    R["w2.m1"] = _n(W2_STEP, lambda s: sx(_w2_step(s, "M1", "delta_lin_vs_M0")))
    R["w2.hplm"] = _n(W2_STEP, lambda s: sx(_w2_step(s, "M-HPLM:M2(f=1)", "delta_lin_vs_M0")))
    R["w2.plasmax"] = _n(W2_STEP, lambda s: sx(_w2_step(s, "M-Plasmax:M2(f=1)", "delta_lin_vs_M0")))
    R["w2.oracle"] = _n(W2_ATTR, lambda s: fx(J(s, W2_ATTR, "human1", "oracle", "ratio"), 3))
    R["w2.netflip"] = _n(W2_ATTR, lambda s: pct(J(s, W2_ATTR, "human1", "net_flip", "L_prob",
                                                  "net_flip_fraction"), 1))

    # --- W3 biomass, tasks, objectives -----------------------------------------------------------------------
    R["w3.tr.bloo"] = _n(W3_SL, lambda s: sci(J(s, W3_SL, "candidates", 0, "delta_lin"),
                                              *J(s, W3_SL, "candidates", 0, "ci")))
    R["w3.tr.bloo.folds"] = _n(W3_SL, lambda s: f"{J(s, W3_SL, 'candidates', 0, 'folds_positive')}/12")
    R["w3.tr.bloo.h2"] = _n(W3_MODES, lambda s: sx(_w3_mode(s, "human2", "B-loo:cofactor_pool_biomass",
                                                            "delta_lin")))
    R["w3.tr.bloo.noce"] = _n(W3_MODES, lambda s: sx(_w3_mode(s, "human1", "B-loo:cofactor_pool_biomass",
                                                              "delta_noCE")))
    R["w3.tr.sym"] = _n(W3_SL, lambda s: sx(J(s, W3_SL, "candidates", 0, "vs_mcr_symmetric"), 3))
    R["w3.tr.metpool"] = _n(W3_SL, lambda s: sx(J(s, W3_SL, "withdrawn", 0, "delta_lin")))
    R["w3.wc"] = _n(W3_SL, lambda s: fx(J(s, W3_SL, "winners_curse_sigma_sqrt_2lnm"), 3))
    for k, mode in (("t3", "T3"), ("t5", "T5"), ("moma", "O-MOMA"), ("ocap", "O-cap"), ("dbs", "D-BS")):
        R[f"w3.tr.{k}"] = _n(W3_MODES, lambda s, m=mode: sx(_w3_mode(s, "human1", m, "delta_lin")))
        R[f"w3.tr.{k}.noce"] = _n(W3_MODES, lambda s, m=mode: sx(_w3_mode(s, "human1", m, "delta_noCE")))
    R["w3.val"] = _n(W3_VAL, lambda s: sci(J(s, W3_VAL, "candidates", 0, "delta_lin"),
                                           *J(s, W3_VAL, "candidates", 0, "delta_ci")))
    R["w3.val.tap"] = _n(W3_VALT, lambda s: fx(_val_tap(s, W3_VALT, "W3:B-loo:cofactor_pool_biomass[H1]")))
    R["w3.val.noce"] = _n(W3_VALT, lambda s: _val_delta(s, W3_VALT, "W3:B-loo:cofactor_pool_biomass[H1]", "noCE"))

    # --- W4 Stage A ---------------------------------------------------------------------------------------------
    for fam in ("F0", "F2", "F3", "F4", "F5"):
        R[f"w4.{fam.lower()}"] = _n(W4_SEL, lambda s, f=fam: sci(_w4_fam(s, f)["outer_delta_lin"],
                                                                  _w4_fam(s, f)["lo"], _w4_fam(s, f)["hi"]))
        R[f"w4.{fam.lower()}.folds"] = _n(W4_SEL, lambda s, f=fam: f"{_w4_fam(s, f)['folds_positive']}/10")
        R[f"w4.{fam.lower()}.h2"] = _n(W4_SEL, lambda s, f=fam: sx(_w4_fam(s, f, "vs_H2")["outer_delta_lin"]))
    R["w4.corda"] = _n(W4_CORDA, lambda s: sx(T(s, W4_CORDA, "delta_vs_H1", rule="R1", family="F3b")))
    for lvl, k in (("reactions", "rxn"), ("tasks57", "task"), ("ko_binary", "ko")):
        R[f"w4.h41.{k}.method"] = _n(W4_H4, lambda s, lv=lvl: fx(J(s, W4_H4, "levels", lv, "H4.1", "r2_method")))
        R[f"w4.h41.{k}.lineage"] = _n(W4_H4, lambda s, lv=lvl: fx(J(s, W4_H4, "levels", lv, "H4.1", "r2_lineage")))
    R["w4.nbuilds"] = _n(W4_MAN, lambda s: ints(json.loads(_manifest_rows(s, W4_MAN)[0]["extra"])["n_builds"]))
    R["w4.buildcpu"] = _n(W4_MAN, lambda s: fx(json.loads(_manifest_rows(s, W4_MAN)[0]["extra"])["build_cpu_h"], 1))

    # --- W5 GPR curation ----------------------------------------------------------------------------------------
    for k, m, ref in (("h1", "H1_W5a", "H1"), ("h2", "H2_W5a", "H2"), ("pid30.h1", "H1_W5a_pid30", "H1"),
                      ("pid30.h2", "H2_W5a_pid30", "H2"), ("tyms", "H1_TYMS", "H1")):
        R[f"w5.{k}"] = _n(W5_SUM, lambda s, m=m, r=ref: sci(_w5(s, m, r)["delta_lin"], _w5(s, m, r)["delta_lo"],
                                                             _w5(s, m, r)["delta_hi"]))
    R["w5.tyms.folds"] = _n(W5_SUM, lambda s: f"{_w5(s, 'H1_TYMS', 'H1')['folds_positive']}/12")
    R["w5.neg.h1"] = _n(W5_SUM, lambda s: "{} to {}".format(*(sx(v) for v in _w5_neg(s))))

    # --- M5 context families on top of B-loo ------------------------------------------------------------------
    for fam in ("F2", "F4"):
        R[f"m5.{fam.lower()}"] = _n(M5_DEC, lambda s, f=fam: sci(J(s, M5_DEC, "families", f, "delta_lin_vs_H1_Bloo"),
                                                                  *J(s, M5_DEC, "families", f, "ci")))
        R[f"m5.{fam.lower()}.folds"] = _n(M5_DEC, lambda s, f=fam: f"{J(s, M5_DEC, 'families', f, 'folds_positive')}/10")

    # --- M6 factorial (training) ---------------------------------------------------------------------------------
    for k, cell in (("c1", "H1+Bloo/M0"), ("c2", "H1/M2"), ("c3", "H1+Bloo/M2")):
        R[f"m6.tr.{k}"] = _n(M6_SUM, lambda s, c=cell: sci(J(s, M6_SUM, "tg", c, "delta_lin"),
                                                            *J(s, M6_SUM, "tg", c, "ci")))
    for k, cell in (("gstar", "H1"), ("c1", "H1+Bloo/M0"), ("c2", "H1/M2"), ("c3", "H1+Bloo/M2"),
                    ("h2c3", "H2+Bloo/M2")):
        R[f"m6.tr.tap.{k}"] = _n(M6_CELLS, lambda s, c=cell: fx(T(s, M6_CELLS, "tap_lin", cell=c)))
    R["m6.tr.inter"] = _n(M6_INT, lambda s: sci(*(T(s, M6_INT, c, parent="H1", label_cell="all")
                                                  for c in ("interaction_lin", "lo", "hi"))))
    R["m6.tr.additive"] = _n(M6_INT, lambda s: sx(T(s, M6_INT, "additive_prediction", parent="H1", label_cell="all")))
    R["m6.tr.mcc.gstar.fix"] = _n(M6_BIN, lambda s: fx(_m6_bin(s, "H1", "fixed-cutoff", "FBA<0.01")))
    R["m6.tr.mcc.p1.fix"] = _n(M6_BIN, lambda s: fx(_m6_bin(s, "H1+Bloo/M2", "fixed-cutoff", "FBA<0.01")))
    R["m6.tr.mcc.p1.p57"] = _n(M6_BIN, lambda s: fx(_m6_bin(s, "H1+Bloo/M2", "Human1/2-paper", "P57")))
    R["m6.tr.mcc.d.fix"] = _n(M6_BINP, lambda s: _m6_binp(s, "H1+Bloo/M2", "fixed-cutoff", "FBA<0.01"))
    R["m6.tr.mcc.d.p57"] = _n(M6_BINP, lambda s: _m6_binp(s, "H1+Bloo/M2", "Human1/2-paper", "P57"))

    # --- M6 factorial (validation look) -----------------------------------------------------------------------
    for k, m in (("c1", "M6:H1+Bloo/M0"), ("c2", "M6:H1/M2(f=1)"), ("c3", "M6:H1+Bloo/M2(f=1)")):
        R[f"m6.val.{k}"] = _n(M6_DEC, lambda s, m=m: sci(_cand(s, M6_DEC, m)["delta_lin"],
                                                          *_cand(s, M6_DEC, m)["delta_ci"]))
    R["m6.val.p1.short"] = _n(M6_DEC, lambda s: sci(_cand(s, M6_DEC, "M6:H1+Bloo/M2(f=1)")["delta_lin"],
                                                    *_cand(s, M6_DEC, "M6:H1+Bloo/M2(f=1)")["delta_ci"], 3))
    R["m6.val.p1.pt3"] = _n(M6_DEC, lambda s: sx(_cand(s, M6_DEC, "M6:H1+Bloo/M2(f=1)")["delta_lin"], 3))
    R["m6.val.selected"] = _n(M6_DEC, lambda s: J(s, M6_DEC, "selected"))
    R["m6.val.wc"] = _n(M6_DEC, lambda s: fx(J(s, M6_DEC, "winners_curse", "term"), 4))
    R["m6.val.adj"] = _n(M6_DEC, lambda s: sx(J(s, M6_DEC, "winners_curse", "best_delta_lin_minus_term")))
    for k, m in (("gstar", "H1"), ("h2", "H2"), ("c1", "M6:H1+Bloo/M0"), ("c2", "M6:H1/M2(f=1)"),
                 ("c3", "M6:H1+Bloo/M2(f=1)"), ("h2c3", "M6:H2+Bloo/M2(f=1)")):
        R[f"m6.val.tap.{k}"] = _n(M6_VAL, lambda s, m=m: fx(_val_tap(s, M6_VAL, m)))
    R["m6.val.p1.noce"] = _n(M6_VAL, lambda s: _val_delta(s, M6_VAL, "M6:H1+Bloo/M2(f=1)", "noCE"))
    R["m6.val.p1.expr"] = _n(M6_VAL, lambda s: _val_delta(s, M6_VAL, "M6:H1+Bloo/M2(f=1)", stratum="expressed_pairs"))
    R["m6.val.h2c3"] = _n(M6_VAL, lambda s: _val_delta(s, M6_VAL, "M6:H2+Bloo/M2(f=1)"))
    R["m6.val.p1.mdap"] = _n(M6_VAL, lambda s: sx(_val_tap(s, M6_VAL, "M6:H1+Bloo/M2(f=1)", col="m_delta_ap_gene")))

    # --- W6 flux information ----------------------------------------------------------------------------------
    for arm, ka in (("H1", "gstar"), ("H1+B-loo", "p1")):
        for an in ("subtypes", "drugs", "metabolomics"):
            R[f"w6.{ka}.flux.{an}.med"] = _n(W6_DEC, lambda s, a=arm, x=an: fx(_w6(s, a, "H_flux", x, "median_delta"), 3))
            R[f"w6.{ka}.flux.{an}.p"] = _n(W6_DEC, lambda s, a=arm, x=an: fx(_w6(s, a, "H_flux", x, "p_holm"), 2))
            R[f"w6.{ka}.sim.{an}.med"] = _n(W6_DEC, lambda s, a=arm, x=an: sx(
                _w6(s, a, "H_sim", x, "median_delta", contrast="FLUX - RXN_RNA"), 3))
        R[f"w6.{ka}.synth"] = _n(W6_SEC, lambda s, a=arm: sx(T(s, W6_SEC, "median_delta", arm=a,
                                                               analysis="synthetic_positive"), 2))
    R["w6.drugs.zero"] = _n(W6_DEC, lambda s: f"{int(_w6(s, 'H1', 'H_flux', 'drugs', 'n_zero'))}/"
                                              f"{int(_w6(s, 'H1', 'H_flux', 'drugs', 'n_units'))}")
    for an in ("subtypes", "drugs", "metabolomics"):
        R[f"w6.n.{an}"] = _n(W6_DEC, lambda s, x=an: ints(_w6(s, "H1", "H_flux", x, "n_units")))

    # --- W7 auto-research ---------------------------------------------------------------------------------------
    R["w7.tr.b6"] = _n(W7_RES, lambda s: fx(_w7(s, "W7-01", "tap_lin"), 3))
    R["w7.tr.05"] = _n(W7_RES, lambda s: fx(_w7(s, "W7-05", "tap_lin"), 3))
    R["w7.tr.05.vsb6"] = _n(W7_RES, lambda s: sci(_w7_num(_w7(s, "W7-05", "d_vs_comparator")),
                                                  _w7_num(_w7(s, "W7-05", "ci_lo")), _w7_num(_w7(s, "W7-05", "ci_hi"))))
    R["w7.tr.05.vsb2"] = _n(W7_RES, lambda s: sx(_w7_num(_w7(s, "W7-05", "d_vs_B2"))))
    for k, idea in (("02", "W7-02"), ("04", "W7-04"), ("12", "W7-12"), ("14", "W7-14"), ("08", "W7-08"),
                    ("15", "W7-15"), ("17", "W7-17")):
        R[f"w7.tr.{k}.d"] = _n(W7_RES, lambda s, i=idea: sx(_w7_num(_w7(s, i, "d_vs_comparator"))))
    R["w7.tr.03.tap"] = _n(W7_RES, lambda s: fx(_w7(s, "W7-03", "tap_lin"), 3))
    R["w7.tr.03.vsp1"] = _n(W7_RES, lambda s: sx(_w7_num(_w7(s, "W7-03", "d_vs_P1"))))
    R["w7.n.logged"] = _n(W7_IDEAS, lambda s: ints(len(s.tsv(W7_IDEAS))))
    R["w7.n.scored"] = _n(W7_RES, lambda s: ints(len(s.tsv(W7_RES))))
    R["w7.n.promoted"] = _n(W7_RES, lambda s: ints(s.tsv(W7_RES)["promoted"].astype(str).str.startswith("YES").sum()))
    R["w7.gcv.05"] = _n(W7_GCV, lambda s: fx(T(s, W7_GCV, "tap_lin", method="W7-05|genecv"), 3))
    R["w7.gcv.b6"] = _n(W7_GCV, lambda s: fx(T(s, W7_GCV, "tap_lin", method="W7-01|genecv"), 3))
    R["w7.val.b6"] = _n(W7_VAL, lambda s: fx(_val_tap(s, W7_VAL, "W7:B6")))
    R["w7.val.05"] = _n(W7_VAL, lambda s: fx(_val_tap(s, W7_VAL, "W7:W7-05")))
    R["w7.val.05.vsg"] = _n(W7_DEC, lambda s: sci(_cand(s, W7_DEC, "W7:W7-05")["delta_lin"],
                                                  *_cand(s, W7_DEC, "W7:W7-05")["delta_ci"], 3))
    R["w7.val.diff"] = _n(W7_DEC, lambda s: sx(_cand(s, W7_DEC, "W7:B6")["diff_vs_best"]))
    R["w7.val.se"] = _n(W7_DEC, lambda s: fx(_cand(s, W7_DEC, "W7:B6")["se_diff"], 4))
    R["w7.val.selected"] = _n(W7_DEC, lambda s: J(s, W7_DEC, "selected"))

    # --- W8 resource -------------------------------------------------------------------------------------------
    R["w8.models"] = _n(W8_SUM, lambda s: ints(J(s, W8_SUM, "n_models")))
    R["w8.pass"] = _n(W8_SUM, lambda s: ints(J(s, W8_SUM, "n_models_pass")))
    R["w8.lines"] = _n(W8_SUM, lambda s: ints(J(s, W8_SUM, "n_lines")))
    R["w8.rows"] = _n(W8_SUM, lambda s: ints(J(s, W8_SUM, "n_line_rows")))
    R["w8.mb"] = _n(W8_SUM, lambda s: ints(J(s, W8_SUM, "model_bytes") / 1e6))
    R["w8.cpu"] = _n(W8_SUM, lambda s: fx(J(s, W8_SUM, "job_cpu_s") / 3600, 1))
    R["w8.wall"] = _n(W8_SUM, lambda s: ints(J(s, W8_SUM, "wall_s") / 60))
    R["w8.h2pass"] = _n(W8_SUM, lambda s: f"{J(s, W8_SUM, 'per_parent', 'human2', 'pass')}/"
                                          f"{J(s, W8_SUM, 'per_parent', 'human2', 'models')}")
    R["w8.h1pass"] = _n(W8_SUM, lambda s: f"{J(s, W8_SUM, 'per_parent', 'human1', 'pass')}/"
                                          f"{J(s, W8_SUM, 'per_parent', 'human1', 'models')}")
    R["w8.ko.pass"] = _n(W8_SUM, lambda s: ints(J(s, W8_SUM, "gate_pass", "gate_random_ko")))
    R["w8.signatures"] = _n(RES_MODELS, lambda s: ints((s.tsv(RES_MODELS)["parent"] == "human2").sum()))
    R["w8.mev.models"] = _n(RES_MODELS, lambda s: ints((~s.tsv(RES_MODELS)["mevalonate_core_lethal"].astype(bool)).sum()))
    R["w8.memote.h2"] = _n(W8_MEMOTE, lambda s: fx(next(m["total_score"] for m in J(s, W8_MEMOTE, "models")
                                                        if "human2" in m["model"]), 2))
    R["w8.memote.h1"] = _n(W8_MEMOTE, lambda s: fx(next(m["total_score"] for m in J(s, W8_MEMOTE, "models")
                                                        if "human1" in m["model"]), 2))

    # --- T1 / T2 / unseal (test tier; committed after the looks) ------------------------------------------------
    def t1ci(s, name, d=4):
        d_ = s.tsv(T1_TSV)
        row = d_[(d_.name == name) & d_.section.isin(["confirmatory", "estimates"])]
        if len(row) != 1:
            raise LookupError(f"{T1_TSV}: {len(row)} rows for {name}")
        r = row.iloc[0]
        return sci(r.estimate, r.ci_lo, r.ci_hi, d)

    for c in ("C1", "C1b"):
        R[f"t1.{c.lower()}"] = _n(T1_GATES, lambda s, c=c: sci(J(s, T1_GATES, "contrasts", c, "estimate"),
                                                               J(s, T1_GATES, "contrasts", c, "ci_low"),
                                                               J(s, T1_GATES, "contrasts", c, "ci_high")))
        R[f"t1.{c.lower()}.claim"] = _n(T1_GATES, lambda s, c=c: next(g["claim"] for g in J(s, T1_GATES, "gates")
                                                                      if g["contrast"] == c))
        R[f"t1.{c.lower()}.p"] = _n(T1_GATES, lambda s, c=c: _pfmt(J(s, T1_GATES, "contrasts", c, "p")))
    R["t1.c1.short"] = _n(T1_GATES, lambda s: sci(J(s, T1_GATES, "contrasts", "C1", "estimate"),
                                                  J(s, T1_GATES, "contrasts", "C1", "ci_low"),
                                                  J(s, T1_GATES, "contrasts", "C1", "ci_high"), 3))
    R["t1.c1.pt3"] = _n(T1_GATES, lambda s: sx(J(s, T1_GATES, "contrasts", "C1", "estimate"), 3))
    R["t1.c1b.status"] = _n(T1_GATES, lambda s: next(g["status"] for g in J(s, T1_GATES, "gates")
                                                     if g["contrast"] == "C1b"))
    R["t1.c3"] = _n(T1_GATES, lambda s: sci(J(s, T1_GATES, "contrasts", "C3", "estimate"),
                                            J(s, T1_GATES, "contrasts", "C3", "ci_low"),
                                            J(s, T1_GATES, "contrasts", "C3", "ci_high")))
    R["t1.c3.p"] = _n(T1_GATES, lambda s: fx(J(s, T1_GATES, "contrasts", "C3", "p"), 2))
    R["t1.c3.claim"] = _n(T1_GATES, lambda s: next(g["claim"] for g in J(s, T1_GATES, "gates") if g["contrast"] == "C3"))
    R["t1.c3.mdap"] = _n(T1_TSV, lambda s: sx(T(s, T1_TSV, "m_delta_ap_gene_p1", name="C3")))
    R["t1.c3.null"] = _n(T1_TSV, lambda s: fx(T(s, T1_TSV, "null_threshold", name="C3"), 4))
    R["t1.c3.ngenes"] = _n(T1_TSV, lambda s: ints(T(s, T1_TSV, "n_genes_eligible", name="C3")))
    R["t1.c3.expr"] = _n(T1_TSV, lambda s: sci(*(T(s, T1_TSV, c, name="C3_expressed_pairs")
                                                 for c in ("estimate", "ci_lo", "ci_hi"))))
    tr = ("transfer_C1",)
    R["t1.tr.t"] = _n(T1_GATES, lambda s: ci(J(s, T1_GATES, *tr, "t_low"), J(s, T1_GATES, *tr, "t_high")))
    R["t1.tr.p"] = _n(T1_GATES, lambda s: fx(J(s, T1_GATES, *tr, "exact_sign_flip_p"), 4))
    R["t1.tr.k"] = _n(T1_GATES, lambda s: f"{J(s, T1_GATES, *tr, 'k_positive')}/{J(s, T1_GATES, *tr, 'n_lineages')}")
    R["t1.n.lines"] = _n(T1_TSV, lambda s: ints(T(s, T1_TSV, "n_lines", section="counts")))
    R["t1.n.noexpr"] = _n(T1_TSV, lambda s: ints(T(s, T1_TSV, "n_lines_no_expression", section="counts")))
    for k, name in (("b2", "C2:B2"), ("b2m", "C2:B2m"), ("b2c", "C2:B2c"), ("b8p", "B8+"), ("b8m", "B8-"),
                    ("b8q", "B8'"), ("wp.h2", "within_parent:H2"), ("p1h2", "vs_G*:P1[H2]")):
        R[f"t1.est.{k}"] = _n(T1_TSV, lambda s, n=name: t1ci(s, n))
    R["t1.headroom"] = _n(T1_TSV, lambda s: "{} {}".format(
        pct(T(s, T1_TSV, "estimate", name="headroom_closed"), 1),
        "[{}, {}]".format(pct(T(s, T1_TSV, "ci_lo", name="headroom_closed"), 1),
                          pct(T(s, T1_TSV, "ci_hi", name="headroom_closed"), 1))))
    R["t1.headroom.den"] = _n(T1_TSV, lambda s: sci(*(T(s, T1_TSV, c, name="headroom_closed")
                                                      for c in ("denominator", "den_ci_lo", "den_ci_hi"))))
    for k, m in (("p1", "P1[H1]"), ("gstar", "H1"), ("p1h2", "P1[H2]")):
        R[f"t1.tap.{k}"] = _n(T1_TSV, lambda s, m=m: fx(T(s, T1_TSV, "tap_lin", name=f"tap:{m}")))
    R["t1.failed"] = _n(T1_TSV, lambda s: ints(s.tsv(T1_TSV).query("section == 'failed_builds' and lineage == 'all'")
                                               ["n_failed"].sum()))
    R["t1.adds"] = _n(T1_GATES, lambda s: "not claimed" if not J(s, T1_GATES, "adds_beyond_generic_expression", "claim")
                      else "claimed")
    R["t1.generalises"] = _n(T1_GATES, lambda s: "holds" if J(s, T1_GATES, "generalises_to_new_tissues") else "fails")
    for proto, call, lab, kk in (("Human1-paper", "P57", "L_ge06", "p57"), ("Human1-paper", "PGR", "L_ge06", "pgr"),
                                 ("fixed-cutoff", "FBA<0.01", "L_ge06", "fix06"),
                                 ("fixed-cutoff", "FBA<0.01", "L_prob", "fixp")):
        for mk, m in (("h1", "H1"), ("h2", "H2"), ("p1", "P1[H1]"), ("p1h2", "P1[H2]")):
            R[f"t1.mcc.{kk}.{mk}"] = _n(T1_TSV, lambda s, p=proto, c=call, lb=lab, m=m: fx(TW(
                s, T1_TSV, "lineage_equal_mean", dict(section="binary_6.8", protocol=p, call_rule=c, label=lb,
                                                      method=m, metric="mcc", nan_policy="exclude",
                                                      gene_set="parent"))))
        if kk != "fix06":
            R[f"t1.mcc.{kk}.d"] = _n(T1_TSV, lambda s, p=proto, c=call, lb=lab: (lambda w: sci(
                TW(s, T1_TSV, "lineage_equal_mean", w), TW(s, T1_TSV, "ci_lo", w), TW(s, T1_TSV, "ci_hi", w)))(
                dict(section="binary_6.8", protocol=p, call_rule=c, label=lb, method="P1[H1] - H1", metric="mcc",
                     nan_policy="exclude", gene_set="parent")))
    for metric in ("f1", "precision", "recall"):
        for mk, m in (("h1", "H1"), ("p1", "P1[H1]")):
            R[f"t1.p57.{metric}.{mk}"] = _n(T1_TSV, lambda s, mt=metric, m=m: fx(TW(
                s, T1_TSV, "lineage_equal_mean", dict(section="binary_6.8", protocol="Human1-paper", call_rule="P57",
                                                      label="L_ge06", method=m, metric=mt, nan_policy="exclude",
                                                      gene_set="parent"))))
    # T2
    for c in ("C4", "C5"):
        R[f"t2.{c.lower()}"] = _n(T2_GATES, lambda s, c=c: sci(J(s, T2_GATES, "contrasts", c, "estimate"),
                                                               J(s, T2_GATES, "contrasts", c, "ci_low"),
                                                               J(s, T2_GATES, "contrasts", c, "ci_high"),
                                                               3 if c == "C4" else 4))
        R[f"t2.{c.lower()}.status"] = _n(T2_GATES, lambda s, c=c: next(g["status"] for g in J(s, T2_GATES, "gates")
                                                                       if g["contrast"] == c))
        R[f"t2.{c.lower()}.p"] = _n(T2_GATES, lambda s, c=c: _pfmt(J(s, T2_GATES, "contrasts", c, "p")))
        R[f"t2.{c.lower()}.t"] = _n(T2_GATES, lambda s, c=c: ci(J(s, T2_GATES, "transfer", c, "t_low"),
                                                                J(s, T2_GATES, "transfer", c, "t_high"), 3 if c == "C4" else 4))
        R[f"t2.{c.lower()}.k"] = _n(T2_GATES, lambda s, c=c: f"{J(s, T2_GATES, 'transfer', c, 'k_positive')}/"
                                                            f"{J(s, T2_GATES, 'transfer', c, 'n_lineages')}")
    for k, m, gs in (("p2", "P2", "all"), ("b6", "B6[L_prob]", "all"), ("p2.noce", "P2", "noCE"),
                     ("b6.noce", "B6[L_prob]", "noCE")):
        R[f"t2.tap.{k}"] = _n(T2_GRID, lambda s, m=m, g=gs: fx(TW(s, T2_GRID, "tap_lin", dict(
            label="L_prob", geneset=g, universe="U", stratum="none", method=m))))
    R["t2.p2.mdap"] = _n(T2_GRID, lambda s: sx(TW(s, T2_GRID, "m_delta_ap_gene", dict(
        label="L_prob", geneset="all", universe="U", stratum="none", method="P2"))))
    # unsealed T1 tables
    for k, row in (("ii", "ii"),):
        R[f"un.attr.{k}"] = _n(UN_ATTR, lambda s, r=row: sci(*(T(s, UN_ATTR, c, contrast=r)
                                                               for c in ("estimate", "ci_lo", "ci_hi"))))
    for lin in ("Breast", "CNS/Brain", "Eye", "Liver", "Ovary/Fallopian Tube", "Pancreas", "Pleura"):
        key = lin.split("/")[0].split()[0].lower()
        R[f"un.lin.{key}.n"] = _n(UN_LIN, lambda s, l=lin: ints(T(s, UN_LIN, "n_lines", lineage=l, method="H1")))
        for mk, m in (("gstar", "H1"), ("p1", "P1[H1]"), ("b2", "B2[L_prob]")):
            R[f"un.lin.{key}.{mk}"] = _n(UN_LIN, lambda s, l=lin, m=m: fx(T(s, UN_LIN, "mean_tap", lineage=l, method=m)))
        for mk, m in (("d.p1", "P1[H1]"), ("d.p1h2", "P1[H2]"), ("d.h2", "H2")):
            R[f"un.lin.{key}.{mk}"] = _n(UN_LIN, lambda s, l=lin, m=m: sx(T(s, UN_LIN, "delta_vs_ref", lineage=l, method=m)))
    R["un.lin.range"] = _n(UN_LIN, lambda s: "{} to {}".format(
        *(sx(v) for v in (lambda x: (x.min(), x.max()))(s.tsv(UN_LIN).query("method == 'P1[H1]'")["delta_vs_ref"]))))
    R["un.gene.n"] = _n(UN_GENE, lambda s: ints(s.tsv(UN_GENE).query("method == 'P1[H1]' and estimable")["entrez"].nunique()))
    for mk, m in (("p1", "P1[H1]"), ("b3", "B3"), ("b2m", "B2m[L_prob]"), ("b8p", "B8+_H1")):
        R[f"un.gene.{mk}.mean"] = _n(UN_GENE, lambda s, m=m: sx(_gene_dn(s, m).mean()))
    R["un.gene.p1.nz"] = _n(UN_GENE, lambda s: ints((_gene_dn(s, "P1[H1]").abs() > 1e-9).sum()))
    R["un.gene.p1.pos"] = _n(UN_GENE, lambda s: ints((_gene_dn(s, "P1[H1]") > 1e-9).sum()))
    R["un.gene.b3.pos"] = _n(UN_GENE, lambda s: ints((_gene_dn(s, "B3") > 1e-9).sum()))
    R["un.task.adding"] = _n(UN_TASK, lambda s: ints((s.tsv(UN_TASK)["n_failing_not_base_lethal"] > 0).sum()))
    R["un.task.n"] = _n(UN_TASK, lambda s: ints(len(s.tsv(UN_TASK))))
    R["un.task.ci_pos"] = _n(UN_TASK, lambda s: ints((s.tsv(UN_TASK)["delta_lo"] > 0).sum()))
    R["un.task.ci_neg"] = _n(UN_TASK, lambda s: ints((s.tsv(UN_TASK)["delta_hi"] < 0).sum()))
    R["un.task.worst"] = _n(UN_TASK, lambda s: sx(s.tsv(UN_TASK)["delta_lin_vs_base"].min()))
    for cls in ("all_errors", "common_essential", "complex_subunit", "paralog_buffered", "housekeeping_adjacent",
                "expected_fba_failure", "medium_sensitive", "missing_from_model", "unexpressed",
                "biomass_artefact_lethal", "cn_amplified"):
        R[f"un.fail.{cls}.fn"] = _n(UN_FAIL, lambda s, c=cls: ints(T(s, UN_FAIL, "n_fn", **{"class": c})))
        R[f"un.fail.{cls}.fp"] = _n(UN_FAIL, lambda s, c=cls: ints(T(s, UN_FAIL, "n_fp", **{"class": c})))
        R[f"un.fail.{cls}.ffn"] = _n(UN_FAIL, lambda s, c=cls: pct(T(s, UN_FAIL, "frac_of_fn", **{"class": c}), 1))
        R[f"un.fail.{cls}.ffp"] = _n(UN_FAIL, lambda s, c=cls: pct(T(s, UN_FAIL, "frac_of_fp", **{"class": c}), 1))
    R["un.fail.upper"] = _n(UN_FAIL, lambda s: pct(T(s, UN_FAIL, "upper_bound", **{"class": "upper_bound_addressable_positives"}), 1))
    R["un.worst.n"] = _n(UN_WORST, lambda s: ints(len(s.tsv(UN_WORST))))
    R["un.worst.lineages"] = _n(UN_WORST, lambda s: ", ".join(sorted(s.tsv(UN_WORST)["lineage"].unique())))
    R["un.worst.fp"] = _n(UN_WORST, lambda s: ints(s.tsv(UN_WORST)["n_fp"].sum()))
    R["un.worst.genes"] = _n((UN_WORST, H1_MAP, H2_MAP), lambda s: _worst_symbols(s))
    for k, st in (("expr", "expressed_pairs"), ("addr1", "addressable_H1"), ("avana", "lib_Avana"),
                  ("humagne", "lib_Humagne"), ("ky", "lib_KY"), ("organoid", "organoid"), ("cn4", "cn_amp4_in"),
                  ("multi", "multicopy_in"), ("hk", "housekeeping_in"), ("moon", "moonlighting_in")):
        R[f"un.str.{k}"] = _n(UN_STRATA, lambda s, st=st: sci(*(T(s, UN_STRATA, c, stratum=st)
                                                                for c in ("estimate", "ci_lo", "ci_hi"))))
        R[f"un.str.{k}.n"] = _n(UN_STRATA, lambda s, st=st: ints(T(s, UN_STRATA, "n_lines", stratum=st)))
    for k, nm in (("ge05", "C1_L_ge05"), ("ge10", "C1_L_ge10"), ("ucap", "C1_Ucap")):
        R[f"un.sens.{k}"] = _n(UN_SENS, lambda s, nm=nm: sci(*(T(s, UN_SENS, c, name=nm)
                                                               for c in ("estimate", "ci_lo", "ci_hi"))))
    for k, m in (("h1", "H1"), ("h2", "H2"), ("h1task", "H1_task"), ("b0", "B0"), ("b3", "B3"),
                 ("b2", "B2[L_prob]"), ("b2m", "B2m[L_prob]"), ("b8p", "B8+_H1"), ("b8m", "B8-_H1"),
                 ("b8q", "B8'_H1"), ("m2", "M2[H1]"), ("p1", "P1[H1]"), ("p1h2", "P1[H2]")):
        for g in ("all", "noCE"):
            R[f"un.w1.{k}.{g.lower()}"] = _n(UN_W1, lambda s, m=m, g=g: fx(TW(s, UN_W1, "tap_lin", dict(
                label="L_prob", geneset=g, universe="U", stratum="none", method=m))))
    R["un.w1.h2h1"] = _n(UN_W1, lambda s: (lambda w: sci(TW(s, UN_W1, "delta_lin_vs_H1", w),
                                                         TW(s, UN_W1, "delta_lo_vs_H1", w),
                                                         TW(s, UN_W1, "delta_hi_vs_H1", w)))(
        dict(label="L_prob", geneset="all", universe="U", stratum="none", method="H2")))
    R["un.all.h1"] = _n(UN_ALL, lambda s: fx(T(s, UN_ALL, "tap_lin", method="H1")))
    R["un.all.b2"] = _n(UN_ALL, lambda s: fx(T(s, UN_ALL, "tap_lin", method="B2[L_prob]")))
    R["un.all.n"] = _n(UN_ALL, lambda s: ints(T(s, UN_ALL, "n_lines", method="H1")))

    # --- W8 full resource (screened + unscreened) ----------------------------------------------------------------
    R["res.models"] = _n(RES_MODELS, lambda s: ints(len(s.tsv(RES_MODELS))))
    R["res.pass"] = _n(RES_MODELS, lambda s: ints((s.tsv(RES_MODELS)["validation_status"] == "pass").sum()))
    R["res.perparent"] = _n(RES_MODELS, lambda s: ints((s.tsv(RES_MODELS)["parent"] == "human2").sum()))
    R["res.mev"] = _n(RES_MODELS, lambda s: ints((~s.tsv(RES_MODELS)["mevalonate_core_lethal"].astype(bool)).sum()))
    R["res.mb"] = _n(RES_MODELS, lambda s: fx(s.tsv(RES_MODELS)["size_bytes"].sum() / 1e6, 1))
    R["res.lines"] = _n(W8_UNS, lambda s: ints(J(s, W8_UNS, "n_lines_total")))
    R["res.rows"] = _n(W8_UNS, lambda s: ints(J(s, W8_UNS, "n_line_rows_total")))
    R["res.uns.lines"] = _n(W8_UNS, lambda s: ints(J(s, W8_UNS, "n_unscreened_lines")))
    R["res.uns.new"] = _n(W8_UNS, lambda s: f"{J(s, W8_UNS, 'n_new_pass')}/{J(s, W8_UNS, 'n_new_models')}")
    R["res.memote.n"] = _n((W8_MEMOTE, W8_MEMOTE6), lambda s: ints(len(J(s, W8_MEMOTE, "models"))
                                                                    + len(J(s, W8_MEMOTE6, "models"))))
    R["res.memote.h2"] = _n((W8_MEMOTE, W8_MEMOTE6), lambda s: fx(_memote_scores(s, "human2"), 3))
    R["res.memote.h1"] = _n((W8_MEMOTE, W8_MEMOTE6), lambda s: fx(_memote_scores(s, "human1"), 3))
    R["res.frog.n"] = _n(W8_FROG, lambda s: ints(len(J(s, W8_FROG, "models"))))
    R["res.frog.fva"] = _n(W8_FROG, lambda s: ints(J(s, W8_FROG, "models", 0, "status_counts", "fva", "optimal")))
    R["res.frog.gd"] = _n(W8_FROG, lambda s: ints(_frog_count(s, "gene_deletion")))
    R["res.frog.rd"] = _n(W8_FROG, lambda s: ints(_frog_count(s, "reaction_deletion")))
    R["res.frog.cpu"] = _n(W8_FROG, lambda s: fx(sum(v["cpu_s"] for v in J(s, W8_FROG, "models", 0, "timing").values())
                                                 / 3600, 2))
    R["res.noise.default"] = _n(W8_NOISE, lambda s: _sci_e(J(s, W8_NOISE, "per_solver", "glpk_default",
                                                             "reloaded:ENSG00000099817", "value")))
    R["res.noise.exact"] = _n(W8_NOISE, lambda s: fx(J(s, W8_NOISE, "per_solver", "glpk_exact",
                                                       "reloaded:ENSG00000099817", "value"), 1))

    # --- runtime ------------------------------------------------------------------------------------------------
    for g in RUNTIME_GROUPS + ("total", "other"):
        R[f"rt.{g}.cpu"] = _n(RUNTIME_SRC, lambda s, g=g: fx(_rt(s, g), 1))
        R[f"rt.{g}.wall"] = _n(RUNTIME_SRC, lambda s, g=g: fx(_rt(s, g, "wall_s"), 1))
    for g in RUNTIME_GROUPS + ("total",):
        R[f"rt.{g}.runs"] = _n(RUNTIME_SRC, lambda s, g=g: ints(_rt_runs(s, g)))
    R["rt.budget.w7"] = _n(BUDGET, lambda s: fx(s.tsv(BUDGET).query("workstream == 'W7'")["cpu_h"].sum(), 1))
    R["rt.elapsed"] = _n(GIT_LOG, lambda s: fx(_elapsed_h(s), 1))
    return R


NUMBERS: dict[str, Number] = _registry()

# =========================================================================================================
# Marker check / update
# =========================================================================================================
MARKER = re.compile(r"<!--n:([A-Za-z0-9_.\-]+)-->")
PLACEHOLDER = re.compile(r"\[(T1|T2)(?: unseal)?:[^\]]*\]")
_NUM = r"[+−\-]?\d+(?:,\d{3})*(?:\.\d+)?"  # commas only as thousands separators


def shape_regex(rendered: str) -> re.Pattern:
    """A regex that matches ``rendered`` with every number free to change (sign included)."""
    parts, last = [], 0
    for m in re.finditer(r"(?<![A-Za-z_])" + _NUM, rendered):  # digits glued to words (Human1, W7-05) stay literal
        parts.append(re.escape(rendered[last:m.start()]))
        parts.append(_NUM)
        last = m.end()
    parts.append(re.escape(rendered[last:]))
    return re.compile("".join(parts) + r"(?!\d|\.\d)")


def _continues_number(rest: str) -> bool:
    """True if the text right after a value would extend it (a digit, or a decimal point and a digit)."""
    return bool(re.match(r"\d|\.\d", rest))


def _block_of(text: str, pos: int) -> str:
    """The paragraph, list or table (text between blank lines) that contains ``pos``. A table also takes in the
    paragraph just before it (its caption) and just after it (its note), where its sources are cited."""
    start = text.rfind("\n\n", 0, pos)
    end = text.find("\n\n", pos)
    start = 0 if start < 0 else start
    end = len(text) if end < 0 else end
    if text[start:end].lstrip("\n").lstrip(" ").startswith("|"):  # tables, also when indented in a list
        prev = text.rfind("\n\n", 0, start)
        nxt = text.find("\n\n", end + 2)
        start = 0 if prev < 0 else prev
        end = len(text) if nxt < 0 else nxt
    return text[start:end]


@dataclass
class Issue:
    key: str
    line: int
    kind: str
    detail: str

    def __str__(self) -> str:
        return f"REPORT.md:{self.line}: {self.kind} [{self.key}] {self.detail}"


def check_text(text: str, src: Source, numbers: dict[str, Number] | None = None) -> tuple[list[Issue], int]:
    """Verify every marker in ``text``. Returns (issues, n_markers)."""
    numbers = NUMBERS if numbers is None else numbers
    issues: list[Issue] = []
    n = 0
    for m in MARKER.finditer(text):
        n += 1
        key = m.group(1)
        line = text.count("\n", 0, m.start()) + 1
        if key not in numbers:
            issues.append(Issue(key, line, "unknown-key", "no entry in NUMBERS"))
            continue
        num = numbers[key]
        try:
            want = num.fn(src)
        except Exception as exc:  # a missing file or row is a failure, not a skip
            issues.append(Issue(key, line, "source-error", f"{type(exc).__name__}: {exc}"))
            continue
        after = text[m.end():]
        if after.startswith(want) and not _continues_number(after[len(want):]):
            pass
        else:
            got = shape_regex(want).match(after)
            shown = got.group(0) if got else after[:40].split("\n")[0]
            issues.append(Issue(key, line, "mismatch", f"REPORT has {shown!r}, source gives {want!r}"))
        block = _block_of(text, m.start())
        missing = [p for p in num.sources if p not in block]
        if missing:
            issues.append(Issue(key, line, "uncited", f"cite {missing} in the same paragraph or table"))
    return issues, n


def update_text(text: str, src: Source, numbers: dict[str, Number] | None = None) -> tuple[str, list[Issue]]:
    """Rewrite each marked value from its source. A value written as ``?`` is filled; otherwise the old value must
    have the same shape as the new one (same words, numbers free), so prose is never overwritten."""
    numbers = NUMBERS if numbers is None else numbers
    out, last, issues = [], 0, []
    for m in MARKER.finditer(text):
        key = m.group(1)
        line = text.count("\n", 0, m.start()) + 1
        out.append(text[last:m.end()])
        last = m.end()
        if key not in numbers:
            issues.append(Issue(key, line, "unknown-key", "no entry in NUMBERS"))
            continue
        try:
            want = numbers[key].fn(src)
        except Exception as exc:
            issues.append(Issue(key, line, "source-error", f"{type(exc).__name__}: {exc}"))
            continue
        after = text[m.end():]
        if after.startswith("?"):
            out.append(want)
            last = m.end() + 1
            continue
        got = shape_regex(want).match(after)
        if not got:
            issues.append(Issue(key, line, "unmatched", f"cannot locate the old value near {after[:30]!r}"))
            continue
        out.append(want)
        last = m.end() + got.end()
    out.append(text[last:])
    return "".join(out), issues


def placeholders(text: str) -> list[str]:
    return [m.group(0) for m in PLACEHOLDER.finditer(text)]


def unmarked_decimals(text: str) -> list[tuple[int, str]]:
    """Lint: decimals with >= 3 places that carry no marker (prose quoted from READMEs, design constants)."""
    hits = []
    for i, line in enumerate(text.splitlines(), 1):
        clean = MARKER.sub("", line)
        marked_spans = []
        for mm in MARKER.finditer(line):
            num = NUMBERS.get(mm.group(1))
            if num is None:
                continue
            g = re.match(r"\S[^<]*", line[mm.end():])
            if g:
                marked_spans.append(g.group(0))
        rest = clean
        for sp in marked_spans:
            rest = rest.replace(sp, "", 1)
        rest = re.sub(r"doi:\S+|10\.\d{4,}/\S+", "", rest)
        for d in re.findall(r"(?<![\w.])[+−\-]?\d+\.\d{3,}", rest):
            hits.append((i, d))
    return hits


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--report", default=str(ROOT / "REPORT.md"))
    ap.add_argument("--check", action="store_true", help="fail if any marked number disagrees with its source")
    ap.add_argument("--update", action="store_true", help="rewrite the marked numbers from their sources")
    ap.add_argument("--sources", action="store_true", help="list the T1/T2 sources and whether they exist")
    ap.add_argument("--worktree", action="store_true", help="read sources from the working tree, not HEAD")
    ap.add_argument("--ref", default="HEAD", help="git ref to read sources from (default HEAD)")
    ap.add_argument("--require-filled", action="store_true", help="also fail while [T1: ...]/[T2: ...] remain")
    ap.add_argument("--lint", action="store_true", help="list unmarked decimals with >= 3 places")
    args = ap.parse_args(argv)

    if args.sources:
        for k, ok in collect_numbers(ROOT / "results").items():
            print(f"{k}\t{'present' if ok else 'missing'}\t{NUMBER_SOURCES[k][1]}")
        return 0
    src = Source(ROOT, None if args.worktree else args.ref)
    path = Path(args.report)
    text = path.read_text()
    if args.update:
        new, issues = update_text(text, src)
        for i in issues:
            print(i, file=sys.stderr)
        if new != text:
            path.write_text(new)
        text = new
        if issues:
            return 1
    issues, n = check_text(text, src)
    ph = placeholders(text)
    for i in issues:
        print(i, file=sys.stderr)
    print(f"markers checked: {n}; issues: {len(issues)}; placeholders open: {len(ph)} "
          f"(sources: {'working tree' if args.worktree else args.ref})")
    if args.lint:
        for line, d in unmarked_decimals(text):
            print(f"lint: REPORT.md:{line}: unmarked decimal {d}")
    if issues or (args.require_filled and ph):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
