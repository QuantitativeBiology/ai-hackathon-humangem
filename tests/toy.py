"""Synthetic (toy) fixtures with known answers for the M1 infrastructure tests. No real label data.

The toy DepMap files carry the real label file *names* (so the tripwire treats them as label files) but
contain synthetic numbers generated here from a fixed seed.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
PY = str(PROJECT / ".venv" / "bin" / "python")
SCRIPTS = PROJECT / "scripts"

# 30 lines: 3 eligible training lineages (6 lines each), 2 rare training lines, validation 5, test 5
LINEAGES = [("LinA", "train", "S1", 6), ("LinB", "train", "S1", 6), ("LinC", "train", "S1", 6),
            ("LinR", "train", "rare", 2), ("LinV", "validation", "S2", 5), ("LinT", "test", "S3", 5)]
N_GENES = 80  # CRISPR columns G0..G79 = Entrez 1000..1079
U_ENTREZ = list(range(1000, 1070))  # U = Entrez 1000..1069 (via the toy baseline maps)
RAW_GE = "CRISPRGeneEffect.csv"
RAW_DP = "CRISPRGeneDependency.csv"
# planted NaNs (row index, column index) -> known per-gene counts
GE_NAN = [(2, 10), (6, 10), (11, 20), (25, 75)]  # gene 1010: 2, 1020: 1 (1075 is outside U)
DP_NAN = [(2, 30), (25, 30), (27, 31)]  # gene 1030: 2, 1031: 1


def git(repo: Path, *args, check=True, env=None) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, check=check, env=env)


def init_repo(repo: Path) -> None:
    repo.mkdir(parents=True, exist_ok=True)
    git(repo, "init", "-q", "-b", "main")
    git(repo, "config", "user.name", "toy")
    git(repo, "config", "user.email", "toy@example.org")
    git(repo, "config", "commit.gpgsign", "false")


def commit_all(repo: Path, msg: str, env=None) -> str:
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", msg, "--allow-empty", env=env)
    return git(repo, "rev-parse", "HEAD").stdout.strip()


DESIGN_TOY = """# toy DESIGN

## 17. Amendment log and revision history

### Amendment log (after the M1 freeze)
| Date | Section | Change | Reason | Before/after the affected analysis |
|---|---|---|---|---|
| 2026-10-27 | 4.1 | A1: toy amendment | test | before |

### Revision history
nothing
"""


@dataclass
class Toy:
    root: Path
    data: Path
    state: Path
    model_ids: list[str]
    tiers: dict[str, str]
    lineage: dict[str, str]
    ge: pd.DataFrame  # lines x Entrez (all 80)
    dp: pd.DataFrame
    env: dict = field(default_factory=dict)

    def ids_in(self, tier: str) -> list[str]:
        return [m for m in self.model_ids if self.tiers[m] == tier]

    def run(self, script: str, *args: str, extra_env: dict | None = None) -> subprocess.CompletedProcess:
        env = {**os.environ, **self.env, **(extra_env or {})}
        return subprocess.run([PY, str(SCRIPTS / script), *args], capture_output=True, text=True, env=env)


def make_toy(root: Path, kappa: str = "high", seed: int = 20) -> Toy:
    """Build toy data + state roots. ``kappa='high'``: dependency tied to gene effect (kappa >= 0.5);
    ``'low'``: dependency independent of gene effect (kappa ~ 0 -> the rule switches to L_ge05)."""
    rng = np.random.default_rng(seed)
    data, state = root / "data", root / "state"
    dm = data / "depmap"
    dm.mkdir(parents=True)
    init_repo(state)

    ids, tiers, lineage, strata = [], {}, {}, {}
    for lin, tier, stratum, n in LINEAGES:
        for _ in range(n):
            mid = f"ACH-{len(ids) + 1:06d}"
            ids.append(mid)
            tiers[mid], lineage[mid], strata[mid] = tier, lin, stratum
    entrez = [1000 + i for i in range(N_GENES)]
    ge = rng.normal(-0.15, 0.35, (len(ids), N_GENES))
    ge[:, :5] = rng.normal(-1.5, 0.2, (len(ids), 5))  # common-essential-like genes
    if kappa == "high":
        dp = 1 / (1 + np.exp(8 * (ge + 0.5)))
    else:
        dp = rng.uniform(0, 1, ge.shape)
    ge[0, 40], dp[0, 40] = -0.55, 0.5  # L_ge05=1, L_ge06=0, L_ge10=0, L_prob=0 (strict >)
    ge[0, 41], dp[0, 41] = -1.2, 0.9  # all four labels 1
    ge[0, 42], dp[0, 42] = -0.6, 0.51  # L_ge06=0 (strict <), L_ge05=1, L_prob=1
    for r, c in GE_NAN:
        ge[r, c] = np.nan
    for r, c in DP_NAN:
        dp[r, c] = np.nan
    cols = [f"G{i} ({1000 + i})" for i in range(N_GENES)]
    ge_df = pd.DataFrame(ge, index=ids, columns=cols)
    dp_df = pd.DataFrame(dp, index=ids, columns=cols)
    ge_df.to_csv(dm / RAW_GE, index_label="")
    dp_df.to_csv(dm / RAW_DP, index_label="ModelID")
    pd.DataFrame({"Essentials": ["G0 (1000)", "G1 (1001)", "G75 (1075)"]}).to_csv(dm / "CRISPRInferredCommonEssentials.csv", index=False)
    pd.DataFrame({"Gene": ["G2 (1002)", "G3 (1003)"]}).to_csv(dm / "AchillesCommonEssentialControls.csv", index=False)
    pd.DataFrame({"Gene": ["G60 (1060)", "G61 (1061)", "G77 (1077)"]}).to_csv(dm / "AchillesNonessentialControls.csv", index=False)

    # toy baseline gene maps: U = Entrez 1000..1069. G5's model symbol differs (joined by Entrez);
    # a model gene "G75" with another Entrez and no CRISPR match is not in U.
    bd = state / "results" / "baseline"
    bd.mkdir(parents=True)
    h1 = pd.DataFrame({"gene": [f"ENSG{e}" for e in range(1000, 1060)], "symbol": [f"G{e - 1000}" for e in range(1000, 1060)],
                       "entrez": [str(e) for e in range(1000, 1060)], "depmap_column": [f"G{e - 1000} ({e})" for e in range(1000, 1060)]})
    h1.loc[5, "symbol"] = "G5_ALIAS"
    h2 = pd.DataFrame({"gene": [f"ENSG{e}" for e in range(1010, 1070)] + ["ENSG_X"],
                       "symbol": [f"G{e - 1000}" for e in range(1010, 1070)] + ["G75"],
                       "entrez": [str(e) for e in range(1010, 1070)] + ["5555"],
                       "depmap_column": [f"G{e - 1000} ({e})" for e in range(1010, 1070)] + [None]})
    h1.to_csv(bd / "human1_gene_mapping.tsv", sep="\t", index=False)
    h2.to_csv(bd / "human2_gene_mapping.tsv", sep="\t", index=False)

    sp = state / "results" / "splits"
    sp.mkdir(parents=True)
    split = pd.DataFrame({"ModelID": ids, "tier": [tiers[m] for m in ids], "stratum": [strata[m] for m in ids],
                          "lineage": [lineage[m] for m in ids],
                          "library_class": ["KY-only" if i < 4 else "Avana-containing" for i in range(len(ids))]})
    split.to_csv(sp / "lineage_split.tsv", sep="\t", index=False)
    sha = hashlib.sha256((sp / "lineage_split.tsv").read_bytes()).hexdigest()
    (sp / "SPLIT_SHA256").write_text(f"{sha}  lineage_split.tsv\n")
    (state / "DESIGN.md").write_text(DESIGN_TOY)
    commit_all(state, "toy init")

    ge_e = pd.DataFrame(ge, index=ids, columns=entrez)
    dp_e = pd.DataFrame(dp, index=ids, columns=entrez)
    toy = Toy(root, data, state, ids, tiers, lineage, ge_e, dp_e)
    toy.env = {"NEXTGEM_DATA_ROOT": str(data), "NEXTGEM_STATE_ROOT": str(state)}
    return toy


def build(toy: Toy) -> subprocess.CompletedProcess:
    return toy.run("build_labels.py", "--data-root", str(toy.data), "--state-root", str(toy.state))


def tap(scores, y) -> float:
    """Tie-aware expected AP (DESIGN §6.1, O(n_g) form). ``y`` binary."""
    s = np.asarray(scores, dtype=float)
    y = np.asarray(y, dtype=int)
    total = y.sum()
    if total == 0:
        return float("nan")
    order = np.argsort(-s, kind="stable")
    s, y = s[order], y[order]
    N = P = 0
    acc = 0.0
    i = 0
    while i < len(s):
        j = i
        while j < len(s) and s[j] == s[i]:
            j += 1
        n, p = j - i, int(y[i:j].sum())
        if p:
            for r in range(1, n + 1):
                frac = (r - 1) * (p - 1) / (n - 1) if n > 1 else 0.0
                acc += (p / n) * (P + 1 + frac) / (N + r)
        N += n
        P += p
        i = j
    return acc / total


def null_ap(y) -> float:
    return tap(np.zeros(len(y)), y)
