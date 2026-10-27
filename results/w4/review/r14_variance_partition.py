"""W4 Stage-A review: own recomputation of the method-vs-tissue partition (label-free).

Reactions level (open masks) and binary KO level (ratio_r < 0.01), A1 builds F0/F2/F3/F4 x {R1, R3} x {real, perm}:
- marginal PERMANOVA R^2 (Anderson 2001) of method, threshold and lineage on real builds (Jaccard dissimilarity);
- H4.1 diff with three intervals: line bootstrap within lineage (as the author did), and a duplicate-free
  half-sampling check (2 of 3 lines per lineage without replacement; 2,000 draws);
- H4.2 for F0/F2/F3: lineage R^2 on real builds vs permuted builds, with a null that relabels lineages across lines.
"""

import json

import numpy as np
import pandas as pd

from nextgem import paths
from nextgem.config import SEEDS

W4 = paths.results_dir() / "w4"
b = pd.read_csv(W4 / "stage_a" / "builds.tsv", sep="\t", usecols=["build_id", "row", "line", "family", "rule", "variant", "genomic", "open_mask"], low_memory=False)
b = b[(b.row == "A1") & (b.genomic == "none") & b.family.isin(["F0", "F2", "F3", "F4"])].reset_index(drop=True)
split = pd.read_csv(paths.splits_dir() / "lineage_split.tsv", sep="\t").set_index("ModelID")
b["lineage"] = split.lineage.reindex(b.line).to_numpy()
n_rxn = 12748
X_rx = np.vstack([np.unpackbits(np.frombuffer(bytes.fromhex(m), dtype=np.uint8))[:n_rxn] for m in b.open_mask]).astype(bool)
ko = pd.read_parquet(W4 / "stage_a" / "ko_A1.parquet", columns=["build_id", "gene", "ratio_r"])
ko["build_id"] = ko.build_id.astype(str)
K = ko.pivot(index="build_id", columns="gene", values="ratio_r")
X_ko = (K.reindex(b.build_id.astype(str)) < 0.01).to_numpy()
assert not np.isnan(K.reindex(b.build_id.astype(str)).to_numpy()).all(axis=1).any()


def D2(X):
    X = X.astype(float)
    inter = X @ X.T
    s = X.sum(1)
    un = s[:, None] + s[None, :] - inter
    J = np.divide(inter, un, out=np.ones_like(inter), where=un > 0)
    return 1 - J


def R2(D, g):
    g = np.asarray(g)
    sst = D.sum() / (2 * len(g))
    ssw = sum(D[np.ix_(i, i)].sum() / (2 * len(i)) for i in (np.flatnonzero(g == k) for k in np.unique(g)))
    return (sst - ssw) / sst if sst > 0 else np.nan


rng = np.random.default_rng(SEEDS[0])
res = {}
for level, X in (("reactions", X_rx), ("ko_binary", X_ko)):
    real = (b.variant == "real").to_numpy()
    D = D2(X[real])
    m = b[real].reset_index(drop=True)
    r_meth, r_lin, r_thr = R2(D, m.family), R2(D, m.lineage), R2(D, m.rule)
    diff = r_meth - r_lin
    lines_by_lin = m.groupby("lineage").line.unique().to_dict()
    idx_of = m.groupby("line").indices
    half = []
    for _ in range(2000):
        pick = [l for ls in lines_by_lin.values() for l in rng.choice(ls, 2, replace=False)]
        ii = np.concatenate([idx_of[l] for l in pick])
        mm = m.iloc[ii]
        half.append(R2(D[np.ix_(ii, ii)], mm.family) - R2(D[np.ix_(ii, ii)], mm.lineage))
    boot = []
    for _ in range(2000):
        pick = [l for ls in lines_by_lin.values() for l in rng.choice(ls, len(ls), replace=True)]
        ii = np.concatenate([idx_of[l] for l in pick])
        mm = m.iloc[ii]
        boot.append(R2(D[np.ix_(ii, ii)], mm.family) - R2(D[np.ix_(ii, ii)], mm.lineage))
    h42 = {}
    for fam in ("F0", "F2", "F3"):
        out = {}
        for var in ("real", "perm"):
            sel = ((b.family == fam) & (b.variant == var)).to_numpy()
            Dm = D2(X[sel])
            mm = b[sel].reset_index(drop=True)
            out[var] = R2(Dm, mm.lineage)
            if var == "perm":
                lines = mm.line.unique()
                lin_of = mm.groupby("line").lineage.first()
                null = []
                for _ in range(999):
                    perm = dict(zip(lines, rng.permutation(lin_of.reindex(lines).to_numpy())))
                    null.append(R2(Dm, mm.line.map(perm)))
                out["perm_null_q95"] = float(np.quantile(null, 0.95))
        out["tissue_informed"] = bool(out["real"] > max(out["perm"], out["perm_null_q95"]))
        h42[fam] = out
    res[level] = {"r2_method": r_meth, "r2_lineage": r_lin, "r2_threshold": r_thr, "diff": diff,
                  "boot_pct": list(np.quantile(boot, [0.025, 0.975])), "boot_mean": float(np.mean(boot)),
                  "half_sampling_pct": list(np.quantile(half, [0.025, 0.975])), "half_mean": float(np.mean(half)),
                  "H4.2": h42}
(W4 / "review" / "r14_variance_partition.json").write_text(json.dumps(res, indent=1, default=float))
print(json.dumps(res, indent=1, default=float))
