"""Reviewer's independent recompute of the M5 A12(b) panel pre-check (training tier only, L_prob, universe U)."""
import json, sys
import numpy as np, pandas as pd
from nextgem import labels, metrics, paths
from nextgem.config import SEEDS

R = paths.results_dir()
panel = pd.read_csv(R / "w4" / "stage_a_panel.tsv", sep="\t")
lines = list(panel.ModelID)
split = labels.load_split().set_index("ModelID")
lin = split.lineage.reindex(lines)
U = labels.gene_universe()
genes = [int(x) for x in U.entrez]

# labels: training tier, panel lines, L_prob only
tl = labels.train(columns=["L_prob"])
tl = tl[tl.ModelID.isin(lines)]
Y = tl.pivot(index="ModelID", columns="entrez", values="L_prob").reindex(index=lines, columns=genes).astype(float)

# Human1 Entrez -> gene map (W1 convention: sort by gene, drop duplicate Entrez, keep first)
b1 = pd.read_parquet(paths.baseline_dir() / "human1_ko.parquet", columns=["gene", "entrez", "ratio_r"])
m = b1.dropna(subset=["entrez"]).sort_values("gene")
m = m.assign(entrez=m.entrez.astype(int)).drop_duplicates("entrez").set_index("entrez").gene
g_of = m.reindex(genes)
have = g_of.notna().to_numpy()

def escore(ratio: pd.Series) -> np.ndarray:
    """e = round(1 - ratio, 6); NaN ratio -> ratio 1 -> e 0; U genes with no Human1 gene -> 0."""
    r = ratio.reindex(g_of[have].to_numpy()).to_numpy(float)
    e = np.zeros(len(genes))
    e[have] = np.round(1.0 - np.nan_to_num(r, nan=1.0), 6)
    return e

ko = pd.read_parquet(R / "m5" / "builds" / "ko.parquet", columns=["build_id", "gene", "ratio_r"])
ko["build_id"] = ko.build_id.astype(str)
tabs = {b: g.set_index("gene").ratio_r for b, g in ko.groupby("build_id")}
bt = pd.read_csv(R / "m5" / "builds" / "builds.tsv", sep="\t", low_memory=False)
print("failed builds in builds.tsv:", int(bt.failed.fillna(False).astype(bool).sum()), "| missing KO tables:",
      sum(f"M5|{l}|{f}|R1|rna|real|none" not in tabs for l in lines for f in ("F2", "F4")))
failed = set()  # none failed (checked above); a failed line would take the all-tied score 0 in both arms of a contrast

def frame(fn):
    return pd.DataFrame(np.vstack([escore(tabs[fn(l)]) for l in lines]), index=lines, columns=genes)

S = {"genericH1": pd.DataFrame(np.tile(escore(b1.set_index("gene").ratio_r), (len(lines), 1)), index=lines, columns=genes),
     "H1+Bloo": frame(lambda l: "generic|none"),
     "F2+Bloo": frame(lambda l: f"M5|{l}|F2|R1|rna|real|none"),
     "F4+Bloo": frame(lambda l: f"M5|{l}|F4|R1|rna|real|none")}
tap = {}
for k, s in S.items():
    pl = metrics.per_line(s, Y)
    tap[k] = pl.tap.where(pl.estimable)
    print(k, "estimable lines", int(pl.estimable.sum()), "lineage-equal tAP %.5f" % metrics.lineage_equal_mean(tap[k].dropna(), lin.reindex(tap[k].dropna().index)))

def own_boot(d, lineage, B, seed):
    """Own lineage-stratified paired bootstrap of the lineage-equal mean (resample lines within lineage)."""
    rng = np.random.default_rng(seed)
    groups = [d[lineage.reindex(d.index) == k].to_numpy() for k in sorted(lineage.reindex(d.index).unique())]
    est = np.mean([g.mean() for g in groups])
    reps = np.empty(B)
    for b in range(B):
        reps[b] = np.mean([rng.choice(g, size=len(g), replace=True).mean() for g in groups])
    return est, np.percentile(reps, 2.5), np.percentile(reps, 97.5)

out = {}
for name, a, b in [("F2 vs H1+Bloo", "F2+Bloo", "H1+Bloo"), ("F4 vs H1+Bloo", "F4+Bloo", "H1+Bloo"),
                   ("H1+Bloo vs genericH1 (sanity)", "H1+Bloo", "genericH1")]:
    d = (tap[a] - tap[b]).dropna()
    l = lin.reindex(d.index)
    bm = metrics.bootstrap_mean(d, l, B=2000, seed=SEEDS[0])
    oe, olo, ohi = own_boot(d, l, 2000, SEEDS[0])
    fm = d.groupby(l).mean()
    cnt = d.groupby(l).size()
    pos = int(((fm > 0) & (cnt >= 2)).sum())
    out[name] = dict(delta_lin=bm.estimate, ci_metrics=[bm.ci_low, bm.ci_high], ci_own_boot=[olo, ohi], own_est=oe,
                     folds_positive=pos, n_folds=int(fm.size), n_lines=int(d.size), delta_line_mean=float(d.mean()))
    print(f"{name}: {bm.estimate:+.5f} CI(metrics) [{bm.ci_low:+.5f},{bm.ci_high:+.5f}] CI(own) [{olo:+.5f},{ohi:+.5f}] "
          f"folds+ {pos}/{fm.size}  n={d.size}")
for f in ("F2", "F4"):
    v = out[f"{f} vs H1+Bloo"]["delta_lin"]
    print(f, "enters Stage B (>= +0.005):", bool(v >= 0.005 - 1e-12))
(R / "m5" / "review" / "rv_m5_result.json").write_text(json.dumps(out, indent=1, default=float))
