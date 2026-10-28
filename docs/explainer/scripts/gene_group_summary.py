"""Which genes does each P1 step move, and how often is each group essential in the 380 test lines?

Label-free model calls (results/m6/scores/sweep_e.parquet; Human1, Ham's signature and the classic-RPMI + serum
signature shared by 1,050 of 1,361 screened lines) joined to the committed, unsealed per-gene test summary
results/m6/test/unsealed/T1/per_gene.tsv (column `prevalence` = fraction of test lines in which the gene is essential).
No label file is read; output is per-group aggregates only.
Usage: /data/tmpA/nextgem/.venv/bin/python gene_group_summary.py --repo /data/tmpA/nextgem --out <dir>
"""
import argparse
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument("--repo", default="/data/tmpA/nextgem"); ap.add_argument("--out", required=True)
a = ap.parse_args(); R = Path(a.repo) / "results"; out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
HAMS, RPMI = "human1:a35fcfdb65d63819", "human1:db9ebc6b123dcb99"
BLOO = "B-loo:cofactor_pool_biomass"
keys = {"s0": f"T0|{HAMS}", "s1": f"{BLOO}|{HAMS}", "sm": f"T0|{RPMI}", "s2": f"{BLOO}|{RPMI}"}
e = pd.read_parquet(R / "m6/scores/sweep_e.parquet")
w = e[e.key.isin(keys.values())].pivot(index="entrez", columns="key", values="e").rename(columns={v: k for k, v in keys.items()})
g = pd.read_csv(R / "m6/test/unsealed/T1/per_gene.tsv", sep="\t"); g = g[g.method == "P1[H1]"].set_index("entrez")
w["frac"] = g.prevalence
mp = pd.read_csv(R / "baseline/human1_gene_mapping.tsv", sep="\t").dropna(subset=["entrez"]); mp["entrez"] = mp.entrez.astype(int)
w["sym"] = mp.drop_duplicates("entrez").set_index("entrez").symbol.reindex(w.index)
L = {k: w[k] >= 0.99 for k in keys}
real = {
    "ALAD", "CPOX", "HMBS", "UROD", "UROS", "ABCB6", "COQ2", "COQ3", "COQ6", "COQ7", "GCH1", "PTS", "SPR", "PDXK",
    "PNPO", "TPK1", "RFK", "MMAA", "SLC5A6", "SLC25A26"}  # the 20 genes F2 keeps lethal (results/f2/README.md §4)
pool = L["s0"] & ~L["s1"]
groups = {
    "G* lethal (generic Human1, Ham's)": L["s0"],
    "P1 lethal (classic RPMI + serum)": L["s2"],
    "P1 partial (0 < e < 0.99)": (w.s2 > 0) & (w.s2 < 0.99),
    "P1 no effect (e = 0)": w.s2 == 0,
    "un-called by pool removal (35)": pool,
    "  of which real-cofactor genes kept lethal by F2 (20)": pool & w.sym.isin(real),
    "  of which catabolite/carnitine artefacts (15)": pool & ~w.sym.isin(real),
    "un-called by medium + serum alone (13)": L["s0"] & ~L["sm"],
    "un-called only when both fixes combine (10)": L["s1"] & L["sm"] & ~L["s2"],
    "newly called by dropping Ham's nucleosides (8)": ~L["s0"] & L["s2"],
    "lethal at every step": L["s0"] & L["s1"] & L["sm"] & L["s2"],
    "all 2,866 scored genes": w.s2 == w.s2,
}
rows = [{"group": k, "n_genes": int(m.sum()), "mean_frac_test_lines_essential": round(float(w.frac[m].mean()), 3),
         "median_frac": round(float(w.frac[m].median()), 3),
         "genes": ";".join(sorted(w.sym[m].astype(str))) if m.sum() <= 40 else ""} for k, m in groups.items()]
pd.DataFrame(rows).to_csv(out / "gene_groups.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).drop(columns="genes").to_string(index=False))
