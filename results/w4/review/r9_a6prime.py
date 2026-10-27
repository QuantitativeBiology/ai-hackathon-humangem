"""W4 re-check: independent reproduction of the A6' line choice (amendment A7), label-free.
Per panel lineage: median log2(TPM+1) over ALL training lines of the lineage with expression (model genes with a
DepMap column, own Entrez join), Euclidean distance of each of the 3 panel lines, minimum (ties: ModelID)."""
import json, re
import cobra
cobra.Configuration().solver = "glpk"
import numpy as np, pandas as pd, pyarrow.csv as pacsv
from nextgem import paths
R = paths.results_dir() / "w4"
m = cobra.io.read_sbml_model(str(paths.data_root() / "Human-GEM-main/model/Human-GEM.xml"))
gt = pd.read_csv(paths.data_root() / "Human-GEM-main/model/genes.tsv", sep="\t", dtype=str).set_index("genes")
ent = gt.geneEntrezID.reindex(sorted(g.id for g in m.genes)).str.split(";").str[0]
f = paths.depmap_dir() / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
col = {}
for c in pd.read_csv(f, nrows=0).columns:
    mm = re.match(r"^(.*) \((\d+)\)$", c)
    if mm: col.setdefault(mm.group(2), c)
want = {g: col[e] for g, e in ent.items() if isinstance(e, str) and e in col}
t = pacsv.read_csv(f, convert_options=pacsv.ConvertOptions(include_columns=["ModelID", "IsDefaultEntryForModel"] + sorted(set(want.values())))).to_pandas()
t = t[t.IsDefaultEntryForModel == "Yes"].set_index("ModelID")
X = pd.DataFrame({g: t[c] for g, c in want.items()}).astype(float)
split = pd.read_csv(paths.splits_dir() / "lineage_split.tsv", sep="\t")
tr = split[(split.tier == "train") & split.has_expression.astype(bool)]
panel = pd.read_csv(R / "stage_a_panel.tsv", sep="\t")
auth = pd.read_csv(R / "a6prime_lines.tsv", sep="\t").set_index("lineage")
rows = []
for lin, g in panel.groupby("lineage"):
    mem = [x for x in tr[tr.lineage == lin].ModelID if x in X.index]
    med = X.loc[mem].median()
    d = {x: float(np.sqrt(((X.loc[x] - med) ** 2).sum())) for x in g.ModelID}
    best = min(sorted(d), key=lambda x: (d[x], x))
    gap = sorted(d.values())[1] - sorted(d.values())[0]
    panel_only = X.loc[list(g.ModelID)].median()
    d2 = {x: float(np.sqrt(((X.loc[x] - panel_only) ** 2).sum())) for x in g.ModelID}
    rows.append({"lineage": lin, "own": best, "author": auth.loc[lin, "ModelID"], "equal": best == auth.loc[lin, "ModelID"],
                 "n_lineage_lines": len(mem), "author_n": int(auth.loc[lin, "n_lineage_lines"]),
                 "dist_own": round(d[best], 4), "dist_author": round(float(auth.loc[lin, "distance"]), 4), "margin_to_2nd": round(gap, 3),
                 "choice_if_panel_only_median": min(sorted(d2), key=lambda x: (d2[x], x)),
                 "all_members_train": bool(split.set_index("ModelID").loc[mem, "tier"].eq("train").all()),
                 "tag_ok": auth.loc[lin, "status"] == "descriptive, not eligible for P1 (amendment A7)"})
df = pd.DataFrame(rows); df.to_csv(R / "review" / "r9_a6prime.tsv", sep="\t", index=False); print(df.to_string())
