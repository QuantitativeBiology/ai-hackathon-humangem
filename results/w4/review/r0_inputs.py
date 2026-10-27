"""W4 review (independent reviewer), part 0: label-free inputs, thresholds, z-statistics, panel draw, A5 grid.

Reads only omics inputs (expression, ProCan via the author's reader) and the canonical split file. No label
file, no labels.* call. GLPK only (no Gurobi session). Writes small JSON/TSV outputs to results/w4/review/ and
regenerable intermediates to results/w4/review/tmp/.
"""

import hashlib
import json
import re
import time

import cobra
import numpy as np
import pandas as pd
import pyarrow.csv as pacsv

cobra.Configuration().solver = "glpk"  # structure only; never opens a Gurobi session

from nextgem import context_inputs as ci  # noqa: E402  (author's code, only for side-by-side comparison)
from nextgem import paths  # noqa: E402
from nextgem.config import SEEDS  # noqa: E402

R = paths.results_dir() / "w4"
OUT = R / "review"
TMP = OUT / "tmp"
TMP.mkdir(parents=True, exist_ok=True)
t0 = time.time()
report = {}

# ---------------------------------------------------------------- model genes and own Entrez mapping
m = cobra.io.read_sbml_model(str(paths.data_root() / "Human-GEM-main/model/Human-GEM.xml"))
gtsv = pd.read_csv(paths.data_root() / "Human-GEM-main/model/genes.tsv", sep="\t", dtype=str).set_index("genes")
model_genes = sorted(g.id for g in m.genes)
entrez = gtsv.geneEntrezID.reindex(model_genes).str.split(";").str[0]

split = pd.read_csv(paths.splits_dir() / "lineage_split.tsv", sep="\t")
train = split[split.tier == "train"]
train_expr_lines = sorted(train[train.has_expression.astype(bool)].ModelID)

# own reader: default rows, SYMBOL (ENTREZ) columns
f = paths.depmap_dir() / ci.EXPR_FILE
hdr = pd.read_csv(f, nrows=0).columns
col_of = {}
for c in hdr:
    mm = re.match(r"^(.*) \((\d+)\)$", c)
    if mm:
        col_of.setdefault(mm.group(2), c)
want = {g: col_of[e] for g, e in entrez.items() if isinstance(e, str) and e in col_of}
tab = pacsv.read_csv(f, convert_options=pacsv.ConvertOptions(
    include_columns=["ModelID", "IsDefaultEntryForModel"] + sorted(set(want.values())))).to_pandas()
tab = tab[tab.IsDefaultEntryForModel == "Yes"].set_index("ModelID")
expr_all = pd.DataFrame({g: tab[c] for g, c in want.items()}).astype("float64")
expr_tr = expr_all.reindex([x for x in train_expr_lines if x in expr_all.index])
report["own_expr"] = {"n_model_genes_with_column": len(want), "n_training_lines_with_expr": len(expr_tr),
                      "n_train_has_expression_flag": len(train_expr_lines)}

# ---------------------------------------------------------------- thresholds (R2/R3) and z statistics
flat = expr_tr.to_numpy().ravel()
flat = flat[~np.isnan(flat)]
p25, p75 = np.percentile(flat, [25, 75])
st = json.loads((R / "threshold_stats.json").read_text())
gm = pd.read_parquet(R / "threshold_gene_means.parquet").set_index("gene")
own_mean, own_sd = expr_tr.mean(), expr_tr.std(ddof=1)
report["thresholds"] = {
    "own_p25": p25, "author_p25": st["rna"]["p25"], "own_p75": p75, "author_p75": st["rna"]["p75"],
    "author_n_lines": st["rna"]["n_lines"],
    "max_abs_diff_gene_mean": float((own_mean - gm.rna_mean.reindex(own_mean.index)).abs().max()),
    "max_abs_diff_gene_sd": float((own_sd - gm.rna_sd.reindex(own_sd.index)).abs().max()),
    "n_genes_mean_author": int(gm.rna_mean.notna().sum()), "n_genes_mean_own": int(own_mean.notna().sum()),
}
# would including validation/test lines change the statistics? (leakage sensitivity: they must NOT be included)
nontrain = split[(split.tier != "train") & split.has_expression.astype(bool)].ModelID
expr_all_screened = expr_all.reindex([x for x in split[split.has_expression.astype(bool)].ModelID if x in expr_all.index])
report["thresholds"]["gene_mean_diff_if_all_tiers_used_max"] = float(
    (expr_all_screened.mean() - gm.rna_mean.reindex(expr_all_screened.columns)).abs().max())

# protein statistics (author's ProCan reader; own statistics)
genes_tab = ci.gene_table(m)
prot_tr = ci.read_procan(genes_tab, train.ModelID)
assert set(prot_tr.index) <= set(train.ModelID)
pf = prot_tr.to_numpy().ravel()
pf = pf[~np.isnan(pf)]
floor = float(pf.min()) - 1.0
pz = prot_tr.where(prot_tr.notna(), floor)
report["protein"] = {
    "n_train_lines_procan": len(prot_tr), "author_n_lines": st["protein"]["n_lines"],
    "own_floor": floor, "author_floor": st["protein"]["floor_value"],
    "own_p25": float(np.percentile(pf, 25)), "author_p25": st["protein"]["p25"],
    "max_abs_diff_mean_floor": float((pz.mean() - gm.protein_mean_floor.reindex(pz.columns)).abs().max()),
    "max_abs_diff_sd_floor": float((pz.std(ddof=1) - gm.protein_sd_floor.reindex(pz.columns)).abs().max()),
}

# ---------------------------------------------------------------- panel draw (own implementation of README §5 item 8)
sizes = train.lineage.value_counts()
eligible = sorted(L for L, n in sizes.items() if n >= 15)
rng = np.random.default_rng(SEEDS[0])
chosen = [eligible[i] for i in rng.permutation(len(eligible))[:10]]
own_panel = []
for L in chosen:
    pool = sorted(train[(train.lineage == L) & train.has_expression.astype(bool)].ModelID)
    idx = sorted(rng.choice(len(pool), 3, replace=False))
    own_panel += [(pool[i], L) for i in idx]
panel = pd.read_csv(R / "stage_a_panel.tsv", sep="\t")
report["panel"] = {
    "eligible_training_lineages": eligible, "chosen": chosen,
    "own_equals_author": [x for x, _ in own_panel] == list(panel.ModelID),
    "all_train": bool(split.set_index("ModelID").loc[panel.ModelID, "tier"].eq("train").all()),
    "all_have_expression": bool(split.set_index("ModelID").loc[panel.ModelID, "has_expression"].astype(bool).all()),
    "n_procan": int(panel.has_procan.sum()), "n_cn": int(panel.has_cn.sum()),
}

# ---------------------------------------------------------------- A5 grid: own z and dedupe vs the author's mapping
acc = json.loads((R / "tmp" / "accepted_losses_panel.json").read_text())
prot_panel = ci.read_procan(genes_tab, list(panel.ModelID))
soluble = genes_tab.soluble
zs = ci.ZStats(gm.rna_mean, gm.rna_sd, gm.protein_mean_floor, gm.protein_sd_floor, st["protein"]["floor_value"])
grid = [(1.0, wp, wc, tau) for wp in (0, 0.5, 1) for wc in (0, 1) for tau in (-1.5, -1, -0.5, 0)]
assert sorted(grid) == sorted(ci.LCP_GRID)
genes_idx = pd.Index(model_genes)
mism, n_builds_own, n_builds_author, n_points, rows = 0, 0, 0, 0, []
for line in panel.ModelID:
    accepted = acc[line]["accepted"] if isinstance(acc[line], dict) else acc[line]
    # own z_RNA
    zr = ((expr_all.loc[line].reindex(genes_idx) - own_mean.reindex(genes_idx)) / own_sd.reindex(genes_idx).where(lambda s: s > 0)).fillna(0.0)
    zp = pd.Series(0.0, index=genes_idx)
    if line in prot_panel.index:
        v = prot_panel.loc[line].reindex(genes_idx).fillna(floor)
        q = genes_idx.isin(prot_panel.columns) & soluble.reindex(genes_idx).to_numpy()
        zz = ((v - pz.mean().reindex(genes_idx)) / pz.std(ddof=1).reindex(genes_idx).where(lambda s: s > 0)).fillna(0.0)
        zp = zz.where(q, 0.0)
    loss = pd.Series(np.where(genes_idx.isin(accepted), -1.0, 0.0), index=genes_idx)
    own_off, auth_hash = {}, {}
    for (wr, wp, wc, tau) in grid:
        s_own = wr * zr + wp * zp + wc * loss - tau
        off = frozenset(genes_idx[(s_own < 0).to_numpy()])
        own_off.setdefault(off, []).append((wp, wc, tau))
        s_auth = ci.lcp_scores(line, wr, wp, wc, tau, genes_tab, expr_all.loc[[line]], prot_panel, zs, accepted)
        h = hashlib.sha1(np.round(s_auth.reindex(genes_tab.index).to_numpy(), 9).tobytes()).hexdigest()
        auth_hash.setdefault(h, []).append((wp, wc, tau))
        # signed-score equality own vs author (to 1e-9)
        mism += int((s_own.reindex(genes_tab.index) - s_auth).abs().max() > 1e-9)
        n_points += 1
    # author dedup partition must be at least as fine as the off-set partition, and identical score vectors only
    part_auth = sorted(sorted(v) for v in auth_hash.values())
    part_off = sorted(sorted(v) for v in own_off.values())
    rows.append({"line": line, "has_procan": line in prot_panel.index, "n_accepted": len(accepted),
                 "author_distinct_builds": len(auth_hash), "distinct_off_sets": len(own_off),
                 "author_partition_refines_offset_partition": all(any(set(a) <= set(b) for b in part_off) for a in part_auth)})
    n_builds_author += len(auth_hash)
    n_builds_own += len(own_off)
a5 = pd.DataFrame(rows)
a5.to_csv(OUT / "a5_dedupe_check.tsv", sep="\t", index=False)
report["a5"] = {"grid_points_total": n_points, "score_vector_mismatches_own_vs_author": mism,
                "author_distinct_builds": n_builds_author, "distinct_off_sets": n_builds_own,
                "all_refine": bool(a5.author_partition_refines_offset_partition.all()),
                "procan_lines_with_no_accepted_loss": int(((a5.has_procan) & (a5.n_accepted == 0)).sum()),
                "rna_lines_with_no_accepted_loss": int(((~a5.has_procan) & (a5.n_accepted == 0)).sum())}
# smoke map consistency (2 lines)
sm = pd.read_csv(R / "stage_a_smoke_a5" / "a5_grid_map.tsv", sep="\t")
report["a5"]["smoke_map_rows"] = len(sm)
report["a5"]["smoke_distinct_builds_by_line"] = sm.groupby("line").build_id.nunique().to_dict()

# save 3-line + panel expression for the rebuild scripts (regenerable, git-ignored)
expr_all.reindex(list(panel.ModelID)).to_parquet(TMP / "panel_expr.parquet")
report["wall_s"] = round(time.time() - t0, 1)
(OUT / "r0_inputs.json").write_text(json.dumps(report, indent=1, default=float))
print(json.dumps(report, indent=1, default=float))
