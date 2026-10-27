"""W4 preparation (label-free): Stage-A panel, threshold statistics, FASTCC consistent set.

Outputs (results/w4/):
  stage_a_panel.tsv            30 training lines (10 lineages x 3), SEEDS[0]; ModelID, lineage, has_cn, has_procan
  threshold_stats.json         R2/R3 statistics per modality (global P25/P75, n training lines)
  threshold_gene_means.parquet gene, rna_mean, protein_mean (R3 inputs); rna_sd, protein_mean_floor,
                               protein_sd_floor (A5 z-scores); training lines only
  fastcc_human2_M0.json        cobra fastcc consistent reactions of Human2 on M0 (FASTCORE input)
  a6prime_lines.tsv            amendment A7: one representative panel line per panel lineage for CORDA (A6')

Reads: the split file (tier and label-free covariates), DepMap expression, ProCan proteomics,
Human-GEM. Opens no label file.

Usage: .venv/bin/python scripts/w4_prepare.py [--skip-fastcc] [--only-a6prime]
"""

import argparse
import json
import resource
import time

import pandas as pd

from nextgem import paths
from nextgem.context import consistent_reactions, make_parent
from nextgem.context_inputs import (A6P_TAG, a6prime_lines, gene_table, read_expression, read_procan, stage_a_panel,
                                    threshold_stats)
from nextgem.model import load_model

OUT = paths.results_dir() / "w4"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-fastcc", action="store_true")
    ap.add_argument("--only-a6prime", action="store_true", help="only (re)write a6prime_lines.tsv")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0, c0 = time.time(), time.process_time()
    split = pd.read_csv(paths.splits_dir() / "lineage_split.tsv", sep="\t")
    panel = stage_a_panel(split)
    if a.only_a6prime:
        import cobra

        cobra.Configuration().solver = "glpk"
        genes = gene_table(load_model("human2"))
        rna = read_expression(genes, split[(split.tier == "train") & split.has_expression.astype(bool)].ModelID)
        a6p = a6prime_lines(panel, rna, split)
        a6p["status"] = A6P_TAG
        a6p.to_csv(OUT / "a6prime_lines.tsv", sep="\t", index=False)
        print(a6p[["lineage", "ModelID", "distance", "n_lineage_lines"]].to_string())
        return
    panel.to_csv(OUT / "stage_a_panel.tsv", sep="\t", index=False)
    print(panel.groupby("lineage").size().to_dict(), "procan:", int(panel.has_procan.sum()), "cn:", int(panel.has_cn.sum()))

    import cobra

    cobra.Configuration().solver = "glpk"  # structure only: no WLS session outside the pool (review B1)
    model = load_model("human2")
    genes = gene_table(model)
    train = split[(split.tier == "train")]
    rna = read_expression(genes, train[train.has_expression.astype(bool)].ModelID)
    prot = read_procan(genes, train.ModelID)
    st_rna, st_prot = threshold_stats(rna, "rna"), threshold_stats(prot, "protein")
    stats = {m: {"p25": s.p25, "p75": s.p75, "n_lines": s.n_lines, "floor_value": s.floor_value,
                 "n_genes": int(s.gene_mean.notna().sum())} for m, s in (("rna", st_rna), ("protein", st_prot))}
    stats["rna"]["R1_log2"] = 1.0
    stats["rna"]["share_values_below_R1"] = float((rna.to_numpy() < 1.0).mean())
    (OUT / "threshold_stats.json").write_text(json.dumps(stats, indent=1))
    prot_z = prot.where(prot.notna(), st_prot.floor_value)  # A5 z_prot: not detected = floor value (README §5 item 7)
    gm = pd.DataFrame({"rna_mean": st_rna.gene_mean, "protein_mean": st_prot.gene_mean,
                       "rna_sd": rna.std(axis=0, ddof=1), "protein_sd_floor": prot_z.std(axis=0, ddof=1),
                       "protein_mean_floor": prot_z.mean(axis=0)}).rename_axis("gene").reset_index()
    gm.to_parquet(OUT / "threshold_gene_means.parquet", index=False)
    print(json.dumps(stats, indent=1))
    a6p = a6prime_lines(panel, rna, split)
    a6p["status"] = A6P_TAG
    a6p.to_csv(OUT / "a6prime_lines.tsv", sep="\t", index=False)
    print(a6p[["lineage", "ModelID", "distance", "n_lineage_lines"]].to_string())

    if not a.skip_fastcc:  # FASTCC on Gurobi (the F3 solver), inside a slot of the shared pool
        from nextgem import gurobi_pool

        _slot = gurobi_pool.hold_slot()
        cobra.Configuration().solver = "gurobi"
        model = load_model("human2")
        parent = make_parent(model)
        tf = time.process_time()
        cons = consistent_reactions(parent, OUT / "fastcc_human2_M0.json")
        print(f"fastcc: {len(cons)} / {len(parent.rxn_ids)} consistent; {time.process_time() - tf:.0f} CPU-s")
    print(f"wall {time.time() - t0:.0f} s, CPU {time.process_time() - c0:.0f} s, peak RSS "
          f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB")


if __name__ == "__main__":
    main()
