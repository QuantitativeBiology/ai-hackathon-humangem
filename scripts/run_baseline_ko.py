"""Baseline single-gene knockouts of all genes in generic Human1 (v1.19.0) and Human2 (v2.1.1).

Medium: Ham's (inputs of the GR task in each release's own metabolicTasks_Essential.txt), uptake
bound 1000; objective: biomass MAR13082. Predictions only: nothing here reads DepMap values (only
CRISPRGeneEffect.csv column labels, to report gene-ID coverage).

Outputs (results/baseline/):
  {model}_gene_mapping.tsv   Ensembl -> Entrez/symbol for every model gene (nextgem.genes.gene_mapping)
  {model}_ko.parquet         gene, entrez, symbol, wt_growth, ko_growth, ratio, ratio_r, status, solver,
                             runtime, map_status (runtime = seconds per KO solve). ratio is raw (solver noise
                             ~1e-13 varies between runs); ratio_r = ratio rounded to 6 decimals, clipped to
                             [0, 1]: ALL downstream scoring uses ratio_r.
  summary.tsv                per-model counts and wall-clock runtimes
  model_differences.txt      Human1 vs Human2: genes, medium exchanges, biomass reaction

The run fails (no outputs for that model) if any KO status is not optimal/infeasible or any ratio is NaN.

Usage: .venv/bin/python scripts/run_baseline_ko.py [--processes 15] [--models human1 human2] [--out DIR]
(--processes counts worker processes; the parent adds one more, so 15 keeps the run at 16
processes / 16 Gurobi WLS sessions.)
"""

import argparse
import time
from pathlib import Path

import pandas as pd

from nextgem import BIOMASS_REACTION, GEMS, hams_medium, load_model, single_gene_knockouts
from nextgem.genes import gene_mapping, mapping_summary
from nextgem.model import PROJECT_ROOT

OUT_DEFAULT = PROJECT_ROOT / "results" / "baseline"
DEPMAP_CRISPR = PROJECT_ROOT / "data" / "depmap" / "CRISPRGeneEffect.csv"
THRESHOLDS = (0.01, 0.5)


def genes_closing_no_reaction(model) -> int:
    """Genes whose KO closes no (not already closed) reaction: ratio 1 by construction."""
    n = 0
    for gene in model.genes:
        before = {r.id: r.bounds for r in gene.reactions}
        with model:
            gene.knock_out()
            n += not any(r.bounds == (0, 0) != before[r.id] for r in gene.reactions)
    return n


def run_model(name: str, processes: int, depmap_columns, out: Path) -> tuple[dict, dict]:
    gem = GEMS[name]
    t0 = time.perf_counter()
    model = load_model(gem.sbml)
    t_load = time.perf_counter() - t0

    mapping = gene_mapping(model, gem.genes_tsv, depmap_columns)
    mapping.to_csv(out / f"{name}_gene_mapping.tsv", sep="\t", index=False)

    t1 = time.perf_counter()
    ko = single_gene_knockouts(model, processes=processes)  # strict: raises on non-optimal statuses
    t_ko = time.perf_counter() - t1
    if not ko.status.isin(["optimal", "infeasible"]).all() or ko.ratio_r.isna().any():
        raise RuntimeError(f"{name}: unexpected statuses {ko.status.value_counts().to_dict()} or NaN ratios")

    ko = ko.merge(mapping[["gene", "entrez", "symbol", "map_status"]], on="gene", how="left", validate="1:1")
    cols = ["gene", "entrez", "symbol", "wt_growth", "ko_growth", "ratio", "ratio_r", "status", "solver", "runtime",
            "map_status"]
    ko = ko[cols].sort_values("gene").reset_index(drop=True)
    ko.to_parquet(out / f"{name}_ko.parquet", index=False)

    summary = {
        "model": name,
        "version": gem.version,
        "genes": len(ko),
        "wt_growth": ko.wt_growth.iloc[0],
        "solver": ko.solver.iloc[0],
        **{f"lethal_ratio_r_lt_{t}": int((ko.ratio_r < t).sum()) for t in THRESHOLDS},
        "partial_ratio_r_0.01_to_1": int(((ko.ratio_r >= 0.01) & (ko.ratio_r < 1)).sum()),
        "ratio_r_eq_1": int((ko.ratio_r == 1).sum()),
        "ratio_r_distinct": int(ko.ratio_r.nunique()),
        "ratio_nan": int(ko.ratio.isna().sum()),
        "genes_closing_no_reaction": genes_closing_no_reaction(model),
        "status_counts": ";".join(f"{k}={v}" for k, v in ko.status.value_counts().items()),
        "medium_exchanges": len(hams_medium(model)),
        **{f"map_{k}": v for k, v in mapping_summary(mapping).items()},
        "processes": processes,
        "load_s": round(t_load, 1),
        "knockouts_wall_s": round(t_ko, 1),
        "ko_solve_cpu_s": round(ko.runtime.sum(), 1),
    }
    print(summary, flush=True)
    biomass = model.reactions.get_by_id(BIOMASS_REACTION)
    facts = {
        "genes": {g.id for g in model.genes},
        "medium": set(hams_medium(model)),
        "medium_names": {r: next(iter(model.reactions.get_by_id(r).metabolites)).name for r in hams_medium(model)},
        "biomass": {m.id: c for m, c in biomass.metabolites.items()},
        "biomass_name": biomass.name,
    }
    return summary, facts


def describe_differences(facts: dict[str, dict]) -> str:
    h1, h2 = facts["human1"], facts["human2"]
    lines = [
        f"Biomass {BIOMASS_REACTION}: Human1 '{h1['biomass_name']}' ({len(h1['biomass'])} species), "
        f"Human2 '{h2['biomass_name']}' ({len(h2['biomass'])} species)",
        f"  species only in Human1: {sorted(set(h1['biomass']) - set(h2['biomass']))}",
        f"  species only in Human2: {sorted(set(h2['biomass']) - set(h1['biomass']))}",
        "  coefficient differences (Human1 vs Human2): "
        + str({m: (h1["biomass"][m], h2["biomass"][m]) for m in sorted(set(h1["biomass"]) & set(h2["biomass"]))
               if abs(h1["biomass"][m] - h2["biomass"][m]) > 1e-9}),
        f"Ham's medium exchanges: Human1 {len(h1['medium'])}, Human2 {len(h2['medium'])}, shared ids {len(h1['medium'] & h2['medium'])}",
        f"  only in Human1: {sorted((r, h1['medium_names'][r]) for r in h1['medium'] - h2['medium'])}",
        f"  only in Human2: {sorted((r, h2['medium_names'][r]) for r in h2['medium'] - h1['medium'])}",
        f"Genes: Human1 {len(h1['genes'])}, Human2 {len(h2['genes'])}, shared {len(h1['genes'] & h2['genes'])}, "
        f"only Human1 {len(h1['genes'] - h2['genes'])}, only Human2 {len(h2['genes'] - h1['genes'])}",
    ]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--processes", type=int, default=15)
    parser.add_argument("--models", nargs="+", default=["human1", "human2"], choices=sorted(GEMS))
    parser.add_argument("--out", type=Path, default=OUT_DEFAULT, help="output directory")
    args = parser.parse_args()

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    depmap_columns = pd.read_csv(DEPMAP_CRISPR, nrows=0).columns[1:]  # labels only, no values
    t0 = time.perf_counter()
    summaries, facts = [], {}
    for name in args.models:
        summary, facts[name] = run_model(name, args.processes, depmap_columns, out)
        summaries.append(summary)
    df = pd.DataFrame(summaries)
    df["total_wall_s"] = round(time.perf_counter() - t0, 1)
    df.to_csv(out / "summary.tsv", sep="\t", index=False)
    if {"human1", "human2"} <= set(facts):
        text = describe_differences(facts)
        (out / "model_differences.txt").write_text(text)
        print(text)
    print(df.T.to_string())


if __name__ == "__main__":
    main()
