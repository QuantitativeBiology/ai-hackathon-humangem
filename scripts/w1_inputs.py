"""W1 label-free inputs (DESIGN.md §4.3, §5.4, §8). Reads no label-bearing file.

Writes (git-ignored, results/w1/tmp/):
- ``expr.parquet``    log2(TPM+1), default entry per model, screened lines x (U ∪ model genes ∪ their paralogs),
                      columns = Entrez IDs (float32)
- ``cn.parquet``      relative WGS CN, default entry, screened lines x U genes (float32)
- ``mut_dam.parquet`` damaging-allele count (0/1/2), default entry, screened lines x U genes (float32)
- ``mut_hot.parquet`` hotspot flag, screened lines x U genes present in the hotspot matrix (float32)
Writes (committed, results/w1/):
- ``gene_features.tsv``  per U gene: in_human1/2, addressable_h1/h2 (KO closes >= 1 reaction, §5.4),
  n_isozyme_genes_h2 (genes in GPR clauses without the gene, over its H2 reactions), neighbour and
  paralog lists are in tmp/gene_neighbours.json (H2 shared non-currency metabolite; Ensembl 116 paralogs,
  either %id >= 30)
- ``line_features.tsv``  per screened line: tier, lineage, base_medium_class, serum flag (formulation text),
  has_expression, has_cn, library_class, model_type

Usage: .venv/bin/python scripts/w1_inputs.py
"""

from __future__ import annotations

import ast
import json
import re
import sys
import time
from collections import defaultdict
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w1_common import BASELINE, DEPMAP, OUT_DIR, PROJECT, TMP_DIR, atomic_write_parquet, atomic_write_tsv, current, run_main  # noqa: E402

DM = PROJECT / "data" / "depmap"
GENE_COL = re.compile(r"^(?P<symbol>.+) \((?P<entrez>\d+)\)$")
# Currency metabolites (compartment-merged names) never define an enzyme neighbourhood. Fixed list plus
# a degree rule: any metabolite name in >= 100 gene-associated reactions of Human2.
CURRENCY = {"H2O", "H+", "ATP", "ADP", "AMP", "NAD+", "NADH", "NADP+", "NADPH", "Pi", "PPi", "CO2", "O2",
            "CoA", "FAD", "FADH2", "NH3", "H2O2", "GTP", "GDP", "UTP", "UDP", "CTP", "CDP", "HCO3-",
            "Na+", "K+", "Cl-", "acetyl-CoA", "glutamate", "2-oxoglutarate", "ubiquinol", "ubiquinone",
            "Fe2+", "Fe3+", "SAM", "SAH", "sodium", "electron"}
NEIGHBOUR_DEGREE_CAP = 100
PARALOG_MIN_ID = 30.0


def omics_matrix(fname: str, entrez_keep: set[int], lines: list[str]) -> pd.DataFrame:
    path = DM / fname
    head = list(pd.read_csv(path, nrows=0).columns)
    gcols = {c: int(m["entrez"]) for c in head if (m := GENE_COL.match(c))}
    keep = [c for c, e in gcols.items() if e in entrez_keep]
    meta = [c for c in ("ModelID", "IsDefaultEntryForModel") if c in head]
    df = pd.read_csv(path, usecols=meta + keep, dtype={c: "float32" for c in keep})
    df = df[df.IsDefaultEntryForModel.astype(str).str.lower().isin(["yes", "true"])]
    if df.ModelID.duplicated().any():
        raise ValueError(f"{fname}: several default rows for one model")
    df = df.set_index("ModelID").drop(columns="IsDefaultEntryForModel")
    df.columns = [gcols[c] for c in df.columns]
    df = df.loc[:, ~pd.Index(df.columns).duplicated()]
    df = df[df.index.isin(lines)].sort_index()
    df.columns = df.columns.astype(str)
    return df


def dnf(node) -> list[frozenset[str]]:
    """GPR ast -> list of AND-clauses (sets of gene ids)."""
    if isinstance(node, ast.Name):
        return [frozenset([node.id])]
    if isinstance(node, ast.BoolOp):
        parts = [dnf(v) for v in node.values]
        if isinstance(node.op, ast.Or):
            return [c for p in parts for c in p]
        out = [frozenset()]
        for p in parts:
            out = [a | b for a, b in product(out, p)]
            if len(out) > 20000:
                raise ValueError("DNF too large")
        return out
    raise TypeError(type(node))


def model_gene_features(name: str):
    import cobra

    from nextgem.model import load_model

    cobra.Configuration().solver = "glpk"  # structure only: no Gurobi session (cap shared with w1_p57.py)
    m = load_model(name)
    disabled = defaultdict(int)
    isozymes: dict[str, set[str]] = defaultdict(set)
    for r in m.reactions:
        gids = [g.id for g in r.genes]
        if not gids:
            continue
        clauses = dnf(r.gpr.body) if r.gpr.body is not None else []
        for g in gids:
            if not r.gpr.eval([g]):
                disabled[g] += 1
            isozymes[g] |= set().union(*[c for c in clauses if g not in c]) if clauses else set()
    return m, disabled, isozymes


def main() -> int:
    t0, c0 = time.time(), time.process_time()
    from nextgem import labels

    split = labels.load_split()
    lines = sorted(split.ModelID)
    U = labels.gene_universe()
    u_set = set(U.entrez.astype(int))

    # ---- model genes, addressable sets, isozymes, neighbours
    gmap = {}
    for name in ("human1", "human2"):
        mp = pd.read_csv(PROJECT / "results" / "baseline" / f"{name}_gene_mapping.tsv", sep="\t", dtype=str)
        gmap[name] = mp.dropna(subset=["entrez"]).assign(entrez=lambda d: d.entrez.astype(int)).set_index("gene").entrez
    feats = U.copy()
    neighbours: dict[int, set[int]] = defaultdict(set)
    h2_model = None
    for name in ("human1", "human2"):
        m, disabled, iso = model_gene_features(name)
        ens2ez = gmap[name]
        addr = {ens2ez[g] for g, n in disabled.items() if n > 0 and g in ens2ez.index}
        feats[f"addressable_{name}"] = feats.entrez.isin(addr)
        if name == "human2":
            h2_model = m
            iso_ez = {ens2ez[g]: {ens2ez[x] for x in s if x in ens2ez.index} - {ens2ez[g]}
                      for g, s in iso.items() if g in ens2ez.index}
            feats["n_isozyme_genes_h2"] = feats.entrez.map(lambda e: len(iso_ez.get(e, ())) if e in iso_ez else np.nan)
    # neighbours (H2): genes sharing a non-currency metabolite (compartment-merged name)
    met_rxn = defaultdict(set)
    for r in h2_model.reactions:
        if r.genes and not r.boundary:
            for mt in r.metabolites:
                met_rxn[mt.name].add(r.id)
    currency = CURRENCY | {k for k, v in met_rxn.items() if len(v) >= NEIGHBOUR_DEGREE_CAP}
    ens2ez = gmap["human2"]
    met_genes = defaultdict(set)
    for r in h2_model.reactions:
        if not r.genes or r.boundary:
            continue
        ez = {ens2ez[g.id] for g in r.genes if g.id in ens2ez.index}
        for mt in r.metabolites:
            if mt.name not in currency:
                met_genes[mt.name] |= ez
    for genes in met_genes.values():
        for e in genes:
            neighbours[int(e)] |= {int(x) for x in genes} - {int(e)}
    # paralogs (Ensembl 116 within-species, either %id >= 30)
    emap = pd.read_csv(PROJECT / "data" / "paralogs" / "ensembl_gene_map.tsv", sep="\t", dtype=str)
    emap = emap.dropna(subset=["NCBI gene (formerly Entrezgene) ID"])
    ens2ez_all = emap.drop_duplicates("Gene stable ID").set_index("Gene stable ID")["NCBI gene (formerly Entrezgene) ID"].astype(int)
    ez2ens = emap.assign(e=lambda d: d["NCBI gene (formerly Entrezgene) ID"].astype(int)).groupby("e")["Gene stable ID"].agg(set)
    u_ens = set().union(*[ez2ens[e] for e in u_set if e in ez2ens.index])
    par = pd.read_csv(PROJECT / "data" / "paralogs" / "ensembl_human_paralogs.tsv", sep="\t",
                      usecols=[0, 1, 3, 4], names=["a", "b", "id_ab", "id_ba"], header=0,
                      dtype={"a": str, "b": str, "id_ab": "float32", "id_ba": "float32"})
    par = par[par.a.isin(u_ens) & ((par.id_ab >= PARALOG_MIN_ID) | (par.id_ba >= PARALOG_MIN_ID))]
    par = par[par.a.isin(ens2ez_all.index) & par.b.isin(ens2ez_all.index)]
    paralogs: dict[int, set[int]] = defaultdict(set)
    for a, b in zip(par.a.map(ens2ez_all), par.b.map(ens2ez_all)):
        if a != b:
            paralogs[int(a)].add(int(b))
    feats["n_paralogs"] = feats.entrez.map(lambda e: len(paralogs.get(e, ())))
    feats["n_neighbours_h2"] = feats.entrez.map(lambda e: len(neighbours.get(e, ())))
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    (TMP_DIR / "gene_neighbours.json").write_text(json.dumps({
        "neighbours_h2": {str(k): sorted(int(x) for x in v) for k, v in neighbours.items() if k in u_set},
        "paralogs": {str(k): sorted(int(x) for x in v) for k, v in paralogs.items() if k in u_set},
        "currency_n": len(currency)}))

    # ---- omics (default entries; screened lines)
    model_genes = set(gmap["human1"].values) | set(gmap["human2"].values)
    expr_keep = u_set | model_genes | {e for v in paralogs.values() for e in v} | {e for k, v in neighbours.items() for e in v}
    expr = omics_matrix("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", expr_keep, lines)
    cn = omics_matrix("OmicsCNGeneWGS.csv", u_set, lines)
    dam = omics_matrix("OmicsSomaticMutationsMatrixDamaging.csv", u_set, lines)
    hot = omics_matrix("OmicsSomaticMutationsMatrixHotspot.csv", u_set, lines)
    for nm, df in (("expr", expr), ("cn", cn), ("mut_dam", dam), ("mut_hot", hot)):
        atomic_write_parquet(df.reset_index(), TMP_DIR / f"{nm}.parquet")
    feats["has_expression_column"] = feats.entrez.astype(str).isin(expr.columns)

    # ---- line features (medium class from the split; serum flag from the formulation text)
    ssm = pd.read_csv(DM / "ScreenSequenceMap.csv")
    smap = pd.read_csv(DM / "CRISPRScreenMap.csv")
    ssm = ssm[(ssm.PassesQC == True) & (ssm.ExcludeFromCRISPRCombined == False) & (ssm.ScreenType != "pDNA")  # noqa: E712
              & ssm.ScreenID.isin(smap.ScreenID)]
    mc = pd.read_csv(DM / "ModelCondition.csv", usecols=["ModelConditionID", "FormulationID"])
    model = pd.read_csv(DM / "Model.csv", usecols=["ModelID", "FormulationID"]).set_index("ModelID")
    txt = ssm.merge(mc, on="ModelConditionID").dropna(subset=["FormulationID"]).groupby("ModelID").FormulationID.agg(lambda s: s.value_counts().index[0])
    text = pd.Series({m: txt.get(m, model.FormulationID.get(m, np.nan)) for m in lines})

    def serum(t):
        if not isinstance(t, str):
            return np.nan
        tl = t.lower()
        if "serum-free" in tl or "serum free" in tl:
            return 0.0
        return 1.0 if re.search(r"\bfbs\b|serum|\bfcs\b|\bhs\b", tl) else 0.0

    lf = split[["ModelID", "tier", "stratum", "lineage", "library_class", "model_type", "medium_source",
                "base_medium_class", "has_expression", "has_cn"]].copy()
    lf["serum"] = lf.ModelID.map(text.map(serum))
    lf["has_expression_matrix"] = lf.ModelID.isin(expr.index)
    atomic_write_tsv(lf, OUT_DIR / "line_features.tsv")
    atomic_write_tsv(feats, OUT_DIR / "gene_features.tsv")

    summary = {"n_lines": len(lines), "expr_lines": int(len(expr)), "expr_genes": int(expr.shape[1]),
               "cn_lines": int(len(cn)), "U": len(u_set), "U_with_expression": int(feats.has_expression_column.sum()),
               "addressable_human1": int(feats.addressable_human1.sum()), "addressable_human2": int(feats.addressable_human2.sum()),
               "serum_flag_counts": lf.serum.value_counts(dropna=False).to_dict(),
               "currency_metabolites": len(currency), "runtime_s": round(time.time() - t0, 1)}
    summary["serum_flag_counts"] = {str(k): int(v) for k, v in summary["serum_flag_counts"].items()}
    (OUT_DIR / "inputs_summary.json").write_text(json.dumps(summary, indent=1))
    current().output(OUT_DIR / "line_features.tsv", OUT_DIR / "gene_features.tsv", OUT_DIR / "inputs_summary.json",
                     TMP_DIR / "gene_neighbours.json", *[TMP_DIR / f"{k}.parquet" for k in ("expr", "cn", "mut_dam", "mut_hot")])
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(run_main(main, inputs=[DEPMAP / f for f in (
        "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", "OmicsCNGeneWGS.csv", "OmicsSomaticMutationsMatrixDamaging.csv",
        "OmicsSomaticMutationsMatrixHotspot.csv", "ScreenSequenceMap.csv", "CRISPRScreenMap.csv", "ModelCondition.csv",
        "Model.csv")] + [BASELINE / "human1_gene_mapping.tsv", BASELINE / "human2_gene_mapping.tsv",
        PROJECT / "data" / "paralogs" / "ensembl_gene_map.tsv", PROJECT / "data" / "paralogs" / "ensembl_human_paralogs.tsv",
        PROJECT / "data" / "Human1" / "model" / "Human-GEM.xml", PROJECT / "data" / "Human-GEM-main" / "model" / "Human-GEM.xml",
        PROJECT / "results" / "splits" / "lineage_split.tsv"], solver="glpk (structure only)"))
