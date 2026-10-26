# Data provenance

Contents of `data/` are not versioned (only this file is). Re-download from the sources below to reproduce.

## Human-GEM — `Human-GEM-main/`
- Source: https://github.com/SysBioChalmers/Human-GEM (main branch snapshot)
- Version: **2.1.1** (`version.txt`; the bundled `pyproject.toml` still says 2.1.0)
- Obtained: 2026-10-26
- Model used: `model/Human-GEM.xml` (12,748 reactions, 8,374 metabolites, 2,870 genes; Ensembl gene IDs, Entrez/symbols in `model/genes.tsv`)
- The repo's Python package `human_gem` is installed editable into `.venv` (creates an ignored `human_gem.egg-info/`).

## DepMap — `depmap/`
- Source: DepMap Portal, https://depmap.org/portal/data_page/
- Release: **DepMap Public 26Q3**
- Obtained: 2026-10-26

| File | Shape (rows × cols, excl. header/index) | Notes |
|---|---|---|
| `CRISPRGeneEffect.csv` | 1,361 models × 18,434 genes | Chronos gene effect, one row per `ModelID` |
| `ScreenGeneEffect.csv` | 1,611 screens × 18,434 genes | Screen-level (map via `CRISPRScreenMap.csv`) |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 2,314 rows × 19,218 genes | Multiple rows per model; filter `IsDefaultEntryForModel` |
| `OmicsCNGeneWGS.csv` | 1,846 rows × 19,344 genes | |
| `OmicsSomaticMutations.csv` | 1,392,076 variants | Long format |
| `OmicsSomaticMutationsMatrixDamaging.csv` / `...Hotspot.csv` | 3,807 rows | |
| `Model.csv` | 2,554 models | Cell line annotation (lineage, OncoTree, etc.) |
| `AchillesCommonEssentialControls.csv`, `CRISPRInferredCommonEssentials.csv` | gene lists | |

Gene columns are labelled `SYMBOL (EntrezID)`; join to Human-GEM on Entrez ID (screened-gene counts per model: `results/design/design_numbers.json` → `genes`).

## Additional datasets (fetched 2026-10-26 by `scripts/download_<name>.py`)
- `proteomics_nusinow/`: `protein_quant_current_normalized.csv.gz` (12,755 protein rows x ~420 TMT channels; columns `Protein_Id, Gene_Symbol, ..., <CellLine_Tissue_TMTplex>`; includes bridge/replicate channels), `Table_S1_Sample_Information.xlsx` (channel -> CCLE name), `Table_S4_...xlsx`. Names are CCLE names, map to DepMap `Model.csv:CCLEName`.
- `metabolomics_ccle/`: `CCLE_metabolomics_20190502.csv` (928 rows = 927 unique lines x 225 metabolites, has `DepMap_ID`), `Cell_lines_annotations_20181226.txt`.
- `prism/repurposing_24Q2/`: single-dose PRISM Repurposing Public 24Q2 (LFC collapsed per compound, `Extended_Primary_Data_Matrix`, treatment and cell-line metadata, README). `prism/repurposing_20Q2_secondary/`: multi-dose curve parameters (AUC, IC50) for the secondary screen. Cell lines are keyed by DepMap ID (24Q2 `Cell_Line_Meta_Data` uses `depmap_id`).
- `gdsc/`: GDSC1 and GDSC2 fitted dose-response (IC50, AUC, Z-score; xlsx), screened compounds, cell line details (COSMIC IDs).
- `complexes/`: `complexportal_9606.tsv` (human complexes, UniProt members), `corum_humanComplexes.txt` + `corum_uniprotCorumMapping.txt` (CORUM 5.1).
- `paralogs/`: `ensembl_human_paralogs.tsv` (gene pairs with %identity both directions, last common ancestor `subtype`, homology type; 3.55 M rows, all biotypes incl. miRNA, each pair appears in both directions), `ensembl_gene_map.tsv` (Ensembl ID, symbol, Entrez ID, biotype, chromosome), `RELEASE.txt`.
- `cellmodelpassports/`: Sanger model list (ids COSMIC/BROAD/CCLE/RRID, growth properties; no medium field).
