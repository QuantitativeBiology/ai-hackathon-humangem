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

Gene columns are labelled `SYMBOL (EntrezID)`; join to Human-GEM on Entrez ID (2,798 / 2,868 GEM genes with Entrez IDs are in `CRISPRGeneEffect.csv`).
