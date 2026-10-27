# Files the user needs to download manually

These are blocked for scripts (Cloudflare, JS pages, paywalls). Put each file in the folder shown, keeping the original file name. Agents will then checksum them and add them to DATA.md. Papers go in `papers/`, which git ignores, so nothing copyrighted is committed.

## Priority 1: needed for the core design

| # | What | Where to get it | Save to |
|---|---|---|---|
| 1 | DepMap 26Q3 `CRISPRGeneDependency.csv` | depmap.org/portal/data_page → DepMap Public 26Q3 | `data/depmap/` |
| 2 | DepMap 26Q3 `AchillesNonessentialControls.csv` (or `CRISPRInferredNonessentials`, whichever exists) | same | `data/depmap/` |
| 3 | DepMap 26Q3 proteomics, if the release has any (e.g. harmonised MS proteomics) | same | `data/depmap/` |
| 4 | DepMap 26Q3 media formulation table, if one exists (maps `MF-xxx-xxx` IDs to components) | same | `data/depmap/` |
| 5 | Gonçalves et al. 2022 ProCan-DepMapSanger protein matrix (949 lines), processed and averaged | cellmodelpassports.sanger.ac.uk/downloads → Proteomics | `data/proteomics_goncalves/` |
| 6 | HPLM composition: Cantor et al. 2017 *Cell*, Table S1 | PMC5421364, supplementary files | `data/media/` |
| 7 | Plasmax composition: Vande Voorde et al. 2019 *Sci Adv*, supplementary materials | doi:10.1126/sciadv.aau7314 | `data/media/` |
| 8 | Formulation sheets (PDF or CSV) for RPMI 1640, DMEM high glucose, DMEM/F12, MEM/EMEM, IMDM, Ham's F12, McCoy's 5A, Leibovitz L-15 | Thermo Fisher (Gibco) product pages → "Formulation" | `data/media/` |

## Priority 2: useful, not blocking

| # | What | Where to get it | Save to |
|---|---|---|---|
| 9 | PRISM Repurposing, latest on the DepMap portal (newer than 24Q2, if it exists) | depmap.org data page | `data/prism/` |
| 10 | CTRPv2 AUC and compound metadata | DepMap portal or the Broad CTRP portal | `data/ctrp/` |
| 11 | De Kegel 2021 *Cell Syst* and Ito 2021 *Nat Genet* supplementary paralog tables | doi:10.1016/j.cels.2021.08.006; doi:10.1038/s41588-021-00967-z | `data/paralogs/` |

## Papers (PDF plus supplements) → `papers/`
- Luo et al. 2026 *PNAS* (Human2): SI Appendix and Datasets S01–S09, doi:10.1073/pnas.2516511123
- Robinson et al. 2020 *Sci Signal* (Human1): supplementary Figs S8–S11 and Table S3, doi:10.1126/scisignal.aaz1482
- Behan et al. 2019 *Nature* (Project Score), doi:10.1038/s41586-019-1103-9
- Hart et al. 2015 *Cell*, doi:10.1016/j.cell.2015.11.015
- Tao et al. 2025 *Cell Rep* (DeepMeta), doi:10.1016/j.celrep.2025.115945
- Ito et al. 2021 *Nat Genet*; De Kegel et al. 2021 *Cell Syst* (main text)
- Optional: GECKO 3.0 (doi:10.1038/s41596-023-00931-7), tINIT (doi:10.1002/msb.145122), iMAT (doi:10.1038/nbt.1487)

## W5 (GPR curation), added 2026-10-27 by the W5 engineer — not downloaded by any agent

| # | What | Where to get it | Save to |
|---|---|---|---|
| 12 | UniProtKB/Swiss-Prot human reviewed proteome as a versioned flat file or TSV with columns: Entry, Entry version, Gene Names, Function [CC], Catalytic activity, EC number, Subcellular location [CC], Transmembrane, Intramembrane, Rhea ID (e.g. the REST "stream" TSV for `reviewed:true AND organism_id:9606`, plus the release number) | https://www.uniprot.org/uniprotkb?query=reviewed:true+AND+organism_id:9606 (Download → TSV with the columns above), or ftp.uniprot.org `knowledgebase/complete/uniprot_sprot.dat.gz` filtered to human | `data/uniprot/` (record the release in DATA.md) |

Why: DESIGN W5 "Remove OR branch" has two paths. Path 1 ("UniProt function unrelated to the reaction", with accession
and entry version cited) needs this file; without it W5 implements only path 2 ("pan-cancer absent"). The same file
would replace W5's label-free *soluble-enzyme proxy* (Human-GEM genes.tsv location + no transport reaction) by the
UniProt Transmembrane/Subcellular-location annotation.
