# Design review r1: data correctness (DepMap/CCLE/Sanger angle)

Reviewer: cancer-genomics data reviewer. All numbers below were computed from files under `data/` (read only). Only label *coverage/prevalence* was inspected (ASNS prevalence by medium, per-line positive rate by library); no model-vs-label performance was computed.

## Verdict
The label join (ModelID, Entrez), the split inputs and the lineage counts are sound. The medium and omics-coverage assumptions are not: 118 Sanger-only lines have no screen medium, and the Sanger-only lines are mostly missing expression and CN. Proteomics and metabolomics overlap with CRISPR is smaller than the headline numbers. Several of these change what W2, W4 and W6 can claim. The must-fix items are cheap to put into the design.

## Verified claims (OK)
- 1,361 ModelIDs in `CRISPRGeneEffect`; 18,434 gene columns, all `SYMBOL (Entrez)`, no duplicate or missing Entrez. `ScreenGeneEffect` has identical columns and 1,611 screens; 1,604 map one-to-one into `CRISPRScreenMap` (all 1,361 models covered). 7 screens (`*.BN*`, Brunello) are in ScreenGeneEffect but not in the map. Using the ModelID-level matrix is right.
- 23 lineages have >=15 screened lines = 1,331 lines; rare lines = 30 (11+9+5+2+1+1+1). Correct. Strata S1/S2/S3 = 8/8/7 lineages.
- Expression, CN, damaging and hotspot matrices all carry `IsDefaultEntryForModel`; one default row per model (expr 2,130; CN 1,763; mutation matrices 2,395).
- Model.csv and ModelCondition.csv: the swap is real (`OnboardedMedia` = MF-xxx id, `FormulationID` = text). ModelCondition has the same swap (`GrowthMedia` = MF id, `FormulationID` = text). OmicsProfiles.GrowthMedia is an MF id.
- Resource scope ~2,100: 2,130 default-expression models minus 42 Fibroblast and a few non-cancer = ~2,085. OK.
- Human-GEM v1.19.0 and v2.1.1 are CC BY 4.0 (LICENSE.md present in both trees).

## Numeric claims that are off (fix the text)
| DESIGN.md | Claim | Actual |
|---|---|---|
| §2.3 | H1 2,825, H2 2,807, overlap ~2,780, union ~2,850 | 2,825 / 2,807 / **2,766** / **2,866** (primary Entrez). DATA.md/README still quote 2,798/2,868 (raw Entrez string, includes the 11 one-to-many genes); make them agree with the baseline README. |
| §2.2/§4 | "Nusinow 375 lines" | 375 map to ModelID, but only **309** are CRISPR-screened (308 with default expression; 293 with CN). |
| §4.1 (c) | test tier >=25 Nusinow lines | Feasible (acceptance in my dry-run of the draw, seeds not the real one: 89% of attempts pass all of a-d; (b)+(c) alone 98%). Proteomics lines are concentrated: Lung 62, Eso/Stomach 25, Bowel 24, Skin 24, Breast 24; **0** in Cervix, Eye; 1 in PNS and Biliary; 4 in Soft Tissue. |
| DATA.md | metabolomics 928 lines | 927 non-null `DepMap_ID`; **652** screened. |

## Coverage of the screened lines (n = 1,361)
| Input | Lines | Notes |
|---|---|---|
| Default expression | 1,292 (69 missing) | 62 of the 69 missing are Sanger-only (KY) lines |
| CN (WGS) | 1,221 (140 missing) | |
| Damaging/hotspot mutation matrices | 1,361 | |
| expr & CN & mutation | 1,214 | |
| Nusinow | 309 | 23% of lines |
| CCLE metabolomics | 652 | 651 with expression |
| PRISM (24Q2 or 20Q2) | 735 | 733 with expression |
| GDSC1/2 drug response (COSMIC join) | 600 | 538 with expression |
| ProCan (949) | pending | no overlap computable |

Per lineage (screened / expr / CN / Nusinow / metabolomics / PRISM / GDSC), the weak spots: Biliary 36/36/34/1/5/28/2; Soft Tissue 50/49/46/4/7/28/8; PNS 47/40/40/1/16/17/22; Bone 59/49/51/7/13/25/23; Cervix 18/18/18/0/0/12/9; Eye 15/13/13/0/0/3/0; Prostate 19/19/15/5/5/5/6. Full table is reproducible from the scratch script logic (join on ModelID; GDSC via `Model.COSMICID`; Nusinow via `Model.CCLEName` after stripping `_TenPxNN`).

## Library and source structure of the labels (not mentioned in DESIGN)
`CRISPRGeneEffect` in 26Q3 merges five libraries: Avana (1,051 screens, Broad), KY (322, Sanger Project Score), Humagne-CD (221, Sanger library run via Broad condition), Brunello (6), TKOv3 (4). By model: Avana-containing 1,048; Humagne-only 194; KY-only 118; other 1. 213 models have 2 screens, 14 have >=3.
- Per-line positive rates (Chronos < -0.5, on U) are similar for Avana-containing (0.113) and KY-only (0.115), higher for Humagne-only (0.128, max 0.21). Chronos < -1: 0.056 / 0.055 / 0.061. So library shift is small but not zero.
- KY-only lines are lineage-skewed (Head and Neck 19 of 77, Ovary 10, Pancreas 8 are KY/Humagne-only) and are the lines lacking omics (below). Humagne-only is concentrated in CNS (44), Bowel (23), Pancreas (23), Eso/Stomach (27).
- 88 gene columns of U have NaN somewhere; 1,355 of 1,361 lines have at least one NaN in U (mean NaN fraction 1.4%). The NaN-exclusion rule in §3 is therefore active for nearly every line; it should be reported per gene, since it is probably library-specific.
- 108 of the screened models are organoids (Bowel 22, Pancreas 23, Eso/Stomach 27, Breast 14, Prostate 9, ...), screened mostly as 3DO (94 screens in total). Organoid media (OPAC:Ad+++, "Organoid Media Formulation #n", "Breast Organoid Media") are not covered by any of the 8 requested formulation sheets.

## Medium annotation (W2) — findings
1. **Screen-level medium is missing for the Sanger lines.** 1,231 of 1,361 models have a screen medium via ScreenSequenceMap.ModelConditionID -> ModelCondition.FormulationID. The remaining 130 are almost all Sanger-sourced: KY screens have a medium in only 8 of 322; 115 of the 118 KY-only lines have none. §4.2 would impute RPMI + 10% FBS for them, although `Model.FormulationID` (the onboarded medium) exists for 56 of the 130 and is a better guess than the modal medium.
2. **Screen medium vs onboarded medium differ** for 66 of 1,231 models (5.4%; 1,165 agree on the full text). Disagreements are real culture changes in the Broad screens (e.g. DMEM:F12 + 10% FBS used for lines onboarded as RPMI or EMEM). So `Model.csv` alone is not the screen medium. The design already prefers the screen medium; keep that, and use Model's as the second tier.
3. Six models have >1 distinct screen medium (`mixed`): ACH-000219, -000472, -000619, -000711, -002029, -002045. Handled correctly by the design.
4. **Parsability by the 8 requested sheets.** Single-base media: 1,225 models; 1,018 have base in {RPMI, DMEM, DMEM:F12, MEM/EMEM, IMDM, F12, McCoy's 5A, L-15}. With the Model.csv fallback, 1,074 of 1,361 lines (79%) are mappable to a requested sheet. Unmapped to a sheet: OPAC organoid (52), Tumor Stem Media (38), other organoid formulations (~47), AlphaMEM (15), DMEM:KSFM, William's E, Waymouth's, RPMI:F12 mixtures, etc. Mixed media ("DMEM:F12", "RPMI:F12 (1:1)") need an explicit rule (average of two bound vectors?). Neither the parser spec nor the download list covers them.
5. **Medium is confounded with lineage.** Screen base medium: Lymphoid 85/97 RPMI, Lung 89/126 RPMI, Biliary 28/36 RPMI, Esophagus/Stomach 46/96 RPMI; CNS 36 DMEM vs 17 RPMI; Head and Neck, Soft Tissue, Skin, Ovary, Pancreas have both. Overall RPMI 584 vs DMEM 187 among mapped lines. A shuffle *within lineage* (the design's negative control) is the right test, but it has no power in lineages with one medium. State that power is limited to the ~10 lineages with >=8 lines in each of two media.
6. **The medium effect on labels is large, so W2 is worth running, and so the non-GEM baselines must see it.** ASNS Chronos < -0.5: DMEM 49.7% (n = 181) vs RPMI 5.7% (n = 560). Within lineage: Head and Neck 67% vs 0%, CNS 40% vs 0%, Skin 69% vs 4%, Lung 50% vs 1%, Ovary 62% vs 0%. This is a clean positive control. It also means a gene-prior baseline that ignores medium is a weak comparator for exactly the genes where the GEM can win.
7. No MF-id -> composition table exists in `data/depmap/` (user download #4 not arrived). The MF family prefix (MF-001 = RPMI, MF-002 = DMEM, MF-005 = IMDM, MF-015 = EMEM, ...) is a usable cross-check on the text parser.

## Genomic backgrounds (W4 level 3): lines available (screened lines)
Definitions I could implement from the files; counts as the data are today.
| Group | Rule used | Lines | Main lineages |
|---|---|---|---|
| KRAS hotspot | hotspot matrix > 0 | 209 | Pancreas 63, Bowel 43, Lung 38 |
| BRAF hotspot | hotspot matrix > 0 (any hotspot) | 95 (V600* only: 81 incl. non-screened, ~76 screened V600E) | Skin 50, Bowel 13 |
| PIK3CA hotspot | hotspot matrix > 0 | 133 | Breast 22, Uterus 20, Bowel 18 |
| STK11 loss | damaging >0 or relCN < 0.3 | 59 (damaging 48, deletion 11); only **9** if TPM < 1 is also required | Lung 27, Cervix 9 |
| KEAP1 loss | same | 27 (damaging 25, deletion 2); only **1** with TPM < 1 | Lung 11 |
| MTAP homdel | relCN < 0.3 | 233 (200 with TPM < 1) | CNS 47, Pancreas 24, Lung 22 |
| ASS1 silenced | log2(TPM+1) < 1 | 161 | Skin 39, Lymphoid 35, CNS 13 |
| IDH1/2 hotspot | hotspot matrix > 0 | **8** | Biliary 2, CNS 2 |
- IDH1/2 cannot meet the "10 training and 5 test lines" rule in any split; drop it or merge it into a pooled "metabolic-driver" group.
- Applying the genomic-constraint rule as written ("biallelic damaging and TPM < 1") leaves STK11 at 9 and KEAP1 at 1 line: these groups are not usable as defined. Use "damaging (matrix >= 1) or deep deletion", with expression reported as a supporting stratum.
- Dry run of the split rule (not the real seed; 2,670 accepted of 3,000 attempts): group feasibility (>=10 training and >=5 test lines) holds in KRAS 97%, PIK3CA 98%, MTAP 100%, ASS1 99%, BRAF 74%, STK11 66%, KEAP1 22%, IDH 0% of draws. KEAP1 and STK11 depend on whether Lung lands in test (it is S1 rank 2; the draw decides). This means which backgrounds can be reported will be a by-product of the split. The design should list the groups per split once the split is drawn, not before.
- Every group is lineage-confounded (BRAF = Skin, STK11/KEAP1 = Lung, KRAS = Pancreas/Bowel/Lung). A group-level model has to be compared with the same-lineage model on the same lines, not only with the generic one.
- **CN scale.** `OmicsCNGeneWGS` is relative CN (diploid ~1, e.g. KRAS max 121). The 0.3 cut is not a deep-deletion call by DepMap, and ploidy is not corrected here. MTAP < 0.3: 336 of 1,763 CN lines (19%), < 0.1: 288, < 0.5: 441. MTAP/CDKN2A deletions are common, so the cut is fairly stable for those two but not generally. No DepMap deletion/gene-level call file is in `data/depmap/`. Say so in §W4 and give the ploidy handling (Model has no ploidy; Cell Model Passports has `ploidy`).
- Mutation matrices: values 0/1/2 (STK11: 4,685 / 25 / 80 over default rows). The meaning of 2 (two alleles/compound hets?) is not documented in DATA.md. "Biallelic" cannot be called from this matrix without verifying that. Alternatively use the long file (`Hotspot`, `LikelyLoF`, `AF`).
- Expression is log2(TPM+1): "TPM >= 1" is the same as value >= 1 (log2(2)). Fine, but state the transform wherever a TPM threshold is given (§3, §W4, §W5).
- 140 screened lines lack WGS CN (28% of KY-only lines have CN, 47% have expression). The genomic constraint, the CN stratum in W1, and `CNloss` in F5 need a stated fallback (mutation-only, or exclusion with the count).

## Proteomics, metabolomics, drug data: what W5/W6 can rest on
- Nusinow 309 screened lines. In W5 the co-abundance rule (r >= 0.5 over >=30 lines where both proteins are quantified) works for common enzymes; lineage-specific enzymes will fail the 30-line rule (Biliary 1, PNS 1, Cervix 0).
- W4 protein arm: the proteomic lines are ~23% of the cohort, lung/GI/skin/breast heavy; LOLO folds on lineages with <5 proteomic lines (8 of 14 plausible training lineages) will have no protein arm. State that the protein arm is evaluated on a subset (pooled across lineages, with lineage as a covariate) and cannot test transfer to unseen tissue except for Lung, Eso/Stomach, Bowel, Skin, Breast, Pancreas.
- ProCan (949 lines): per user list #5, still pending; there is no way to know overlap until it arrives. The design should not assume it.
- W6 drug/metabolomic comparisons: 652 (metabolomics), 735 (PRISM), 600 (GDSC) screened lines. The flux features are computed under the CRISPR-screen medium; PRISM pools are cultured in RPMI-type conditions and GDSC lines in the "Screen Medium" column of `Cell_Lines_Details` (D/F12 or R). The medium used for W6 flux states should be that of the assay, or the medium effect becomes a nuisance covariate. Add that choice to §W6.
- PRISM 24Q2 meta has no medium; 20Q2 secondary has `depmap_id`. Fine.

## Licences / redistribution of the SBML resource (W8)
- Human1/Human2: CC BY 4.0, derivative models are fine under CC BY 4.0 with attribution. Keep the Human-GEM LICENSE.md and a NOTICE in the resource.
- The models embed only reaction presence/bounds derived from DepMap (CC BY 4.0), Nusinow (no licence stated), and media tables. Reaction presence derived from proteomics is a derived product; the risk is low but Nusinow's missing licence should be recorded as "no restrictions stated; cite" and the resource should not ship per-protein values.
- Media tables (user downloads #6-8): Plasmax (Sci Adv, CC BY-NC) and HPLM (Cell, copyright, supplementary) are not CC BY. If `data/media/compositions.tsv` or bounds derived from them go into the resource, the effective licence is not clean CC BY 4.0. Options: ship only the pipeline and bounds (concentrations are facts, but the table layout is not safe to assume), or release the HPLM/Plasmax-variant models separately with a NOTICE. Gibco formulation sheets are product data, fine as facts but cite.
- DATA.md says "CC BY 4.0 per site (not re-verified)" for Cell Model Passports and GDSC. Re-verify before the release, since the whole resource's claim "CC BY 4.0" depends on the weakest input.
- DepMap: lines under the CC BY 4.0 licence; some CCLE/Broad cell-line-derived data are under "DepMap terms of use"; record the 26Q3 terms page version.
- CORUM 5.1 is read through a third-party Zenodo mirror. Used for curation evidence only; no redistribution.
- RAVEN-python / Gurobi: the models do not embed solver code; no restriction on SBML. The WLS academic licence restricts *use*, not output.

## Changes to DESIGN.md
### Must-fix
1. **§W2 / §4.2 medium tiers.** Replace the single "imputed RPMI" rule with a 3-tier `medium_source`: `screen` (ModelCondition via ScreenSequenceMap; 1,225 single-base models), `onboarded` (Model.FormulationID where no screen medium; +56 lines, mainly Sanger), `imputed` (the rest; ~74). Attribution analyses use `screen` only, with `onboarded` as a sensitivity set. Report the 5.4% screen-vs-onboarded disagreement. State that 115 of the 118 Sanger-only lines have no screen medium.
2. **§4.1 / §3 library covariate.** Add `library_set` (Avana / KY / Humagne-CD / mixed), `n_screens`, `model_type` (cell line vs organoid) and `omics_available` flags to `lineage_split.tsv`; report headline metrics stratified by library_set and by organoid vs cell line (reporting only). Chronos is one merged matrix, but KY-only and Humagne-only lines are 23% of the cohort and non-random in lineage and in omics coverage.
3. **§4.2 / §6.3 evaluation set.** The primary set is "test lines with a default expression profile". 62 of the 69 lines without expression are Sanger-only; say how many test lines this drops (count after the split is drawn) and that CN is missing for 140 lines. For W4's genomic constraint specify the fallback when CN is missing.
4. **§W4 definition of backgrounds.** Remove "biallelic damaging and TPM < 1" as the rule: it yields STK11 9 and KEAP1 1 lines. Use hotspot/damaging matrix >= 1 or relCN < 0.3, with expression as a separate stratum; document the 0/1/2 coding or switch to the long mutation file; drop IDH1/2 (8 lines); state that groups are decided *after* the split from the counts; require within-lineage comparison for groups. Specify the CN scale (relative, not ploidy-adjusted) and the sensitivity of the 0.3 cut (0.1/0.5: 288/441 for MTAP).
5. **§2.2-2.3, DATA.md numbers.** Correct |U| = 2,866, H1 2,825, H2 2,807, overlap 2,766 (the script will report exact), and say Nusinow = 375 lines of which 309 are screened. Fix the DATA.md "2,798/2,868" and metabolomics 927/652.
6. **§W2 baseline fairness.** Add a medium-aware non-GEM baseline (B2m = gene prior stratified by base medium class within fitting lines, and put medium class in B6 features). The ASNS DMEM-vs-RPMI difference (50% vs 6%) shows a medium-blind prior understates the non-GEM competitor on the genes that W2 targets; C5 would otherwise overstate the GEM's value.

### Should-fix
7. **§W2 parser and downloads.** Add rules for mixed media (DMEM:F12, RPMI:F12 1:1, DMEM:KSFM 2:1, IMDM:RPMI 4:1), AlphaMEM, William's E, Waymouth's, McCoy's, and organoid/stem-cell media. The 8 requested sheets cover 1,074 of 1,361 lines (79%); say that organoid/proprietary media (~150 lines) get `unmappable` and are excluded from attribution. Use the MF-id family as a parser cross-check.
8. **§W2 shuffle control power.** State that within-lineage shuffling is only informative in lineages with >=8 lines in each of two media (about 10 lineages); report the shuffle at the pooled level plus per eligible lineage. Use the ASNS result as an explicit label-based positive control (already planned as "Lagziel replication"; the data support it).
9. **§W6 assay medium.** For PRISM, GDSC and CCLE metabolomics, flux states should be computed in the assay's medium (or the medium is a covariate). Give the overlap counts (652 / 735 / 600) and a minimum per-drug line count.
10. **§4.2 / §W4 proteomics arm.** State that Nusinow lines are 309, 8 of the 14 plausible training lineages have <5 proteomic lines, and the protein arm is pooled with a lineage covariate. Do not claim protein-based transfer to unseen tissues beyond six lineages. Do not plan on ProCan until it is on disk.
11. **§W8 licence.** Add a NOTICE rule: the resource is CC BY 4.0 only if no CC BY-NC (Plasmax) or restricted (HPLM table) content is redistributed; otherwise the HPLM/Plasmax variants ship separately. Re-verify the GDSC and Cell Model Passports licences before M8, and record the Nusinow "none stated" status. Also state that 838 resource lines (2,130 - 1,292) have no CRISPR label and cannot be benchmarked; they carry `benchmarked = false` in the manifest.
12. **§3 NaN labels.** Report NaN per gene and per library (88 gene columns of U have NaN; 1,355 lines affected) and require the same gene x line set for all methods (already stated); add the per-library NaN table to `dropcounts.tsv`.

### Optional
13. Record screen-type (2DS/3DO/2DN, etc.) in the split file; 3D screens (94) have a different growth regime.
14. Cross-check lineages with `OncotreeLineage` vs `OncotreePrimaryDisease` for the small S3 lineages (Pleura 21, Prostate 19, Cervix 18, Eye 15); lineage-level bootstraps for Eye/Prostate/Cervix have little proteomic, metabolomic or drug support, so drop the omics-dependent per-lineage tables for them rather than reporting empty cells.
15. `Model.csv` has `PediatricModelType` (245 screened lines); lineage-neutral, so it can be a reporting stratum.
16. Cell Model Passports `ploidy` could normalise CN; not needed unless the 0.3 cut is shown to be unstable.

## Not an issue
- Entrez joins: 11 one-to-many genes per model are handled in the baseline mapping (first Entrez = primary); no many-to-one; DepMap columns have unique Entrez. The 2 unmapped model genes are counted. No further change.
- ScreenGeneEffect vs CRISPRGeneEffect: only the ModelID-level matrix is needed; screen-level is not required for any workstream.
- Filter `PassesQC == True, ExcludeFromCRISPRCombined == False, ScreenType != pDNA`: all 1,604 mapped screens already pass at screen level (26 screens have a failing replicate row); the filter does not remove models. It is harmless; apply it at the ScreenSequenceMap row level.
