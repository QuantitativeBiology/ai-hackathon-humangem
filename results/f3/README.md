# F3: pFBA traces of P1's 20 worst T1 cells

**post hoc / descriptive (A15).** This is the pre-registered DESIGN §10.3 item ("The 20 worst lineage × gene cells are traced by pFBA") that was not produced at T1. It was run in the A15 follow-up phase (§17), after the T2 unseal. It fits no rule, and it changes no model, score, call or headline. It does not replace T1/T2.

**Author:** F3 engineer agent (Claude Opus). **Date:** 2026-10-28.
**Code:** `scripts/f3_trace.py` (one run, 6 workers, GLPK, BLAS pinned to 1 thread, 114 s wall time; row in `MANIFEST.tsv`). **Test:** `tests/test_f3.py` (11 label-free toy tests).

**Label hygiene:**
- The session-start probe (`head -c 50` of the raw CRISPR gene-effect file) was **blocked** by the label guard, as required.
- No label file and nothing in `data/locked/` was opened.
- The only test-derived inputs are the unsealed, committed `results/m6/test/unsealed/T1/worst_cells.tsv` and `failure_account.tsv`.
- Every other input is label-free:
  - line assignment;
  - `sweep_e`;
  - the W3 B-loo spec;
  - the W2 media library;
  - W1 line and gene features and omics;
  - the P1 gene annotations.
- The failure-account class `common_essential` is label-derived, so it is **not reproduced** here.

## What was traced

**The cells.** All 20 cells are **CNS/Brain false negatives**: essential in DepMap and called viable by P1 in 126 of 126 lines (n_fp = 0).
- CNS/Brain is the largest test lineage, so its cells top the `n_errors` ranking.
- Many CNS/Brain cells tie at 126/126. `evaluate_test.py` breaks ties by ascending Entrez ID, so these 20 are the lowest-Entrez genes among the tied all-lines FNs, not a "worst of the worst".
- They are all pan-essential housekeeping genes.

**The model.** P1 = Human1 + W3 `B-loo:cofactor_pool_biomass` (MAM10012c removed from MAR13082) on documented medium M2(f=1), gate-free.
- Each of the 126 lines has one M2(f=1) bounds signature (`results/m6/scores/line_assignment.tsv`).
- The 126 lines use 6 signatures, all with serum: RPMI 81 lines, DMEM 35, DMEM:F12 4, IMDM 4, DMEM+pyruvate 1, RPMI+AAHT 1.
- Every cell was traced on all 6 signatures, with bounds rebuilt by `m6_sweeps.bounds_by_id` and the KO engine settings of `w3.ko_sweep`.
- The mechanism class is the same on all 6 media for every cell (`class_constant_across_media`).

**Check.** The recomputed gene-KO growth ratios equal P1's scores bitwise (`1 − e` in `sweep_e`, max deviation 0). All 20 genes have P1 ratio_r = 1.000 on every medium: the KO has no effect on growth.

**Per (cell, medium):**
- **(i) Gene KO:** `Gene.knock_out()`.
- **(ii) All-reactions KO:** a counterfactual that closes every reaction the gene appears in, as if every isozyme and alternative complex were off.
- **(iii) pFBA:** WT pFBA, and KO pFBA where the gene's reactions carry WT flux and the KO stays viable. The top 10 rerouted reactions and top 8 exchanges by |v_KO − v_WT| go to `reroutes.tsv`.
- **(iv) FVA** of the gene's reactions:
  - at ≥ 99.9 % of WT growth (`fva_opt`);
  - with growth free, which flags reactions as `blocked`.
- **(v) Bypass test:** each of the top 5 exchanges whose uptake rises in the KO pFBA is closed in the KO and in the WT (`bypass_tests.tsv`). It counts as a medium bypass if KO growth without it is < 1 % of WT growth without it.

**Classes (first match):**
1. `outside_model`: the gene is not in Human1.
2. `isozyme_in_gpr`: the all-reactions KO is lethal (< 0.01), so only the GPR OR hides the lethality.
3. `function_not_represented`: the gene's reactions carry no WT pFBA flux, because they are blocked or unused. Nothing the biomass objective needs depends on the modelled function.
4. `medium_uptake_bypass`: the reactions are used, the all-reactions KO is viable, and the KO depends on one extra uptake.
5. `alternative_pathway`: as 4, but with no single-exchange bypass.

The FP branch (which biomass precursor the KO cannot make, `blocked_biomass_components`) is implemented and tested on a toy, but was not needed: there are no FP cells.

## Per-cell table

The full table is `cells.tsv`, with media, GPR detail, FVA, alternative-gene expression and the auto note. `cells_by_medium.tsv` has one row per medium.

| Gene | Human1 reactions (GPR as seen from the gene) | Gene KO ratio | All-rxn KO ratio | WT pFBA flux on its rxns | Class | Failure-account classes |
|---|---|---|---|---|---|---|
| ALDOA | 4 aldolase rxns, OR (ALDOB, ALDOC) | 1.00 | 0.80 (1.00 on DMEM:F12, RPMI+AAHT) | 506 | alternative_pathway | paralog_buffered |
| ATP2A2 | MAR07629 Ca²⁺-ATPase c→e, OR of 9 ATP2A/B/C | 1.00 | 1.00 | 0 | function_not_represented | paralog_buffered, complex_subunit |
| ATP6V0C | MAR07799 (OR of 10), MAR00080 lysosomal V-ATPase (OR-of-AND, closed by KO) | 1.00 | 1.00 | 0 | function_not_represented | complex_subunit, housekeeping_adjacent, expected_fba_failure |
| COPB1 | 4 glucitol / ATP transport rxns, OR of 6 coatomer genes | 1.00 | 1.00 | 0 | function_not_represented | complex_subunit |
| DAD1 | MAR07285 OST N-glycan transfer, OR-of-AND (closed by KO) | 1.00 | 1.00 | 0 | function_not_represented | complex_subunit |
| DHPS | 3 spermidine / deoxyhypusine rxns, single gene | 1.00 | 1.00 | 0 | function_not_represented | — |
| DPAGT1 | MAR07264 first dolichol N-glycan step, single gene | 1.00 | 1.00 | 0 | function_not_represented | housekeeping_adjacent |
| DUT | 3 dUTP/dUDP hydrolysis rxns, 2 single-gene + 1 OR (ITPA) | 1.00 | 1.00 | 0 | function_not_represented | — |
| EEF2 | MAR04136 *extracellular* GTP hydrolysis, OR of 19 GTPases (blocked) | 1.00 | 1.00 | 0 | function_not_represented | paralog_buffered, complex_subunit |
| GAPDH | MAR04373 GAPDH, OR (GAPDHS) | 1.00 | 0.60 (0.66 DMEM:F12, 0.71 RPMI+AAHT) | 1000 (at bound) | alternative_pathway | paralog_buffered, complex_subunit |
| HSPA9 | not in Human1 (addressable in Human2) | — (e = 0 by convention) | — | — | outside_model | missing_from_model, complex_subunit |
| POLD1 | MAR07160 DNA synthesis, OR-of-AND of 23 genes | 1.00 | **0.00** | 1.25 to 2.25 | isozyme_in_gpr | paralog_buffered, complex_subunit |
| POLE | MAR07160, as above | 1.00 | **0.00** | 1.25 to 2.25 | isozyme_in_gpr | complex_subunit |
| POLE2 | MAR07160, as above | 1.00 | **0.00** | 1.25 to 2.25 | isozyme_in_gpr | complex_subunit |
| RFC5 | MAR07160, as above (RFC5 is itself a single-gene OR branch) | 1.00 | **0.00** | 1.25 to 2.25 | isozyme_in_gpr | complex_subunit |
| SOD1 | 3 SOD rxns (n, x single gene; c OR SOD2/SOD3), all blocked | 1.00 | 1.00 | 0 | function_not_represented | — |
| UBA1 | MAR09491 "isolated" protein ubiquitination, OR-of-AND of 176 genes (blocked) | 1.00 | 1.00 | 0 | function_not_represented | paralog_buffered |
| UBE2I | MAR09491, as above | 1.00 | 1.00 | 0 | function_not_represented | paralog_buffered, complex_subunit |
| OGT | 4 galactosyl-ceramide transferase rxns, OR of 48 glycosyltransferases | 1.00 | 1.00 | 0 | function_not_represented | complex_subunit |
| RNGTT | MAR09507/09508 mRNA capping, "isolated", single gene (blocked) | 1.00 | 1.00 | 0 | function_not_represented | — |

Ratios are KO / WT biomass flux. "1.00" means unchanged on all 6 media. The flux column is the maximum over media of Σ|v| in WT pFBA, in a.u. (uptake bounds of 1000).

## Mechanism note per cell

- **ALDOA (alternative pathway, behind an isozyme).**
  - The gene KO is hidden by ALDOC (median log2 TPM+1 3.7 in these lines) and ALDOB (0.06) in the OR.
  - Even with every aldolase reaction closed, growth stays at 0.80. pFBA reroutes glucose-6-P through the oxidative pentose-phosphate pathway (G6PD, PGLS, PGD to their 1000 bounds; PGI reverses), and transketolase/transaldolase supply GAP without aldolase.
  - No single uptake is required (closing glutamine, serine, etc. scales both arms alike).
  - On DMEM:F12 and RPMI+AAHT the aldolase flux is not needed at all (ratio 1.00).
- **GAPDH (alternative pathway, behind an isozyme).**
  - The OR partner is sperm-specific GAPDHS, unexpressed in all these lines (median 0.00). An expression-aware GPR would close the reaction.
  - Even then P1 would not call it lethal: the reaction KO keeps 0.60 to 0.71 of growth.
  - pFBA abandons lower glycolysis: glucose uptake falls from 533 to 20, and the mitochondrial H⁺ transport (Htmi) and malate export collapse. ATP and precursors come from glutamine, aspartate, alanine and histidine catabolism (aspartate and alanine aminotransferases reverse).
  - In WT, GAPDH runs at its 1000 bound. The arbitrary uptake scale makes the WT flux bound-limited, which is part of why a 40 % loss still reads as "viable" at the 0.01 cutoff.
- **POLD1, POLE, POLE2, RFC5 (isozyme in GPR; the only cells that GPR curation could flip).**
  - DNA synthesis is one lumped reaction (MAR07160) whose rule ORs 4 replicative complexes with 13 single-gene branches, among them POLB, POLG, POLK, REV3L, TENT4A and RFC5 itself.
  - Closing the reaction is lethal (0.00), because the biomass needs DNA. Any one expressed single-gene branch (POLB 4.9, POLG 4.6, POLA2 4.8, ...) keeps it open.
  - Leading/lagging-strand replication by Pol δ/ε (with RFC loading PCNA) is non-interchangeable in vivo. The OR encodes "some polymerase" and not "the replisome".
- **HSPA9 (outside model).** This mitochondrial Hsp70 (mortalin) is not in Human1 (it is in Human2), so P1 scores it e = 0 by the W1 convention. Its essential role, mitochondrial protein import and Fe-S cluster assembly, has no Human1 reaction.
- **EEF2 (not represented).**
  - EEF2 is one of 19 OR genes on an *extracellular* GTP hydrolysis reaction (MAR04136, compartment e), which is blocked on every medium because there is no extracellular GTP.
  - Translation elongation is costed only inside the lumped protein-pool biomass term, with no gene attached.
- **UBA1, UBE2I (not represented).**
  - Both sit on MAR09491, an "isolated" protein-ubiquitination reaction with 176 genes in the OR. It is blocked because ubiquitin has no source.
  - UBE2I is in fact the SUMO E2. Ubiquitin/SUMO conjugation is not a biomass requirement in the model.
- **RNGTT (not represented).** mRNA capping (MAR09507/09508) is "isolated" and blocked. RNA in the biomass is made without a cap step.
- **SOD1 (not represented).**
  - All three SOD reactions are blocked: no reaction can produce superoxide under these media (FVA with growth free).
  - Superoxide toxicity, the reason SOD1 is essential, is not a constraint of FBA.
- **DUT (not represented / alternative route in WT).**
  - dTMP is made from dUDP → dUMP (reverse dUMP kinase, MAR03968) and then TYMS, so the dUTPase reactions carry no flux.
  - The real essentiality (keeping dUTP out of DNA) is not modelled, because the DNA reaction cannot take dUTP.
  - The third DUT reaction also has ITPA in its OR.
- **DPAGT1, DAD1 (not represented).**
  - These are the first dolichol N-glycan step and an OST subunit. DPAGT1's reaction is feasible but unused; DAD1's is closed by the KO but carries no flux.
  - Protein N-glycosylation (ER folding and quality control) is not part of the biomass.
- **DHPS (not represented).** In Human1 its reactions are spermidine-cleavage side reactions. Hypusination of eIF5A, its essential function, has no biomass demand.
- **ATP6V0C (not represented; a pre-registered expected FBA failure).**
  - Both V-ATPase reactions (lysosomal H⁺ pumping) carry no flux. Lysosomal acidification has no use in FBA, because there is no pH or lysosomal-degradation demand.
  - The KO closes the lysosomal MAR00080, but nothing needs it.
- **ATP2A2 (not represented).** SERCA2 is modelled as a cytosol → extracellular Ca²⁺ ATPase (not ER), with 9 OR genes and zero flux. ER Ca²⁺ homeostasis is not a model function.
- **COPB1 (not represented).** Coatomer genes are attached to glucitol and ATP transport reactions (OR of 6), all with zero flux. COPI vesicle trafficking is not represented.
- **OGT (not represented, misannotated).** In Human1, OGT is one of 48 OR genes on galactosyl-ceramide transferase reactions (zero flux). O-GlcNAc transferase activity on proteins is absent.

## Summary of classes

| Mechanism class | Cells | Genes |
|---|---|---|
| function_not_represented | 13 | ATP2A2, ATP6V0C, COPB1, DAD1, DHPS, DPAGT1, DUT, EEF2, OGT, RNGTT, SOD1, UBA1, UBE2I |
| isozyme_in_gpr | 4 | POLD1, POLE, POLE2, RFC5 (all through the lumped DNA-synthesis OR) |
| alternative_pathway | 2 | ALDOA, GAPDH (behind an OR isozyme, and the reaction KO keeps 0.60 to 0.80 of growth) |
| outside_model | 1 | HSPA9 |
| medium_uptake_bypass | 0 | none: no single uptake carries any of the reroutes |

- 12 of the 20 genes are in Human1 but non-addressable (an OR alternative survives on every reaction), 7 are addressable (the KO closes ≥ 1 reaction: ATP6V0C, DAD1, DHPS, DPAGT1, DUT, SOD1, RNGTT), and 1 is absent (HSPA9).
- For all 7 addressable genes the closed reactions carry no flux, so the KO changes nothing.
- Only 4/20 (the DNA polymerase subunits) would flip to lethal under a GPR correction. These are AND/OR curation of a lumped reaction, not a paralog rule.
- For GAPDH and ALDOA even the reaction KO stays above P1's 0.01 cutoff.
- For the other 14 cells, the gene's essential function is not represented in Human1 or is outside the model. No GPR, medium or task change in P1's lever set would make it lethal.
- This is consistent with the §10.3 upper bound (54.7 % of positives addressable): these top errors are dominated by non-metabolic housekeeping biology.

## Cross-tab against the failure account (`crosstab.tsv`, `crosstab_long.tsv`)

A cell counts in every failure-account class it carries. Gene classes come from `results/m6/t1/p1_gene_annotations.parquet`. A cell counts as `unexpressed` / `cn_amplified` when ≥ 50 % of its lines are (expr < 1, CN ≥ 4). `common_essential` is not reproduced (label-derived).

| Failure-account class (share of all T1 FNs) | Worst cells | Their mechanism classes |
|---|---|---|
| complex_subunit (48.8 %) | 13 | not represented 7, isozyme 4, alternative 1, outside 1 |
| paralog_buffered (37.0 %) | 7 | not represented 4 (ATP2A2, EEF2, UBA1, UBE2I), alternative 2 (ALDOA, GAPDH), isozyme 1 (POLD1) |
| housekeeping_adjacent (11.6 %) | 2 | not represented 2 (ATP6V0C, DPAGT1) |
| expected_fba_failure (4.8 %) | 1 | not represented 1 (ATP6V0C) |
| missing_from_model (3.1 %) | 1 | outside 1 (HSPA9) |
| medium_sensitive (3.9 %), unexpressed (2.7 %), CN-amplified (0.2 %), biomass_artefact_lethal (0.5 %) | 0 | — |
| common_essential (39.1 %) | not reproduced | — |

**What the cross-tab shows:**
- **`paralog_buffered` overstates fixable redundancy.** Of its 7 worst cells, only POLD1 is truly hidden by the OR alone. For 4, the reactions are blocked or unused, so removing the paralog would change nothing. For 2, an alternative pathway keeps ≥ 60 % of growth.
- **`paralog_buffered` also misses the clearest GPR cases.** POLE, POLE2 and RFC5 are isozyme cases, but they are not flagged, because their OR partners are other polymerases and not Ensembl paralogs.
- **`complex_subunit` mostly marks function that is not represented.** It is not a model-fixable class.
- **The medium never mattered here.** No cell is medium-sensitive, no bypass exchange was found, and the classes are identical on all 6 media.

## Caveats

- pFBA optima are not unique, so the rerouted-reaction lists are one parsimonious solution. The class rule uses only growth ratios and Σ|flux| on the gene's own reactions, which are robust to that. FVA at ≥ 99.9 % of optimum confirms no reaction of a not-represented gene is required.
- Growth is in arbitrary units (uptake bounds 1000). GAPDH and the PPP enzymes hit their 1000 bounds, so the size of the partial losses (0.60 to 0.80) depends on that scale.
- The 20 cells come from a tie-break by Entrez ID among all-lines FNs in one lineage. They illustrate P1's dominant FN mode, but they are not a random or stratified sample of its errors.

## Files

| File | Content |
|---|---|
| `cells.tsv` | One row per cell: GPR, KO ratios, flux, FVA, class (modal over media, line-weighted), failure-account classes, expression of the gene and its OR alternatives, auto note |
| `cells_by_medium.tsv` | One row per cell × medium (120 rows): gene-KO ratio vs P1's ratio_r, reaction-KO ratio, WT pFBA flux, FVA, blocked/required, class |
| `reroutes.tsv` | Top 10 rerouted reactions and 8 exchanges (pFBA KO − WT) per traced KO and medium (ALDOA, GAPDH reaction KOs) |
| `bypass_tests.tsv` | Exchange-closure tests behind the medium-bypass rule |
| `crosstab.tsv`, `crosstab_long.tsv` | Mechanism class × failure-account class |
| `f3_summary.json` | Counts and checks |
| `MANIFEST.tsv` | Provenance row (inputs and outputs SHA-256, solver, versions) |
