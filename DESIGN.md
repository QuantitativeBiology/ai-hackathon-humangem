# DESIGN.md — nextgem pre-registered design

**Version v3 — FROZEN at M1** (review rounds 1–3 and user decision D2 applied), 2026-10-26. Author: design agent. Owner: PI. Freeze approved by the PI after the round-3 verifier's v3 final check; all its edits applied.
**Status:** frozen. The design freezes at milestone M1, after at most three review rounds. After M1, changing labels, splits, metrics, controls, candidate sets or decision rules counts as an *amendment*. An amendment is dated and logged in §17 before the analysis it affects is run. A change made after its analysis has run is reported as post hoc.

**Numbers.** Every numeric fact about the data comes from `scripts/design_numbers.py`, which writes `results/design/design_numbers.json`. The JSON key is given in brackets, e.g. **[dn:cohort.n_screened_lines]**.
- That script opens **no label-bearing file** (the screened-line set comes from `CRISPRScreenMap.csv`, the gene universe from the baseline mapping tables) and computes no performance or prevalence (`_meta.label_values_read = false`).
- Runtime and cost figures are not data facts. They cite their source: results/baseline/README.md, or the compute review probe in notes/design_review_r1/algorithms_compute.md.

**Read with:** PROMPT.md; LITERATURE.md (§8.1, §9); DATA.md; results/baseline/README.md; notes/design_review_r1/ and notes/design_review_r2/ (seven reviews + PI_decisions.md each); notes/design_review_r3/verification.md; **notes/apriori_lists_approved.md** (the frozen a priori lists, used verbatim); notes/label_exposure.md; RESEARCH_TREE.md decision D2. Appendix D holds the review dispositions and §17 the revision history.

---

## 0. Summary

- **Question.** Do optimised genome-scale models of cancer cell lines predict metabolic-gene essentiality better than generic Human1/Human2? Whole lineages are held out. **ftINIT is excluded by user decision D2**; the only GEM comparator is generic Human1/Human2 FBA, together with label-free and gene-prior baselines. The optimisation levers are:
  - medium;
  - tasks and objectives;
  - context-specific extraction;
  - GPR curation.
- **Metric.** Tie-aware expected average precision (tAP; §6.1). The null is the expected AP of a random ranking.
- **Primary endpoint (C1).**
  - *Measure:* the paired per-line ΔtAP between the selected GEM pipeline P1 and the stronger generic model G\*, on locked test lines, summarised as the **lineage-equal-weighted mean**. That is the mean, over the 7 test lineages, of each lineage's mean line ΔtAP.
  - *Inference:* a line bootstrap stratified by lineage, which is conditional on the 7 lineages, plus a sign-flip permutation p-value.
  - *Claim tiers (§6.4):* "detectable" if the one-sided α = 0.025 test passes; "improves" if, in addition, the point estimate is at least SESOI = 0.01; "meaningful" if the CI lower bound is at least SESOI. "Transfers to unseen tissue" requires the lineage-level t-interval together with the descriptive ≥5-of-7 rule.
- **Confirmatory sequence.** Fixed-sequence gatekeeping C1 → C1b → C3 → C4 → C5, each at one-sided α = 0.025. Everything else is estimation.
  - **C1b** = C1 with common essentials removed, using a list derived from training ∪ validation lines only.
  - **C3** = context specificity: per-gene ΔtAP of P1 versus expression only, above a within-lineage permutation null.
- **Required secondary endpoints (D2).** MCC, F1, precision, recall and accuracy of binary essentiality calls, computed under the essentiality protocols of the Human1 paper (Robinson 2020) and the Human2 paper (Luo 2026), on our splits (§6.8). They serve comparability only; tAP stays primary.
- **Leakage.**
  - A single label door, `src/nextgem/labels.py`, with fold-scoped access and an audit-hook tripwire.
  - Three tiers of whole lineages from a label-blind, size-stratified draw: 12 training lineages plus rare lines; 4 validation; 7 locked test. The split seed is **supplied by the user after the freeze commit**; the old seed 20261026 is retired because its draw was published in review round 1 (§4.1).
  - Hash-chained ledgers with look caps enforced in code.
  - Test detail beyond the pre-declared scalars stays sealed until P2 is frozen.

---

## 1. Goal → design map (PROMPT.md)

| PROMPT requirement | Where |
|---|---|
| FBA baseline on Human1 and Human2, all single KOs, AUPRC per gene, per lineage and overall | W1, §6 |
| Context-specific models from genomics, transcriptomics and proteomics: by tissue, by genomic background, and combined | W4 |
| Methods beyond CORDA/ftINIT; explicit test of the "cluster by method" failure mode | W4 (F0, F2–F5, F3b CORDA; variance partition with ≥ 3 established non-ftINIT families; permuted control). ftINIT is excluded by D2. |
| Exchange bounds from the documented medium; physiological variants; error attributable to medium alone | W2, §10.1–10.2 (attribution) |
| Tasks as objectives or constraints, from an a priori set; which tasks help in which contexts | W3 |
| GPR curation from proteomics and paralogs; curation vs constraints | W5, §10.1–10.2 |
| Flux states vs input omics (subtypes, metabolic drugs, metabolomics) | W6 |
| DepMap data used for fitting held out of evaluation; whole lineages withheld | §4, §5, §11 |
| Best predictor through an auto-research loop; runtime alongside accuracy | W7, §6.7 |
| SBML resource; benchmark vs generic; per-task and per-lineage breakdown; failure account | W8, §10.3 |
| Reproducibility: seeds, toy tests, controls, non-author review | §7, §11, §12 |

---

## 2. Prior art and reuse (PI decision 34)

Code that already exists on disk is reused, not re-implemented. Each reuse has an equality test against our wrapper.

| Existing code | What it does | How we use it |
|---|---|---|
| Human-GEM `code/DepMapGeneEss/*.m` (MATLAB): `generate_DepMap_models_ftINIT.m`, `getTaskEssentialGenes.m`, `evaluateDepMapEssentialityPredictions.m` | The Robinson/Gustafsson DepMap protocol: per-line ftINIT, then task-essential genes compared with Achilles. Default series is `1+0`. | **Read, not run.** MATLAB RAVEN is not available, and ftINIT is excluded (D2). Read for gene-ID handling and task-essential logic only. |
| Human-GEM `code/test/gradedEssentiality.py`, `estimateEssentialGenes.py`, `taskEssentialGenes.py`, `evaluateHart2015Essentiality.py` (Python) | ftINIT `1+1`, removeGenes = true, `resolve_ties`, `prove_abs_gap = 1`, time limit 900 s. Viability vs capability task scoring. Growth ratio on Ham's. Sharding and checkpointing. | **No ftINIT is run** (D2). W3's task-augmented scores and the paper-protocol binary calls (§6.8) call `taskEssentialGenes.py`/`raven_toolbox.tasks.check`. Equality test: our KO ratio = `growth_ratios()` on the generic parents. |
| `raven_toolbox` 3.0.0b1 (installed) | `tasks.check`, `tasks.tasklist`, `gapfilling` (`init.ftinit` is not used; D2) | task checks (W3, §6.8) |
| cobra 0.32.1 | `add_moma`, `fastcc`, `single_gene_deletion` | O-MOMA; F3 core consistency; cross-check of our KO engine on a sample |
| Repository Hart-2015 results (`data/testResults/gene-essential_summary.md`) | growth AUPRC per Hart line, pooled AUPRC and base rate [dn:hart_reference] | Context only. They come from ftINIT models, so they are not a control under D2. |

**What is new relative to this code and the literature** (LITERATURE §8.1 lists what we searched for and did not find):
1. A lineage-held-out, tie-aware AUPRC benchmark of generic H1/H2 on DepMap 26Q3, with binary metrics under the Human1- and Human2-paper protocols for comparability.
2. A medium ladder with screen-level medium annotation and an attribution estimate.
3. An a priori task set with negative controls.
4. Evidence-based GPR curation with a curation × constraints factorial.
5. A method-vs-tissue variance partition with a permuted-expression control.
6. Flux-information tests against capacity-matched omics.
7. A versioned per-line SBML resource for DepMap 26Q3.

---

## 3. Fixed inputs

### 3.1 Models
| Name | Version | Role |
|---|---|---|
| Human1 | Human-GEM **v1.19.0** (commit 022ed5c2) | Generic baseline B4 and a **full arm** at T1: P1 is re-applied to the H1 parent (PI decision 12), with LP families only (D2) |
| Human2 | Human-GEM **v2.1.1** | Generic baseline B5; default parent for optimisation |
| Human1-published | **v1.3.0** (fetched by tag with a pinned SHA-256) | W1 sensitivity only; first item on the cut list (§13) |

- The generic baselines use the models as released.
- Human1's reversible TYMS (MAR04644) bypasses UMPS/DHODH on thymidine-containing media (results/baseline/README.md). In **optimised** H1 pipelines it is set irreversible as a label-free W5 curation item, counted inside Δ_curation (domain S8, accepted as curation).
- Biomass is MAR13082 in both models.
- **M0 (Ham's)** = the 44 GR-task inputs [dn:tasks.hams_gr_inputs_n] at bound 1000. It includes hypoxanthine and thymidine and supplies cysteine, not cystine [dn:tasks.hams_gr_inputs]. As a result, the hypoxanthine-rescued de novo purine genes (PPAT, GART, PFAS, PAICS, ATIC) have ratio_r = 1 in both parents [dn:baseline_ko.human2.panel_ratio_r], a systematic, medium-driven false-negative block that W2 targets. ADSL and IMPDH2 also have ratio_r = 1, but hypoxanthine cannot rescue them biochemically (IMPDH1 OR-logic; ADSL to be traced); these are model errors, not medium errors (approved lists §1.2).

### 3.2 Data (all in DATA.md with SHA-256 and a `label_derived` column; §11.2)
- **DepMap 26Q3:**
  - CRISPRGeneEffect: 1,361 screened lines [dn:cohort.n_screened_lines] (counted from `CRISPRScreenMap.csv`, which is label-free).
  - CRISPRGeneDependency: present since 2026-10-26. Its ID sets are checked against CRISPRGeneEffect by `labels.py` at M1. (A v1 run of `design_numbers.py` found them equal; the script no longer opens label files.)
  - Expression (log2(TPM+1), default rows); WGS CN (relative); damaging and hotspot mutation matrices; Model; ModelCondition; ScreenSequenceMap; CRISPRScreenMap.
  - Common-essential and control lists, used **for reporting only**.
- **Proteomics:**
  - Nusinow (CCLE TMT): 375 lines mapped, 309 of them screened [dn:coverage_screened.nusinow_raw_mapped_total, .nusinow_screened];
  - ProCan-DepMapSanger: 947 mapped, 598 screened [dn:coverage_screened.procan_mapped_total, .procan_screened];
  - either source: 659 screened [dn:coverage_screened.any_proteomics_screened].
- **Metabolomics and drug screens:**
  - CCLE metabolomics: 927 lines, 652 screened [dn:coverage_screened.metabolomics_total, .metabolomics_screened];
  - PRISM 24Q2/20Q2: 735 screened [dn:coverage_screened.prism_24q2_or_20q2_screened];
  - CTRP: 612 screened [dn:coverage_screened.ctrp_screened];
  - GDSC1/2 (via `openpyxl`, pinned per PI decision r2-17): 610 screened lines, 546 of them with expression [dn:coverage_screened.gdsc_any_screened, .gdsc_any_screened_with_expr]. GDSC screen media are RPMI ("R") or DMEM/F12 ("D/F12") [dn:coverage_screened.gdsc_screen_medium_values]; both sheets are on disk.
- **Other:**
  - media composition files as they arrive [dn:medium.media_files_present];
  - CORUM 5.1 and Complex Portal;
  - Ensembl 116 paralogs;
  - the Hart 2015 inputs shipped with Human-GEM (not used under D2);
  - the Human1- and Human2-paper PDFs and supplements in `papers/`, which are the sources of the §6.8 protocols.

### 3.3 Gene universe
- **U**: genes of H1 or H2 whose primary Entrez ID is a CRISPR column. |U| = 2,866 [dn:genes.U_union]. Of these, H1 contributes 2,825 and H2 2,807 [dn:genes.human1_screened, .human2_screened].
- **U∩** = 2,766 genes present in both parents [dn:genes.U_intersection]. It is used as a sensitivity set.
- Joins use Entrez IDs only.
- **Missing genes.** A model that lacks a gene in U scores it e = 0. The gene joins the bottom tied block, so under tAP it is "uninformative", not "confidently non-essential". The same holds when extraction removed all of a gene's reactions.
- The DepMap inferred common-essential list is used only for a reporting sensitivity analysis; C1b uses its own list (§5.3). Its overlap with U is computed by `labels.py` at M1, because `design_numbers.py` no longer opens label-bearing files.

---

## 4. Splits (PI decision 3)

### 4.1 Procedure: `scripts/make_splits.py`, run exactly once at M1, with a seed supplied after the freeze (PI decision r2-1)
- **Why the seed changed.** The v0/v1 seed 20261026 is **retired**. The round-1 statistics review evaluated the draw for that seed under the v0 rules, named the resulting test and validation lineages in a committed file, and computed label aggregates on them. The draw function takes test and validation lineages as prefixes of the same per-stratum permutations, so that published draw also reveals most of a v1 draw with the same seed. v1's claim that "the real draw is unseen" was therefore wrong.
- **Seed procedure (commit-reveal). In this order:**
  1. **Freeze commit.** The PI commits the frozen DESIGN.md, FROZEN.tsv, `make_splits.py` and the draw function, which moves to `src/nextgem/splits.py` (importing `scripts/design_numbers.py` would run the whole script); both scripts import it from there. The commit SHA is recorded in the test ledger. **No draw with the production seed is computed before this commit**, by anyone; the production seed does not exist yet.
  2. **The user supplies an integer.** After the freeze commit, the user sends a positive integer S of their choosing by message. The PI records it verbatim in `results/splits/SEED_SOURCE.md` with the message timestamp and the freeze-commit SHA. A commit SHA, a timestamp or any value an agent can choose or grind is not acceptable as S.
  3. **`make_splits.py` runs once:** `make_splits.py --seed S --freeze-commit <SHA>`. It refuses to run if:
     - any attempt seed S + attempt (attempt 0…99) falls in 20261026…20261125, i.e. S lies in 20260927…20261125; the check is repeated per attempt inside the draw loop (V-2);
     - the freeze commit is not an ancestor of HEAD, or `git diff --quiet <freeze SHA> HEAD -- src scripts tests DESIGN.md FROZEN.tsv` fails, or those paths are dirty in the working tree (only ledger rows and `results/splits/SEED_SOURCE.md` may be committed between the freeze commit and the draw) (V-1);
     - `results/splits/SPLIT_SHA256` already exists (PI decision r1-17).
  4. **Draw.** `numpy.random.default_rng(S + attempt)`, starting at attempt 0. Every attempt, its acceptance-check results and every rejected draw are written to `results/splits/attempt_log.tsv`. The output records S, the freeze-commit SHA, the script SHA and the number of redraws.
  5. **Ledger.** A ledger row records the split SHA. The split is never redrawn.
- **Cross-check.** Before the production draw, `make_splits.py` must reproduce the reference draws for seeds 1–20 recorded in `results/design/design_numbers.json` [dn:split_simulation.reference_draws_seeds_1_to_20] and the acceptance statistics for seeds 1–1000.
- **Residual exposure.** Lineages named in the round-1 review may land in the new test tier by chance. What is known about them is aggregate label statistics only (§11.5). This is stated as a limitation (§15).

Steps:
1. **Universe.** The 1,361 screened lines. The lineage of each line is `OncotreeLineage`.
2. **Eligible lineages.** A lineage is eligible if it has ≥15 screened lines: 23 lineages, 1,331 lines [dn:cohort.n_eligible_lineages, .n_lines_eligible_lineages]. Lines of rarer lineages (30 lines [dn:cohort.n_lines_rare]) always go to training and are never held out.
3. **Strata.** Rank eligible lineages by size, descending, with ties broken alphabetically. S1 = ranks 1–8, S2 = 9–16, S3 = 17–23 [dn:split_simulation.strata].
4. **Draw.** Use `numpy.random.default_rng(S + attempt)`, where S is the user-supplied seed, and permute each stratum in the order S1, S2, S3.
   - **Test:** the first 2 of S1, the first 2 of S2 and the first 3 of S3 (7 lineages).
   - **Validation:** the next 1 of S1, the next 1 of S2 and the next 2 of S3 (4 lineages).
   - **Training:** the remaining 12 lineages plus the rare lines.
5. **Label-blind acceptance checks only.** If any check fails, increment `attempt` (maximum 100). Every rejected draw and the reason it failed are logged.
   - (a) every test and validation lineage has ≥10 lines with expression;
   - (b) test has ≥250 lines with expression;
   - (c) test has ≥40 lines with any proteomics;
   - (d) Lymphoid and Myeloid are not both in test ∪ validation. This keeps a haematopoietic, suspension-culture lineage in training for medium and serum effects.

   Simulated over seeds 1–1000 (never the production seed, which does not exist until after the freeze):
   - acceptance is 0.854, and only (d) ever fails [dn:split_simulation.acceptance_rate, .constraint_fail_rate];
   - accepted test tiers have a median of 345 expression lines (5–95%: 297–396) [dn:split_simulation.accepted_test_expr_lines_median_p05_p95];
   - the largest test lineage holds a median share of 0.30 [dn:split_simulation.accepted_largest_lineage_share_median_p95];
   - validation has a median of 181 expression lines [dn:split_simulation.accepted_val_expr_lines_median_p05_p95].
6. **Output.** `results/splits/lineage_split.tsv`, with one row per screened line. Covariates (PI decision 24; data MF2, O13):
   - tier, stratum and lineage;
   - `library_class`: Avana-containing / KY-only / Humagne-only / other = 1,048 / 118 / 194 / 1 [dn:screens.library_class_counts];
   - `n_screens`, `screen_type` and `model_type`, which flags the 108 organoids [dn:cohort.model_type_screened];
   - `medium_source` and `base_medium_class` (§8 W2);
   - `has_expression`, `has_cn`, `has_nusinow`, `has_procan`, `has_metabolomics`, `has_prism`, `has_ctrp`.

   The script also writes `split_summary.tsv`, the attempt log and the SHA.
7. **Background groups.** The genomic-background groups that are evaluable (W4) are listed at M1 from counts on this file, before any scoring.

### 4.2 Tier use
| Tier | Labels may be used for | Never for |
|---|---|---|
| Training: 12 lineages + rare lines; LOLO over the 12 eligible lineages, with nested inner LOLO for any choice | fitting thresholds, calibrations, ML, gene priors, training-fold co-essentiality, the C1b common-essential list; screening candidates | — |
| Validation: 4 lineages | choosing among ≤3 shortlisted candidates per decision. A look reveals **only** scalar endpoints with CIs, for every workstream (PI decision 21). There are at most **40** looks in total, ≤20 of them for W7 (ext-PI O2). | fitting |
| Test: 7 lineages | T1 and T2 only (§10.2) | anything else |

- **Inputs of test lines may be used.** Expression, CN, proteomics and media of test lines are label-free and are used to build their models.
- **Fitted statistics come from fitting lines only.** These include gene-specific thresholds, PCA loadings and the CellFie-inspired percentiles. For test evaluations, every fitted component of candidates *and* baselines is refit once on training ∪ validation, with hyperparameters frozen.

### 4.3 Evaluation sets and fallbacks (PI decision 25)
- **Primary evaluation set:** test lines with a default expression profile. Lines without expression (69 screened lines [dn:coverage_screened.no_expression]) are scored for generic models and baselines only, and their count is reported. They are excluded from every context-specific comparison.
- **Lines without WGS CN** (140 [dn:coverage_screened.no_cn]): the CN-based genomic constraint is omitted (flagged). B6 handles the missing CN natively (the gradient-boosting learner accepts NaN).
- **Lines without proteomics:** excluded from protein-arm comparisons. On test, those comparisons are **descriptive**, with their lineage composition shown (stats SF8).
- **Missing medium annotation:** see W2 tiers. The split itself ignores medium.

---

## 5. Labels (frozen at M1)

### 5.1 Definitions
| Label | Definition | Role |
|---|---|---|
| **L-prob** | `CRISPRGeneDependency` > 0.5 | **Primary**. The file is present; `labels.py` verifies at M1 that its ID sets equal those of CRISPRGeneEffect. |
| **L-ge05** | Chronos gene effect < −0.5 | Sensitivity. **It becomes primary if κ(L-prob, L-ge05) < 0.5 on training lines.** This rule is fixed now (stats r1 SF6). |
| **L-ge10** | Chronos gene effect < −1 | Sensitivity |

- **How labels are produced and accessed.** `src/nextgem/labels.py` builds the labels at M1 (§11.1). Every gene × line pair whose label source is NaN is excluded for all methods. Exclusions are counted per gene and per library in `dropcounts.tsv` (data SF12).
- **Known dependency on all lines.** L-prob, and DepMap's own common-essential list, are produced by DepMap from fits over **all** lines, test lines included. This is an accepted **label-generating** dependency and not a fitted one (stats SF6). It is noted in §11.
- **κ rule mechanics** (stats r2, R2-4).
  - κ is computed **once**, inside `labels.py` at M1, on training lines, pooled over all U pairs.
  - It is logged in the validation ledger as a non-look event before any metric is computed, and the outcome is written to FROZEN.tsv in the same run.
  - The rule is mechanical: no agent or PI decides after seeing κ; any deviation is an amendment.
  - κ without CE_fit is reported for information only and does not enter the rule.
- **No label statistics on non-training lines.** From now on, no reviewer or agent computes any label statistic on validation or test lines outside `evaluate_val.py`/`evaluate_test.py` (AI r2, R2-11).
- **Gene priors use the same label as the evaluation.** B2, B2m and C1b's common-essential list are always computed from the label being evaluated.

### 5.2 Estimability (PI decision 4)
- **Per line.** AP_line requires ≥5 positives in that line within the cell's gene set. Lines with fewer are NA, are counted, and are excluded identically for every method.
- **Per gene.** AP_gene requires ≥5 positive and ≥5 negative lines.
- **Per cell.** A cell of the reporting grid (§6.6) is reported as **NA** if fewer than 50% of its lines, or fewer than 30 of its genes, are estimable. The statistics review expects the L-ge10 × no-common-essentials cells to be NA; the rule decides, not the expectation.

### 5.3 Common essentials (PI decision 4)
- **CE_fit**: genes positive in ≥90% of the *fitting* lines for the evaluation at hand, under the label being evaluated:
  - training-minus-fold for LOLO;
  - training for validation looks;
  - training ∪ validation for test.

  C1b removes CE_fit.
- The DepMap 26Q3 list is a reporting sensitivity analysis and is **never used for fitting**.

### 5.4 Reporting strata (fixed now; files frozen at M1)
- **Expressed pairs:** log2(TPM+1) ≥ 1, i.e. TPM ≥ 1. This stratum is a required secondary analysis for every context-specific claim (PI decision 8).
- **Selective genes:** prevalence on fitting lines in [0.05, 0.95].
- **Addressable genes:** the KO closes ≥1 reaction in the parent. That is 1,120 genes in H2 [dn:baseline_ko.human2.addressable_genes_ko_closes_reaction] (domain/ext-PI SF7). The recall ceiling, i.e. the fraction of positives outside this set, is computed on training at M2.
- **Library class** and **organoid** status (PI decision 24).
- **CRISPR-artefact strata** (domain S5, accepted in part):
  - CN ≥ 4 (relative) amplicons;
  - the multi-copy loci, including their completion rule;
  - the housekeeping-adjacent families;
  - the moonlighting list.

  These three lists are those of `notes/apriori_lists_approved.md` §1.6, used verbatim; `data/controls/gene_panels.tsv` is generated from that file at M1 (Appendix A). A guide-multiplicity stratum from the Avana/KY guide maps is approved at low priority (PI decision r2-17). It is added only if the maps arrive, and **no claim depends on it**.

---

## 6. Metrics and statistics

### 6.1 Scores and AP
- **Score.**
  - ratio_r = clip(round(KO/WT, 6), 0, 1) and e = 1 − ratio_r.
  - Task-augmented scores: e = 1 if the KO fails a required task, otherwise 1 − ratio_r.
  - MOMA: ratio_r of the biomass flux.
- **Tie-aware AP (tAP)** (PI decision 1; McSherry & Najork 2008, LITERATURE ref 92). This is the exact expectation of AP over uniformly random orderings within each tied score block.
  - For a block g of size n_g, with p_g positives, preceded by N items of which P are positive, the j-th positive of the block falls at in-block rank r with hypergeometric order-statistic probability C(r−1, j−1)·C(n_g−r, p_g−j)/C(n_g, p_g).
  - Its expected precision contribution is (P + j)/(N + r).
  - tAP is the sum of these contributions over blocks, positives and ranks, divided by the total number of positives.
  - Implemented in `src/nextgem/metrics.py`. The binomials overflow for large blocks, so the code uses the equivalent O(n_g) form (stats r2, verified identical): a block contributes Σ_{r=1..n_g} (p_g/n_g)·[P + 1 + (r−1)(p_g−1)/(n_g−1)]/(N + r), with the bracketed fraction taken as 0 when n_g = 1. Both forms are implemented and asserted equal on test cases. Blocks with p_g = 0 contribute nothing but still advance N.
  - Verified (§12) against brute-force enumeration on toy cases and against Monte-Carlo averages of score + 1e-9 noise (10^5 permutations).
  - sklearn `average_precision_score` is kept as a **sensitivity row only**.
- **Null.** E[AP_random](n, n_pos) = the tAP of a constant scorer, which equals the expected AP of a random ranking. Then:
  - **Lift** = tAP / E[AP_random];
  - **ΔAP_null** = tAP − E[AP_random];
  - prevalence π is shown alongside for readability.
- **Binary calls** (PI decision 13). MCC, precision and recall at the fixed cutoff ratio_r < 0.01 are reported as secondary estimates for every FBA-type method, including P1 vs G\* at T1. Non-FBA scorers are dichotomised at their training-tier positive rate, which is fixed on fitting lines.

### 6.2 Units
| Unit | Definition | Role |
|---|---|---|
| tAP_line(l) | over genes, for line l | primary unit |
| **Δ̄_lin** | mean over lineages k of [mean over lines l ∈ k of ΔtAP_line(l)] | **primary summary** (equal weight per lineage) |
| Δ̄_line | mean over lines of ΔtAP_line | secondary (line-weighted) |
| tAP_gene(g) | over lines; scored if ≥5 positive and ≥5 negative lines. For a line-constant scorer it equals E[AP_random] exactly. | context specificity (C3) |
| mΔAP_gene | mean over eligible genes of tAP_gene − E[AP_random]. In each bootstrap replicate, genes with 0 positives or 0 negatives are dropped, and the mean count dropped is reported. | C3 |
| per-lineage | headline = mean tAP_line within the lineage, with Δ vs G\* and lift. Pooled AP within the lineage is secondary. Lineages are never ranked by raw AP (stats SF7). | breakdown |
| pooled | tAP over all gene × line pairs | overall (gene-identity-inflated) |
| AUROC | same units | secondary |

### 6.3 Confirmatory contrasts (PI decisions 2, 4, 5, 7)
Gatekeeping runs in the fixed order **C1 → C1b → C3 → C4 → C5** (PI decision r2-2 as amended by D2: C6 is removed with ftINIT). Each contrast is tested one-sided at α = 0.025, and testing stops at the first failure. Contrasts after a failure are reported as estimates only.

**Consequence of the simple chain, stated in advance.** A C3 failure (P1 not beating expression at per-gene ranking, a high bar) removes confirmatory status from C4 and C5, which are then reported as estimates. A graphical α split (Bretz 2009) was proposed in review and rejected by PI decision r2-2.

| ID | Statistic (test tier; parent and pipeline as stated) | Test | Look |
|---|---|---|---|
| **C1** | Δ̄_lin(P1 − G\*), all genes in U | sign-flip permutation: each line's ΔtAP sign is flipped at random, within lineage, and the **lineage-equal-weighted** Δ̄_lin is recomputed; 10,000 flips from SEEDS[0]; p = (1 + #{≥ observed})/(1 + 10,000). 95% CI from a lineage-stratified line bootstrap (B = 10,000) | T1 |
| **C1b** | as C1, on U minus CE_fit (from training ∪ validation) | same | T1 |
| **C3** | mΔAP_gene(P1) − mΔAP_gene(B3), all genes | the line-bootstrap CI lower bound > 0, **and** mΔAP_gene(P1) above the **97.5th** percentile of a within-lineage label-permutation null (consistent with one-sided α = 0.025). That null permutes each gene's labels across lines within each lineage: **2,000** permutations spawned from SEEDS[0]. | T1 |
| **C4** | Δ̄_lin(P2 − G\*) | as C1 | T2 |
| **C5** | Δ̄_lin(P2 − B6) | as C1 | T2 |

- **G\*** = the better of generic H1 and H2 on validation Δ̄_lin, fixed at M2.
- **Estimates only (no test; point estimate and CI):**
  - C2 = P1 vs B2 (and vs B2m, B2c);
  - P1 vs **B8+ and B8−** (both tie-break directions, PI decision r2-4) and B8′. "Adds information beyond generic + expression" requires P1 to exceed **all three** (stats r2);
  - within-parent P1(H1) − B4 and P1(H2) − B5 (PI decision 12);
  - MCC contrasts;
  - Δ̄_line;
  - library-class and organoid strata;
  - the expressed-pairs versions of C1 and C3.
- **C3 claim condition.** The claim "P1 captures line-specific dependency beyond expression" needs C3 to pass **and** the expressed-pairs version to have the same sign with CI lower bound > 0 (PI decision 8).

### 6.4 Effect size, claim language and transfer (PI decisions 2, 6)
- **SESOI for C1, C4 and C5 = 0.01 absolute in Δ̄_lin (all genes).** Basis (PI decision r2-5):
  - It is a **convention fixed by expert judgement** before any model was scored, proposed independently by two reviewers. Whether 0.01 is *meaningful* rests on that judgement, not on data.
  - The power simulation answers a different question. It shows that 0.01 is *detectable* (statistics r2: on a ~350-line, 7-lineage tier the minimum detectable effect is ≈ 0.19·s, where s is the within-lineage SD of line ΔtAP). It is not evidence that 0.01 matters.
  - Under tAP, float-tie and tie-breaking artefacts are zero by construction, so 0.01 is not within metric noise.
  - REPORT.md also translates Δ into an interpretable unit: the change in the mean number of true positives among each line's top-k genes, with k = that line's number of positives.
- **SESOI for C1b (no common essentials), lift-scaled.** Without common essentials the null AP is several-fold lower, so an absolute 0.01 would demand a much larger relative gain. Therefore SESOI_C1b = 0.01 × [mean E[AP_random] on U minus CE_fit] ÷ [mean E[AP_random] on U]. Both means are computed on **training lines at M2, before any model is scored**, and frozen. Basis: the same gain relative to the random-ranking null as 0.01 represents on all genes.
- **C3 has no SESOI;** only "detectable" wording applies to it.
- **Claim language** (C1, C1b, C4, C5, each with its own SESOI):
  - "detectable improvement": the contrast passes its test;
  - "improves": it passes and its point estimate ≥ its SESOI (0.01; SESOI_C1b for C1b);
  - "meaningfully improves": its CI lower bound ≥ that SESOI.
- **Headroom closed** = Δ̄_lin(P1 − G\*) / Δ̄_lin(B2 − G\*), with a bootstrap CI. It is reported only if the denominator's CI excludes 0; otherwise both raw contrasts are reported.
- The abstract-level sentence must give C1 and C1b together, plus C4/C5 once T2 has run.
- **Transfer to unseen tissue:**
  - (1) the C1 CI is **conditional on the 7 test lineages**, and the headline says so;
  - (2) a lineage-level t-interval with df = 6 on the 7 lineage means is reported, together with the exact sign-flip p over 2^7 configurations;
  - (3) **descriptive criterion:** Δ_k > 0 in ≥ ⌈0.7·7⌉ = 5 of 7 test lineages.

  "Generalises to new tissues" is written only if both (2) has a lower bound > 0 and (3) holds.
- **Optimistic supplementary estimate** (stats O3): the frozen P1 protocol's LOLO Δ over the 16 non-test eligible lineages, labelled selection-contaminated.

### 6.5 Uncertainty
- **Bootstrap.** Lines are resampled within lineage, using SEEDS[0]. B = 10,000 for C1–C5; B = 2,000 otherwise. Intervals are percentiles; BCa is used only if the training-tier distribution of line ΔtAP has |skew| > 1 (checked at M2; stats O6).
- **Pairing.** Methods are paired on identical gene × line sets, and an assertion fails otherwise.
- **Conditionality.** The CI is conditional on the fitted models: predictions are fixed and lines are resampled.
- **Multiplicity.** Benjamini–Hochberg (q = 0.05) is used, labelled exploratory, wherever many units are tested.
- **W6 tests.** W6 runs one Wilcoxon test per analysis, Holm-corrected across its three analyses (stats O5).

### 6.6 Reporting grid
For every evaluated method, every cell of this grid is reported, with value, E[AP_random], π, lift and CI, or NA under §5.2:
- labels: L-prob, L-ge05, L-ge10;
- gene sets: all genes, minus CE_fit, minus the DepMap common-essential list;
- universes: U, U∩;
- units: line, gene, lineage, pooled;
- strata: §5.4.

### 6.7 Runtime
Every result table carries wall and CPU time per line and in total (build, KO sweep, ML fit and predict), peak RSS, solver and worker count. These come from the provenance manifest (§11.4).

---

### 6.8 Required secondary endpoints: binary metrics under the Human1- and Human2-paper protocols (user decision D2)

**Purpose and status.**
- These metrics exist for comparability with the published essentiality evaluations of Human1 (Robinson et al. 2020, *Sci Signal* 13:eaaz1482) and Human2 (Luo et al. 2026, *PNAS* 123:e2516511123).
- They are **estimates**: outside the gatekeeping chain, never used for selection or for any gate, and reported for G\*, generic H1, generic H2, P1(H1) and P1(H2), each with lineage-stratified bootstrap CIs (§6.5) and paired differences vs G\*.
- tAP stays the primary metric.

**What "under the paper's protocol" means here.** Both papers built per-line tINIT/ftINIT models. Those are excluded by D2, so the protocol is applied as **the papers' KO call rule, label definition and metric definitions, used on our models and our splits**, not as their extraction step.

**Protocol details extracted from the papers** (sources in `papers/`; page numbers are PDF pages):

| Element | Human1 paper (Robinson 2020) | Human2 paper (Luo 2026) |
|---|---|---|
| Model per line | tINIT model per line from RNA-seq; 621 DepMap lines (DepMap 19Q3) plus 5 Hart 2015 lines (main text pp. 5–6; SM `aaz1482_sm.pdf` pp. 13–14) | ftINIT, series "1+1", per line (SI Appendix `pnas.2516511123.sapp.pdf` p. 2, lines 25–33); 5 Hart lines plus 625 DepMap 19Q3 lines (SI p. 3, lines 62–65); 908 Project Score lines (SI Fig. S3A, p. 11) |
| KO call rule | the gene is deleted; reactions with no compensating isozyme get bounds 0; **"essential" if the deletion impairs at least one of the 57 metabolic tasks** (data file S3, sheet "HumanGEM") (SM p. 14) | **identical**: "a gene was considered 'essential' … if its deletion impaired at least one of these 57 metabolic tasks" (SI pp. 2–3, lines 58–62) |
| Biomass-only variant | essential only if the growth (biomass) task fails (SM p. 14; Fig. S11, p. 29) | not reported |
| Label (DepMap) | gene effect **< −0.6** = essential, taken from the dataset's own publication; robustness grid "from −1.0 … to 0.0" (SM p. 14; Fig. S10, p. 28) | **not stated** in the main text, the SI or Datasets S1–S9 (S1–S9 contain no essentiality data) (repository script: < −0.6) |
| Metrics | confusion matrix per line; accuracy, sensitivity, specificity, MCC; hypergeometric enrichment with BH (main text p. 6; SM Fig. S9, p. 27). MCC is the headline, because accuracy is uninformative under imbalance (main text p. 6) | MCC (main text Fig. 2D legend, p. 3; SI Fig. S2D, p. 10; Fig. S3A, p. 11); "accuracy" (main-text Methods, p. 9) |
| Aggregation over lines | per-line values, shown as distributions over 621 lines, plus the mean MCC (main text p. 6; SM Fig. S10, p. 28) | per-line values over 625 / 908 lines (SI Figs. S2D, S3A); summary statistic not stated |
| Gene set scored | not stated in the papers; repository script: parent-model genes ∩ screened, non-NaN | not stated in the papers; repository script: parent-model genes ∩ screened, non-NaN |
| Medium for task checks | each task's own inputs (data file S3, column "IN"); GR = Ham's | the same 57-task list (SI p. 2, line 30: "growth on Ham's media") |

**Our pre-specified implementation** (fixed now; every choice for a detail the papers do not state is marked as ours):
- **Call rule P57.** A gene is essential if its KO makes ≥ 1 of the 57 tasks of `metabolicTasks_Essential.txt` infeasible (`raven_toolbox.tasks.check`, each task's own inputs). This is the same 57-task list as Robinson data file S3 (to be confirmed by a script diff at M1). Applied to every GEM method.
- **Call rule PGR (biomass only)**, from the Human1-paper protocol only: essential if the GR task becomes infeasible.
- **Label L-ge06** = Chronos gene effect < −0.6. It is a new, frozen, secondary-only label, built by `labels.py` at M1.
  - Robinson applied −0.6 to CERES (19Q3); we apply it to Chronos (26Q3). The scale difference is stated as a limitation.
  - The Human2 paper does not state a DepMap threshold. Not stated in the paper; the Human-GEM evaluation script uses gene effect < −0.6 (line 117), i.e. L-ge06.
  - Human1-paper robustness grid: 21 thresholds, −1.0 to 0.0 in steps of 0.05, as in SM Fig. S10 (p. 28) and `linspace(-1,0,21)` in `evaluateDepMapEssentialityPredictions.m` (line 106).
- **Gene set (not stated in either paper; taken from the repository script, lines 121–146).** The parent model's genes (H1 for generic H1 and P1(H1); H2 for generic H2 and P1(H2)) that have a CRISPR column and a non-NaN label in that line. A parent gene absent from a line model is predicted non-essential. U-wide scoring is a sensitivity row.
- **Metrics.**
  - Per line: TP, FP, TN, FN, then MCC, F1, precision, recall (= sensitivity), specificity and accuracy. The papers report MCC, accuracy, sensitivity and specificity; we add precision and F1 from the same confusion matrix.
  - MCC is the headline binary metric. Accuracy is reported, but it is dominated by negatives.
  - A metric with a zero denominator is NaN, as in the repository script; NaN lines are excluded from the means and counted per method and protocol. A sensitivity row sets them to 0. F1 = 2TP/(2TP + FP + FN), as in the script.
- **Summaries.** Three summaries:
  - the lineage-equal-weighted mean (our convention, §6.2);
  - the line-weighted mean (the papers' convention; Robinson 2020 reports the mean MCC);
  - the per-line median (descriptive; Luo 2026 shows per-line violins).
- **What each protocol reports:**
  - *Human1-paper protocol:* P57 and PGR × L-ge06 (+ grid); MCC, accuracy, sensitivity, specificity, plus precision and F1.
  - *Human2-paper protocol:* P57 × L-ge06; MCC and accuracy, plus F1, precision and recall.
- **Not reproduced:**
  - the papers' per-line extraction (D2);
  - Hart 2015 and Project Score labels (not DepMap 26Q3 and not on our splits; Project Score is label-derived per §11.2);
  - the hypergeometric enrichment test (it duplicates what tAP measures).
- **Non-GEM baselines.** They have no task-based call. For context only, they are dichotomised at their fitting-line positive rate, flagged "not a paper protocol".
- **When reported:**
  - M2: generic H1/H2 on training and validation;
  - validation looks: scalars only;
  - T1: revealed scalars for G\* and P1 under both protocols;
  - per-lineage detail sealed until T2 (§10.2).

## 7. Seeds
- **SPLIT_SEED = S, supplied by the user after the M1 freeze commit** (§4.1). It is used for the split only. **20261026 is retired** and refused by `make_splits.py`.
- **SEEDS = [11, 23, 37, 41, 53, 67, 79, 83, 97, 101]**, in `src/nextgem/config.py`. Uses:
  - bootstraps and permutation streams: SEEDS[0], expanded with `SeedSequence.spawn`;
  - random baselines, pseudo-tasks, random GPR edits and permuted-expression controls: SEEDS[0:5];
  - Gurobi `Seed` with Threads = 1: SEEDS[0];
  - LP-family stability study (Stage-A row A7): SEEDS[0:3];
  - ML and CV shuffles: SEEDS[0:3].
- No unseeded randomness. A test greps for `default_rng()` and `random_state=None`, and flags sklearn estimators constructed without `random_state` (AI-review note).

---

## 8. Baselines (PI decisions 9–11)

| ID | Baseline | Fitting |
|---|---|---|
| B0 | Constant score (tAP = E[AP_random]) | none |
| B1 | Uniform random scores, SEEDS[0:5] | none |
| **B2** | Gene prior: the fraction of fitting lines in which the gene is positive (same label) | within tier |
| **B2m** | Medium-aware gene prior: the prior within the line's base-medium class, falling back to B2 if the class has <10 fitting lines (data MF6) | within tier |
| B2c | Mean Chronos effect over fitting lines (continuous; sensitivity analysis) | within tier |
| B3 | Expression only: log2(TPM+1) of the gene | none |
| B4 / B5 | Generic H1 / H2 FBA on M0 | none |
| **B6** | Non-GEM ML: HistGradientBoosting on B2, B2m, own expression, CN, mutation, neighbour-enzyme expression (shared non-currency metabolite in H2; our re-implementation of the Lagziel-style feature, not their code), max paralog expression, isozyme count, **medium class + serum flag** from the **same tiered medium annotation as W2**, imputed lines included, so B6 and the GEM arms see identical medium information (PI decision 11; prior-work r2). No lineage features. | training LOLO; grid in W7 |
| **B8+ / B8−** | Generic + expression tie-break: rank by e_{G\*}, and break ties by **higher** (B8+) or **lower** (B8−) expression. **Both directions are reported** (PI decision r2-4); neither is chosen on test. B8′ (sensitivity) = e_{G\*} with genes at TPM < 1 set to 0. | none |

**Why B8 is needed.** Generic FBA scores are near-binary (§3.1, W1). Any pipeline that orders the large tied block by an expression-correlated quantity would gain without metabolic insight. B8 measures that gain directly, and "adds information beyond generic + expression" requires P1 > B8+, B8− and B8′ (estimates; ext-PI r1 MF1, stats r2).

**Gate.** An optimised candidate is eligible for selection only if its paired validation Δ̄_lin vs G\* has CI lower bound > 0. It is always shown next to B2, B2m, B6 and B8± (D2: the comparator set is generic H1/H2 as G\*, plus B0–B3, B6 and B8±/B8′).

---

## 9. Workstreams

**Template.** Each workstream gives its question, inputs, a priori candidates, method, metric, controls, decision rule and compute.

**Training gate (TG), harmonised with W7 (stats SF3):** a candidate passes TG on training (nested LOLO) if all three hold:
1. mean Δ̄_lin vs the current reference ≥ 0.005;
2. the paired line-bootstrap CI lower bound is > 0;
3. Δ > 0 in ≥ 9 of the 12 LOLO folds. A fold with < 5 estimable lines counts as "not positive"; it is never dropped, so the denominator stays 12 (stats r2). For gates computed on the 10-lineage Stage-A panel (3 lines per lineage), the bar is **≥ 8 of 10** folds, and a panel fold counts as "not positive" if fewer than 2 of its 3 lines are estimable; the < 5-line rule applies to the 12 training-LOLO folds only (ext-PI r2 S3; r3 N-1).

**Shortlist rule.** At most 3 candidates go to validation. Choice there follows the **paired-difference one-SE rule** (stats SF2): among candidates within one paired-difference SE of the best, pick the simplest in this pre-declared order: generic < medium < tasks < curation < extraction < combination. Fewer fitted parameters wins first.

**Reporting.** The number of candidates screened per workstream is logged. The expected winner's-curse optimism σ·√(2 ln m) is reported for each validation decision.

### W1 — Baseline benchmark
- **Hypotheses.**
  - H1.1: lift > 1 for generic H1 and H2 (mAP_line and pooled).
  - H1.2: H2 > H1. Luo 2026 claim this for DepMap 19Q3 in SI Fig. S2D (cited in the main text, p. 3); main-text Fig. 2D is the Hart-2015 version series.
  - H1.3: B2 > generic FBA.
- **Inputs.** The baseline KO parquets (149 / 154 lethal [dn:baseline_ko.*.lethal_lt_0.01]; 24 / 32 distinct ratio_r [dn:baseline_ko.*.distinct_ratio_r]); v1.3.0 KO; labels; split.
- **Method.**
  - At M2 the full grid is computed on training and validation.
  - The **full-cohort benchmark appears only at T1**, and its per-lineage and per-gene parts stay sealed until T2.
  - **All error inspection before T1 is on training lines only** (AI SF9).
- **Outputs.**
  - The §6.8 binary endpoints (Human1- and Human2-paper protocols) for generic H1 and H2, on training and validation at M2.
  - A per-gene table (ext-PI SF1): generic call, fraction of training lines positive, positive predictive value of called genes, miss rate of uncalled genes, and error class.
  - The addressable-set mAP and recall ceiling (§5.4).
- **Failure-mode tests (label-free unless stated).**
  - H1 TYMS diagnostic.
  - A pFBA trace of H2 β-oxidation lethality (ACADL, HADHA, ECHS1 have ratio 0 in H2 [dn:baseline_ko.human2.panel_ratio_r]).
  - A trace of the biomass-artefact lethals CYP4F2, BBOX1 and TMLHE.
  - On training only: errors stratified by CN, expression, library and CE.
  - v1.3.0 vs v1.19.0 on U∩.
- **Decision.** W1 fixes G\*. If the H2 lift on training is ≤ 1, a **label-free audit** follows (Entrez join, NaN counts, label orientation). It is run by an agent that has not seen the metric, and the pre- and post-audit numbers are both reported (AI SF2).

### W2 — Medium and exchange bounds (PI decisions 23, 27)

**Medium source tiers.**

| Tier | Source | Lines |
|---|---|---|
| `screen` | QC-passing combined screens → ScreenSequenceMap.ModelConditionID → ModelCondition text (the column named `FormulationID`; the column names are swapped, which a content assertion checks). `PassesQC` is applied at replicate-row level, so a screen with one failed replicate still contributes through its passing rows. If several texts share one base medium, the modal text is used (0 such models today [dn:medium.screen_multi_supplement]). | 1,225 [dn:medium.tier_counts.screen] |
| `screen_mixed` | more than one base medium across screens | 6 [dn:medium.tier_counts.screen_mixed] → unmappable stratum |
| `onboarded` | Model.csv text, used when there is no screen medium | 56 [dn:medium.tier_counts.onboarded] |
| `imputed` | RPMI + 10% FBS, flagged | 74 [dn:medium.tier_counts.imputed] |

- **Where the medium information is missing.** Of the 118 KY-only lines, 115 lack a screen medium [dn:medium.kyonly_without_screen_medium].
- **Screen vs onboarded.** Of the 1,220 lines that have both a screen medium and an onboarded medium, 57 disagree [dn:medium.screen_vs_onboarded_text_disagree, .screen_and_onboarded_both_present].
- **Base classes** [dn:medium.base_class_counts]:
  - **sheet-mappable in principle:** RPMI, DMEM, DMEM:F12, MEM/EMEM, IMDM, F12, McCoy's 5A, L-15. In the screen tier: 1,018 lines, or 1,034 including two-sheet mixtures [dn:medium.screen_tier_sheet_mappable, .screen_tier_sheet_or_mixture].
  - **Sheets actually on disk** (PI decision r2-9): RPMI, DMEM, DMEM/F12 and IMDM [dn:medium.sheet_bases_on_disk]. MEM/EMEM, F12, McCoy's 5A and L-15 are missing [dn:medium.sheet_bases_missing]; the user has been asked for them.
  - **Medium-attribution set** = screen-tier lines whose base-medium sheet is on disk **at M3**. Today that is 835 lines: RPMI 557, DMEM 180, DMEM/F12 60, IMDM 38 [dn:medium.attribution_set_screen_tier_on_disk, .attribution_set_by_base]. It is recomputed by `design_numbers.py` at M3 and frozen.
  - **Fallback lines.** Lines whose base is sheet-mappable but whose sheet is missing at M3 (198 today, all tiers [dn:medium.fallback_lines_sheet_missing_any_tier]) are flagged `fallback_sheet_missing`. They get the imputed medium (RPMI + 10% FBS) in the integrated pipeline, with M0 as a sensitivity, and are excluded from the attribution set.
  - **unmappable:** organoid, stem-cell and other defined media. 195 lines [dn:medium.unmappable_organoid_stem_other] form their own stratum (PI decision 23).
- **Parser rules** (data SF7, domain S2):
  - mixtures are linear blends of concentrations at the stated ratio;
  - supplements are parsed into components: NEAA, pyruvate, GlutaMAX → Gln, 2-mercaptoethanol → cysteine allowed, ITS/B27/N2, HEPES, hydrocortisone;
  - the serum % is parsed from the text, and "no FBS" means serum-free (160 lines [dn:medium.no_fbs_in_text]);
  - DMEM with unstated glucose or pyruvate gets the default high glucose with no pyruvate, plus a with-pyruvate sensitivity variant, flagged;
  - every distinct string has a unit test, and the MF-ID family cross-checks the text.
- **Compositions.** Formulation sheets go into `data/media/compositions.tsv`, and `exchange_map.tsv` maps components to exchange IDs for both parents. Both are produced by **two independent agents**, script-diffed, and every discrepancy is resolved against the source page or row (AI SF3). They are frozen in FROZEN.tsv before the first W2 scoring.

**Variant ladder (a priori; one factor changes per step; same genes and lines throughout).**

| ID | Bounds | Factor |
|---|---|---|
| M-open | all uptakes 1000 | negative control |
| M0 | Ham's GR inputs at 1000 | baseline |
| M1 | documented base medium and supplements, availability only (present = 1000, absent = 0); **cystine instead of cysteine** (cysteine closed unless 2-ME is present); plus the **biomass micronutrient floor** | base medium |
| M2(f) | M1 plus the serum model: each **synthesisable** serum component bounded at **f × req_c**, f ∈ {0.001, 0.01, 0.1, 1} (f chosen by TG); serum-free media get no serum but do get their parsed supplements | serum |
| M3 | a calibrated quantitative regime: glucose, glutamine and O2 bounded at literature-anchored rates (values and citations fixed in `data/media/calibration.tsv` before M3), lactate/NH4 secretion free, the rest as in the best M2 | quantitative bounds |
| **+g (toggle)** | the **transporter gate**, applied on top of any of M1, M2(f), M3 (giving M1g, M2g, M3g): Asp and Glu uptake = 0 unless SLC1A1, SLC1A2 or SLC1A3 has log2(TPM+1) ≥ 1 in the line (domain r1 M2). One boolean per line, so 2 states per medium. **The gate uses line expression, so it counts as a *context* lever, not medium** (PI decision r2-8). | context (expression-based) |
| M-HPLM, M-Plasmax | the HPLM / Plasmax composition plus the M2 serum model, applied to all lines | physiological counterfactual |

- **Biomass micronutrient floor** (compute r1 MF1; PI decision r2-7). This is the set of biomass-required inputs the parent **cannot synthesise** and that standard media lack.
  - It is derived at M3 by `scripts/medium_floor.py` (drop-one-input on each parent, with the result frozen). The compute probe expects linoleate, linolenate, α/γ-tocopherol, retinoate, lipoate and cobalamin.
  - Floor components are bounded at a **fixed, non-limiting level that is not scaled by f**: 10 × req_c. The precheck asserts that no floor component is binding: WT growth with the floor at 10 × req equals growth with the floor at 1000, within 1e-6 relative. If it does not, the bound is raised to 100 × req and logged.
  - The floor is supplied in every medium, including serum-free ones (domain and compute r2, R2-5). Where a supplement documents a floor component (B27/N2: tocopherol, retinyl acetate, linoleic and linolenic acid; DMEM/F12: lipoate, B12), its source is that supplement; otherwise it is flagged `floor_source = modelling`.
  - **Cholesterol is never in the floor**, because the serum-free sterol control depends on it.
  - Floor-free sensitivity: B-loo of the floor components (W3), which removes the requirement rather than supplying it.
- **Serum model.**
  - Synthesisable components scaled by f: palmitate, stearate, oleate, cholesterol, carnitine, ethanolamine. Linoleate and linolenate are floor components (not f-scaled).
  - Selenium and trace metals are at availability bound 1000 (not limiting for biomass).
  - **No nucleosides and no hypoxanthine in serum** (domain M1). Nucleosides enter only through documented compositions (F12-based media, HPLM).
  - The list and the FBS reference concentrations, with citations, are in `serum_model.tsv`, reviewed by a domain agent and frozen.
- **Requirement req_c** (PI decision r2-7; domain r2 N1, compute r2 NMF1). It is computed once per parent on a stated **reference medium**: RPMI-1640 (the modal base, sheet on disk), plus the floor open at 1000, every other input as in that sheet.
  - **True auxotrophs** (floor components): req_c = FVA-min uptake of c at WT growth.
  - **Synthesisable serum components:** req_c = the uptake of c in a pFBA solution at ≥ 99% of WT growth **with c open at 1000**, i.e. the amount biomass takes up when c is free. Then f is applied.
  - Each component's rule is recorded in `serum_model.tsv`. A unit test fails if any serum component has req_c = 0 unless it is flagged.
- **Prechecks (gates).**
  - WT growth > 1e-3 × the M0 growth on every (medium × parent) before any sweep. A table is produced; a failure blocks that medium and is logged.
  - M3 growth calibration: WT growth must fall within 0.01–0.06 h⁻¹ (domain M1d). If it does not, M3 is reported as uncalibrated and excluded from selection.
- **Computation.** Generic sweeps run per distinct (medium × variant × gate state × parent), not per line. The gate has 2 states (compute r2 SF-f).
- **Metrics.**
  - Paired Δ̄_lin, mΔAP_gene and pooled AP for each ladder step, on the **medium-attribution set** (above). The onboarded tier is a sensitivity analysis; imputed, fallback and unmappable lines are a separate stratum.
  - **Attribution**, all computed on identical sets:
    - (i) **"Medium alone"** = ΔtAP between M0 and the best **gate-free** documented variant (M1, M2(f) or M3) on the generic parent (PI decision r2-8). The gate never enters this number;
    - (ii) the fraction of generic binary FP/FN calls flipped to correct, net of correct calls flipped to wrong;
    - (iii) a **training-only medium oracle** (ext-PI SF5): per line, the best ladder variant chosen with its own labels. Reported as (doc − M0)/(oracle − M0).
  - **Label-free flip analysis** (domain O1): which generic calls any single nutrient addition or removal can flip. Expected to be tens of genes; a small attributable share is pre-registered as a plausible outcome, not a failure.
- **Controls.**
  - *Positive, label-free:* the medium-dependent panel of `notes/apriori_lists_approved.md` §1.2, with the expected conditions stated there (e.g. ASNS, the sterol branch, GOT1/GOT2 under the gate, SLC7A11 under cystine).
  - *Positive, training labels:* ASNS by medium (Lagziel replication).
  - *Mevalonate core* (the genes listed in `notes/apriori_lists_approved.md` §1.1): expected essential in every medium. It is reported as a **diagnostic panel at every W2 step**. Ubiquinol is already a biomass component in H2, so HMGCR, MVK, PMVK and MVD should stay essential when serum cholesterol is supplied; a loss of essentiality under serum flags a modelling error (HMGCS1, IDI1, FDPS and GGPS1 are non-essential already on M0, for isozyme/GPR reasons; approved lists §1.1, §5). The W2 choice is re-checked in the M6 factorial, which crosses W2 with W3 B-add (domain r2 S-a).
  - *Note on ASNS:* asparagine uptake is not gated, so the ASNS control tests medium only (domain r2 S-c).
  - *Negative:* **shuffled medium** within lineage × serum-free status, SEEDS[0:5]. It needs 0 CPU-h, because it is a relabelling of cached sweeps. It has power only in lineages with ≥ 8 attribution-set lines in each of two media: **9 lineages on the four sheets on disk** [dn:medium.n_lineages_with_2_media_ge8_on_disk_sheets, .lineages_with_2_media_ge8_on_disk_sheets] (13 if all eight sheets arrive [dn:medium.n_lineages_with_2_media_ge8]). It is reported pooled and per eligible lineage. M-open is also a negative control.
- **Physiological variants.** They are not expected to improve agreement with labels measured in culture. Reported:
  - Δ vs the documented medium;
  - the list of physiology-sensitive genes;
  - the expected-direction checks: purines lose essentiality in HPLM; pyruvate relaxes GOT/MDH and OXPHOS essentiality;
  - an audit of exchange coverage for HPLM/Plasmax components (domain S3).
- **Decision.** The best TG-passing ladder variant goes to the shortlist. If none passes, W2 is negative and M0 is kept. The attribution is reported regardless.

### W3 — Objectives and metabolic tasks (PI decisions 30, 31)
- **A priori sets.**
  - The Essential task list: 57 tasks (BS 18, SU 17, IC 12, ER 9, GR 1) in 11 functional groups [dn:tasks.essential_n_tasks, .essential_by_category, .essential_n_groups].
  - The CellFie consensus list: 195 tasks [dn:tasks.cellfie_n_tasks].
- **Capacity matrix (compute SF3).** For each (parent, selected medium), compute once the gene × task matrix of max output capacity after KO, using deduplicated KO sets. Every criterion mode is then post-processing.
- **Task scoring.** It reuses `taskEssentialGenes.py` and `raven_toolbox.tasks.check`, with an equality test on the generic parents.
- **Candidate sets and controls** are those of `notes/apriori_lists_approved.md` §5, used verbatim.

| ID | Mode | Definition |
|---|---|---|
| T0 | criterion | biomass only |
| T1 | criterion | biomass + ER (viability) |
| T2 | criterion | biomass + ER + IC + BS |
| T3 | criterion | biomass + all 56 non-GR tasks |
| T4_k (k = 1…11) | criterion | biomass + group k, one at a time |
| T5 | criterion, line-specific | biomass + tasks called active by a **CellFie-inspired** rule (local min-max-mean thresholds at the 25th/75th percentiles from fitting-line expression). Renamed because it is not validated against CellFie output (prior-work SF8). |
| D-BS, D-ER | constraint | **product demands only** (compute MF6). Each BS product, and ER as ATP and NAD(P)H demand, gets a demand reaction with lb = 1% of its max under the same medium with biomass free. Demands are added jointly; any that is infeasible with WT is dropped and logged. Only the continuous score is used, because the binary calls nearly coincide with the criterion modes. |
| O-cap | objective | e = 1 − min(ratio_biomass, min_t capacity ratio), over the T3 tasks |
| O-NGAM_k, k ∈ {0, 0.1, 0.3} | objective | ATP-maintenance flux fixed at k × the **maximum** MAR03964 flux on the medium (computed once; the reaction ID is verified in both parents) |
| O-MOMA | objective | QP MOMA (unique optimum; needs Gurobi, so it takes Gurobi slots) built once with `add_moma` per (model, pFBA-WT reference), with KO bounds toggled. Linear MOMA with a lexicographic biomass step is the fallback. **Dropped if the M3 pilot exceeds 1 CPU-h per model.** |
| B-loo_j | biomass | leave one component out. The components are the direct reactants of MAR13082 except ATP and H2O, plus the cofactor-pool groups linked to the floor or to artefacts: retinoid and vitamin-D catabolites, tocopherol catabolites, carnitine, heme, ubiquinol, lipoate and cobamide. They are enumerated per parent by script at M3 and frozen (V-13; approved lists §5). |
| B-add_c | biomass | add one demand at a time, and then all together, at a fixed coefficient of 1e-4 mmol/gDW per unit biomass (1e-3 as sensitivity). The demand list is `notes/apriori_lists_approved.md` §5 (B-add row), verbatim. **Any component already present in the parent's biomass is skipped for that parent and logged** (V-14). Metabolite IDs per parent are frozen at M3. |

- **Composition scaling** (lipid ±50% etc.) is **rejected**: it changes coefficients, which LITERATURE §5 says barely affect essentiality, whereas B-loo and B-add cover presence.
- **Metrics.**
  - Paired Δ vs T0.
  - A **per-task table** for all 57 tasks (ext-PI SF8): the KO-failing genes (label-free) and the precision and AP of those calls. Reported on training at M4 and on test, sealed, at T1. Exploratory.
- **Context claims (stats SF4).** The contexts are line-level: base medium, serum-free status and background group. A claim requires BH q < 0.05 on training LOLO **and** the same sign on validation lines of the same context. Lineage contexts are reported only as within-lineage exploratory results.
- **Controls.**
  - *Negative:* random pseudo-tasks, where outputs are swapped for random producible metabolites of the same compartment (SEEDS[0:5]).
  - *Negative:* matched call rate. The same number of extra calls is assigned to the genes with the highest B2.
  - *Positive:* ETF capability-only false positives under T3 (label-free).
  - *Positive:* the mevalonate-core control as restated in `notes/apriori_lists_approved.md` §5: "**no loss of mevalonate-core essentiality** when serum cholesterol is added under B-add".
- **Decision.** At most 2 candidates that pass TG and beat both negative controls go to the shortlist.

### W4 — Context-specific models (PI decisions 28, 29, 32)
- **Genomic constraint** (applies to every family; PI decision r2-11):
  - **Loss rule.** A gene is lost if it is:
    - homozygously damaged (damaging-matrix code 2); or
    - heterozygously damaged (code 1) **and** shows LOH (proxied by relative CN < 0.5 until an LOH file exists) or TPM < 1; or
    - homozygously deleted (relative CN < 0.3, with 0.1 and 0.5 as sensitivities).
  - **Codes.** The damaging matrix holds DepMap's per-gene count of damaging alleles: 0 none, 1 heterozygous, 2 homozygous or two hits [dn:mutation_matrix_values] (data r2 N2). This is re-verified from the DepMap documentation at M1.
  - **LOH proxy.** No LOH file is on disk, so LOH is proxied by relative CN < 0.5 at the gene. A DepMap LOH call replaces the proxy if one is added to DATA.md before M1.
  - **Old rule as sensitivity.** The v1 rule (any damaging allele, code ≥ 1) is kept as a sensitivity analysis. Losses in screened lines, new rule vs v1 rule: ME2 23 vs 27, KEAP1 20 vs 27, STK11 52 vs 59, FH 0 vs 1 [dn:loss_rule_primary_counts, dn:loss_rule_sensitivity_v1_rule_counts].
  - **Order-independent drop.**
    1. Test each loss **alone** on the parent with the line's medium.
    2. Drop every loss that alone kills WT growth, and log it as a *model-contradicted loss*. These form an error-analysis output (domain r1 M4).
    3. Apply the remaining losses jointly.
    4. If the joint set still kills growth, drop losses one at a time, in descending order of the lost gene's expression (ties broken by ascending Entrez ID), until growth is restored, and log them.
  - Lines without CN use mutations only (flagged).
  - Ploidy sensitivity: absolute CN ≈ relative CN × Cell Model Passports ploidy, with a homozygous deletion if < 0.5.
- **Repair rule** (compute MF2; fixed now; applies to F0, F2, F3, F3b and F5 and after the genomic constraint):
  1. Solve one L1 LP that minimises Σ_r w_r·|v_r| over the closed reactions, with w_r = 1 + the expression deficit of r, subject to WT biomass ≥ 1% of the generic WT on the line's medium.
  2. Viability tasks are closed-system tests with their own inputs, so they cannot share that LP (compute r2 SF-a). Solve **one L1 LP per viability task** (GR and the 9 ER tasks) with the same weights, subject to that task being feasible.
  3. Re-open the **union** of reactions with |v_r| > 1e-6 across these LPs. The compute probe re-opened 1–2% of closed reactions in < 0.5 s per line.

  Reactions closed by an accepted loss are never re-opened. The re-opened sets are logged per line, and the repair size is reported per family.
- **GPR → score convention** (pre-registered; prior-work SF6):
  - presence scoring (F0, F2, F3, F5): AND = min, OR = max;
  - E-Flux capacity (F4): AND = min, OR = sum;
  - CORDA (F3b) uses the same presence convention for its confidence levels.

**Families.**

| ID | Family | Implementation and specification |
|---|---|---|
| F0 | GPR pruning plus repair (naive reference) | own |
| (F1) | ftINIT: **excluded by user decision D2** (in the user's own analysis it underperformed the generic model). The label F1 is retired, not reused. | — |
| F2 | GIMME: required function 0.9 × max biomass; cutoff per threshold rule | `troppo` if it installs and passes toy tests at M3 (pinned); otherwise own |
| F3 | FASTCORE: core = reactions scoring above threshold; consistency via cobra `fastcc`. F3b = CORDA (`corda` 0.5.1), the only PROMPT-named established method still in scope. Its panel scope is set by its M3 pilot cost (row A6) | own / `corda` |
| F4 | **E-Flux** (Colijn 2009): internal bounds scaled by GPR capacity, with a floor of 10; GPR-less reactions and exchanges keep their default or medium bounds | own |
| F5 | **Label-calibrated presence (LCP)**, tagged "label-calibrated" in every table: gene off if w_r·z_RNA + w_p·z_prot + w_c·loss < τ. A 24-point grid (Stage-A row A5) is fitted by **nested** LOLO on the panel. Excluded from the variance panel. | own |

- **Thresholds:**
  - R1: global, TPM 1;
  - R2: global 25th percentile;
  - R3: local, the gene's mean over **training lines**, bounded to the 25th–75th percentiles. R3 is a label-free statistic, so it needs no LOLO.
- **Modalities:**
  - RNA;
  - protein: ProCan primary, Nusinow as sensitivity;
  - RNA + protein: protein overrides RNA where quantified;
  - each with and without the genomic constraint.

  Protein **absence** counts as evidence only for soluble enzymes; transporters use RNA only (domain S6). The protein arm is pooled across lineages with a lineage covariate, and it is descriptive on test (data SF10).

**Staged run matrix (PI decisions r1-29, r2-12; LP-only under D2).**
- **Stage A: 30-line training panel.** The panel takes 10 training lineages × 3 lines, sampled with SEEDS[0]; if a lineage has fewer lines, the gap is filled from the largest training lineages. n_p = panel lines with proteomics, counted at M1; the costs below assume n_p ≤ 20.
  - **Stage-A build table (frozen).** Unit cost: LP build + repair + KO sweep ≈ 1.5 CPU-min per line (compute r2 probe, with the vectorised closure map).

| Row | Configuration | Line set | Seeds | Builds | Expected CPU-h |
|---|---|---|---|---|---|
| A1 | F0, F2, F3, F4 × {R1, R3} × {real, permuted expression}, RNA (variance partition) | panel (30) | SEEDS[0] | 480 | 12 |
| A2 | F0, F2, F3, F4 × R2 × real RNA | panel (30) | — | 120 | 3 |
| A3 | F0, F2, F3, F4 × {R1, R2, R3} × {protein, RNA+protein} | panel lines with proteomics (≤ 20) | — | ≤ 480 | 12 |
| A4 | ± genomic constraint, applied **post-extraction** (gene KO + repair + sweep), never as a rebuild, to every real-expression build of A1–A3 | as A1–A3 | — | ≤ 840 re-sweeps | 21 |
| A5 | F5 LCP: 24-point grid (coarse 12 + local refinement 12); nested panel LOLO reuses these label-free builds | panel (30) | — | 720 | 18 |
| A6 | F3b CORDA × {R1, R3} × {real, permuted}, RNA, if its M3 pilot is ≤ 0.25 CPU-h per model. If the pilot is in (0.25, 1] CPU-h per model: CORDA × R1 × {real, permuted} on 15 panel lines. Above 1 CPU-h per model it is dropped and the drop is reported (V-9; verification Part 5 item 9). | panel (30), or 15 lines | SEEDS[0] | ≤ 120 | ≤ 30 |
| A7 | Stability: F2 (GIMME, alternate optima via random objective perturbation) and F3 (FASTCORE, expression-noise perturbation) | 20 panel lines | SEEDS[0:3] | 120 | 3 |
| | **Total Stage A** | | | **≈ 2,880 LP builds** | **≈ 100 (0 MILP)** |

  - Selection keeps ≤ 2 thresholds per family and the modality. Label-blind gates come first (growth, viability tasks, repair size ≤ 5% of reactions), then nested panel LOLO with the ≥ 8-of-10 bar.
- **Stage B: all lines.** At most 3 LP configurations per line. Builds for validation and test lines are made only for configurations that passed TG (ext-PI SF6).
- **Stage C: group models.** About 150 models × ≤ 4 LP configurations:
  - per lineage, from the median expression of fitting lines; a test lineage's model uses its own lines' inputs;
  - per genomic background;
  - lineage × background, where ≥ 8 lines exist.

**Genomic backgrounds (PI decision 32).** The groups, their expected dependencies, their status (approve, modify, drop, in-silico only) and the negative controls (the **primary** control is size-matched random within-lineage member sets, SEEDS[0:5]; KRAS_hotspot is a **secondary**, non-metabolic-driver control) are those of `notes/apriori_lists_approved.md` §2, used verbatim. The counts below are label-free [dn:backgrounds].

| Group | Rule | Screened lines | Lineages eligible for within-lineage contrast |
|---|---|---|---|
| MTAP_del | relCN < 0.3 | 233 | 19 |
| ASS1_silenced | log2(TPM+1) < 1 | 161 | 14 |
| NRF2_active | KEAP1 loss (r2-11 rule) or NFE2L2 hotspot | 43 | 5 |
| STK11_loss | loss (r2-11 rule) | 52 | 4 |
| NAPRT_silenced | log2(TPM+1) < 1 | 30 | 4 |
| FH_SDH_loss | loss (r2-11 rule) | 2 | 0 |
| PHGDH_amp | relCN ≥ 2 | 30 | 5 |
| ME2_loss | loss (r2-11 rule) | 23 | 3 |
| ENO1_del | relCN < 0.3 | 0 | 0 |
| IDH_mut | IDH1/2 hotspot | 8 | 0 |
| KRAS_hotspot | hotspot | 209 | 11 |

Notes on the table:
- Counts are from [dn:backgrounds.\*.n_screened_members] and [dn:backgrounds.\*.within_lineage_eligible_lineages]. "Eligible" means ≥3 members and ≥3 non-members in the lineage.
- **Minimum n (a priori).** A group is evaluated in a tier only if ≥2 lineages in that tier are eligible and those lineages contain ≥10 members.
- **Evaluation.** The contrast is pooled within lineages (lineage as stratum): the group model vs the same-lineage model on member lines. Two metrics are used: ΔtAP on the pre-registered expected-dependency genes, and Δ̄_lin.
- **Feasibility on test.** The share of simulated accepted draws whose test tier can evaluate the group is given in [dn:split_simulation.background_test_feasibility_share]: MTAP 1.0, ASS1 0.995, KRAS 0.98, NRF2 0.70, STK11 0.61, NAPRT 0.49, PHGDH 0.47, ME2 0.29, FH/SDH 0. The evaluable list is fixed at M1 from the real split.
- **All background contrasts are estimates (exploratory)** and are never part of the confirmatory family. Groups below 50 lines (NRF2, ME2, NAPRT, PHGDH) are flagged as small (data r2 N4).

**Method-vs-tissue test (pre-registered).**
- *Builds:* F0, F2, F3, F4 (+ F3b CORDA per row A6) × {R1, R3} × {real expression, permuted expression}, where each gene's values are permuted across lines, on the 30-line panel (rows A1, A6). These are **≥ 3 established non-ftINIT families** (GIMME, FASTCORE, E-Flux, and CORDA where kept) plus F0 as the naive reference, which satisfies D2.
- *Power:* with 10 lineages × 3 lines, the lineage R² is estimated with limited power. A null H4.2 result is reported as "not detected at this panel size", not as absence.
- *Analysis:* PERMANOVA (999 permutations) with R² for method, threshold and lineage, at three levels:
  - reaction-content Jaccard;
  - the 57-task capability profile;
  - the binary KO profile.

  The same partition on the input expression is the reference.
- *Hypotheses:*
  - **H4.1:** R²(method) > R²(lineage), tested by a line bootstrap of the difference.
  - **H4.2:** within each method, R²(lineage | real) > the 95th percentile of R²(lineage | permuted). A method that fails is labelled "method-dominated".
  - **H4.3:** H4.1 and H4.2 on KO profiles.
- *Positive controls:* TYR/DCT in Skin models; urea cycle absent from most lines.
- *Stability:* row A7 (GIMME alternate optima; FASTCORE under expression noise; 20 lines × SEEDS[0:3]).

**Metric and decision.**
- TG is computed against the best shortlisted generic pipeline. TG reports must include the **expressed-pairs** Δ (ext-PI MF4) and mΔAP_gene vs B3.
- At most 2 configurations go to the shortlist. The best configuration is carried to M6 even if it fails TG, and is flagged.

### W5 — GPR curation
- **Targets (label-free).** For each reaction in H2 (and H1), enumerate every single edit: OR→AND, AND→OR, or removal of one OR branch.
  - Rank reactions by the number of genes whose binary call flips.
  - Curate the **top 200**.
  - Re-evaluate only the genes in the edited GPR, deduplicated (compute SF4).
- **Edits are generated by script** from local CORUM, Complex Portal, Ensembl and proteomics files. A validator checks that every evidence ID exists locally (AI SF3).

| Edit | Required evidence (all label-free) |
|---|---|
| OR → AND | the branch genes are **not** Ensembl paralogs of each other **and** are co-members of a CORUM or Complex Portal complex **and** have protein co-abundance r ≥ 0.5 over ≥ 30 lines (ProCan, then Nusinow) |
| AND → OR | the genes are paralogs **and** are not co-members of any complex |
| Remove OR branch | the gene has a UniProt function unrelated to the reaction (accession and entry version cited from a downloaded UniProt file added to DATA.md); or the gene is a **soluble enzyme** never detected in either proteomics set **and** has TPM < 1 in ≥ 95% of lines ("pan-cancer absent", tagged) |
| Never edited automatically | the isozyme-exception list of `notes/apriori_lists_approved.md` §3, used verbatim, including its notes on RRM2/RRM2B and on tissue-restricted isoforms |
| Optimised H1 only | MAR04644 made irreversible (TYMS) |

- **Arms.**
  - W5a uses label-free evidence only.
  - W5b adds training-fold co-essentiality: Chronos r ≥ 0.4 between OR-branch genes supports OR → AND. It is recomputed inside every fold.
  - The De Kegel 2021 and Ito 2021 tables are **denylisted** as evidence and serve only as post-hoc annotation (PI decision 16).
- **"Paralog-buffered"** (§10.3) requires the partner to be in the same compartment and catalyse the same reaction.
- **Ablation.** A 2 × 2 factorial on identical sets: {generic, curated GPR} × {generic, best constraint pipeline}.
  - Curated GPRs are applied to the **already-extracted** reaction sets.
  - Re-extracting with curated scores is a P1-only sensitivity analysis, to avoid alternate-optima noise.
  - The analysis is broken down by edit class.
- **Controls.**
  - *Negative:* random edits of matched number and type among reactions without evidence (SEEDS[0:5]).
  - *Positive:* recall of the v1.19 → v2.1.1 GPR changes on the targeted reactions. This control is **possibly label-informed**, because Luo et al. benchmarked on DepMap and Project Score and whether their GPR edits were tuned on DepMap is unknown (the SI has not been read). It is reported with that caveat (PI decision 16).
- **Review.** A non-author agent **without label access** checks a random 20% of edits and every edit that flips ≥ 5 genes.
- **Decision.** Curated GPRs enter the shortlist if they pass TG on the generic parent. Δ_curation is reported regardless.

### W6 — Information in flux states beyond the inputs
- **Features from P1** (frozen at M6):
  - **X_KO**: the e-vector;
  - **X_KOflux** (compute SF6): the pFBA flux of KO models for the curated drug-target genes, normalised by biomass flux, then PCA;
  - **X_flux**: WT pFBA flux, normalised by biomass flux;
  - **X_rxn**: reaction presence.
- **Comparators:**
  - X_RNA and X_prot;
  - **X_rxnRNA / X_rxnProt**: GPR-mapped reaction expression without simulation (ext-PI SF4);
  - the combinations omics + flux.
- **Assay medium (data SF9, domain S7).** Flux features for each dataset are computed in that assay's medium:
  - PRISM: the assay's own medium (recorded from the PRISM documentation at W6 start);
  - GDSC: its "Screen Medium" column (RPMI or DMEM/F12, both on disk [dn:coverage_screened.gdsc_screen_medium_values]);
  - metabolomics: the line's culture medium.
- **Matched capacity.** Every block uses PCA-20 (PCA-50 as sensitivity), ridge regression, and identical folds (LOLO, or 5-fold within lineage for subtypes) with SEEDS[0:3]. Lineage is a covariate.
- **Analyses:**
  1. Subtypes. Within-lineage OncotreeSubtype and `OmicsInferredMolecularSubtypes`, for lineages with ≥ 2 subtypes of ≥ 5 lines (list generated by script before fitting). Metrics: silhouette and kNN macro-F1.
  2. Metabolic drugs. The compound list is **curated a priori**: compounds whose primary mechanism targets a Human-GEM enzyme or transporter; the rules, classes and compounds are those of `notes/apriori_lists_approved.md` §4, used verbatim (kinase inhibitors excluded; prodrug flag and activating genes; ≥ 100 lines with features, checked by script).
     - Response data: PRISM, GDSC2 and CTRP; a compound needs ≥ 100 lines with features.
     - Metric: held-out Spearman ρ. A fit-free test is also run: within-lineage ρ between the target's line-specific e and the response, compared with the target's expression.
     - Positive-control pairs: those of `notes/apriori_lists_approved.md` §4.
  3. Metabolomics. The 225 CCLE metabolites, scored by held-out ρ.
- **Controls:**
  - permuted-flux features (shuffled within lineage);
  - a synthetic positive control (injected nonlinear information is recovered);
  - a synthetic negative control (flux as a linear function of expression).
- **Decision.**
  - "Flux adds information": median Δρ(omics + flux − omics) > 0 by one Wilcoxon test per analysis (Holm across the 3 analyses), above the permuted control.
  - "Simulation adds beyond network mapping": the same criteria for X_flux/X_KO vs X_rxnRNA.
  - Per-target BH results are exploratory. A null result is reported as such.

### W7 — Auto-research loop
- **Feature whitelist (PI decision 16):**
  - expression, CN, mutation, proteomics, medium class / serum;
  - GEM outputs (from any pipeline);
  - within-fold B2 and B2m;
  - neighbour and paralog features.

  Not allowed: drug screens, any external essentiality set (Project Score, CoRe, Hart CEG/NEG, DEMETER2, older DepMap releases, DepMap common-essential lists), and full-matrix co-essentiality.
- **Seed ideas:**
  - H-a: B6 + GEM features;
  - H-b = B6 (ablation);
  - H-c: flux-cone features, only on generic and group models, capped at 100 CPU-h (compute SF7);
  - H-d: stacking P1 + B6.
- **Protocol.**
  1. Each idea is logged in `autoresearch/ideas.tsv` before it runs.
  2. Screen by nested 12-fold LOLO with SEEDS[0:3]. An idea is promoted under TG against the current best.
  3. A promoted idea gets one validation look, which reveals a scalar only.
  4. Gene-grouped CV is reported for promoted ideas.
  5. Ideas that name a lineage's biology must name a **training** lineage (AI O1).
  6. **Stop** when any of these happens:
     - 10 consecutive ideas fail promotion;
     - 300 CPU-h have been used (read from `results/budget.tsv`);
     - 24 h remain before the deadline;
     - a fresh agent finds no unlogged promising idea.
  7. **P2** is the best validated idea. If it does not beat the runner-up by more than one paired SE, the simpler of the two is P2 (stats O4).
  8. Ideas tried, validation looks used and the shrinkage from LOLO → validation → test are all reported.

### W8 — Deliverable resource
- **Scope.** P1, the GEM pipeline, is applied to every 26Q3 model that has a default expression profile and a cancer lineage. That is 2,030 models, after excluding 100 that are non-cancer or Fibroblast [dn:resource_scope.resource_models, .excluded_noncancer].
  - 1,277 of these models are screened [dn:resource_scope.resource_models_screened].
  - 753 are not, and they carry `benchmarked = false` in the manifest [dn:resource_scope.resource_models_not_benchmarkable].
- **What is released.**
  - Parameters are those fitted on training ∪ validation, exactly as tested at T1.
  - The build starts **right after T1**, in parallel with W6/W7 (ext-PI SF6).
  - Builds cached by (config hash, input hash) are reused (compute SF8).
  - **Human2-parent models are released first. Human1-parent models follow if the budget allows; otherwise REPORT.md says they were not released and why** (ext-PI r1 SF8).
  - **H1 resource scope.** Under D2, P1 is LP-only, so P1(H1) for all 2,030 models costs ≈ 50 CPU-h (§13). It is built after the H2 release.
- **Format and naming.**
  - SBML L3V1 FBC v2, gzipped, at `resource/v{X.Y.Z}/{ModelID}_{parent}-{version}_{pipelineID}.xml.gz`.
  - Notes record the medium source tier, pipeline ID, git commit, input SHAs and the repaired reactions.
  - `resource/manifest.tsv` lists sha256, size, build time, counts, flags and validation status.
- **Validation gates (every model).**
  - libSBML reports no errors beyond those of the parent.
  - A cobra round-trip gives identical sets and bounds.
  - WT growth matches within 1e-6 relative.
  - 20 random KOs match within 1e-6.
  - The viability tasks pass.
  - MEMOTE (separate venv) and FROG are run on 10 models and reported, not gated.
- **Licence (data SF11).**
  - CC BY 4.0, inherited from Human-GEM.
  - HPLM/Plasmax-variant models ship separately with a NOTICE, because those composition sources are not CC BY.
  - The resource contains no per-protein values.
  - The GDSC and Cell Model Passports licences are re-verified before M8.
- **Size.** The compute probe measured 2.1 MB gzipped for the full Human2, so the resource is ≤ 4.5 GB.

---

## 10. Integration, attribution and test protocol

### 10.1 Selection
1. Shortlists come out of W2, W3, W4 and W5.
2. The **combination factorial** (≤ 24 combinations) is first run on the 30-line panel. The top 3 are run on all training lines with nested LOLO (compute SF5) and then on validation.
3. The paired one-SE rule picks P1.
4. P1 must pass the validation gate (§8). If no candidate passes, P1 is the best point estimate, labelled "gate not passed", and is still tested.
5. **Attribution (scope reduced by the user's compute decision, v3).** The **full 2^4 factorial of P1's levers runs on validation only**: main effects, the curation × constraints interaction, and Shapley shares. The levers are:
   - medium (gate-free);
   - context (extraction, the genomic constraint **and the transporter gate**, PI decision r2-8);
   - tasks;
   - curation.

   On **test**, only these **pre-declared contrasts** are built and evaluated (all sealed at T1):
   - (i) P1 vs G\*;
   - (ii) "medium alone" = M0 → the best gate-free documented variant, on the generic parent;
   - (iii) curation alone = curated vs generic GPRs on the generic parent;
   - (iv) constraints alone = P1 with generic GPRs vs G\*;
   - (v) curation given constraints = P1 vs P1 with generic GPRs.

   Shapley shares are therefore validation-only estimates. On test, the curation-vs-constraints answer rests on (iii)–(v).

### 10.2 Test looks
| Look | When | Revealed at the look | Computed and **sealed** (hash in the ledger; revealed after P2 is frozen) |
|---|---|---|---|
| Rehearsal | before T1 | `evaluate_test.py --rehearsal` runs on synthetic labels of the same shape, then on validation. The validation run **counts as a validation look** and gets a ledger row (AI r2 S6). | — |
| **T1** | M6 | scalars and CIs of C1, C1b, C3 and C2; the B8± and B8′ contrasts; within-parent contrasts; the §6.8 binary endpoints (both paper protocols) for G\* and P1; our fixed-cutoff MCC; Δ̄_line; the lineage t-interval; the k-of-7 count; and **failed-build counts per method and lineage** | per-lineage and per-gene tables; the full W1 benchmark; the per-task table; library strata; the pre-declared test contrasts (i)–(v) of §10.1, including "medium alone" (M0 → best gate-free documented variant on the generic parent; PD r2-8); the failure account |
| **T2** | after the W7 stop | C4, C5 and P2's scalar grid; then all sealed T1 outputs are unsealed | — |

- **Where sealed outputs go** (PI decision r2-16). `evaluate_test.py` writes sealed outputs only to `data/locked/sealed/<look>/`, which is git-ignored and protected by the hook and the tripwire. Only their SHA-256 values go into the ledger and git. The pre-commit guard against test-keyed tables stays on until the **T2 unseal** ledger row exists.
- **Single first number.** The first T1 number stays the headline in REPORT.md's summary (AI SF2).
- **Third look.** One more look is allowed only after a bug fix that is demonstrated by a label-free failing regression test and applied identically to P1 and every baseline. It is reported as post hoc.

### 10.3 Failure account (sealed at T1, reported at the end)
P1's test errors are classified by these classes:
- common essential;
- paralog-buffered (same compartment);
- medium-sensitive;
- CN-amplified;
- unexpressed;
- gene-less bypass;
- missing from the model;
- complex subunit;
- housekeeping-adjacent;
- biomass-artefact lethal;
- model-contradicted loss;
- a pre-registered expected FBA failure (MTAP→MAT2A/PRMT5, V-ATPase, GPX4, PGD).

The 20 worst lineage × gene cells are traced by pFBA. The report states an upper bound on GEM performance: the fraction of positives in U that have a metabolic mechanism.

---

## 11. Leakage control and agent guards (PI decisions 14–20)

### 11.1 Label access
- **Single door.** `src/nextgem/labels.py` is the only module that may open label-bearing files: `CRISPRGeneEffect`, `CRISPRGeneDependency`, `ScreenGeneEffect`, `CRISPRInferredCommonEssentials`, `AchillesCommonEssentialControls` and `AchillesNonessentialControls`.
- **Files written at M1.** It writes `labels_train.parquet`, `labels_val.parquet` and `labels_test.parquet` (the last two in `data/locked/`), and records their hashes in FROZEN.tsv.
- **API.**
  - `labels.fold(k)` gives the training labels minus fold k, and `labels.train()` gives all training labels. These are the only label handles fitting code receives.
  - All CV goes through the `nextgem.cv.lolo` harness, with nested inner LOLO for every choice.
  - Validation labels are reachable only from `scripts/evaluate_val.py`; test labels only from `scripts/evaluate_test.py`.
  - `labels.ids()` exposes ModelIDs and the gene header without values. `design_numbers.py` is refactored to use it once `labels.py` exists.
- **Tripwire** (PI decision r2-15).
  - A `sitecustomize.py` in `.venv` registers it for **every interpreter launch of `.venv`**, whether or not `nextgem` is imported.
  - `sys.addaudithook` covers plain `open` and pandas readers.
  - Wrappers cover the pyarrow readers, which raise no audit event: `pyarrow.parquet.read_table`, `ParquetFile`, `pyarrow.dataset.dataset` and `pyarrow.csv.read_csv`.
  - The `subprocess.Popen` audit event is checked for protected paths in its arguments.
  - Any read of a label file or of `data/locked/` raises `PermissionError` unless **both** hold:
    - the calling frame is in `labels.py`;
    - the entry point (`__main__.__file__`) is `build_labels.py`, `evaluate_val.py` or `evaluate_test.py`, or the call goes through the fold-scoped training API.

    This closes the "call a private reader of `labels.py`" bypass (AI r2 S2). Every denied or allowed event is logged to `logs/label_access.log`.
  - All data paths come from `nextgem.paths`, which honours `NEXTGEM_DATA_ROOT`.
  - **Residual risk, documented.** Interpreters outside `.venv` (system `python3`, `Rscript`, `sqlite3`) do not load `sitecustomize`. CLAUDE.md mandates the venv.
- **PreToolUse hook** (user-approved; now active and fail-closed). It denies direct agent reads of raw label files, `labels_val`, `labels_train`, `data/locked/` and recursive searches over `data/`. **Every new agent session first runs a liveness probe:** a Read of the sentinel `data/locked/HOOK_PROBE` must be *denied*. The result is recorded in RESEARCH_TREE.md, and a session whose probe is allowed stops.
  - Requested hardening of the PI-owned hook config (AI r2 §3, S2): deny `cd` into, and `cp -r`/`rsync`/`tar` of, `data/depmap` or `data/locked`; deny non-venv interpreters on commands that mention `data/` or `.parquet`.
- **Content-based label-copy scanner** (PI decision r2-14): `scripts/scan_label_copies.py`.
  - **Where it looks:** the project tree (`results/`, `notes/`, everything outside `data/locked/`), the session scratchpad roots and `/tmp/claude-*`.
  - **What it flags**, in two ways:
    - (a) **row fingerprints.** At M1, `labels.py` writes salted hashes of rounded label rows to `data/locked/fingerprints`. Any table row matching a fingerprint is a hit.
    - (b) **shape heuristic.** An `ACH-` index or column, ≥ 50 columns matching `SYMBOL (Entrez)`, and either values < −0.3 (Chronos) or a 0/1 matrix keyed by validation or test ModelIDs (AI r2 N2).
  - **When it runs:** at every milestone, from `status.py`, in the pre-commit hook and before every push.
  - **On a hit:** the commit is blocked; a scratch copy is deleted and the deletion logged in the register.
- **Not adopted.**
  - Encrypting the raw label files with a user-held key (PI decision r1-15).
  - Fernet encryption with an agent-readable key (AI r2 S1). The hook, the tripwire and the content scanner (PI decisions r2-14 to r2-16) cover the careless-access paths, and the PI chose that guard set.
- **Leakage canary** (at M4 and M6, on every shortlisted pipeline): `NEXTGEM_DATA_ROOT` points to a copy in which the **validation- and test-line rows of every raw label file** are replaced by noise. All training-tier outputs must stay bitwise identical.
- **Pre-commit hook** (AI r1 SF6, r2 N4) rejects:
  - files under `data/` or `papers/`;
  - files > 50 MB;
  - any hit of the content scanner;
  - any table keyed by test ModelIDs × genes, until the **T2 unseal** ledger row exists.
- **PI action (not editable by the design author).** Add to CLAUDE.md: "Never write label values outside `data/locked/` or `labels_train.parquet`; cache derived statistics, not label matrices."

### 11.2 Label-derived sources (denylist; PI decision 16)
- **Where the denylist is recorded.** DATA.md gains a `label_derived` column, and `config.LABEL_DERIVED_SOURCES` lists the same sources.
- **Restricted sources.** The accessors for these sources require `purpose="reporting"` or `purpose="annotation"`:
  - the De Kegel 2021 and Ito 2021 paralog tables (now on disk);
  - Project Score (Behan 2019 supplement, now in `papers/`);
  - CoRe, Hart CEG/NEG, DEMETER2 and older DepMap releases;
  - every DepMap common-essential or control list;
  - DepMap co-essentiality outside training folds.
- **Drug screens** are W6 targets only, never features of the essentiality predictor.
- **Label-informed by construction** (stated limitations):
  - L-prob and DepMap's common-essential list, which DepMap derives from all lines;
  - Human2's GPR curation, which may have been DepMap-informed (unknown).
- **Hart 2015.** No longer used: the ftINIT acceptance test was removed with D2.

### 11.3 Ledgers, caps and frozen artefacts (PI decisions 17, 18)
- **Ledgers.** `nextgem.ledger.append()` maintains `results/{test_ledger,validation_ledger}.tsv`. The ledgers are hash-chained (each row holds the SHA of the previous row), file-locked, **committed on every append**, and checked by `scripts/verify_ledgers.py` (chain intact; additions only in git history; every output has a row).
  - `ledger.append()` is the **only non-PI commit path**. It uses `git commit --only <ledger>` under a lock file shared with the PI's commits.
  - The `started` row of any look is committed **before** any number is computed. If that commit fails, the look aborts, so no unlogged reveal is possible (AI r2 S3).
- **Look caps in code.**
  - Validation: 40 looks, ≤20 of them for W7. Raising the cap needs a committed `results/authorisations/val_cap_raise_<n>.md`.
  - Test looks are idempotent. Each writes a `started` row, then the output SHA, then a `completed` row. A look ID that already exists is refused. A crashed look resumes only with the same config hash and still counts as a look.
- **Frozen artefacts.**
  - `FROZEN.tsv` covers `metrics.py`, `labels.py`, `splits`, `paths.py`, `cv.py`, the evaluate scripts, the known-answer tests and tolerances, the split and label hashes, and every a priori list (media tables, serum model, floor, calibration, gene panels, task sets, background groups, compound list).
  - `scripts/check_frozen.py` runs as a git pre-commit hook and at the start of every evaluate script. A frozen file may change only through a §17 amendment row that names the new SHA and is committed before first use.
  - `results/test_runs.tsv` logs the counts of passing, failing, xfail and skip tests and the hash of the test-ID list. A known-answer test that disappears or is downgraded fails the check.
  - **Rule:** a failing test is never edited in the same commit as the code fix.

### 11.4 Provenance and execution (PI decisions 19, 20)
- **Provenance wrapper.** Every result-producing script runs under `nextgem.run.provenance()`.
  - It refuses to run when `src/`, `scripts/` or `tests/` are dirty (AI r2 S4). Uncommitted outputs under `results/` do not block it.
  - In pinned-worktree runs it asserts that `nextgem.__file__` lies inside the worktree, so the editable install cannot shadow it.
  - It writes a `MANIFEST.tsv` row with: argv, commit, input SHAs, seeds, solver and version, workers, wall and CPU time, RSS, output SHA and author agent.
- **Long jobs.** A job longer than 1 CPU-h runs from a pinned git worktree (`/data/tmpA/nextgem_runs/<commit>`).
- **Writes.** Writes are atomic (tmp → fsync → rename). A free-space check of ≥ 20 GB runs before any job that writes > 1 GB.
- **Solver infrastructure.**
  - Gurobi slots: a `flock` pool on `/dev/shm` sized to the licence session limit verified at M3. A worker waits for a slot and never falls back silently to another solver.
  - **Workers** (PI decision r2-13 as amended by D2 and r3): there is no MILP lane. Start at **24 LP workers + ≤ 4 Gurobi sessions** (≈ 41 GB; ≈ 1.45 GB per LP worker, measured). All 28 run on Gurobi only if an M3 check confirms that the WLS licence admits 28 concurrent sessions, in which case the `flock` pool is set to 28; otherwise LP workers beyond the licence pool run GLPK, only for families that pass the GLPK = Gurobi equality check and whose M3 GLPK per-line cost is ≤ 2× Gurobi. Raising to 28 LP workers at M3 requires measured RSS to keep ≥ 20 GB free. A job aborts below 10 GB.
  - Licence, environment and network errors raise `SolverInfraError`, which **aborts the job with no output**.
  - Every KO parquet must have `solver` equal to its declared solver. LP sweeps may be declared GLPK (compute SF2), with a per-family check that GLPK and Gurobi give equal ratio_r. QP MOMA always uses Gurobi.
- **Time limits** (PI decision r2-6, reduced to its remaining object under D2). No MILP arm remains. QP MOMA solves that hit their time limit are failed solves for that gene (NaN, counted). The MOMA arm is dropped if more than 0.1% of its pairs fail.
- **Failed build.** Only these count as a failed build:
  - no feasible model after the method's build (infeasibility);
  - WT growth ≤ threshold after the method's own repair rule;
  - a numerical non-optimal LP status that persists after one retry with NumericFocus = 2.

  Licence and solver errors are not failures: they still abort the job (r1 decision 20).
- **No-information score.** A failed build gets the no-information score, applied **identically in every arm**, including the H1 arm.
  - Per line: all genes are tied, so tAP_line = E[AP_random].
  - Per gene (for tAP_gene across lines): the failed line takes, for each gene, the method's **median score for that gene over its non-failed lines in the same tier**. That value is placed as a tie and carries no line-specific information (stats r2 SF-d).
  - It never gets e = 0 for single genes, and never the generic parent's scores.
- **Sensitivities and flags.**
  - Two sensitivities are reported: (a) lines failed by **any** compared method excluded for all methods; (b) generic fallback.
  - Failed-build counts per method and lineage are revealed as T1 scalars.
  - A method with > 5% failed lines is flagged as failing.

### 11.5 Process
- **Single writer.** Only the PI writes RESEARCH_TREE.md, LESSONS.md, DESIGN.md, DATA.md and the checksums, and only the PI commits, never with `git add -A`. Sub-agents write under `results/<ws>/`, and their summaries go to `results/<ws>/agent_summaries/`.
- **Job tracking.** A job registry (`results/jobs.tsv`) and `scripts/status.py` are run after any context reset. Each job declares a CPU budget and is killed at 2× it. Any step estimated above 10 CPU-h first runs a pilot of ≤ 10 lines. `results/budget.tsv` tracks spending. At most 4 sub-agents run at once (AI SF4, SF5).
- **Label-exposure register** (created; see `notes/label_exposure.md`, which is authoritative; AI r1 SF1, r2 S5, V-16). It logs every viewing of label *values* or label-derived summaries. Entries as of v3:
  - the baseline review read the CRISPR header only;
  - the literature agent saw the all-line DepMap common-essential count;
  - the **round-1 statistics review**:
    - evaluated the retired seed 20261026 under the v0 rules;
    - named the old test tier (Pancreas, Soft Tissue, PNS, Prostate, Eye) and validation tier (CNS/Brain, Ovary/Fallopian Tube, Bladder/Urinary Tract, Cervix);
    - computed label summaries on the old tier: positives per line, eligible-gene counts and a prevalence/SD remark;
    - also computed all-line prevalences and the common-essential share of positives;
  - the **round-1 data reviewer** computed ASNS prevalence by medium and per-library positive rates on all lines;
  - the **design author**, while inspecting a file header with `head`, printed about 10 L-prob values for one line (ACH-002730);
  - the **all-line label copy** `ge_U.parquet` (CRISPRGeneEffect restricted to U, all 1,361 lines) and `splits.pkl`, written by a round-1 reviewer in the shared scratchpad and **deleted by the PI on 2026-10-26**. The remaining label-reading scratch scripts (`testtier.py`, `stats1.py`, `prev.py`, `splits_sim.py`, `apnull.py`, `genes2.py`) were deleted later the same day (V-10);
  - the **PI's header-only read** of a label file while testing the hook;
  - the **PI** read reviewer summaries containing all-line aggregate label statistics (no gene × line values);
  - the **round-3 verifier** read the register and the round-2 statistics review, which contain all-line aggregates (no gene × line values, nothing on the new test tier).

  **Consequences.** The retired seed's tiers are known, which is why the seed is replaced by a user-supplied one (§4.1). All other exposures are aggregate (no gene × line information on the new test tier) and method-symmetric.

  **Authorship rule** (PI decision recorded in the register; resolves F-4). Aggregate all-line label statistics do **not** disqualify an agent from authoring or approving a priori lists. **Disqualifying exposure** is gene × line label values, or any label summary on the post-freeze validation or test tiers.
  - The design author (about 10 values for one line) therefore does not author lists.
  - The round-3 verifier (aggregates only) authored `notes/apriori_lists_approved.md`, which DESIGN.md uses verbatim.
  - Whoever writes an a priori list attests their register status in the commit message.
- **Reviews.**
  - Each `REVIEW.md` has front matter (reviewed output SHA, commit, author, reviewer, verdict, independent code path).
  - `scripts/check_reviews.py` blocks a milestone if a cited output lacks a matching ACCEPT, or if the author and reviewer are the same agent.
  - Reviewers receive code and outputs, not conclusions. Opus and Sonnet review each other where possible (AI SF9, O3).
- **Reports and citations.**
  - `results/claims.tsv` and `scripts/verify_report.py` check every number in REPORT.md against the manifest outputs.
  - `scripts/verify_citations.py` resolves every DOI in LITERATURE.md and REPORT.md (AI SF8).
  - Wording rules (§6.4) are applied by the report generator, which also carries the flags "gate not passed", "method-dominated" and "label-calibrated".
- **Audits.**
  - Red-team leakage audit by a fresh agent at M4 and M6 (AI O4).
  - `scripts/reproduce.py --sample 10%` at each milestone (AI O5).

---

## 12. Verification plan (non-author review for every result the PI relies on)

| Component | Known-answer and toy tests | Positive control | Negative control |
|---|---|---|---|
| tAP and metrics | brute-force enumeration of tAP on toy rankings with ties; the O(n_g) form equals the hypergeometric form; constant → E[AP_random] exactly; failed-line per-gene median rule; mean over 1,000 seeds of score + 1e-9·noise equals tAP within MC error; perfect and reversed rankings; the ≥5-positive NA rule; set-mismatch assertion | synthetic known gap recovered | **bootstrap calibration**: 500 synthetic Δ = 0 datasets with a lineage random effect, near-binary score blocks and training-like per-line prevalence; coverage of the conditional CI and of the lineage t-interval reported (stats SF5) |
| Splits | determinism; reproduces the simulated draws in design_numbers; no lineage in two tiers; all lines assigned; refuses to overwrite; **refuses any S whose attempt range overlaps the retired range 20261026…20261125; refuses code, DESIGN.md or FROZEN.tsv differing from the freeze commit, and dirty code** | — | **raw-level leakage canary** |
| Agent guards | scanner flags a planted synthetic label copy (fingerprint and shape) and does not flag expression or CN matrices; hook liveness probe is denied; tripwire raises on `pyarrow.parquet.read_table` of a label file and on a private `labels.py` reader called from a non-entry-point script | — | the same reads from `build_labels.py` succeed |
| §6.8 binary endpoints | toy confusion matrices with known MCC, F1, precision, recall and accuracy; undefined-MCC rule; the P57 call equals `taskEssentialGenes.py` on the generic parents; script diff confirming that `metabolicTasks_Essential.txt` = Robinson data file S3 (57 tasks) | — | constant-call method gives MCC 0 |
| Labels | Entrez join on a toy matrix; NaN counting; fold scoping (a memorising method gets null AP on the held-out fold) | κ(L-prob, L-ge05) on training | tripwire raises on a raw open from a test module |
| KO engine | the 36/37 existing checks; zero-flux shortcut is exact (kept for correctness only: it saves ≈ 5 s per 209 s, compute probe); GLPK = Gurobi per family; nextgem KO = `growth_ratios()` of `gradedEssentiality.py` | CAD/UMPS lethal (H2) | LAP3/TGDS ratio 1; M-open collapses |
| W2 | toy nutrient removal; every formulation string parsed; per-medium WT-growth precheck table; floor non-binding assertion; req_c ≠ 0 for every unflagged serum component; growth calibration | ASNS DMEM vs RPMI; sterol branch serum-free; GOT2 with the gate | shuffled medium; M-open |
| W3 | toy task; `add_moma` vs `moma` on 5 genes; D-mode feasibility log | ETF capability FP; no loss of mevalonate-core essentiality under serum + B-add | pseudo-tasks; matched call rate |
| W4 | toy network per family; repair rule restores growth on the 12 lines the compute probe killed under F0; loss-rule drop procedure gives the same result under shuffled input order (ties broken by Entrez ID); GIMME/FASTCORE/E-Flux/CORDA toy networks with known answers (label-free positive controls replacing the removed ftINIT Hart test) | TYR/DCT in melanoma models; urea cycle absent from most lines | permuted expression; F0 |
| W5 | toy edit flips the KO; evidence validator | v1 → v2 recovery (caveated) | random edits |
| W6 | synthetic information injection | drug positive-control pairs | linear-of-expression null; permuted flux |
| W7 | label permutation → AP ≈ E[AP_random] | B6 ≥ B2 on training | gene-grouped CV |
| W8 | toy SBML round-trip | GR passes | a corrupted SBML is rejected |

---

## 13. Compute plan (v3: LP-only under D2; user target ≈ 1–1.5 days wall-clock)

**Unit costs** (measured or cited; the M3 pilots replace the estimates):
- generic LP sweep: 3.5 CPU-min (baseline README);
- per-line LP build + repair + KO sweep: ≈ 1.5 CPU-min with the vectorised closure map (compute r2 probe);
- 57-task (P57) calls for one line model: ≈ 0.1 CPU-h with a per-task zero-flux shortcut (**estimate; checked by the M3 pilot**);
- QP MOMA: Gurobi.

Every job declares a budget and is killed at 2× it (§11.5).

| Block | CPU-h | of which Gurobi |
|---|---|---|
| W1 + W2 (generic per-medium × variant × gate sweeps, prechecks, floor) | 40 | 0 |
| W3 (capacity matrix, QP MOMA, NGAM, B-loo/B-add, pseudo-tasks) | 25 | 5 |
| W5 (enumeration, local re-evaluation, 2×2 on extracted sets) | 60 | 0 |
| M3 pilots (floor, repair, CORDA, P57 cost, RSS) | 10 | 0 |
| **Stage A** (frozen LP-only table, 30-line panel; W4) | 100 | 0 |
| Stage B (≤ 3 LP configurations × TG-passing lines) | 100 | 0 |
| Stage C (≈ 150 group models × ≤ 4 LP configurations) | 15 | 0 |
| H1 arm: P1(H1) on validation + test lines (LP) | 15 | 0 |
| §6.8 P57 calls for P1(H1) and P1(H2) on validation + test lines (≈ 526 lines × 2 parents) | 110 | 0 |
| Factorial: 24 combinations on the panel → top 3 on training → **2^4 on validation only**; pre-declared test contrasts (i)–(v) | 160 | 0 |
| W6 (including the three assay media) | 40 | 0 |
| **W7 (cap)** | **100** | ≤ 10 |
| W8: H2 and H1 resource builds (2 × 2,030 LP models) and validation | 150 | 0 |
| Subtotal | 925 | ≈ 15 |
| **Re-run buffer (+10%)** | 95 | — |
| **Total** | **≈ 1,020** | **≈ 15** |

**Workers.** 24–28 LP workers plus ≤ 4 Gurobi sessions (§11.4).

**Critical path** (compute-bound; hours from the M1 freeze, H0):

| Step | Work | CPU-h | Wall |
|---|---|---|---|
| 1 | W1 scoring; W2/W3/W5 generic sweeps ‖ Stage A builds | ≈ 235 | ≈ 9 h (H0 → H9) |
| 2 | Shortlists → Stage B + Stage C | ≈ 115 | ≈ 4 h (H9 → H13) |
| 3 | Combination factorial (panel → training) and the validation 2^4 factorial | ≈ 150 | ≈ 6 h (H13 → H19) |
| 4 | P1 frozen → test contrasts, P1(H1) and P57 calls on test; T1 | ≈ 75 | ≈ 3 h (H19 → H22) |
| 5 | W8 ‖ W6 ‖ W7 (W7 time-boxed at **8 h wall**, ≤ 100 CPU-h); T2 at the W7 stop, W8/W6 finish after | ≈ 290 | ≈ 11 h (H22 → H33) |
| 6 | Resource validation complete; REPORT.md | — | ≈ 4 h (H33 → H37) |

- **Expected wall-clock: ≈ 36–42 h (≈ 1.5–1.75 days) at 24–28 LP workers** from the freeze to REPORT.md; cut-list items 2 and 5 recover ≈ 4–6 h. That includes the 10% buffer and assumes the pipeline code (labels, guards, families, metrics) is implemented and verified **before** H0.
- Implementation time is not included in this estimate (open item F-6).
- **Budget.** ≈ 1,020 CPU-h, well under the earlier ≈ 2,650 CPU-h envelope. Disk ≤ 10 GB.
- **Cut list, in order:**
  1. v1.3.0 sensitivity;
  2. M3 (quantitative bounds);
  3. O-NGAM reduced to one level;
  4. CORDA reduced to its 15-line fallback (row A6);
  5. the W7 time box cut to 4 h.

  The lineage × background combinations are never cut (PROMPT).

---

## 14. Milestones and go/no-go (hours from the freeze, H0)

| Milestone | Content | Go/no-go |
|---|---|---|
| **M1** freeze (H0) | `labels.py`, tripwire (`sitecustomize`), scanner, ledgers, FROZEN.tsv, pre-commit hooks; hook liveness probe denied; **freeze commit → user-supplied seed → split drawn once** (§4.1); labels built (including L-ge06); κ rule applied; tAP and §6.8 metric libraries with toy tests; `data/controls/gene_panels.tsv` generated from `notes/apriori_lists_approved.md`; evaluable background groups listed | all §12 metric, split, label and guard tests pass; scanner clean |
| **M2** baseline (H+1) | W1 on training and validation, including the §6.8 binary endpoints for generic H1/H2; G\* fixed; B0–B6, B8±; addressable set; SESOI_C1b frozen from training; label-free Stage-A and W2/W3/W5 sweeps already running | H2 lift on training > 1, otherwise label-free audit (AI r1 SF2) |
| **M3** pilots (H+2) | medium floor and growth prechecks; requirement table; attribution set recomputed from the sheets on disk; calibration; `add_moma`; repair rule; troppo and CORDA decisions (row A6 cost bands); P57 per-line cost; GLPK = Gurobi; RSS → worker counts frozen | a family is dropped if it fails its toy tests; MOMA is dropped above 1 CPU-h per model; if P57 > 0.2 CPU-h per line model, the P57 calls run on the generic parents and on P1(H2) only (stated) |
| **M4** generic levers (H+9) | W2, W3 and W5 shortlists and attribution on training/validation; canary; red-team audit | only TG-passing candidates go forward; otherwise the workstream is marked [-] |
| **M5** context models (H+9 → H+13) | Stage A verdict (H4.1–H4.3); Stage B; Stage C | method-vs-tissue verdict recorded; ≤ 2 configurations shortlisted |
| **M6** integration + **T1** (H+13 → H+22) | combination factorial; validation 2^4 factorial; P1 frozen; test contrasts; rehearsal; canary; **T1** (scalars only) | C1 decides the headline |
| **M7** (H+22 → H+33) | W8 builds; W6; W7 until its stop rule or its 8 h time box; P2 frozen; **T2**; unseal | C4/C5; stop rule documented |
| **M8** resource (H+33) | validation gates pass for 100% of models; manifest | — |
| **M9** REPORT.md (H+37) | claims register; `verify_report`, `verify_citations`; non-author review | all checks pass |

At each milestone: the RESEARCH_TREE status and decision note are updated, LESSONS.md is pruned, and a fresh agent is asked "what should we stop doing / what worked well?".

---

## 15. Limitations stated a priori
- **The test CI is conditional on 7 lineages.** Transfer to new tissues is weakly powered.
- **Residual exposure of the retired seed.** Lineages named by the round-1 review may fall in the new test tier by chance. Only aggregate label statistics on them were ever computed (§11.5).
- **ftINIT is excluded by user decision D2**, because it underperformed the generic model in the user's own analysis.
  - The method-vs-tissue test covers LP families only (GIMME, FASTCORE, E-Flux, CORDA where kept, plus F0).
  - The comparison with the Human1/Human2 papers (§6.8) uses their call rules, labels and metrics on our models, not their per-line extraction.
- **The method-vs-tissue panel has 30 lines (10 lineages × 3),** so lineage R² is estimated with limited power.
- **Attribution is split by tier.** Shapley shares are estimated on validation only; test carries only the pre-declared contrasts (i)–(v) of §10.1.
- **L-ge06 applies the Human1 paper's −0.6 cutoff to Chronos.** That cutoff was used on CERES (19Q3), and the two scales differ; the Human2 paper states no threshold.
- **Sanger-only (KY-only) lines are under-represented** among the lines evaluated for context-specific models, because they most often lack expression. The library stratum cannot be evaluated for context models on them (data r2 N3).
- **LOH is proxied by low relative CN**, because no LOH file is on disk.
- **Some labels and lists are all-line by construction.** L-prob and DepMap's common-essential lists are derived from all lines, and Human2's curation may have been DepMap-informed.
- **Uptake bounds are not rates.** No uptake-rate data exist for DepMap lines, so most bounds stay availability-based; quantitative bounds are confined to the M3 calibration (Huijer 2026 [LIT 3] shows what quantitative exchange data would add).
- **Enzyme constraints are out of scope** (ext-PI O1). We found no evidence that they improve binary essentiality (LITERATURE §3.1). E-Flux with protein input is a crude capacity stand-in.
- **Some essentiality mechanisms are outside FBA** (MTA/PRMT5, ferroptosis, V-ATPase, non-metabolic roles). This bounds what any GEM can achieve.

---

## 16. Open items for the PI freeze check (v3)

Round-2 questions R2-1 to R2-11 were closed in v2. In v3, D2 makes R2-1 (C6) and R2-3 (B7 lines) moot. **F-1 (compute) is resolved:** there is no per-line ftINIT, and the total is ≈ 1,020 CPU-h. **F-4 is resolved by the PI rule in `notes/label_exposure.md`:** aggregate-only exposure does not disqualify an agent from authoring lists. The a priori lists are `notes/apriori_lists_approved.md`, used verbatim.

Remaining items:
- **F-2. PI housekeeping.** Done per `notes/label_exposure.md` and CLAUDE.md:
  - the register is created;
  - the round-1 statistics note is in place;
  - the CLAUDE.md label rules are in place;
  - the label-reading scratch scripts are deleted (V-10);
  - the DATA.md metabolomics count is corrected.

  - V-17 done (register redacted; PI note on the round-2 statistics review);
  - hook hardening done (fail-closed; validation-label and `data/locked/` paths, recursive data searches).
- **F-3. LOH.** LOH is proxied by relCN < 0.5. Should a DepMap LOH file be fetched before M1?
- **F-5. Missing medium sheets.** MEM/EMEM, F12, McCoy's 5A and L-15. If they arrive after M3, their lines stay `fallback_sheet_missing`.
- **F-6. Implementation time.** The ≈ 30–36 h estimate covers compute and agent iteration from the freeze. It assumes the pipeline code and its §12 tests are ready at H0. Implementing them is a separate, unbudgeted phase.
- **F-7. v3 choices by the design author, for PI ratification.** The design author has gene × line exposure (one line), so these protocol and threshold choices made in v3 should be ratified by the PI or by the unexposed verifier:
  - the L-ge06 robustness grid steps (0.2);
  - "undefined MCC = 0";
  - the gene set and summaries in §6.8;
  - the CORDA cost bands (row A6);
  - the P57 fallback at M3;
  - the 30-line panel composition (10 × 3);
  - the W7 8-hour time box;
  - the 28-worker count.

  They are protocol parameters, not a priori lists. **Resolved:** ratified or amended by the unexposed round-3 verifier (`notes/design_review_r3/verification.md`, "v3 final check", §D): grid → 21 steps of 0.05; undefined metrics → NaN, excluded (0 as sensitivity); gene set → parent-model genes; workers → 24 + ≤ 4; CORDA bands, P57 fallback, 30-line panel (with the N-1 fold rule), W7 time box and per-line means ratified.

---

## 17. Amendment log and revision history

### Amendment log (after the M1 freeze)
| Date | Section | Change | Reason | Before/after the affected analysis |
|---|---|---|---|---|
| 2026-10-26 | 4.1 (steps 1-3), 12 | A1: The split code did not exist at the design freeze 61b0beb, so V-1's diff check runs against the reviewed M1 code commit C (`make_splits.py --code-commit C`), with freeze ⊑ C ⊑ HEAD, a clean tree, and C..HEAD touching no src/scripts/tests/DESIGN.md/FROZEN.tsv. In addition: (i) DESIGN.md may change after 61b0beb only by rows added to this table; (ii) before drawing, splits.draw/checks and the accepted attempt must equal the functions extracted verbatim from scripts/design_numbers.py at 61b0beb (seeds 1-1000 and S..S+99), together with the 20 reference draws and the 1-1000 acceptance statistics; (iii) the draw inputs are pinned: input files against data/checksums.tsv at 61b0beb; harmonized_MS_CCLE_Gygi.csv and CTRPAUCMatrix.csv, which that file lacks, through the split.inputs_sha256 row of FROZEN.tsv; per-lineage expression/proteomics counts equal to coverage_by_lineage of the frozen JSON; and FROZEN.tsv rows split.* (inputs, line→lineage map, per-lineage counts, draw parameters) plus the SHA-256 values of make_splits.py, splits.py and cohort.py; (iv) production mode refuses --reference-commit ≠ 61b0beb, --repo and NEXTGEM_* root overrides, and sandbox seeds must be ≥ 500000. Deviation recorded: the user stated S = 38 before the freeze commit (RESEARCH_TREE D4), whereas §4.1 step 2 requires after. The design authors did not know S, and S was kept out of all files until 5a8c218, except one slip: the PI wrote it into an exposure-register entry in an unpushed local commit, amended away before any agent read it (LESSONS.md). Every post-reveal author (M1 infrastructure, metrics, D5 list re-derivation, the M1 reviewers, the PI) attests that they did not evaluate draw() on seeds 38-137 and were not told the tiers. Test IDs retired by the M1 review are listed in results/test_retirements.tsv. | The literal V-1 rule cannot pass, because 61b0beb contains no make_splits.py. The pins make the production tiers a function of pre-reveal artefacts only (review M1-M3, §3 a-e). | Before the split draw; no label or draw has been computed |
| 2026-10-26 | 5.4, W2–W6, App. A (a priori lists) | A2: The a priori lists referenced as `notes/apriori_lists_approved.md` are replaced by `notes/apriori_lists_final.md`. Per user decision D5, a fully unexposed agent derived the lists blind (`notes/apriori_lists_blind.md`, committed 57670ff before it saw the approved lists), then reconciled them with the approved lists item by item on biology, with a difference table and citations re-verified; both attestations are in the file. All DESIGN references to the approved lists now resolve to the final file. | The approved lists were authored by an agent with aggregate-only label exposure after the PI relaxed the authorship rule post hoc (M1 process review); the user chose re-derivation rather than ratification. | Before any analysis that uses the lists; no label or draw has been computed |
| 2026-10-26 | 11.1 (content scanner, pre-commit test-keyed rule) | A3: The canonical split file results/splits/lineage_split.tsv (ModelID, tier, label-free covariates; written once by make_splits.py) is exempt from the scanner, at that path only and only while its SHA-256 equals SPLIT_SHA256 and split_meta.json. Every copy, every modified version and every other table keyed by test ModelIDs (e.g. per-line results) is still flagged. Changed files: src/nextgem/scan.py, scripts/scan_label_copies.py, tests/test_scanner.py; and tests/test_frozen.py, tests/test_splits.py (two assertions that assumed no amendment rows and no drawn split). | After the draw, the pre-commit test-keyed rule blocked committing the split itself: an integer covariate column (n_screens) parses as a gene column, so the file looked like a long ModelID x gene table. | After the split draw (639f3a9 → 62000d3); before any label build or scoring; the split is unchanged |
| 2026-10-27 | 9 (W4), 13 (Stage A rows A5, A6), 11.4 | A4: Stage-A specifications where DESIGN is silent, fixed by the PI after the label-free W4 pilots and before any context-model build is scored. (a) KO gene-state convention for presence families: genes off under expression or the genomic constraint are non-functional; an open reaction whose GPR is still True without the off genes uses that pruned GPR (an off isozyme does not buffer; mirrors ftINIT removeGenes=true); a reaction kept or re-opened despite a False pruned GPR keeps its full GPR; E-Flux keeps full GPRs and a KO zeroes the gene capacity. (b) Items 2-12 of results/w4/README.md §5 (sha256 16f6aac0…a4) are ratified as written (GPR no-evidence rule, FASTCORE ε=1 and core definition, CORDA bins, E-Flux bounds, F4 failed-build rule, thresholds R1-R3, panel draw, permuted control, A7 magnitudes, loss rule, solver tolerances). (c) One declared solver per family, fixed for all reruns: GLPK for F0, F2, F4; Gurobi for F3 (and F3b); Gurobi sessions are taken from a flock slot pool. (d) A5 (F5 LCP) uses a fixed, label-free 24-point grid instead of "coarse 12 + local refinement 12": w_r = 1; w_p ∈ {0, 0.5, 1}; w_c ∈ {0, 1}; τ ∈ {−1.5, −1, −0.5, 0}; z_prot = 0 for lines without proteomics (identical builds deduplicated); loss = −1 for an accepted genomic loss and 0 otherwise, so the term pushes lost genes towards off; F5 is built on GLPK (it reuses the F0 machinery, GLPK-verified). Nested panel LOLO selects among these 24 points only. (e) A6 is dropped under its frozen rule (CORDA pilot 1.5-1.7 CPU-h per model > 1); the drop is reported. (f) Stage A is built on medium M0 (Ham's); per-line media enter at the combination factorial (§10), so Stage A isolates the context lever. Measured unit cost 1.8-3.3 CPU-min per build (above the 1.5 assumed); projected Stage A ≈ 59 CPU-h + A5 ≈ 18 CPU-h, within the ≈100 CPU-h row total. | The frozen A5 row gives no grid values and its local refinement would be label-dependent; the other items are implementation details the design left open, which the pilots had to choose to run. | After label-free pilots; before Stage A and before any context-model scoring |
| 2026-10-27 | 11.3 (test-run check) | A5: Two pre-freeze toy tests of the not-yet-committed scripts/evaluate_val.py were renamed and strengthened during the M2 review fixes (notes/reviews/M2_W1_review.md) and are retired in results/test_retirements.tsv: test_gstar_rule → test_gstar_sign_rule (+ test_gstar_matches_sign_in_look); test_scalar_guard_rejects_unit_keyed_output → test_guard_rejects_unit_keys_in_cells, test_guard_rejects_unit_key_columns, test_guard_rejects_numeric_per_unit_columns. No assertion was weakened. In the same freeze, data/labels/label_derived_lists.parquet (DepMap common-essential/control gene lists; used only for the "minus DepMap common essentials" reporting cells) gets no direct label_file row, because check_frozen.py --add cannot hash protected files; it is pinned transitively by its raw sources in data/checksums.tsv and by the frozen scripts/build_labels.py that wrote it (c54e867). tests/test_tripwire.py::test_entry_points_must_be_frozen asserted the M1 state (no evaluate script exists, so both are denied); it now asserts that a frozen evaluate script is trusted for its own tier only and denied otherwise (same name, no assertion weakened). src/nextgem/frozen.py: in staged (pre-commit) mode, git-ignored frozen inputs (data/, results/**/tmp/), which can never be in the index, are checked in the working tree instead of being reported missing (new test tests/test_frozen.py::test_git_ignored_frozen_file_checked_in_working_tree). | The test-run check counts any test that ever passed in a full run, including runs of uncommitted work in progress. | Before evaluate_val.py was frozen and before the W1 validation look |
| 2026-10-27 | 11.3, 11.5, 17 | A6 (user decision D7): (1) Guard-internal changes to frozen code (test renames, scanner and check-tool fixes, test-harness plumbing) are logged in notes/guard_changelog.md with FROZEN.tsv rows refreshed under this amendment id, instead of a new §17 row each time. This holds only while every frozen label, metric, split and list SHA is unchanged and no reportable number changes; otherwise a §17 row is required. (2) New validation looks reuse the reviewed evaluate_val.py harness. A new look needs its config entry reviewed by a non-author agent, a label-free dry run and the pre-flight check. A full independent code review is required only when evaluate_val.py, metrics.py or another frozen harness file changes. | M2 process review (notes/milestone_reviews/M2_process_review.md): amendments for guard plumbing and per-look code reviews added overhead without leakage benefit. Ratified by the user as D5 requires. | After the W1 look; before any other validation look |
| 2026-10-27 | 13 (Stage A row A6), 9 (W4) | A7 (user decision D8): CORDA (F3b) is not dropped. It runs as a reduced comparison row A6′: CORDA × {R1, R3} × real RNA on 10 representative panel lines, one per panel lineage. In each lineage the line is chosen label-free as the panel line with the smallest Euclidean distance to its lineage's median log2(TPM+1) profile over model genes, with ties broken by ModelID. That is 20 builds, ≈ 33 CPU-h, on Gurobi through the slot pool. Comparison: against F0, F2, F3 and F4 built on the same lines and thresholds in A1, report reaction-set Jaccard, the agreement of KO essential-gene sets and ratio_r correlation, WT growth, repair size, and build cost. A label-based comparison on these 10 training lines is reported as tAP_line through labels.train()/fold() only. CORDA is descriptive and **not eligible for P1 selection**. It enters the method-vs-tissue variance partition only as a sensitivity analysis on the 10 lines. | The user preferred a costed comparison to a silent drop: CORDA is a widely used family, and its 1.5–1.7 CPU-h per-model cost allows 10 lines within the compute budget. | Before Stage A; no context model has been scored |
| 2026-10-27 | 9 (W4/W5), 6 (secondary scores) | A8: The GPR transform (Machado et al. 2016, PLoS Comput Biol, doi:10.1371/journal.pcbi.1005140; user suggestion) was assessed label-free in notes/methods/GPRTransform_assessment.md. With an exact isozyme "arm" encoding, gene-level KO equals the Boolean KO for every gene of both parents (max |Δratio| 8e-11). The paper's own encoding inflates isozyme capacity (Human2 WT +25%). Gene-level E-Flux added no distinct KO values on 3 panel lines, so it is **not adopted** as a family. **Adopted as an exploratory secondary score only**: gene-pFBA KO-cost tie-break, score = round(ratio_r − 10⁻³·min(1, (C_g/μ_g)/(C*/μ*) − 1), 6). Here C* = min Σ w_g·u_g at growth ≥ 0.9μ*, C_g is the same with u_g = 0, and w_g = 1 + max(0, −s_g) from the line's RNA scores (w_g = 1 on generic parents). It runs on the generic parents and the 30 Stage-A panel lines (≈ 3.5 CPU-h, GLPK). It is reported against ratio_r and against B8± on training (LOLO) only. It never replaces ratio_r, enters no gate, selection or confirmatory contrast, and is labelled "exploratory". | User-suggested method; cheap mechanistic tie-break, compared with the expression tie-break baseline B8. | Before any context-model scoring |
| 2026-10-27 | 9 (W2), §1.2 a priori expectations, 13 (M3) | A9: W2 medium decisions after the blind two-agent extraction and review (notes/reviews/M3_W2_review.md; results/w2/review/extraction_diff.tsv is the two-agent record). (a) The media library (compositions, exchange map) lives in results/w2/media_library/, not data/media/, and is frozen by SHA in FROZEN.tsv after the fixes (thiamine HCl → 2 Cl⁻; α-aminobutyrate name; IMDM cystine 292 µM). (b) The cysteine rule is DESIGN-literal: cysteine closed unless 2-ME. Documented free cysteine is a sensitivity analysis. (c) Plain "DMEM"/"DMEM:F12" strings receive standard glutamine; a sensitivity analysis runs without it. (d) The RPMI base is classic RPMI 1640 without pyruvate; the on-disk ATCC sheet (with pyruvate) is a sensitivity analysis (0 genes change). (e) The H2 micronutrient floor at 100× req_c still binds in 20 medium cases (≤ 1.9%; retinoate). This is accepted and reported. (f) Serum % is not used for scaling; M2 bounds are f × req_c. A sensitivity f × pct/10 is computed at the chosen f only. (g) The a priori HPLM-uridine expectation is corrected to Plasmax only (Cantor 2017 Table S1 has no uridine). (h) An explicit ATP supplement maps to adenosine + phosphate (extracellular ATP is hydrolysed); this affects 1 line. (i) M3 is reported as uncalibrated and excluded from selection under its frozen rule. Growth is ≈ 12.3 × the O2 bound because other uptakes stay at 1000 ("the rest as in the best M2"). No calibration.tsv was written. (j) The sterol-branch and CAD/uridine positive controls fail because of model structure: the biomass needs vitamin D3 products reachable only from lathosterol via 7-dehydrocholesterol. They are handed to W3/W5 and are not fixed through the medium. Controls whose expected and observed values are equal are scored "non-informative", not "passed". | Two-agent reconciliation required by DESIGN W2; frozen-text ambiguities resolved on biology (Ishii & Bannai 1985/1987 for cysteine oxidation; Cantor 2017). | Before any medium-variant scoring |

### Revision history → v1 (review round 1)
v0.1 (2026-10-26) → v1 (2026-10-26). Every must-fix is mapped below to the section that addresses it, or to the PI decision (PD) that modifies or rejects it. The disposition of should-fix and optional items is in Appendix D.

| Reviewer | Must-fix | Resolution | Section |
|---|---|---|---|
| External PI | MF1 generic + expression baseline | B8 / B8′; MCC at fixed cutoff (PD10, PD13) | §8, §6.1 |
| External PI | MF2 confirmatory structure, SESOI | fixed sequence C1→C1b→C3→C4→C5; C2 estimate; SESOI 0.01, headroom closed (PD5, PD6) | §6.3, §6.4 |
| External PI | MF3 transfer claim | lineage t-interval + ≥5 of 7; conditional wording. **Modified by PD3**: test raised to 7 lineages, against this reviewer's advice | §4.1, §6.4 |
| External PI | MF4 expressed pairs | required secondary analysis; C3 claim condition; in TG reports (PD8) | §5.4, §6.3, W4 |
| External PI | MF5 attribution on test | full factorial of P1's levers on validation and test, with Shapley shares; medium alone on the generic parent; **sealed until T2** (PD22) | §10.1, §10.2, W2 |
| External PI | MF6 Human1 full arm | P1(H1) at T1; within-parent contrasts (PD12) | §3.1, §6.3, §13 |
| Statistics | MF1 tie-aware AP | tAP; null E[AP_random]; sklearn as sensitivity; toy tests (PD1) | §6.1, §12 |
| Statistics | MF2 C3 passes on noise | ΔAP vs E[AP_random]; paired vs B3; within-lineage permutation null (PD7) | §6.3 |
| Statistics | MF3 common essentials | C1b with CE_fit from training ∪ validation; ≥5 positives; NA cells (PD4) | §5.2, §5.3, §6.3 |
| Statistics | MF4 multiplicity | fixed-sequence gatekeeping, no Holm; sign-flip p (PD5) | §6.3 |
| Statistics | MF5 lineage-level inference | lineage-equal-weighted primary; lineage t-interval; 7 test lineages (PD2, PD3) | §6.2–6.4, §4.1 |
| Statistics | MF6 T1 leaks into P2 | T1 scalars only; detail sealed; validation looks scalar-only (PD21) | §10.2, §4.2 |
| AI failure modes | MF1 raw label files | `labels.py` single door, tripwire, three label files, raw-level canary. **Encryption rejected (PD15)**; PreToolUse hook proposed to the user | §11.1 |
| AI failure modes | MF2 label-derived sources | denylist, `label_derived` column, W7 whitelist, Human2 caveat (PD16) | §11.2, W5, W7 |
| AI failure modes | MF3 fold-scoped access | `labels.fold(k)`, `cv.lolo`, nested choices (PD14) | §11.1, §9 (TG) |
| AI failure modes | MF4 ledgers and caps | hash chain, commit per append, caps in code, idempotent looks, rehearsal, split refuses to overwrite (PD17) | §11.3, §10.2, §4.1 |
| AI failure modes | MF5 frozen manifest | FROZEN.tsv, check_frozen, test-run ledger (PD18) | §11.3 |
| AI failure modes | MF6 provenance | wrapper, clean tree / pinned worktree, atomic writes (PD19) | §11.4 |
| AI failure modes | MF7 solver failures | abort on infrastructure errors; slot pool; solver assertion; MILP gap rule. Mathematical failure gives a no-information line score (open question R2-2) (PD20) | §11.4 |
| Algorithms/compute | MF1 medium infeasible | micronutrient floor, growth precheck, M1 redefined (PD27) | W2 |
| Algorithms/compute | MF2 pruning repair | fixed LP repair rule; per-loss viability (PD28) | W4 |
| Algorithms/compute | MF3 ftINIT spec/cost/acceptance | B7 configuration pinned to the repository Python pipeline; Jaccard ≥ 0.95 one-shot test; measured cost (PD9, PD29) | §8, W4, §12, §13 |
| Algorithms/compute | MF4 run matrix | Stages A/B/C; ≤ 4 configurations per line, ≤ 1 ftINIT (PD29) | W4, §13 |
| Algorithms/compute | MF5 MOMA | `add_moma` once, QP for uniqueness, drop if > 1 CPU-h per model (PD30) | W3 |
| Algorithms/compute | MF6 task modes | product-demand D-modes; NGAM from max flux (PD31) | W3 |
| Data | MF1 medium tiers | screen > onboarded > imputed; unmappable stratum (PD23) | W2 |
| Data | MF2 library covariate | split-file covariates and strata (PD24) | §4.1, §5.4 |
| Data | MF3 evaluation set and fallbacks | expression-line set; CN fallback (PD25) | §4.3 |
| Data | MF4 backgrounds | metabolism-linked groups, minimum n, within-lineage contrast, groups fixed after the split (PD32) | W4 |
| Data | MF5 numbers | all counts from `design_numbers.py` (PD26) | throughout |
| Data | MF6 medium-aware baseline | B2m; medium in B6 (PD11) | §8 |
| Prior work | MF1 published protocol baseline | B7 + C6 (estimate; open question R2-1) (PD9) | §8, §6.3 |
| Prior work | MF2 prior art | new §2 (PD34) | §2 |
| Prior work | MF3 B6 side information | medium features (PD11) | §8 |
| Prior work | MF4 Vieira | LITERATURE corrected (PD33) | LITERATURE §2.1, §3.4 |
| Domain | M1 serum units / nucleosides | log-fraction bounds; nucleosides removed from serum; calibration check (PD27) | W2 |
| Domain | M2 transporter gate, cystine | M1g; cystine (PD27) | W2 |
| Domain | M3 sterol control | sterol branch vs mevalonate core; B-add (PD27) | W2, W3 |
| Domain | M4 gene loss kills a viable model | per-loss drop and log (PD28) | W4 |
| Domain | M5 backgrounds | metabolism-linked groups; KRAS as negative control; IDH dropped (PD32) | W4 |

**Other v1 changes.**
- LITERATURE.md fixes (PD33): Vieira 733 lines and MCF7 parameter selection; the ftINIT DepMap-scale claim reworded; DeepMeta's first author corrected to Wu and its PMID confirmed; the reference gap 3–4 filled with Huijer 2026 and Cheng 2023; Guil 2024, Barrena 2023, Hart 2015 (ref 91) and McSherry & Najork (ref 92) added; a §8.1 "searched, not found" section.
- CRISPRGeneDependency arrived, so L-prob is primary (under the κ rule).
- ProCan arrived and is the primary protein modality.

### Revision history → v2 (review round 2)
v1 (2026-10-26) → v2 (2026-10-26, candidate for freeze). Every round-2 must-fix is mapped below to the section that addresses it; "PD" refers to notes/design_review_r2/PI_decisions.md. Should-fix dispositions are in Appendix D (round 2).

| Reviewer | Must-fix | Resolution | Section |
|---|---|---|---|
| External PI | N1 test tier not unseen; re-seed | seed 20261026 retired; commit-reveal with a user-supplied integer after the freeze commit; refusals in `make_splits.py`; register corrected (PD1) | §4.1, §7, §11.5, §12, §15 |
| External PI | N2 "medium alone" contaminated by the transporter gate | gate made an orthogonal toggle; "medium alone" = M0 → best gate-free variant; gate attributed to context (PD8) | §9 W2, §10.1 |
| External PI | N3 MILP failure rule handicaps B7 | time-limited incumbents accepted for every MILP arm; failure = no incumbent, infeasible or WT ≤ threshold; symmetric no-information score (PD6) | §11.4, §8 |
| Statistics | NM1 split seed compromised | as ext-PI N1 (PD1); `make_splits.py` refuses 20261026 and its attempt range | §4.1, §7, §11.5, §12 |
| AI failure modes | N1 test tier known; un-grindable seed | user-supplied integer after the freeze commit; a commit SHA is not accepted as the seed (PD1) | §4.1 |
| AI failure modes | N2 all-line label caches; content scanner | `ge_U.parquet` and other scratch copies deleted by the PI (register); `scan_label_copies.py` with row fingerprints plus a shape heuristic, run at milestones, pre-commit and pre-push (PD14) | §11.1, §11.5, §12 |
| AI failure modes | N3 hook not live and fails open | hook now user-approved and fail-closed; liveness probe at every session start; tripwire via `sitecustomize` with pyarrow wrappers, Popen audit and an entry-point check; residual risk documented (PD15) | §11.1, §12 |
| AI failure modes | N4 sealed T1 outputs location | `data/locked/sealed/<look>/` (hook-protected, git-ignored); only hashes committed; pre-commit guard until the T2 unseal (PD16) | §10.2, §11.1 |
| Algorithms/compute | NMF1 floor and serum f-scaling cap growth | floor fixed and non-limiting (10 × req, non-binding assertion); f only on synthesisable serum components; req on a stated reference medium (PD7) | §9 W2, §12 |
| Algorithms/compute | NMF2 MILP failure rule vs B7 | as ext-PI N3 (PD6) | §11.4 |
| Algorithms/compute | NMF3 F1 off-M0 growth | the repair rule applies to F1 on any non-M0 medium; B7 unrepaired on M0 (PD10) | §9 W4 |
| Algorithms/compute | NMF4 Stage A not enumerated | frozen Stage-A build table (rows A1–A9) whose total matches §13 (PD12) | §9 W4, §13 |
| Algorithms/compute | NMF5 RAM over-allocation | 12 MILP + 16 LP to start; final counts from ftINIT peak RSS at M3; ≥ 20 GB kept free; MILP lane from M2 (PD13) | §11.4, §13, §14 |
| Data | N1 media sheet coverage overstated | attribution set = screen-tier lines with the sheet on disk at M3 (835 today); others flagged fallback; two-media lineages recomputed (9) (PD9) | §9 W2, App. C |
| Domain | N1 req_c zero for substitutable components | pFBA uptake at ≥ 99% WT with the component open, for synthesisable components; FVA-min for auxotrophs; unit test that req ≠ 0 (PD7) | §9 W2 |
| Domain | N2 heterozygous damaging must not knock out | loss = code 2, or code 1 + LOH proxy / low CN / TPM < 1, or deletion; order-independent drop (PD11) | §9 W4 |
| Prior work | (no must-fix) | B7 renamed (PD3); C6 last (PD2) | §8, §6.3 |

**Other v2 changes.**
- C6 joins the sequence last, with B7g as its only confirmatory comparator (PD2).
- B7 renamed; it runs on the panel, validation and test lines, with training deferred and the F1 training-gate gap stated (PD3).
- Both B8 directions are reported (PD4).
- SESOI basis restated, with a lift-scaled SESOI for C1b (PD5).
- `openpyxl` 3.1.5 installed and pinned in `requirements.txt`, giving GDSC coverage; guide maps approved at low priority, and no claim depends on them (PD17).
- `design_numbers.py` no longer opens any label-bearing file.

### Revision history → v3 (review round 3 and user decision D2; final, pending PI freeze)
v2 (2026-10-26) → v3 (2026-10-26). Sources: `notes/design_review_r3/verification.md` (items V-1 to V-17 and Part 5); RESEARCH_TREE.md decision D2 (user); the PI's v3 instructions (binary endpoints; a priori lists verbatim; compute target). All are pre-freeze changes, made before any analysis.

| Item | Change | Section |
|---|---|---|
| V-1 | The HEAD check is replaced by "the freeze commit is an ancestor of HEAD + no diff in code/DESIGN/FROZEN + clean tree"; only ledger rows and `SEED_SOURCE.md` may be committed in between | §4.1, §12 |
| V-2 | Retired-range check applied to S + attempt (S ∉ 20260927…20261125; re-checked per attempt) | §4.1, §12 |
| V-4 | Claim tiers use each contrast's own SESOI | §6.4 |
| V-5 | SESOI for C1, C4, C5 = 0.01; C5 named in the claim language | §6.4 |
| V-6 | "Medium alone" = M0 → best gate-free documented variant, also in §10.2 | §10.1, §10.2 |
| V-7 | Loss rule: LOH proxied by relCN < 0.5; drop order tie-broken by ascending Entrez ID | §9 W4 |
| V-8 | Small-group flag includes NRF2 | §9 W4 |
| V-9 | CORDA cost bands aligned to row A6 | §9 W4 |
| V-10 | Register updated: remaining scratch scripts deleted by the PI | §11.5 |
| V-11 | Stability study redefined and costed (row A7) | §9 W4 |
| V-12 | Moot (B7 removed) | — |
| V-13 | B-loo index defined from MAR13082's direct reactants plus the named cofactor-pool groups; enumerated by script at M3 | §9 W3 |
| V-14 | B-add skips components already in biomass; list per `notes/apriori_lists_approved.md` §5 | §9 W3 |
| V-15 | Appendix A replaced by a reference to `notes/apriori_lists_approved.md` (verbatim); the v2 addendum marked superseded | App. A, App. D |
| V-16 | Register status: created | §11.5, §16 |
| V-17 | Done by the PI: register redacted; warning note on the round-2 statistics review | §16 F-2 |
| D2 / Part 5 items 1–14 | ftINIT removed everywhere. B7, B7g/B7t, C6, per-line F1, ftINIT Stage-A/C rows, H1 ftINIT pre-builds, the conditional resource builds, the Hart one-shot test and Appendix B are deleted. Chain C1→C1b→C3→C4→C5. Comparators are generic H1/H2 (G\*) plus B0–B3, B6 and B8±/B8′. The method-vs-tissue test uses ≥ 3 non-ftINIT families (GIMME, FASTCORE, E-Flux; CORDA where kept) plus F0. Stage A is LP-only. MILP rules are reduced to QP MOMA. Workers are LP only. | §0–§3, §6.3–6.5, §7, §8, §9 W1/W4/W8, §10, §11.2, §11.4, §12–§16 |
| Part 5 item 15 / PI (c) | New §6.8: MCC, F1, precision, recall, specificity and accuracy under the Human1-paper (Robinson 2020) and Human2-paper (Luo 2026) protocols, extracted from the PDFs with page citations. Undocumented details are marked, with our pre-specified choices. New secondary label L-ge06. | §6.8, §12, §14, §15 |
| PI (d) | All a priori lists are referenced to `notes/apriori_lists_approved.md` verbatim (W2 controls, W3 sets and controls, W4 backgrounds and negative controls, W5 isozyme exceptions, W6 drugs, §5.4 strata, App. A) | as listed |
| PI (e) | Compute: 2^4 factorial on validation only plus pre-declared test contrasts (i)–(v); W7 cap 100 CPU-h with an 8 h time box; re-run buffer 10%; Stage A LP-only on a 30-line panel. ≈ 1,020 CPU-h; expected ≈ 30–36 h wall-clock. | §10.1, §13, §14 |
| PI (f) | F-4 resolved by the PI authorship rule | §11.5, §16 |

---

## Appendix A — Pre-registered lists
The a priori lists are **`notes/apriori_lists_approved.md`, used verbatim**. They were authored by the round-3 verifier, who is not disqualified under the PI rule in `notes/label_exposure.md`. The file is frozen at M1, with its SHA-256 in FROZEN.tsv.
- §1: gene panels (pan-medium, medium-dependent, buffering pairs, false-positive and false-negative classes, CRISPR-artefact strata);
- §2: genomic-background groups and their negative controls;
- §3: isozyme exceptions;
- §4: metabolic-drug rules, classes and positive-control pairs;
- §5: W3 task sets, B-loo/B-add and the W3 controls.

DESIGN.md does not restate, alter or reinterpret these lists. Where earlier versions listed items inline, the text now points to the file. `data/controls/gene_panels.tsv` and the W6 compound table are generated from it by script at M1.

(Appendix B, the ftINIT score, was deleted in v3 under D2.)

## Appendix C — Data facts used (from `results/design/design_numbers.json`; no data fact changed in v3)
- **Cohort:** 1,361 screened lines (from CRISPRScreenMap); 30 lineages; 23 eligible.
- **Coverage:** expression 1,292; CN 1,221; expression and CN 1,214.
- **Proteomics:** Nusinow 309, ProCan 598 (2 ProCan IDs unmapped [dn:coverage_screened.procan_unmapped_ids]), either 659 (screened lines).
- **Other data:** metabolomics 652; PRISM 735; CTRP 612; GDSC 610 (546 with expression).
- **Libraries:** Avana-containing 1,048, Humagne-only 194, KY-only 118, other 1.
- **Screens:** 227 models with ≥ 2 screens [dn:screens.models_with_2plus_screens].
- **Medium tiers:** screen 1,225, mixed 6, onboarded 56, imputed 74.
- **Medium sheets:** on disk: RPMI, DMEM, DMEM/F12, IMDM. Attribution set 835. Sheet-missing fallback 198. Two-media lineages: 9 on the sheets on disk (13 if all eight arrive).
- **Genes:** |U| = 2,866; U∩ = 2,766.
- **Baseline KO:** 149 / 154 lethal; addressable genes in H2 = 1,120.
- **Loss rule:** r2-11 rule vs v1 rule: ME2 23/27, KEAP1 20/27, STK11 52/59, FH 0/1.
- **Resource:** 2,030 models, of which 753 are not benchmarkable.
- **Split simulation:** seeds 1–1000; acceptance 0.854; median test expression lines 345.

## Appendix D — Disposition of should-fix and optional items (review round 1)
Each item is marked A (accepted) or R (rejected), with a one-line reason. Items already covered by a PI decision are marked "PD".

**External PI**
- SF1 A: informative per-gene baseline table (W1).
- SF2 R (in part): library-class strata from the combined matrix answer the batch question (PD24). Relabelling from ScreenGeneEffect would open a second label copy and enlarge the leakage surface.
- SF3 A: known-answer background cases in gene_panels.
- SF4 A: X_rxnRNA comparator.
- SF5 A: training-only medium oracle.
- SF6 A: panel builds start at M2; W8 starts after T1; cut list.
- SF7 A: addressable set and recall ceiling.
- SF8 A: per-task table; H1 resource or a stated reason.
- SF9 A: labels present; training-derived common-essential list (CE_fit).
- O1 A: enzyme constraints out of scope (§15).
- O2 A: 40 validation looks, 20 for W7.
- O3 A: F5 tagged "label-calibrated".
- O4 A (modified): Jaccard against the repository's Python pipeline, since MATLAB is unavailable.
- O5 A: keep the "gate not passed" fallback.
- O6 superseded by PD27.

**Statistics**
- SF1 A (PD6).
- SF2 A: paired-SE one-SE rule.
- SF3 A: harmonised TG with fold consistency.
- SF4 A: line-level contexts.
- SF5 A: calibration test with lineage effects.
- SF6 A: κ outcome pre-specified; B2 from the same label; B2c.
- SF7 A: per-lineage headline is the mean tAP_line.
- SF8 A: proteomic subsets descriptive.
- O1 A (modified): cross-check against the simulated draws, since the real draw is not computed in advance.
- O2 A: U stays primary.
- O3 A: optimistic LOLO transfer estimate.
- O4 A.
- O5 A.
- O6 A.
- O7 A, via the NA rule.

**AI failure modes**
- SF1 A: label-exposure register.
- SF2 A: label-free audit; the first T1 number is the headline.
- SF3 A: double extraction; edits generated by script.
- SF4 A: single writer, job registry, status script.
- SF5 A: budgets, pilots, 1.5× buffer, ≤ 4 agents.
- SF6 A: pre-commit guard.
- SF7 A: one-shot ftINIT test.
- SF8 A: claims register, report and citation verification.
- SF9 A: review front matter and check.
- O1 A.
- O2 A (as C1b).
- O3 A.
- O4 A.
- O5 A.

**Algorithms/compute**
- SF1 A: one parent per worker; deduplicated KOs; reaction sets, not pickles.
- SF2 A: LP sweeps on GLPK under an equality check; MILPs on Gurobi.
- SF3 A: capacity matrix.
- SF4 A: local re-evaluation; curation applied to extracted sets.
- SF5 A: factorial on the panel first.
- SF6 A: X_KOflux.
- SF7 A: H-c capped at 100 CPU-h.
- SF8 A: reuse builds; libSBML relative to the parent.
- SF9 A: family specifications.
- Loopless note A: KO ratios are unaffected by loops, so no loopless step is added.
- ftINIT score formula A (Appendix B).
- H1 ftINIT panel-only: R as a default, because PD12 requires P1(H1) on validation and test.
- Shuffled-medium cost A (0 CPU-h).

**Data**
- SF7 A: parser rules.
- SF8 A: shuffle power stated (13 lineages).
- SF9 A: assay medium in W6.
- SF10 A: protein arm pooled, descriptive on test; ProCan now primary.
- SF11 A: licence NOTICE; re-verify licences.
- SF12 A: NaN counts per gene and per library.
- O13 A: screen type in the split file.
- O14 A: omics-dependent per-lineage tables are dropped for lineages without that omics.
- O15 R: PediatricModelType adds strata unrelated to the metabolic question.
- O16 covered by PD32 (ploidy sensitivity).

**Prior work**
- SF5 A: Huijer and Cheng cited (LITERATURE, §15).
- SF6 A: GPR convention pre-registered.
- SF7 A: reuse `gradedEssentiality`/`taskEssentialGenes` with an equality test.
- SF8 A: renamed "CellFie-inspired".
- SF9 A (modified): troppo is tried at M3 and pinned if it installs and passes the toy tests; otherwise our own implementation.
- SF10 A (done in LITERATURE).
- SF11 A: §8.1 added. The query strings are still to be logged by the literature agent.
- O12 A: check Metabolic Atlas before REPORT.md.
- O13 A: training-only anchors.
- O14 A: gap filled.
- O15 A: full-text reads of Gonçalves 2022 and Gopalakrishnan 2023 before M1.

**Domain metabolism**
- S1 A in part: B-add accepted. Composition scaling R, because coefficients barely change essentiality (LITERATURE §5), and B-loo/B-add cover the presence effects that matter.
- S2 A: supplements, DMEM ambiguity, pyruvate.
- S3 A: HPLM/Plasmax coverage audit.
- S4 A: gene_panels file.
- S5 A in part: CN ≥ 4, multi-copy, housekeeping and moonlighting strata. Guide multiplicity deferred (open question R2-9). Doubling time not available locally.
- S6 A: isozyme exception list; same-compartment buffering; soluble-only protein absence.
- S7 A: curated drug list, prodrug flags, assay medium, controls.
- S8 A as W5 curation: counted in Δ_curation, not as a silent parent patch.
- O1 A: label-free flip analysis.
- O2 A: MTAP as a labelled failure case.
- O3 R: TIFM. No composition file is on disk, and HPLM/Plasmax already cover the physiological counterfactual.
- O4 R: hypoxia variant. No hypoxic labels exist, and it is outside the culture-condition question.

### Round 2 (v2) dispositions
Each round-2 should-fix and optional item is marked A (accepted) or R (rejected), with a one-line reason.

**External PI**
- S1 R: graphical α split. PI decision r2-2 fixes the simple chain with C6 last; the consequence (a C3 failure blocks C4–C6) is stated in §6.3.
- S2 A: SESOI scope; lift-scaled SESOI_C1b, no SESOI for C3.
- S3 A: ≥ 8 of 10 folds on the Stage-A panel.
- R2-2 A: sensitivity that excludes failed lines for all methods.
- R2-8 note: one direction was judged sufficient, but PI decision r2-4 requires both.

**Statistics**
- SF-a A: 97.5th percentile, 2,000 permutations.
- SF-b A: lift-scaled C1b SESOI (PD5).
- SF-c R: PI decision r2-2 puts C6 last and keeps the simple chain; the consequence is stated.
- SF-d A: symmetric failed-line rule, per-gene median value, two sensitivities, counts at T1.
- SF-e A: sign-flip statistic and p-value definition.
- Implementation note A: O(n_g) tAP form, cross-checked.
- R2-4 conditions A (§5.1).
- R2-7 fold-count rule A (§9).

**AI failure modes**
- S1 R: Fernet encryption with an agent-readable key. The PI's guard set (hook, tripwire, scanner) covers the careless paths.
- S2 A: entry-point check, pyarrow wrappers and Popen audit in the tripwire; hook hardening passed to the PI as requested config.
- S3 A: the ledger is the only non-PI commit path; `started` row committed first.
- S4 A: clean-tree check scoped to code; worktree assertion.
- S5 A: exposure register content specified; file creation is a PI action (F-2).
- S6 A: the rehearsal on validation counts as a look.
- CLAUDE.md rule A, as a PI action.

**Algorithms/compute**
- SF-a A: one repair LP per viability task, union re-opened.
- SF-b A: vectorised closure map.
- SF-c A: H1 resource scope fixed now.
- SF-d A: W6 raised to 90 CPU-h.
- SF-e A: H1 pre-builds at M5.
- SF-f A: the gate has 2 states.

**Data**
- N2 A: matrix codes stated; ≥ 2 or het + LOH proxy (PD11).
- N3 A: KY-only under-representation in §15.
- N4 A: small groups are exploratory; all background contrasts are estimates.
- Modal formulation text A.
- Row-level PassesQC note A.
- ProCan unmapped count printed A (2).
- GDSC count cross-check A (610 in our join; the reviewer's COSMIC-based count was 600).

**Domain metabolism**
- S-a A: mevalonate-core diagnostic at every W2 step; W2 choice re-checked in the factorial with B-add.
- S-b A: Appendix A expectations, purines and pyrimidines (below).
- S-c A: note that ASNS tests medium only.

**Prior work**
- SF1 A: B7 renamed (PD3). The optional panel-only 1+0 run is R: budget, and no claim depends on it.
- SF2 A: B7 is H2 only; B7(H1) is an estimate only if the H1 pre-builds exist.
- SF3 A: B6 feature provenance and identical medium information.
- SF4 A: claims about Huijer, Cheng, Guil and Barrena are kept to their abstracts.

**Appendix A addendum (domain r2 S-b). Superseded in v3:** the gene-panel expectations below were replaced by `notes/apriori_lists_approved.md` §1, which corrects them (e.g. ADSL and IMPDH2 are not hypoxanthine-rescued; TYMS and DHFR are medium-dependent). Kept only as a historical record.
- De novo purine genes (PPAT, GART, PFAS, PAICS, ADSL, ATIC, IMPDH2): **expected essential in RPMI/DMEM + serum, relaxed in F12-based media and HPLM** (hypoxanthine).
- CAD/DHODH/UMPS: expected essential even with thymidine (F12). Only uridine or cytidine rescues the pyrimidine ring.
