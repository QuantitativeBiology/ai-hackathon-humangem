# DESIGN.md — nextgem pre-registered design

Version **v0.1 (draft for review round 1)**, 2026-10-26. Author: design agent. Owner: PI.
Status: **not frozen**. The design is frozen at milestone M1, after at most three review rounds. After M1, any change to labels, splits, metrics, controls, candidate sets or decision rules is an *amendment*. Each amendment is dated and logged in §16 together with its reason, and it must be made **before** the affected analysis is run. A change made after the analysis has run is reported as post hoc.

Read with: PROMPT.md (fixed goal), LITERATURE.md (§9 implications), DATA.md, results/baseline/README.md, notes/review_baseline_ko.md.

---

## 0. Summary

- **Question.** Do optimised genome-scale metabolic models (GEMs) of cancer cell lines predict metabolic-gene essentiality better than generic Human1 and Human2? Medium, tasks/objectives, context-specific extraction and GPR curation are the optimisation levers. Each lever is tested under whole-lineage hold-out.
- **Primary endpoint.** On the locked test lineages, the paired difference in **mean per-line AUPRC** (mAP_line) between the selected GEM pipeline P1 and the stronger generic model G\*. The 95% CI comes from a lineage-stratified bootstrap over lines. Success means the CI lower bound is > 0. The justification is in §6.3.
- **Mandatory honesty contrasts.** Every result is shown against prevalence, the gene prior, expression only, generic Human1 and generic Human2. Every result is reported with and without DepMap common essentials.
- **Leakage control.** Three tiers of whole lineages (training / validation / locked test) are produced by `scripts/make_splits.py` with seed 20261026. Test labels sit in a locked file. At most two test looks are planned (T1, T2), and every look is recorded in a ledger.

---

## 1. Mapping of the fixed goal (PROMPT.md) to this design

| PROMPT requirement | Where |
|---|---|
| Baseline FBA on Human1 and Human2, single knockouts of all metabolic genes, AUPRC per gene, per lineage, overall | W1, §6 |
| Context-specific models from genomics, transcriptomics and proteomics; per cancer/tissue type, per genomic background, and combinations | W4 |
| Methods beyond CORDA/ftINIT; explicit test of "cluster by method, not tissue" | W4 (F0–F5; variance partitioning, permuted control) |
| Exchange bounds from the documented medium; physiological variants; error attributable to medium alone | W2 |
| Alternative tasks as objectives or constraints, with the candidate set fixed a priori; which tasks help in which contexts | W3 |
| GPR curation with proteomics and paralogs; curation vs constraints | W5 |
| Flux states vs input omics (subtypes, metabolic drugs, metabolomics) | W6 |
| DepMap data used for fitting held out of evaluation; whole lineages withheld | §4, §9 |
| Auto-research loop for the best predictor; runtime alongside accuracy | W7, §6.6 |
| SBML resource for DepMap 26Q3 lines; benchmark vs generic; per-task and per-lineage breakdown; failure account | W8, §9.4, REPORT.md |
| Reproducibility: fixed seeds, toy tests, controls, non-author review | §5, §10, §11 |

---

## 2. Fixed inputs

### 2.1 Models
| Name | Version | Role |
|---|---|---|
| Human1 | Human-GEM **v1.19.0** (commit 022ed5c2; `data/Human1/`) | Generic baseline; parent for H1 replication |
| Human2 | Human-GEM **v2.1.1** (`data/Human-GEM-main/`) | Generic baseline; **default parent** for all optimisation |
| Human1-published | Human-GEM **v1.3.0** (to fetch by tag with a pinned SHA-256) | Sensitivity check in W1 only (generic KO, same medium) |

All models are used as released. The Human1 TYMS reversibility (MAR04644) is **not** fixed in the baseline. It is reported as a model difference, and as a diagnostic variant in W1. The biomass reaction is MAR13082. The baseline medium ("Ham's", M0) is the 44 inputs of the GR task, with uptake bound 1000 (results/baseline/README.md).

### 2.2 Data (all in DATA.md with SHA-256)
- DepMap 26Q3: CRISPRGeneEffect, (pending) CRISPRGeneDependency, expression (TPM log2(x+1), `IsDefaultEntryForModel == "Yes"`), OmicsCNGeneWGS, damaging/hotspot mutation matrices, Model, ModelCondition, ScreenSequenceMap, CRISPRInferredCommonEssentials.
- Other data: Nusinow 2020 proteomics (375 lines); ProCan-DepMapSanger (pending, user download #5); CCLE metabolomics; PRISM 24Q2 and 20Q2 secondary; GDSC1/2 r8.5; CORUM 5.1 and Complex Portal; Ensembl 116 paralogs; Hart 2015 RNA-seq and labels shipped with Human-GEM (`data/Human-GEM-main/data/datasets/`, used only as a positive control for our ftINIT port).
- Pending media compositions (user downloads #6–8). The fallback if they are missing is in W2.

### 2.3 Gene universe
- **U** = the union of Human1 v1.19.0 and Human2 v2.1.1 genes whose primary Entrez ID is a column of `CRISPRGeneEffect.csv`. Expected size is ≈2,850 (H1 2,825, H2 2,807, overlap ≈2,780); the script reports the exact count.
- Joins use Entrez IDs only, never symbols. The two unmapped genes (GK3P; SLC37A4 alternative-locus ENSG00000281500) and unscreened genes (13 MT- genes, POLR2A, …) are excluded and counted.
- A model that lacks a gene in U, or whose context-specific extraction removed every reaction of that gene, gets **essentiality score 0** for that gene (predicted non-essential). Every method is therefore scored on the identical gene × line set.
- Sensitivity analysis: the same metrics on U∩ = genes present in both models.

---

## 3. Labels (fixed a priori; frozen at M1)

| Label | Definition | Role |
|---|---|---|
| **L-prob** | `CRISPRGeneDependency` (26Q3) probability > 0.5 | **Primary**, if the file is present and verified at M1 |
| **L-ge05** | Chronos `CRISPRGeneEffect` < −0.5 | Pre-specified **fallback primary** if L-prob is absent at M1; otherwise sensitivity |
| **L-ge10** | Chronos `CRISPRGeneEffect` < −1.0 | Sensitivity |

- **Freeze rule.** The primary label is decided once, at M1, by the existence and verification of the file. If `CRISPRGeneDependency.csv` arrives after M1, it becomes a sensitivity label. The primary is never switched after any scoring.
- **Verification of the L-prob file** (label-free, before use):
  - Its ModelID set and gene columns equal those of CRISPRGeneEffect, or the differences are listed.
  - Its SHA-256 is added to DATA.md.
  - On **training lines only**, Cohen's κ between L-prob and L-ge05 is reported. κ < 0.5 triggers a PI review before use.
- **Missing values.** A gene × line pair with NaN in the label source is excluded for all methods, and the count is reported. Genes absent from CRISPRGeneEffect are excluded (§2.3), never counted as non-essential.
- **Common essentials.** The DepMap 26Q3 `CRISPRInferredCommonEssentials` list is used. Every metric is reported (a) on all genes in U and (b) with these genes removed. The list is derived by DepMap from all lines, test lines included, so it is used **only to stratify reporting**. It is **never** used as a feature or for fitting (§11).
- **Additional strata (reporting only):**
  - *selective* genes: training-line prevalence in [0.05, 0.95];
  - *expressed* pairs: TPM ≥ 1 in that line (LITERATURE §2.3 pitfall 7: unexpressed genes are easy negatives).

---

## 4. Splits (leakage control)

### 4.1 Procedure: `scripts/make_splits.py`, run once at M1
Inputs, all label-free:
- Model.csv (`OncotreeLineage`);
- the ModelID index of CRISPRGeneEffect (the values are not read);
- the list of default expression profiles;
- the Nusinow sample map.

Steps:
1. **Universe:** the 1,361 ModelIDs in CRISPRGeneEffect. Lineage = `OncotreeLineage`.
2. **Eligible lineages** have ≥ 15 screened lines. Today that is 23 lineages and 1,331 lines. Lines of smaller lineages always go to **training**, flagged `rare_lineage`, and are never held out or reported per lineage. Today those are Thyroid 11, Ampulla of Vater 9, Testis 5, Vulva/Vagina 2, Fibroblast 1, Adrenal Gland 1 and Other 1, 30 lines in total.
3. **Size strata:** rank the eligible lineages by line count, descending, with ties broken alphabetically. S1 = ranks 1–8, S2 = ranks 9–16, S3 = ranks 17–23.
4. **Draw:** `rng = numpy.random.default_rng(20261026 + attempt)`, with attempt starting at 0. Permute each stratum. The **test** tier takes the first 1 lineage of S1, the first 2 of S2 and the first 2 of S3. The **validation** tier takes the next 1 of S1, the next 1 of S2 and the next 2 of S3. All other lineages go to **training**.
5. **Acceptance constraints** (label-free; if any fails, attempt += 1, up to 100 attempts, and every attempt is logged):
   - (a) every test and validation lineage has ≥ 10 lines with a default expression profile;
   - (b) the test tier has ≥ 150 lines;
   - (c) the test tier has ≥ 25 lines with Nusinow proteomics;
   - (d) Lymphoid and Myeloid are not both in test ∪ validation, so that training keeps at least one haematopoietic lineage.
6. **Outputs:**
   - `results/splits/lineage_split.tsv` (ModelID, lineage, stratum, tier, has_expression, has_proteomics, medium_status);
   - `split_summary.tsv`;
   - `SPLIT_SHA256`;
   - the accepted attempt number.

   The outputs are committed, never re-drawn, and loaded by every analysis through one function (`nextgem.splits.load()`). That function asserts the SHA-256.

Expected size: 14 training lineages plus the rare lines, about 900–1,000 lines; 4 validation lineages, about 150–200 lines; 5 test lineages, about 200–260 lines.

### 4.2 Use of the tiers
| Tier | Labels may be used for | Never used for |
|---|---|---|
| Training (with leave-one-lineage-out CV, 14 folds; rare lines always in the training folds) | Fitting thresholds, calibrations, ML models, the gene prior, training-fold co-essentiality; screening candidates | — |
| Validation | Choosing among ≤ 3 shortlisted candidates per decision (§9.1); at most 60 validation looks in total (ledger) | Fitting parameters |
| Locked test | Pre-declared evaluations T1 and T2 only (§9.3) | Anything else, including error inspection before T1 |

- **Physical separation.** `scripts/build_labels.py` writes `labels_trainval.parquet` and `labels_test.parquet`. The test file goes to `data/locked/`, with mode 0400 and a SHA-256 recorded. Only `scripts/evaluate_test.py` reads it, and that script first appends a row to `results/test_ledger.tsv`.
- **Inputs (X) of test lines** may be used without restriction: their expression, proteomics, CN and media are used to build their models, and that is label-free. Fitted **statistics** are computed from training lines only (validation added at refit). These include gene-specific expression thresholds, normalisation constants and PCA loadings.
- **Refit for test evaluations.** Hyperparameters are frozen at selection. Then every fitted component, of the candidates and of the baselines alike (gene prior, ML models, thresholds), is refit once on training ∪ validation and evaluated on test.
- **Lines without a medium annotation** stay in their tier. In W2 they get the default documented medium (RPMI + 10% FBS, the modal formulation), flagged `medium_source = imputed`. W2's attribution analysis is restricted to lines whose medium could be mapped (§W2).
- **Lines without an expression profile** are scored for the generic models and baselines only. The primary evaluation set is the **test lines that have a default expression profile**, so every method is scored on the same lines. Lines without proteomics are excluded only from the proteomics-specific comparisons, which are reported separately.

---

## 5. Seeds

- **SPLIT_SEED = 20261026** (used only by the split; §4.1).
- **SEEDS = [11, 23, 37, 41, 53, 67, 79, 83, 97, 101]**, stored in `src/nextgem/config.py`.

| Use | Seeds |
|---|---|
| Bootstrap CIs (all) | SEEDS[0] |
| Random-score baseline; random task sets; random GPR edits; permuted-expression controls | SEEDS[0:5]; all 10 at T1/T2 |
| MILP / solver randomness (Gurobi `Seed`, Threads=1): main builds | SEEDS[0] |
| MILP / solver randomness: alternate-optima study on 20 lines | SEEDS[0:3] |
| ML training (W7) and CV shuffles (W6) | SEEDS[0:3]; mean ± SD reported |

- Any further random stream is derived with `numpy.random.SeedSequence(seed).spawn()`. No unseeded randomness is allowed; a test greps for `default_rng()` and `random.` calls without a seed.

---

## 6. Scores, metrics and statistics

### 6.1 Essentiality score
- **Knockout ratio:** ratio_r = clip(round(KO_growth / WT_growth, 6), 0, 1), as in results/baseline/README.md and LESSONS.
- **Essentiality score:** e = 1 − ratio_r, so higher means more essential.
- **Task-augmented variants (W3):** e = 1 if the KO makes any required task infeasible; otherwise e = 1 − ratio_r.
- **MOMA variants:** ratio_r is computed from the biomass flux of the MOMA solution, with the same rounding.
- **Missing genes:** a gene absent from the model gets e = 0 (§2.3).
- **NaN** (solver status not optimal or infeasible): e = 0, counted and reported. A method with more than 0.1% NaN pairs is flagged, and the cause is investigated before its numbers are used.
- **Infeasible KO LP** gives growth 0, so e = 1.
- **Failed context models** (WT growth ≤ 1e-6, or a build error): the line gets the generic parent's scores, flagged. More than 5% failures marks the method as failing. A sensitivity analysis excludes the flagged lines.
- **Ties** are kept as ties. AP is `sklearn.metrics.average_precision_score`: step-wise and uninterpolated, with equal scores forming one threshold and no random tie-breaking. A constant score therefore gives AP = prevalence exactly; this is checked by a toy test.
- **Binary operating point:** lethal if ratio_r < 0.01. Precision, recall, F1 and MCC are reported for FBA-type methods only, as a secondary result.

### 6.2 Metric definitions
For a set of lines L and genes G, with y the label, s the score and prevalence π shown next to every AUPRC:

| Metric | Definition | Measures |
|---|---|---|
| **AP_line(l)** | AP over g ∈ G for line l | Ranking a line's metabolic genes, i.e. what a per-line model is for |
| **mAP_line** | mean of AP_line over l ∈ L | Primary (§6.3) |
| **AP_gene(g)** | AP over l ∈ L for gene g; scored only if g has ≥ 5 positive and ≥ 5 negative lines in L | Ranking lines by their dependency on g (context specificity). A method that is constant across lines (generic FBA, gene prior) has AP_gene = π_g **exactly**, so lift = 1. This metric is therefore blind to generic-model changes (W2 medium applied identically, W3 tasks, W5 curation) and responds only to line-specific information. |
| **mΔAP_gene** | mean over eligible g of (AP_gene − π_g); also median lift and the number of eligible genes | Key secondary: context specificity |
| **AP_lineage(k)** | AP pooled over all gene × line pairs of lineage k | Per-lineage breakdown; reported with π_k and lift |
| **AP_pooled** | AP over all gene × line pairs in L | Overall; inflated by gene identity (LITERATURE §2.3) |
| **Lift** | AP / π (any unit) | Prevalence-adjusted comparison across units |
| AUROC | same units | Secondary |

### 6.3 Primary endpoint and justification
**Primary endpoint:**

> Δ_primary = mAP_line(P1) − mAP_line(G\*)

computed on test lines with an expression profile, genes U, label L-prob (or L-ge05 per §3), and all genes, common essentials included.

- G\* is whichever of generic Human1 or generic Human2 has the higher mAP_line on **validation** (fixed at M2).
- P1 is the GEM pipeline selected on validation (§9.1).
- **Success:** the lower bound of the 95% CI is > 0.

Why this endpoint, and not the alternatives:
1. **It responds to every lever in the goal.** Medium, tasks, GPR curation and extraction all change it. Per-gene AUPRC cannot detect any lever applied identically to all lines, because generic models sit at π_g by construction. If per-gene AUPRC were primary, "beating the generic model" would reduce to "having any across-line signal", and the curation and medium work could not register.
2. **Pairing is clean.** Both methods are scored against identical labels within each line, so differences in prevalence between lines cancel. Pooled AUPRC mixes lines with different prevalence and different screen quality, and it depends on scores being comparable across lines. Per-line AP needs neither.
3. **The statistical unit matches the resampling unit.** Lines are the bootstrap unit (§6.4), and mAP_line is a mean over lines.
4. **It matches the deliverable**, which is one model per line used to rank that line's genes.
5. **Its known weakness is handled.** Gene identity (common essentials) inflates it. The contrast against the gene prior (C2) and the result without common essentials are therefore mandatory. Context-specificity claims rest on mΔAP_gene (C3), not on the primary.

### 6.4 Uncertainty and comparisons
- **Bootstrap.** Lines are resampled with replacement *within lineage* (stratified). The number of resamples is B = 10,000 for the confirmatory contrasts and B = 2,000 otherwise. Intervals are percentile 95% CIs, with seed SEEDS[0].
  - mAP_line, AP_pooled, AP_lineage: resample lines.
  - mΔAP_gene: resample lines and recompute per-gene AP; the set of eligible genes is fixed from the original sample.
  - AP_pooled sensitivity: a two-way bootstrap that resamples lines and genes independently.
  - Two-stage bootstrap: resample test lineages, then lines within them. It is reported for C1, because with five test lineages it measures transfer to new tissues. It is **not** the decision rule (open question Q2).
- **Paired comparisons.** The same resample indices are applied to both methods on the identical gene × line set. An assertion fails if the sets differ.
- **Confirmatory family** (test tier only): Holm correction at α = 0.05 over C1–C5. Everything else is an estimate with a CI. Where many units are tested (genes, lineages, tasks, drugs), Benjamini–Hochberg at q = 0.05 is used and labelled exploratory.

| ID | Contrast (test tier) | When |
|---|---|---|
| C1 | mAP_line(P1) − mAP_line(G\*) > 0 (**primary**) | T1 |
| C2 | mAP_line(P1) − mAP_line(gene prior), two-sided | T1 |
| C3 | mΔAP_gene(P1) > 0 (context specificity) | T1 |
| C4 | mAP_line(P2) − mAP_line(G\*) > 0 (best predictor, may be a hybrid) | T2 |
| C5 | mAP_line(P2) − mAP_line(B6, best non-GEM ML) > 0 (does the GEM add value) | T2 |

### 6.5 Reporting grid (every evaluated method)
Every method is reported in a 3 × 2 × 2 × 4 grid:
- label: L-prob, L-ge05, L-ge10;
- genes: all, minus common essentials;
- universe: U, U∩;
- unit: line, gene, lineage, pooled.

Each cell shows the value, π, the lift and the CI. The primary cell is marked. Per-lineage tables list every lineage of the tier.

### 6.6 Runtime
Every results table has runtime columns:
- wall time and summed CPU time for the model build (per line and total);
- KO sweep time (per line and total);
- ML training and inference time;
- peak RSS;
- hardware: 48-core shared host, Gurobi 13 WLS, number of workers.

Runtime is logged automatically by a `nextgem.timing` context manager into every output parquet.

---

## 7. Baselines (mandatory, on every comparison)

| ID | Baseline | Score | Fitting |
|---|---|---|---|
| B0 | Prevalence | constant (AP = π) | none |
| B1 | Random | uniform random scores, SEEDS[0:5], mean reported | none |
| B2 | **Gene prior** | fraction of fitting lines where the gene is essential (training-minus-fold in LOLO; training for validation; training ∪ validation for test); genes with no data get the overall training prevalence | label-based, within tier |
| B3 | Expression only | log2(TPM + 1) of the gene in the line (higher = more likely essential); missing = minimum | none |
| B4 | Generic Human1 v1.19.0 FBA, M0 Ham's | e from §6.1 | none |
| B5 | Generic Human2 v2.1.1 FBA, M0 Ham's | e | none |
| B6 | Non-GEM ML | gradient-boosted trees (sklearn `HistGradientBoostingClassifier`) on gene prior, own expression, CN, neighbour-enzyme expression (Lagziel-style: mean expression of genes sharing a metabolite in the H2 network, currency metabolites excluded via a fixed list), max paralog expression and isozyme count; no lineage features | training LOLO; hyperparameter grid fixed in W7 |

- **B6** is the reference that every hybrid in W7 must beat (C5).
- **Gate.** Every optimised candidate is reported against B0–B5. It must beat G\* (§9.1) to be eligible for selection, and it is always shown next to B2. "Beats generic but loses to the gene prior" is a reportable result, and it is reported as such.

---

## 8. Workstreams

Each workstream lists its question, inputs, a priori candidates, method, metric, controls, ablations or failure tests, and a decision rule.
- **Training gate (TG):** paired LOLO Δ mAP_line against the current reference on training, with the CI lower bound > 0.
- **Shortlist:** at most 3 candidates go to validation (§9.1).

### W1 — Baseline benchmark (generic Human1 and Human2)
- **Hypotheses:**
  - H1.1: generic H1 and H2 have AP_pooled lift > 1 and mAP_line lift > 1.
  - H1.2: H2 > H1 (Luo 2026 claims gains from GPR fixes).
  - H1.3: B2 gene prior > generic FBA on mAP_line (expected, per LITERATURE §0.3).
- **Inputs:** `results/baseline/{human1,human2}_ko.parquet`; v1.3.0 KO (to be computed); labels; splits.
- **Method:** scoring and bootstraps per §6, on training and validation at M2. The full benchmark on test, and therefore on all 1,361 lines, is produced **only at T1**, together with P1. Generic models have no fitted parameters, but looking at their test-lineage errors would leak into curation targets (W5).
- **Failure-mode tests:**
  1. Human1 TYMS diagnostic variant (MAR04644 lb = 0): report the change in metrics. This is not a candidate.
  2. Mechanism check of β-oxidation lethality in H2 (ACADL, HADHA/B, ECHS1, …) by pFBA flux tracing, before any biological interpretation.
  3. Errors stratified by copy number (amplified vs not), by expression (TPM < 1 vs ≥ 1), by common-essential status and by "KO closes no reaction" (1,750 H2 genes, ratio 1 by construction).
  4. v1.3.0 vs v1.19.0 on U∩.
- **Decision:** none to make; this establishes G\* and the reference rows. If the H2 AP_pooled lift on training is ≤ 1, the M2 go/no-go stops the project for a label/join audit.

### W2 — Medium and exchange bounds
- **Question:** How much baseline error comes from medium misspecification alone? Which physiologically motivated bounds help?
- **Inputs:**
  - Per-line medium: DepMap screens with `PassesQC == True`, `ExcludeFromCRISPRCombined == False` and ScreenType ≠ pDNA, linked through ScreenSequenceMap.ModelConditionID → ModelCondition. Use the **text formulation column** (values like "RPMI + 10% FBS"). In 26Q3 this is the column named `FormulationID`, because the names are swapped; this is checked by a content assertion, not by the column name.
  - A model with more than one distinct base medium across its screens is labelled `mixed` and treated as unmappable.
  - Composition tables: user downloads #6–8, entered with sources into `data/media/compositions.tsv`, plus a reviewed mapping of components to H1/H2 exchange IDs in `data/media/exchange_map.tsv`.
  - Formulation parser: base medium, serum %, serum-free flag and supplements, with unit tests over all distinct strings. Unparsed strings are listed, never dropped silently.
- **Variant ladder (a priori; one factor changes per step; all genes and lines identical):**

| ID | Exchange bounds | Changed factor |
|---|---|---|
| M-open | all 1,654 uptakes open (1000) | negative control (expected collapse of nutrient-synthesis essentiality) |
| M0 | Ham's GR task, 44 inputs, 1000 | baseline |
| M1 | documented base medium + supplements, availability only (present = 1000, absent = 0), plus a fixed inorganic set (O2, H2O, Pi, Fe, ions) | base medium |
| M2a / M2b | M1 + serum model S1 (lipids: palmitate, stearate, oleate, linoleate, linolenate, cholesterol) or S2 (S1 + carnitine, ethanolamine, hypoxanthine, nucleosides, selenium, trace metals). Serum uptake bound = 10 (1% of the default); lines marked serum-free get no serum | serum |
| M3 | M2 (better of a/b on training) with concentration-proportional bounds: bound_i = 1000 · c_i / c_glucose; serum c = %FBS × FBS reference concentration | quantitative bounds |
| M-HPLM, M-Plasmax | the HPLM / Plasmax composition + S2, applied to all lines | physiological counterfactual |

  - The serum component lists and the FBS reference concentrations are fixed in `data/media/serum_model.tsv`, with citations, before M4. A domain-expert agent reviews them.
  - Serum bound sensitivity: ×0.1 and ×10 (reported, not selected).
- **Computation:** media are shared across lines, so KO sweeps run **per distinct medium × variant** (≈ 15 base media × 6 variants × 2 models), not per line.
- **Metrics:**
  1. Paired Δ mAP_line, Δ mΔAP_gene and Δ AP_pooled for each ladder step on identical sets: **mappable-medium lines only**. This is the attribution set; a secondary run uses all lines with imputation.
  2. Attribution:
     - (i) the fraction of generic false positives and false negatives (binary call at ratio < 0.01) that flip to correct when only the medium changes, net of correct calls that flip to wrong;
     - (ii) the ΔAP of M1/M2/M3 vs M0.
- **Controls:**
  - *Positive (label-free):*
    - ASNS predicted essential in DMEM (no Asn) but not in RPMI;
    - de novo purine genes lose essentiality under HPLM (hypoxanthine);
    - cholesterol-synthesis genes gain essentiality under serum-free media.
  - *Positive (training labels):* ASNS label prevalence in DMEM vs RPMI lines (Lagziel replication).
  - *Negative:*
    - **shuffled medium**: line–medium assignments permuted *within lineage* (SEEDS[0:5]). A documented-medium gain must exceed the shuffled gain; otherwise the gain is not line-specific;
    - M-open.
- **Physiological variants:** DepMap labels were measured in culture, so HPLM/Plasmax are **not expected** to improve agreement. They are reported as (a) Δ vs documented medium, which is expected to be ≤ 0, and (b) the set of predictions that change (a "physiology-sensitive gene" list). An external check runs on the K562 HPLM vs RPMI screen (Rossiter/Cantor; user download, conditional).
- **Decision:** the ladder variant with the best training mAP_line that passes TG goes to the shortlist. If none passes, W2 is reported as negative ([-]) and M0 is kept. The attribution numbers are reported whatever the outcome.
- **Fallback if compositions never arrive by M3:** W2 is restricted to M-open / M0 / M0-minus-component ablations (drop each of the 44 Ham's inputs; label-free attribution of which inputs drive calls). The documented-medium claim is reported as not tested.

### W3 — Objectives and metabolic tasks
- **Question:** Which tasks or objectives, added to biomass, improve essentiality prediction, and in which lineages?
- **A priori candidate set** (enumerated now from files on disk; no labels seen):
  - `metabolicTasks_Essential.txt`: 57 tasks (ER 9, IC 12, SU 17, BS 18, GR 1) in 11 functional groups, given by the file's `#` headers (e.g. "De novo synthesis of nucleotides", "Beta oxidation of fatty acids", "Vitamins and co-factors").
  - `metabolicTasks_CellfieConsensus.txt`: 195 tasks.

| ID | Mode | Definition |
|---|---|---|
| T0 | criterion | biomass only (baseline) |
| T1 | criterion | biomass + ER (viability, as in the Human-GEM repo) |
| T2 | criterion | biomass + ER + IC + BS |
| T3 | criterion | biomass + all 56 non-GR Essential tasks (Robinson-style) |
| T4_k (k = 1..11) | criterion | biomass + Essential group k (one at a time) |
| T5 | criterion, line-specific | biomass + CellFie tasks called active in the line (CellFie defaults: local "minmaxmean" thresholds at the 25th/75th percentiles; percentiles computed on training-line expression) |
| C-T1, C-T3, C-T5 | constraint | maximise biomass subject to each required task's output ≥ 1% of its WT max (tasks infeasible on the line's medium are dropped and logged) |
| O-cap | objective | e = 1 − min(ratio_biomass, min_t ratio_t), where ratio_t = task-output capacity after KO / before KO, for the T3 tasks |
| O-NGAM {0, 0.1, 0.3} | objective | biomass with ATP maintenance forced at the given fraction of WT ATP-hydrolysis flux (the reaction ID to be verified in each model; expected MAR03964) |
| O-MOMA | objective | linear MOMA KO solution relative to the WT pFBA flux |
| B-loo_j (j = 1..13) | biomass | leave-one-component-out of the 13 MAR13082 species (attributes calls to biomass components) |

- **Note on medium dependence:** task checks (RAVEN) use their own inputs and are medium-independent. The C-modes add the line's medium (W2 selection).
- **Metrics:** paired Δ vs T0 on mAP_line, AP_lineage and mΔAP_gene (T5, C-T5). Per-task and per-group breakdown: Δ per T4_k per lineage, LOLO on training. "Helps in context k" means BH q < 0.05 across 11 groups × 14 training lineages, and the effect replicates in direction on the validation lineages. This is exploratory and labelled as such.
- **Controls:**
  - *Negative 1 (random pseudo-tasks):* for each group, swap each task's required output for a random metabolite of the same compartment that is producible in WT (SEEDS[0:5]).
  - *Negative 2 (matched call rate):* add the same number of extra "essential" calls as a task set produces, assigned to the genes with the highest B2 gene prior. This tests whether tasks add more than calling frequent essentials.
  - *Positive:* capability-only false positives (ETFA/ETFB/ETFDH under T3) are recovered (label-free check from LITERATURE §1.3).
- **Decision:** at most 2 candidates that pass TG and beat both negative controls go to the shortlist. Otherwise W3 is negative and T0 is kept.

### W4 — Context-specific models
- **Inputs:**
  - expression (TPM);
  - proteomics: Nusinow, and ProCan if it arrives;
  - CN (OmicsCNGeneWGS);
  - damaging and hotspot mutations.
- **Genomic constraint (applies to every family):** a metabolic gene with a homozygous deletion, or with a biallelic damaging mutation and TPM < 1, is knocked out in that line's model. A homozygous deletion is DepMap's deep-deletion call if one is available; otherwise relative CN < 0.3 and TPM < 1. The thresholds are fixed now.
- **Extraction families (≥ 3 standard plus our own; all built on the shortlisted medium):**

| ID | Family | Implementation |
|---|---|---|
| F0 | Direct GPR pruning (naive reference): close reactions whose GPR is False given absent genes; no network repair | own |
| F1 | **ftINIT** (iMAT-like MILP; Human-GEM protocol "1+0"; protects the 57 Essential tasks) | RAVEN-python 3.0.0b1 `raven_toolbox.init.ftinit` |
| F2 | GIMME (LP) | own, toy-tested; or `troppo` if it installs cleanly (pinned) |
| F3 | FASTCORE (MBA-like; core = high-expression reactions) | own or `troppo`; CORDA (`corda` package) as F3b if it installs |
| F4 | **Own: expression-weighted bounds (EWB)**: no pruning; ub/lb of each reaction scaled by GPR-mapped expression (AND = min, OR = sum), normalised to the line's 90th percentile, with floor ε = 0.01 × 1000 | own |
| F5 | **Own: label-calibrated presence rules (LCP)**: gene "off" if w_r·z_RNA + w_p·z_protein + w_c·CNloss < τ; (w, τ) chosen on a fixed 4 × 4 × 3 grid by training LOLO mAP_line; reactions closed by GPR as in F0 | own |

- **Threshold rules (a priori):**
  - R1: global 1 TPM (ftINIT default);
  - R2: global 25th percentile of the line's model-gene expression;
  - R3: local gene-specific threshold = the gene's mean over training lines, bounded to [the 25th, 75th] global percentiles (StanDep/"localT2"-like).
  - Proteomics: protein *detected* means present. Missing protein means **no evidence** (fall back to RNA), never absence.
- **Data modality arms:** RNA, protein (proteomic lines only), RNA + protein (protein overrides RNA where quantified), and each of these with or without the genomic constraint.
- **Levels:**
  1. **per line** (resource unit);
  2. **per lineage**: one model from the median expression of the lineage's training lines; a test lineage's model uses the median expression of its own lines' inputs (label-free);
  3. **per genomic background**: the a priori groups are KRAS hotspot, BRAF V600, PIK3CA hotspot, STK11 loss, KEAP1 loss, MTAP homozygous deletion, ASS1 silenced (TPM < 1) and IDH1/2 hotspot; a group is used if it has ≥ 10 training lines and ≥ 5 test lines;
  4. **lineage × background** combinations with ≥ 8 lines.

  Each line is scored with its model at each level, so the levels are compared on identical lines.
- **Method-vs-tissue test (pre-registered):**
  - *Panel:* 100 training lines, 10 lines from each of 10 training lineages, sampled with SEEDS[0].
  - *Builds:* F0–F4 × {R1, R3} × {real expression, permuted expression}. Permutation shuffles each gene's values across lines, so the per-gene marginals are kept and line coherence is destroyed. This gives 2,000 builds.
  - *Analysis:* PERMANOVA (999 permutations), with R² for method, threshold and lineage, at three levels: Jaccard distance on reaction content, task-capability profiles (57 tasks) and binary KO essentiality profiles. The same partition is run on the input expression as the reference.
  - *Hypotheses:*
    - **H4.1:** R²(method) > R²(lineage) for reaction content (replicating Richelle 2019). Decided by a bootstrap over lines of the difference.
    - **H4.2:** within each method, R²(lineage | real) > the 95th percentile of R²(lineage | permuted). A method that fails H4.2 is reported as "method-dominated".
    - **H4.3:** H4.1 and H4.2 on KO essentiality profiles.
  - *Positive controls* (label-free): TYR/DCT pigment reactions present in Skin/melanoma models; urea cycle absent in most lines; ASS1-silenced lines need arginine import.
  - *Stability:* 20 lines × F1 × SEEDS[0:3] measure alternate-optima Jaccard.
- **Performance evaluation:**
  - Each family × threshold × modality runs at the per-line level on all training lines (LOLO for F5 and R3).
  - TG is computed against the best shortlisted generic pipeline (medium + tasks).
  - Key metric: mΔAP_gene, since context specificity is the point here, together with mAP_line.
  - Ablations: ± genomic constraint; RNA vs protein vs both on the proteomic lines; per line vs lineage vs background vs combination.
- **Decision:** at most 2 configurations go to the shortlist. If no family passes TG, the best one is still carried to M6 and the W4 result is reported as negative.

### W5 — GPR curation
- **Targets (label-free):**
  - A *single GPR edit* is any of: top-level OR → AND; AND → OR; removal of one OR branch.
  - For each reaction, every single edit is enumerated in H2, and in H1 for replication, on the shortlisted medium.
  - Reactions are ranked by how many genes' binary lethal call flips under at least one edit.
  - The top **200** reactions are curated.
- **Evidence rules (pre-specified; every edit is logged in `curation/gpr_edits.tsv` with evidence IDs):**

| Edit | Required evidence |
|---|---|
| OR → AND (complex) | all branch genes are co-members of ≥ 1 CORUM 5.1 or Complex Portal complex **and** have protein co-abundance r ≥ 0.5 (Nusinow; ≥ 30 lines where both are quantified) |
| AND → OR (isozymes) | genes are Ensembl paralogs and are not co-members of any complex |
| Remove OR branch | the branch gene has a UniProt function unrelated to the reaction (non-metabolic); or it is never detected by Nusinow in any of 375 lines **and** has TPM < 1 in ≥ 95% of DepMap lines ("pan-cancer absent", tagged separately) |
| Training-fold co-essentiality support (arm W5b only) | OR-branch genes with Chronos correlation ≥ 0.4 across **training-fold** lines support OR → AND |

  - The Ito 2021 and De Kegel 2021 paralog synthetic-lethality tables (user download #11) are added as AND → OR / OR evidence if they arrive.
- **Arms:**
  - W5a: label-free evidence only.
  - W5b: W5a plus training-fold co-essentiality, recomputed inside each LOLO fold and on training ∪ validation for test.
- **Ablation (curation vs constraints):** 2 × 2 factorial on identical sets: {generic GPR, curated GPR} × {generic model, best shortlisted constraint pipeline (W2 + W3 + W4)}. Reported quantities: Δ_curation, Δ_constraints and the interaction, each with a paired CI. Curation is broken down by edit class.
- **Controls:**
  - *Negative:* random edits of the same number and type, drawn among eligible reactions that lack supporting evidence (SEEDS[0:5]). Curated edits must beat their distribution.
  - *Positive (label-free):* apply the pipeline to v1.19.0 and measure recall of the GPR changes Luo et al. made from v1 to v2 on the targeted reactions (from the diff of v1.19.0 vs v2.1.1).
- **Review:** a non-author agent re-verifies the evidence for a random 20% of edits and for **every** edit that flips ≥ 5 genes.
- **Decision:** curated GPRs enter the shortlist if they pass TG on the generic model. The Δ_curation share is reported regardless.

### W6 — Information in flux states beyond the input omics
- **Features per line**, from the P1 pipeline (frozen at M6):
  - X_KO: the e-score vector across U;
  - X_flux: WT pFBA flux, z-scored per reaction;
  - X_rxn: reaction presence.
- **Comparators:** X_RNA (expression of model genes) and X_prot (Nusinow); also the combinations RNA + flux and prot + flux.
- **Matched capacity:**
  - each block is reduced to k = 20 PCs (k = 50 as sensitivity), with loadings fitted inside the training folds;
  - every block uses the same learner (ridge, α by inner CV) and the same folds;
  - folds are leave-one-lineage-out, or 5-fold within lineage for subtype tasks, with SEEDS[0:3].
- **Analyses:**
  1. **Subtypes:** within-lineage OncotreeSubtype, plus `OmicsInferredMolecularSubtypes`. Subtypes are used for lineages with ≥ 2 subtypes of ≥ 5 lines each; the list is generated by script before any fitting. Metrics are the silhouette in PC space (no fitting) and the kNN (k = 5) macro-F1.
  2. **Metabolic drugs:** PRISM 24Q2 LFC and GDSC2 ln IC50/AUC, for compounds whose annotated target is in U. The compound list is fixed by script before analysis. Metric: Spearman ρ between held-out predictions and observations, with lineage as a covariate. Also a **direct, fit-free test:** the within-lineage Spearman ρ of the drug target's line-specific e-score vs response, compared with the target's expression vs response.
  3. **Metabolomics:** the 225 CCLE metabolites; held-out ρ per metabolite.
- **Controls:**
  - permuted-flux features (rows shuffled within lineage);
  - a synthetic positive control, where nonlinear information injected into a simulated phenotype is recovered;
  - a synthetic negative control, where the flux is a linear function of expression and no gain is expected.
- **Decision:**
  - Claim "flux adds information" if the median across targets of Δρ(omics + flux − omics) > 0, with a Wilcoxon signed-rank test across drugs/metabolites (BH q < 0.05), **and** the gain exceeds the permuted-flux control.
  - Claim "flux better than its input" only if Δρ(flux − omics) satisfies the same criteria.
  - A null result is reported as such.

### W7 — Auto-research loop for the best essentiality predictor
- **Scope:** any predictor, including hybrid ML + FBA, evaluated under lineage hold-out.
- **Seed ideas (a priori):**
  - H-a: B6 features + GEM features (e from P1 and from generic H1/H2, closed-reaction count, isozyme count, flux through the gene's reactions);
  - H-b: H-a without the GEM features (= B6);
  - H-c: flux-cone-learning features (sampling the KO flux cone), only if a 10-line pilot costs ≤ 2 CPU-h per line;
  - H-d: stacking of P1 and B6.
- **Loop protocol:**
  1. **Propose.** Ideas come from agents or the PI. Each is logged in `autoresearch/ideas.tsv` *before it runs*, with: id, hypothesis, proposer, features/model, expected Δ, CPU budget.
  2. **Screen.** 14-fold training LOLO, SEEDS[0:3]. An idea is promoted only if the mean Δ mAP_line vs the current best is ≥ 0.005 and the paired-bootstrap lower bound is > 0.
  3. **Validate.** A promoted idea gets one validation look, which reveals only the scalar mAP_line and its CI; no per-gene output. W7 may use at most 30 of the project's 60 validation looks.
  4. **Secondary check:** gene-grouped CV within training, where the genes in the evaluated fold are unseen, shows whether a model learns from features or from gene identity.
  5. **Stop** when the first of these happens:
     - 10 consecutive ideas fail promotion;
     - a cumulative budget of 300 CPU-h;
     - 24 h before the project deadline;
     - a fresh agent, asked "what promising ideas remain?", returns none that are not already logged.
  6. **P2** = the best validated idea, frozen and evaluated once at T2. The number of ideas tried and the validation shrinkage (training LOLO → validation → test) are reported.

### W8 — Deliverable resource
- **Scope:**
  - The GEM pipeline P1 (not a hybrid ML) is applied to every DepMap 26Q3 model with a default expression profile and a cancer OncotreeLineage. Non-cancer and "Fibroblast" lines are excluded and listed. The count comes from the script.
  - Parameters are those fitted on training ∪ validation, i.e. exactly what was tested at T1, with no all-data refit.
  - Human1-derived models are also released if the H1 replication of P1 is within the CI of H2.
- **Format:**
  - SBML L3V1 FBC v2, gzipped, with medium bounds set and biomass as the objective. The notes record parent, version, pipeline ID, git commit, medium source (documented / imputed) and input checksums.
  - Path: `resource/v{MAJOR.MINOR.PATCH}/{ModelID}_{human2-2.1.1}_{pipelineID}.xml.gz`. The SBML id is `nextgem_{ModelID}_{parent}`.
- **Manifest:** `resource/manifest.tsv`, with sha256, size, build time, reaction/gene counts, medium, flags and validation status.
- **Validation (every model is a gate):**
  - libSBML consistency check with 0 errors;
  - a cobra round-trip that gives identical reaction, metabolite and gene sets and bounds;
  - WT growth equal to the in-memory value within 1e-6 (relative);
  - KO ratios for 20 random genes equal within 1e-6;
  - the GR task passes.
  - MEMOTE snapshot on 10 random models (reported, not a gate).
- **Storage and licensing:** the models are not committed (≈ 5 MB each gz, ≈ 10 GB in total). The manifest and build script are committed. Licence CC BY 4.0, inherited from Human-GEM.

---

## 9. Integration, selection and the test protocol

### 9.1 Selection
1. **Per workstream:** TG on training, then a shortlist of ≤ 3.
2. **Combination at M6:** a full factorial of the shortlisted W2 × W3 × W4 × W5 options (≤ 24 combinations, each including "generic"). It runs first on training LOLO; the top 3 go to validation.
3. **Validation choice:** the highest mAP_line, with a **one-SE parsimony rule**. Among the candidates within one bootstrap SE of the best, pick the simplest, in the pre-declared order generic < medium < tasks < curation < extraction < combinations, and fewer fitted parameters first.
4. **Gate:** P1 must have a validation Δ vs G\* with CI lower bound > 0. If no candidate passes, P1 is the best point estimate and is still tested, labelled "gate not passed".

### 9.2 Validation ledger
`results/validation_ledger.tsv` records: timestamp, git commit, candidate ID, config hash, metric value, and the purpose of the look. There is a cap of 60 looks. Exceeding the cap requires a PI note in RESEARCH_TREE.md.

### 9.3 Test evaluations
| Look | When | What is evaluated (all frozen) |
|---|---|---|
| **T1** | M6 | P1; B0–B5; the ablation ladder (generic → +medium → +tasks → +curation → +extraction); W1 full benchmark; per-lineage and per-task tables; C1–C3 |
| **T2** | end of W7 | P2, B6 (refit), P1 unchanged; C4–C5 |

- `scripts/evaluate_test.py` refuses to run unless the config hash matches a frozen entry in `results/frozen_configs.tsv`.
- Each run appends a ledger row: timestamp, commit, config hash, look ID, invoking agent, PI authorisation line, and the SHA-256 of the output.
- A third look is allowed **only** to re-evaluate after a verified bug fix (§11). Both numbers are reported, and the re-evaluation is labelled post hoc.

### 9.4 Failure account (part of REPORT.md)
For P1 on test, errors are categorised by:
- common essential;
- paralog-buffered (an OR partner is expressed);
- medium-sensitive (prediction changes across W2 variants);
- CN-amplified;
- unexpressed;
- gene-less bypass (the KO is rescued through a reaction without a GPR);
- missing from the model;
- complex subunit.

Lineages and genes where P1 < G\* are listed, with the mechanism traced by pFBA for the top 20.

---

## 10. Verification plan
Every test lives in `tests/` or `scripts/verify_*.py`, and passing is part of each milestone gate. Every result that the PI relies on gets a non-author agent review (`results/<ws>/REVIEW.md`, verdict ACCEPT / ACCEPT WITH FIXES / REJECT), using an independent code path where feasible (LESSONS: the float-tie bug was found this way).

| Component | Toy / known-answer tests | Positive control | Negative control |
|---|---|---|---|
| Metrics | perfect ranking → AP 1; reversed; constant → AP = π exactly; ties as one threshold; per-gene eligibility; set-mismatch assertion | synthetic scores with known AUPRC gap recovered | **bootstrap calibration**: 500 synthetic paired datasets with Δ = 0; CI covers 0 in 93–97% |
| Splits | determinism (same seed, same SHA); no lineage in two tiers; all 1,361 lines assigned; constraints hold | — | **leakage canary**: replace test labels with noise; all training/validation outputs bitwise identical |
| Labels | Entrez join on toy matrix; NaN exclusion counted | κ(L-prob, L-ge05) on training | genes not screened never labelled 0 |
| KO engine | existing 36/37 checks (results/baseline/verification.tsv) + zero-flux shortcut exactness (skip KO if all closed reactions carry 0 WT flux) vs the full sweep on H1/H2 | CAD/UMPS lethal (H2) | LAP3/TGDS ratio 1; M-open collapses |
| W2 medium | toy: nutrient removed → synthesis gene essential; parser on all distinct formulation strings | ASNS in DMEM vs RPMI; purines under HPLM | shuffled medium; M-open |
| W3 tasks | toy task requiring a product; infeasible-task logging | ETF capability-only FP reproduced under T3 | random pseudo-tasks; matched call rate |
| W4 extraction | toy network with known expression → known reaction set per family; GIMME/FASTCORE vs published toy examples | **our ftINIT on the Hart 2015 lines reproduces the repository's AUPRC (0.25–0.33) within ±0.05**; TYR/DCT in melanoma | permuted expression; F0 naive reference |
| W5 GPR | toy edit flips the expected KO; evidence parser | recovery of v1→v2 GPR fixes | random edits |
| W6 | synthetic information injection recovered | — | linear-function-of-expression null; permuted flux |
| W7 | label-permutation run → AP ≈ π | B6 ≥ B2 on training (sanity) | gene-grouped CV |
| W8 | SBML round-trip on 3 toy models | GR task passes | corrupted SBML is rejected |

---

## 11. AI-agent failure modes and guards

| Failure mode | Guard | Detection |
|---|---|---|
| **Leakage** (test labels, test-derived statistics, common-essential list or full-matrix co-essentiality used in fitting) | locked label file; fitted statistics from training only; co-essentiality only inside folds; common-essential list used only for stratification | leakage canary (§10); code review greps for reads of `data/locked/` outside `evaluate_test.py` and for imports of `CRISPRInferredCommonEssentials` outside `reporting/` |
| **Benchmark peeking** (looking at test or full-cohort errors to pick curation targets or tasks) | W5 targets are label-free; W1 full-cohort numbers only at T1 | test ledger audit at every milestone; diff of ledger vs RESEARCH_TREE |
| **Metric drift** (changing the metric, label, threshold, gene set or tie handling after seeing results) | definitions frozen in `src/nextgem/metrics.py` at M1; metric code SHA in every ledger row; amendments only via §16 *before* running | CI test compares the metric module SHA with the frozen SHA |
| **Silent data drops** (joins, NaNs, unmappable media, failed builds) | every join writes `dropcounts.tsv` (in, out, reason); hard assertions on expected counts (1,361 lines; |U|) | reviewer checks drop tables; the report lists all exclusions |
| **Hallucinated citations or data** | every reference must be in LITERATURE.md with a DOI checked against PubMed or the publisher; every number from a paper quoted with its table/figure; every data claim backed by a command output; no hardcoded data (compositions come from files with a source column) | reviewer spot-checks 20% of citations and every composition table |
| **Overfitting through many iterations** | validation cap of 60 looks; W7 promotion threshold; one-SE parsimony rule; at most 2 test looks | report the number of candidates and looks; validation → test shrinkage |
| **Unverified "fixes"** | each bug fix starts with a failing regression test; the full verification suite is re-run; a non-author reviewer confirms; a fix that changes test-tier numbers is reported with before and after values | git log links each fix to its test; verification.tsv diff |
| **Selective reporting** | the full reporting grid (§6.5) for every evaluated method, including negatives ([-] in the tree) | the report generator iterates over `frozen_configs.tsv`, not over a hand-picked list |
| **Solver nondeterminism and float ties** | Threads = 1, Gurobi Seed fixed, ratio_r rounded to 6 dp | serial = parallel test; rerun bitwise check on ratio_r |

---

## 12. Compute plan (16 concurrent Gurobi sessions; other cores for statistics and ML)

The unit costs come from the baseline: a generic KO sweep is ≈ 3.5 CPU-min per model. Other costs are estimates to be replaced by M3 pilot measurements.

| Step | Unit cost (est.) | Count | CPU-h | Wall-clock time (16 workers) |
|---|---|---|---|---|
| W1 v1.3.0 KO + scoring | 5 CPU-min | 1 | < 1 | < 0.5 h |
| W2 media sweeps (per distinct medium) | 4 CPU-min | ≈ 200 (15 media × 6 variants × 2 models + controls) | 15 | 1 h |
| W3 tasks/objectives (task-KO ≈ 1 CPU-h per variant) | 0.1–1 CPU-h | ≈ 45 variants × 2 models | 50 | 3–4 h |
| W4 pilots (10 lines × 6 families) | — | 60 | 5 | 0.5 h |
| W4 variance panel (2,000 builds + KO) | ≈ 3.5 CPU-min (ftINIT ≈ 6) | 2,000 | 130 | 8 h |
| W4 per-line, ≤ 4 configurations × 1,361 lines (build + KO with the zero-flux shortcut) | ≈ 6 CPU-min | 5,400 | 540 | 34 h |
| W4 group models | 6 CPU-min | ≈ 300 | 30 | 2 h |
| W5 edit enumeration + ablation | — | — | 70 | 5 h |
| W6 features + ridge | — | — | 20 | 2 h (non-Gurobi cores) |
| W7 loop (cap) | — | — | 300 | ≤ 20 h |
| W8 resource (≈ 2,100 lines build + validation) | 7 CPU-min | 2,100 | 245 | 15 h |
| **Total** | | | **≈ 1,400** | **≈ 90 h** |

- **Disk.** Context models are stored as reaction-presence bit matrices plus KO parquets until W8 (< 5 GB). The resource is ≈ 10 GB gz. The project ceiling is 30 GB of the ~100 GB free.
- **Overrun rule.** If the ftINIT pilot exceeds 15 CPU-min per model, W4's full per-line build is restricted to F0/F2/F3/F4/F5. F1 then runs only on the variance panel and on the lines with proteomics, and this is flagged.

---

## 13. Milestones and go/no-go

| Milestone | Content | Go/no-go trigger |
|---|---|---|
| **M1** Design freeze (target D+1) | ≤ 3 review rounds done; split script run and committed; labels built and frozen (L-prob or fallback); metric library and its toy tests pass; ledgers created | All §10 metric/split/label tests pass; no open blocking review item. Otherwise fix before any scoring |
| **M2** Baseline (D+1) | W1 on training + validation; G\* fixed; B0–B5 rows | H2 AP_pooled lift on training > 1, or STOP for a label/join audit |
| **M3** Pilots (D+2) | media tables (or fallback); ftINIT Hart reproduction; runtime pilots for F0–F5; zero-flux shortcut verified | ftINIT within ±0.05 of the repo AUPRC (else F1 flagged and demoted); a family is kept only if it passes its toy tests and costs ≤ 15 CPU-min per model |
| **M4** Generic levers (D+2–3) | W2 ladder + attribution; W3 candidates; W5 curation + label-free controls; shortlists | Per workstream: shortlist only TG-passing candidates; otherwise the workstream is reported negative |
| **M5** Context models (D+3–4) | W4 variance partitioning (H4.1–H4.3); per-line/group builds for all lines; TG | Method-vs-tissue verdict recorded; ≤ 2 configurations shortlisted (the best is carried even if it fails TG, flagged) |
| **M6** Integration + **T1** (D+4–5) | factorial; validation selection; P1 frozen; **T1** | C1 decides the headline; REPORT draft of the ablation ladder |
| **M7** W6 + W7 + **T2** (D+5–7) | flux-information analyses; auto-research loop to its stop rule; P2 frozen; **T2** | C4/C5; the stop rule is met and documented |
| **M8** Resource (D+7) | SBML for all eligible lines; manifest; validation | 100% pass the validation gate |
| **M9** REPORT.md (D+8) | approach, best methods, test results, per-task and per-lineage breakdown, failure account, runtime, limitations, remaining ideas | Non-author review of REPORT.md against the ledgers |

At each milestone the RESEARCH_TREE.md status summary and decision note are updated, LESSONS.md is pruned, and a fresh agent is asked "what should we stop doing / what worked well?" (PROMPT.md).

---

## 14. Output layout
Outputs go in `results/{splits,labels,w1_baseline,…,w7_autoresearch}/`, each with a README (command, runtime) and a REVIEW.md. The ledgers are `results/{test_ledger,validation_ledger,frozen_configs}.tsv`. Curation edits go in `curation/gpr_edits.tsv` and media tables in `data/media/` (with sources in DATA.md). The resource has `resource/manifest.tsv`, which is committed; the model files themselves are not committed.

---

## 15. Open questions for reviewers

1. **Primary endpoint (deviation from LITERATURE §9(b)3).** LITERATURE recommends per-gene AUPRC as primary. I chose mAP_line because per-gene AUPRC cannot detect generic-model levers (medium applied to all lines, tasks, GPR curation), and the goal names those levers. Per-gene AUPRC is confirmatory contrast C3 instead. Is this the right trade-off, or should C3 be co-primary for the W4 "context-specific" claim?
2. **Test tier size and the claim of transfer to unseen tissue.** Five test lineages (~17% of lines) give tight line-bootstrap CIs, but the two-stage (lineage) bootstrap will be wide, and "transfer to unseen tissue" rests on n = 5 lineages. Is 1 + 2 + 2 by stratum right? Is constraint (d) (haematopoietic lineages) justified, or is it hand-tuning?
3. **Label freeze at M1, and the common-essential list.** If `CRISPRGeneDependency.csv` does not arrive by M1, L-ge05 becomes primary permanently. Is that acceptable to the PI, or should M1 wait? The DepMap common-essential list is derived from all lines, test included, and is used only to stratify reporting. Should we instead recompute it from training lines with the same rule (sensitivity)?
4. **Missing genes score 0 (union universe U).** This penalises the model that lacks a gene. It also gives extraction methods "easy negatives" on unexpressed genes. The expressed-pairs stratum and U∩ are reported, but which should be the primary universe?
5. **Serum model and bounds.** With uptake bounds of 1000, medium effects are availability effects. The serum bound of 10 and the ×0.1/×10 sensitivity are arbitrary. Is the a priori serum list adequate, given that no uptake-rate data exist for most lines?
6. **ftINIT implementation.** RAVEN-python 3.0.0b1 is a beta port, and MATLAB RAVEN is unavailable. Is the Hart-2015 reproduction within ±0.05 a sufficient acceptance test?
7. **F5 (label-calibrated presence rules)** fits on training labels. Its 48-point grid is small, but it is the one extraction family with label-based fitting. Should it be confined to W7?
8. **Disagreement flag (mild), the gate wording.** The PI requires that every optimised model beat generic FBA. I made this a selection gate on validation (CI lower bound > 0) with a fallback ("gate not passed", still tested). That avoids a project with no test result, but reviewers may prefer a strict stop.

---

## 16. Amendment log
| Date | Section | Change | Reason | Before/after any affected analysis? |
|---|---|---|---|---|
| 2026-10-26 | — | v0.1 draft | — | — |
