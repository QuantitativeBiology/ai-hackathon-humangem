> **PI note (2026-10-26):** this review ran the real split seed 20261026 and reports label summaries on the resulting tiers. That seed is **retired** (DESIGN.md §4.1; PI decision r2-1); the tiers named below are not the project's splits. Logged in `notes/label_exposure.md`. Agents authoring a priori lists must not read this file.

# Design review, round 1: evaluation statistics

Reviewer: biostatistics agent (not the author). Reviewed: DESIGN.md v0.1. Also consulted: PROMPT.md, LITERATURE.md §2.2–2.3 and §8–9, and results/baseline/README.md.

Scope: the primary endpoint, AUPRC mechanics, prevalence, the bootstrap, multiplicity, power, splits, test and validation looks, baselines and labels.

**What I computed (labels only, no model scored).** Label prevalences from `CRISPRGeneEffect.csv` on the U genes (union of the H1 and H2 `depmap_column`), and an independent reimplementation of the §4.1 split procedure (label-free). The AP null and tie behaviour were simulated on **synthetic** labels and scores. No model or baseline score was compared with any DepMap label. Scripts are in my scratchpad (`prev.py`, `stats1.py`, `splits_sim.py`, `testtier.py`, `apnull.py`). They are not project files and are not needed for the design.

---

## A. Facts the design should be built on

| Quantity | Value |
|---|---|
| \|U\| (genes with a CRISPR column) | **2,866** (§2.3 says ≈2,850) |
| Common essentials (DepMap 26Q3 list) inside U | 438 |
| NaN label pairs | 1.4%. 88 genes have some NaN, and 1,355 of the 1,361 lines contain at least one NaN |
| L-ge05, pooled prevalence | 0.116. Per line: median 0.115, 5–95% range 0.099–0.136, range 0.084–0.207 |
| L-ge05, share of positives that are common essentials | **86.5%** |
| L-ge05 without common essentials, per-line prevalence | median 0.018 (5–95% range 0.0075–0.035), median **37–38 positives per line** |
| L-ge10, pooled prevalence | 0.057. Common essentials make up **97.4%** of positives |
| L-ge10 without common essentials | median **3 positives per line**. 11 of 192 predicted test lines have 0 positives, and 136 of 192 have < 5 |
| Lineage ICC of per-line prevalence (L-ge05) | 0.19 with all genes, 0.10 without common essentials. Lineage means range 0.107–0.128 (all genes) and 0.011–0.026 (without common essentials) |
| Genes never positive at L-ge05, all lines | 286. 518 genes are "selective" (prevalence 0.05–0.95) |

**Predicted split.** This is my reimplementation of §4.1 with seed 20261026. The design's script must confirm it.
- Attempt 0 passes.
- **Test** = Pancreas 71, Soft Tissue 50, Peripheral Nervous System 47, Prostate 19, Eye 15. That is 202 lines, 192 of them with expression; Pancreas is 36% of the expression lines. The tier has 28 proteomic lines, 18 of them Pancreas.
- **Validation** = CNS/Brain, Ovary/Fallopian Tube, Bladder/Urinary Tract, Cervix. That is 237 expression lines, and CNS/Brain is 53% of them.

Over 3,000 independent seeds, the accepted test tier had:
- a median of 228 expression lines (5–95% range 194–280);
- a largest-lineage share with median 0.38 (95th percentile 0.49);
- an acceptance rate of 0.89 per draw. Constraint (d) rejects 9.7% of draws and (c) rejects 1.6%; (a) and (b) never bind.

**Per-gene eligibility on the predicted test tier** (≥ 5 positive and ≥ 5 negative lines):

| Label | Genes | Eligible genes | Median positives among eligible |
|---|---|---|---|
| L-ge05 | all | 820 | 24 |
| L-ge05 | without common essentials | 474 | 9 |
| L-ge10 | all | 365 | — |
| L-ge10 | without common essentials | 45 | — |

**AP null and ties** (synthetic):
- **Random ranking.** E[AP] of a random continuous ranking is π + 0.019–0.025 at n = 192 lines with 5–50 positives. In lift terms that is 1.94 at 5 positives, 1.51 at 9 and 1.18 at 24. It is π + 0.003 per line at 37 positives in 2,428 genes.
- **Tie-breaking.** For a near-binary scorer, i.e. a top block of 150–600 genes and everything else tied at 0 (the shape of generic FBA: about 2,650 of 2,870 genes have ratio_r = 1), `sklearn` step-wise AP is **0.025–0.053 lower** than the expected AP under random tie-breaking.
- **No positives.** With no positive in `y_true`, sklearn 1.9.1 returns **0.0** with a warning; it does not return NaN.

---

## B. Must-fix

### MF1. sklearn step-wise AP rewards uninformative jitter. Switch the primary metric to tie-aware AP (§6.1 "Ties", §6.2, §10 Metrics)
- **Problem.** sklearn scores a block of tied scores as a single threshold, with the block's end precision. Breaking ties at random, i.e. adding pure noise, raises the expected AP. On a generic-FBA-shaped line it rises by about 0.03–0.05; on per-gene AP at 5–24 positives it rises by about 0.02.
- That is the same order as any plausible lever effect. It favours, with no information:
  - every continuous scorer over a near-binary one: B2 gene prior, B3, B6 and the W7 hybrids vs FBA, so C2 and C5 are affected;
  - any GEM variant that produces more distinct values: MOMA, O-cap, F4 EWB, partial ratios and MILP alternate-optima noise, so C1 and C3 are affected.
- The auto-research loop (W7) will find this for free. "Ties are kept as ties" is honest bookkeeping, but it makes the metric non-invariant to adding noise.
- **Change:**
  1. Make the primary AP the **expected AP under uniformly random tie-breaking**, computed exactly. This is the closed form of McSherry & Najork 2008 (ECIR) for tied blocks. Implement it in `metrics.py` with an independent Monte-Carlo cross-check.
  2. Keep sklearn step-wise AP as a pre-specified sensitivity row.
  3. Replace the toy test "constant → AP = π exactly" with these:
     - a constant score gives E[AP_random](n, n_pos) exactly;
     - score + 1e-9·noise has a mean over 1,000 seeds equal to tie-aware AP(score), within MC error;
     - perfect and reversed rankings unchanged.
- **Null for lift and ΔAP.** Under tie-aware AP, the expectation over permutations of *any* score vector equals the analytic random-ranking baseline E[AP_random](n, n_pos) (Bestgen 2015 gives the closed form). Use that, not π, as the null in "Lift" and in mΔAP_gene (see MF2). Show π as well, for readability.

### MF2. C3 (mΔAP_gene > 0) passes on noise as written (§6.2, §6.4 C3)
- **Problem.** AP_gene − π_g is about +0.02 for a random continuous score at the test tier's positive counts: median 24 positives (L-ge05, all genes) and 9 positives without common essentials. Generic and line-constant scorers sit at exactly π_g. Any line-varying score therefore yields mΔAP_gene > 0 with a tight CI over 474–820 genes, even when it carries zero information. Examples are F1 alternate optima, F4 scaling, and P1 itself.
- **Change:**
  1. Define ΔAP_gene = tieawareAP_gene − E[AP_random](n_l, n_pos,g).
  2. Add a mandatory **within-gene, within-lineage permutation control**: permute P1's per-gene score vector across lines inside each lineage, 5 seeds. C3 must exceed it.
  3. State C3 as **paired against B3** (expression only), i.e. mΔAP_gene(P1) − mΔAP_gene(B3) > 0, with "> 0 vs random" reported. The scientific claim, "the GEM adds context specificity", needs the comparator that has the same input. Expression alone already carries line-specific signal, and P1 is built from it.
- **Resamples.** In a bootstrap resample, a gene with 0 positives or 0 negatives is dropped from that replicate. Record the mean number dropped. Do not use sklearn's 0.0 for such a gene.

### MF3. Common essentials dominate the primary; add a confirmatory non-common-essential endpoint and a zero-positive rule (§3, §6.3, §6.5)
- **Problem.** 86.5% of L-ge05 positives in U are DepMap common essentials (97% at L-ge10). mAP_line on all genes is therefore mostly a measure of how well each model ranks the same ~438 genes in every line. Reason 5 in §6.3 acknowledges this, but the without-common-essentials result is only "reported".
- **Change:**
  1. Add **C1b** = Δ mAP_line(P1 − G\*) on U minus common essentials. Place it in the confirmatory hierarchy (MF4) directly after C1.
  2. For C1b, define the common-essential set from **training ∪ validation lines only**. Replicate DepMap's rule if feasible; otherwise use "positive in ≥ 90% of fitting lines". Then no test label defines the confirmatory gene set. The DepMap list becomes a sensitivity analysis, which answers open question Q3.
- **Zero-positive rule.** AP_line is computed only for lines with ≥ 5 positives in that cell. Excluded lines are counted, and the same lines are excluded for all methods.
  - This is required because sklearn returns 0.0 for a line with no positives and silently drags means down.
  - With this rule, **the L-ge10 × without-common-essentials per-line cell is not estimable** (median 3 positives per line), and its per-gene cell has only 45 eligible genes. Mark both NA in the §6.5 grid rather than reporting unstable numbers.

### MF4. The multiplicity rule contradicts the success rule; replace Holm with a fixed sequence (§6.3, §6.4)
- **Problem.** §6.3 declares success as "95% CI lower bound > 0", i.e. one-sided α = 0.025, unadjusted. §6.4 applies Holm over C1–C5, which tests the smallest p at 0.01. Two rules can disagree on C1.
- **Further problems:**
  - C2 is two-sided, and H1.3 already expects the gene prior to win on generic FBA, so C2 is not a hypothesis we want to "reject".
  - C4 and C5 are tested at a second look, after C1–C3 are known.
- **Change:**
  - Use a **pre-specified fixed-sequence (gatekeeping)** order: C1 → C1b → C3 → C4 → C5. Each is tested at one-sided α = 0.025, matching the 95% CI lower bound > 0, and testing stops at the first failure. Hypotheses after a failure are reported as estimates only. This controls FWER without splitting α and keeps the stated success rule for C1.
  - Make **C2 estimation-only**: point estimate and 95% CI, no test, outside the family.
  - Define the p-value. Use a paired, lineage-stratified sign-flip permutation test over lines (10,000 flips, SEEDS[0]) alongside the bootstrap CI.

### MF5. The "transfer to unseen tissue" claim needs lineage-level inference that works with 5 lineages (§6.4, Q2)
- **Problem.** The lineage-stratified line bootstrap holds the 5 test lineages fixed. Its CI is **conditional on these 5 tissues**: Pancreas, Soft Tissue, PNS, Prostate and Eye in the predicted draw, with Pancreas at 36% of lines.
- **Two-stage percentile bootstrap at K = 5 is unreliable:**
  - there are only 126 distinct multisets of lineages;
  - the bootstrap variance is low by a factor of (K − 1)/K;
  - the percentile tails are lumpy.
- **Limit of an exact lineage-level test.** An exact sign-flip test cannot reach p < 1/32 = 0.031 one-sided.
- **Change:**
  1. Write the C1 estimand explicitly as "mean over test lines of these five held-out lineages" and phrase the headline accordingly.
  2. Replace the two-stage bootstrap with:
     - per-lineage Δ (mean AP_line within lineage) with its line-bootstrap CI;
     - a **lineage-level t-interval** (df = K − 1) on the K lineage means;
     - the exact sign-flip p (32 configurations).
  3. Pre-register the transfer criterion as **descriptive**: "Δ > 0 in ≥ 4/5 test lineages and lineage-level t-interval reported".
  4. If the PI wants transfer to be confirmatory, raise the test tier to **≥ 7 lineages** (e.g. S1:2, S2:2, S3:3). At K = 7 the minimum sign-flip p is 1/128 and the t half-width is 0.93τ instead of 1.24τ; training drops to 12 lineages. Decide this before M1; it is the only fix here that changes the split.
- **Line-weighting sensitivity.** Add **lineage-equal-weighted mAP_line** as a pre-specified sensitivity row for C1, because one lineage can carry 30–50% of test lines.

### MF6. Information leaks from T1 into P2; make T2 blind to T1 detail (§9.3, §9.4, W7)
- **Problem.** At M6, T1 releases:
  - per-lineage and per-task test tables;
  - the §9.4 failure account, which categorises test errors and pFBA-traces the top-20 lineages and genes where P1 < G\*.
- W7 then runs for days and builds P2 from ideas proposed by agents and the PI who have seen those test errors. C4 and C5 at T2 are then not a clean held-out test.
- **Change (either option is acceptable):**
  - **(a) Seal T1 detail.** At M6, `evaluate_test.py` writes T1 outputs to `data/locked/t1/` and reveals only the C1/C1b/C3 scalars and CIs. The per-gene, per-lineage and failure-account outputs are released only after P2 is frozen. W7 agents get no read access to them. A grep check, analogous to §11, enforces this.
  - **(b) Single unblinding.** Freeze P1 at M6 and evaluate T1 and T2 together after M7.
- **Validation looks.** Specify what a validation look reveals for **all** workstreams, not only W7 (§9.2): the scalar mAP_line and the paired Δ with its CI only, with no per-gene or per-line output. Otherwise W5 and W3 can be tuned on validation errors.

---

## C. Should-fix

### SF1. Pre-register a smallest effect size of interest (SESOI) (§6.3)
- **Line-level CI widths.** These are tight. For the predicted test tier (N = 192 lines) with the lineage-stratified bootstrap, SE ≈ s/√N, where s is the within-lineage SD of per-line ΔAP:

| s (SD of per-line ΔAP) | 95% CI half-width | MDE (80% power, one-sided α = 0.025) |
|---|---|---|
| 0.01 | 0.0014 | 0.002 |
| 0.02 | 0.0028 | 0.004 |
| 0.04 | 0.0057 | 0.008 |
| 0.08 | 0.0113 | 0.016 |

- **Lineage-level CI widths.** With between-lineage SD τ of lineage-mean Δ, the half-width is 1.24τ at K = 5, 0.93τ at K = 7, and 0.43τ for a 23-lineage LOLO.
- **Implication.** The primary will "succeed" for an effect of about 0.003, which is smaller than the tie artefact in MF1.
- **Change.** Declare a SESOI before M1, and require "CI LB > 0 **and** point estimate ≥ SESOI" before using the words "improves" or "better". For example:
  - Δ mAP_line ≥ 0.01 absolute; or
  - ≥ 10% of the B2 − G\* gap measured on training at M2.
- Estimate s from training LOLO at M2, using *generic-model variants only*, to replace the scenario table above.

### SF2. The one-SE parsimony rule must use the paired SE (§9.1 step 3)
- The "bootstrap SE" of a candidate's mAP_line is dominated by between-line AP variance, which all candidates share. It is about 5–20× larger than the SE of the paired difference.
- Used as written, the rule will almost always pick the simplest candidate.
- **Change:** candidates within one **paired-difference** SE of the best.

### SF3. The training gate (TG) is near-automatic and multiplicity-blind (§8 header, W2–W5)
- **Problem.** On about 950 training lines the line-level LB > 0 passes trivially small effects. Dozens of candidates are screened: W3 has about 45, W4 has families × thresholds × modality arms × levels, and the factorial has 24.
- **Change.** Harmonise TG with W7's promotion rule:
  - mean Δ ≥ 0.005, **and** LB > 0, **and** Δ > 0 in ≥ 9 of 14 LOLO folds (lineage sign consistency; binomial p ≈ 0.21 under the null, a weak but meaningful filter).
- **Bootstrap scope.** State that the CI is conditional on the fitted fold models, i.e. resampling lines with predictions fixed.
- **Reporting.** Log the number of candidates screened per workstream. Report the expected winner's-curse optimism at validation: E[max of m noise Δ] ≈ σ·√(2 ln m), where σ is the paired validation SE.

### SF4. W3's "helps in context k" cannot replicate on validation (§W3 Metrics)
- **Problem.** The contexts are training lineages. The validation lineages are *different* lineages, so "replicates in direction on the validation lineages" is not a replication of "helps in lineage k".
- **Change (pick one):**
  - **(a)** Define the context-specific claim on line-level contexts that exist across tiers: genomic background (W4 groups), base medium and serum-free status. Replicate it on validation lines of the same context.
  - **(b)** Keep lineage contexts, but replicate by a leave-one-line-out-within-lineage check, and call it within-lineage exploratory only.

### SF5. The bootstrap calibration test must include lineage heterogeneity and realistic ties (§10 Metrics, negative control)
- **Problem.** The 500 synthetic Δ = 0 datasets should be generated:
  - with a lineage random effect τ > 0 on Δ;
  - with near-binary score blocks and per-line prevalences drawn from **training-tier** labels (median 320 positives per line at L-ge05; 37 without common essentials).
- **Report.** Coverage of:
  - the conditional (stratified) CI for the conditional estimand, where 0.93–0.97 is expected;
  - the lineage-level t-interval for the population estimand.
- This documents what the stratified CI does and does not cover. Monte-Carlo SE of the coverage at 500 sims is 0.0097, so 0.93–0.97 is ±2 SE, which is fine.

### SF6. Specify the outcome of the κ check and the label lineage of the gene prior (§3, §7 B2)
- **κ outcome.** "κ < 0.5 triggers a PI review" leaves a label decision open after training labels are seen. Pre-specify the outcome now: if κ < 0.5, L-ge05 is primary and L-prob is sensitivity. Otherwise L-prob is primary.
- **File status.** `CRISPRGeneDependency.csv` is still absent locally, so the fallback is likely to bind. That is acceptable: the decision depends on file existence only and is label-free.
- **Both lists use test lines.** L-prob is produced by DepMap from a fit over all lines, test included, just like the common-essential list. Note this in §11 as an accepted, label-generating (not fitted) dependency.
- **B2 label.** B2 must be computed **from the same label as the one evaluated** (an L-ge05 prior for L-ge05 cells, and so on). Its ties (286 never-positive genes at fraction 0) are handled by the tie-aware AP of MF1.
- Add **B2c** = the mean Chronos score over fitting lines (continuous, fewer ties) as a sensitivity baseline. It is the stronger "gene identity" reference for C2.

### SF7. Per-lineage breakdown metric (§6.2 AP_lineage, §6.5)
- **Problem.** Pooled AP within a lineage depends on whether scores are comparable across lines. Lineage prevalences differ (0.107–0.128 overall; 0.011–0.026 without common essentials, a 2.3-fold range).
- **Change:**
  - Make the per-lineage headline the **mean AP_line within lineage**, i.e. the same unit as the primary, with Δ vs G\* and lift over E[AP_random].
  - Keep pooled AP_lineage as secondary.
  - Never rank lineages by raw AP. Lineage prevalence tracks screen quality as well as biology: Pancreas and Bowel have higher prevalence and SD.

### SF8. Make proteomics and other small-subset comparisons on test descriptive (§4.2, W4 ablations)
- In the predicted draw, the test tier has 28 proteomic lines, and 18 of them are Pancreas. RNA vs protein vs both on test is effectively a single-lineage comparison.
- **Change:** label these comparisons descriptive, and report them with lineage composition. Constraint (c) (≥ 25 proteomic lines) is not enough to support inference.

---

## D. Optional

- **O1. Split procedure.** I see no label-selection bias.
  - The acceptance constraints are label-free and the seed sequence is deterministic. Attempt 0 passes under my reimplementation, so there is no garden of forking paths if the script runs once and is logged.
  - Constraint (d) shifts Lymphoid's test inclusion probability from 0.125 to 0.085. This is a label-free design choice; keep it, but justify it in one line (training needs a haematopoietic suspension-culture lineage for medium and serum effects).
  - Constraints (a) and (b) never bind in 3,000 draws and could be dropped. Keeping them is harmless.
  - Have the split script assert my predicted tiers as an independent cross-check, or explain any difference (e.g. use of `shuffle` vs `permutation`, or the order of the RNG calls).
- **O2. Universe.** Update "≈2,850" to the computed 2,866 once the script confirms it.
  - Under tie-aware AP, missing genes (e = 0) join the bottom tied block. This is fine, but U∩ as sensitivity remains necessary.
  - I would keep U as the primary universe (Q4). The "easy negatives" concern applies equally to every method on identical pairs, so it is paired out of Δ. The expressed-pairs stratum covers the rest.
- **O3. Supplementary transfer estimate.** At T1, also report the frozen P1 protocol's LOLO Δ over all 18 non-test eligible lineages, as a lineage-level distribution (K = 18 → t half-width ≈ 0.5τ). Label it optimistic, because those lineages informed selection. The gap between this and the test is itself informative about selection bias.
- **O4. Auto-research (W7).** Report the expected validation optimism for the number of promoted ideas (SF3 formula). Consider requiring P2's validation Δ over the runner-up to exceed one paired SE; otherwise take the simpler of the two.
- **O5. W6 decision rule.** "Wilcoxon signed-rank across drugs (BH q < 0.05)" mixes one global test with a per-unit correction. Specify:
  - one Wilcoxon across targets per analysis (subtypes, drugs, metabolomics: 3 tests, Holm);
  - BH applied only to per-target Δρ, reported as exploratory.
- **O6. Bootstrap resamples.** B = 10,000 with percentile intervals is fine for a mean over about 200 lines. Use BCa only if the line-level ΔAP distribution is strongly skewed (check at M2 on training).
- **O7. Label threshold.** L-ge05 (prevalence 11.6%) is the sensible fallback primary. L-ge10 (5.7%, almost all common essentials) is a reasonable sensitivity analysis for all genes, but it is uninformative without common essentials (MF3). No change is needed beyond marking those cells NA.

---

## E. Items needing no change
- Per-line rather than per-gene as the primary unit: §6.3 reasoning 1–4 is sound. Per-gene AP is structurally blind to line-constant levers, which the goal names explicitly. Keep it, with C3 fixed per MF2 and C1b added per MF3.
- Pairing on identical gene × line sets, with an assertion: correct.
- Using the common-essential list only for reporting: correct, but see MF3 for the confirmatory set.
- The label freeze at M1 decided by file existence: correct.
- The test ledger, frozen-config hash and post-hoc labelling of a bug-fix re-look: correct. MF6 concerns *what* is revealed, not how it is logged.
- The two-way (line × gene) bootstrap for AP_pooled sensitivity: correct.
