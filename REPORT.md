# Context-specific genome-scale metabolic models and cancer gene essentiality: a pre-registered benchmark on DepMap 26Q3

> **Read first.** This article was written by the AI pipeline. The human-written disclaimer at the top of [README.md](README.md) applies to it: the project is an exercise in agentic AI for bioinformatics (Beltrao-lab hackathon; INESC-ID & IST), no result has yet been verified by a human, and the text is not a traditional publication.

> **Status, 2026-10-28: draft before the test looks.** Training and validation results are final. The test looks T1 (P1) and T2 (P2) have not run, so every test-tier number is a marked placeholder, `[T1: …]` or `[T2: …]`. Every other number is regenerated from the committed file cited next to it by `scripts/report_build.py --check`, which fails if the text and the file disagree. A few counts quoted from README prose are cited but not machine-checked.

## Abstract

**Background.** Published tests of genome-scale metabolic models (GEMs) against CRISPR screens report binary metrics after choosing among many settings, and none holds out whole tissues (LITERATURE.md §2.1). We asked whether medium-aware, curated, context-specific versions of Human1 (Human-GEM v1.19.0) and Human2 (v2.1.1) predict metabolic gene essentiality in DepMap 26Q3 better than the generic models on lineages they have never seen.

**Design.** The benchmark was pre-registered (DESIGN.md v3, frozen at 61b0beb; amendments A1–A13). <!--n:split.n.total-->1,361 screened lines were split by lineage into 7 test, 4 validation and 12 training lineages (results/splits/split_meta.json). The primary metric is tie-aware average precision (tAP) per line, averaged with equal weight per lineage. Four levers were screened on training: medium, biomass and tasks, context-specific extraction, and GPR curation. At most 40 scalar-only validation looks were allowed; <!--n:looks.used-->5 were used (results/validation_ledger.tsv).

**Results on training and validation.** Generic Human1 scores validation tAP <!--n:w1.val.h1-->0.302 against <!--n:w1.val.h2-->0.295 for Human2, Δ <!--n:w1.val.h2h1-->−0.0071 [−0.0079, −0.0063], so G\* = Human1 (results/w1/val_look_W1_M2_baseline_gstar.json). A label-only gene prior (B2) scores <!--n:w1.val.b2-->0.892, and the generic model with an expression tie-break (B8+) <!--n:w1.val.b8p-->0.430 (results/w1/val_look_W1_M2_baseline.tsv). No context-specific family beat G\* on training, evidence-based GPR curation lowered tAP, and a documented culture medium helped only as much as a randomly shuffled one. Removing the cofactor/vitamin pool from the biomass reaction gained <!--n:w3.val-->+0.0322 [0.0308, 0.0334] on validation (results/w3/val/val_look_W3_M4_shortlist_decision.json). The selected pipeline P1 (Human1, biomass pool removed, documented medium) beats G\* on validation by <!--n:m6.val.p1.short-->+0.052 [0.050, 0.054] (results/m6/val/val_look_M6_factorial_decision.json). This estimate is selection-biased and carries four flags (W3-D1, W3-MEV, W2-SHUF, W2-MOPEN). P1 stays below B8+ and far below B2.

**Test.** [T1: C1, Δ̄_lin(P1 − G\*) on all genes, with 95% CI, sign-flip p and §6.4 claim tier; reported together with C1b, the same contrast without CE_fit, with CI, p and claim tier against SESOI_C1b; transfer wording from the lineage t-interval and the k-of-7 count — filled by scripts/report_build.py after T1]. [T2: C4, Δ̄_lin(P2 − G\*), and C5, Δ̄_lin(P2 − B6), with CIs, p and claim tiers, confirmatory only if C3 passed — filled after T2].

**Other findings.** Simulated flux states carried no information beyond the omics and the medium. An auto-research loop found that GEM features add nothing to gene-prior machine learning, so the best predictor, P2, is a non-GEM model. We release <!--n:w8.models-->40 SBML models that cover <!--n:w8.lines-->1,361 DepMap lines (results/w8/p1/summary.json). The gains a GEM can make on unseen lineages come from removing structural false positives of the generic model and from a realistic medium, not from line-specific modelling. The strongest predictor of CRISPR essentiality here is a gene prior learned from the labels, not a GEM.

## 1. Question and pre-registration

**Question.** PROMPT.md asks for optimised Human1/Human2 models that predict DepMap 26Q3 metabolic gene essentiality better than the generic model, judged by AUPRC per gene, lineage and overall. It also asks for an SBML resource, a flux-information test and an account of failure. DESIGN §1 maps these to eight workstreams: W1 baseline, W2 medium, W3 tasks and biomass, W4 context-specific models, W5 GPR curation, W6 flux information, W7 auto-research, W8 resource.

**Pre-registered hypotheses** (DESIGN §9).
- **W1.** H1.1: the generic models beat a random ranking. H1.2: Human2 beats Human1, as its paper reports [2]. H1.3: the gene prior B2 beats generic FBA.
- **W2–W5.** A documented medium beats Ham's medium and a shuffled medium; a priori tasks or biomass edits beat their negative controls; curated GPRs beat the generic ones.
- **W4 (H4.1–H4.3).** Extraction method explains more variance in model content and KO profiles than lineage. Within each method, lineage explains more than under permuted expression; methods that fail are "method-dominated".
- **W6.** Flux states predict subtypes, drug response and metabolomics beyond the omics they were built from.

**Design.** DESIGN.md v3 was frozen at 61b0beb after three review rounds by seven reviewer agents. Every later change is a dated row in §17 (A1–A13), and every user or PI decision is logged in RESEARCH_TREE.md (D0–D14). The user excluded ftINIT, because it underperformed the generic model in the user's own analysis, and set a compute target of about 1–1.5 days (D2).

**Endpoints and claim language** (DESIGN §6.3–6.4).
- The primary endpoint, C1, is Δ̄_lin(P1 − G\*): the lineage-equal mean of per-line ΔtAP on the 7 locked test lineages.
- Gatekeeping runs C1 → C1b → C3 → C4 → C5, each one-sided at α = 0.025. C1b removes CE_fit (genes positive in ≥ 90% of fitting lines). C3 tests P1's per-gene ranking against expression. C4 and C5 compare P2 with G\* and with B6.
- SESOI is 0.01 (C1, C4, C5). SESOI_C1b is <!--n:sesoi.c1b-->0.0057, fixed on training before any model was scored (results/w1/sesoi_c1b.json).
- Claim tiers: "detectable improvement" if the test passes; "improves" if the point estimate also reaches the SESOI; "meaningfully improves" if the CI lower bound does. "Generalises to new tissues" needs a positive lineage t-interval and Δ > 0 in at least 5 of 7 test lineages.
- Flags that travel with a result: gate not passed, method-dominated, label-calibrated, exploratory, post hoc.

**Selection rules** (DESIGN §8–9).
- **Training gate (TG).** Δ̄_lin ≥ 0.005 vs the reference, a paired-bootstrap CI lower bound > 0, and Δ > 0 in ≥ 9 of 12 leave-one-lineage-out (LOLO) folds (≥ 8 of 10 on the 30-line panel).
- **Validation.** At most 3 candidates go to a look. The validation gate (CI lower bound > 0 vs G\*) comes first; then the paired one-SE rule picks the simplest candidate in a pre-declared order. Each look reports the winner's-curse term σ·√(2 ln m).

## 2. Data and splits

**Inputs** (DATA.md; SHA-256 per file in data/checksums.tsv): Human-GEM v1.19.0 and v2.1.1; DepMap Public 26Q3 (CRISPR, expression, copy number, mutations, models and screen media); ProCan and Nusinow proteomics; CCLE metabolomics; PRISM, GDSC and CTRP drug screens; CORUM, Complex Portal and Ensembl paralogs; and medium formulation sheets, including HPLM [55] and Plasmax [56].

**Split.** `make_splits.py` drew the split once, with seed <!--n:split.seed-->38 supplied by the user after the freeze (commit-reveal). The earlier seed was retired because a round-1 review had published its draw. Lineages with ≥ 15 screened lines were stratified by size, and the <!--n:split.rare.n-->30 lines of rarer lineages always train. Counts are from results/splits/split_summary.tsv and results/splits/split_meta.json.

| Tier | Lineages (screened lines) | Lines | With expression | With proteomics |
|---|---|---|---|---|
| Test (7) | <!--n:split.list.test-->CNS/Brain (133), Pancreas (71), Breast (67), Ovary/Fallopian Tube (62), Liver (25), Pleura (21), Eye (15) | <!--n:split.n.test-->394 | <!--n:split.expr.test-->380 | <!--n:split.prot.test-->170 |
| Validation (4) | <!--n:split.list.val-->Lymphoid (97), Soft Tissue (50), Bladder/Urinary Tract (34), Prostate (19) | <!--n:split.n.val-->200 | <!--n:split.expr.val-->198 | <!--n:split.prot.val-->85 |
| Training (12 + rare) | <!--n:split.list.train-->Lung (126), Esophagus/Stomach (96), Bowel (85), Head and Neck (77), Skin (75), Bone (59), Peripheral Nervous System (47), Myeloid (45), Uterus (38), Biliary Tract (36), Kidney (35), Cervix (18), plus rare lines | <!--n:split.n.train-->767 | <!--n:split.expr.train-->714 | <!--n:split.prot.train-->404 |

**Labels** (DESIGN §5; results/labels/build_report.json). The primary label, L_prob, is DepMap gene dependency > 0.5. Chronos [22] gene effect < −0.5 and < −1 are sensitivities, and < −0.6 (L-ge06) serves the Human1-paper protocol. L_prob stayed primary because its agreement with L-ge05 on training was κ = <!--n:labels.kappa-->0.934, above the pre-registered 0.5 switch (results/validation_ledger.tsv, row 1). Labels are read only through `src/nextgem/labels.py`, fold-scoped. Validation and test labels are sealed.

**Gene universe.** U is the <!--n:genes.U-->2,866 Human1 or Human2 genes that have a CRISPR column. U∩, the <!--n:genes.Ucap-->2,766 genes in both parents, is a sensitivity set (results/design/design_numbers.json).

## 3. Methods as run

**3.1 Scores and statistics.** The knockout (KO) score is e = 1 − ratio_r, with ratio_r = clip(round(KO/WT biomass flux, 6), 0, 1).
- **Metric.** Generic scores are near-binary, so we use tAP, the exact expected AP over random orderings within tied blocks (McSherry & Najork [92]).
- **Uncertainty.** A lineage-stratified line bootstrap, with methods paired on identical gene × line sets. The test looks add a sign-flip permutation (DESIGN §6.3, §6.5).

**3.2 Baselines** (DESIGN §8).
- B0: constant.
- B2: gene prior (the fraction of fitting lines in which the gene is positive). B2m: the same within medium class.
- B3: expression.
- B6: gradient boosting on gene priors, omics, paralog and neighbour features, and medium class.
- B8+ / B8−: the generic score with ties broken by higher or lower expression. B8′: the generic score with unexpressed genes zeroed.

B8 measures what any pipeline gains by ordering the tied block with an expression-correlated quantity.

**3.3 W1.** Single-gene KOs of both parents on Ham's medium (M0). The binary endpoints follow the Human1 [1] and Human2 [2] papers' call rules: P57 (the KO fails ≥ 1 of 57 metabolic tasks) and PGR (biomass only), against L-ge06 (DESIGN §6.8). Our P57 calls equal Human-GEM's `taskEssentialGenes` exactly (results/w1/README.md).

**3.4 W2, medium.** Each line received its documented screen medium, parsed into components. The medium library was built by two blind agents and frozen (A9). The variant ladder:
- M-open: all uptakes open (negative control).
- M0: Ham's medium.
- M1: base medium plus supplements and a non-limiting micronutrient floor.
- M2(f): M1 plus a serum model at f × requirement.
- HPLM and Plasmax: physiological counterfactuals.

The negative control shuffles media within lineage × serum status. The calibrated regime M3 was excluded as uncalibrated, because growth tracked the O2 bound (A9).

**3.5 W3, tasks, objectives, biomass.** The modes:
- Task criterion modes, including T5, a CellFie-inspired [50] line-specific set.
- Objectives, including MOMA [64].
- Demand constraints.
- Biomass leave-one-out (B-loo) of each reactant or cofactor-pool group, and biomass add-ons.

Pseudo-tasks and a matched-call-rate control are the negative controls; the mevalonate core is a positive control. Amendment A11 recorded three rulings made after training aggregates were seen:
- B-loo modes qualify on TG alone, because the call-rate control is defined only for added calls (flag **W3-D1**).
- The cofactor-pool leave-out carries **W3-MEV**: under serum it makes the mevalonate core non-lethal.
- The metabolite-pool leave-out was withdrawn before its look.

**3.6 W4, context-specific models.** Stage A built five families on a 30-line training panel (10 lineages × 3), each at three thresholds, with RNA or protein and with or without a genomic-loss constraint:
- F0: GPR thresholding with repair.
- F2: GIMME [34].
- F3: FASTCORE [40].
- F4: E-Flux [43].
- F5: label-calibrated presence.

CORDA cost 1.5–1.7 CPU-h per build, so it ran as a 10-line descriptive comparison (A7, D8). A fixed repair rule restores growth and the viability tasks. A4 fixed the implementation details before any build was scored, and A10 fixed time limits. Method vs tissue was tested by PERMANOVA against permuted-expression builds, the failure mode reported by [7].

**3.7 W5, GPR curation.** For the 200 GPR groups where single edits flip most calls, evidence edits were scripted:
- OR → AND for co-complex members;
- AND → OR for paralogs that share no complex;
- removal of branches absent pan-cancer.

W5b added training-fold co-essentiality. Human1's reversible thymidylate synthase reaction was made irreversible as a separate item.

**3.8 M5 and M6, integration.** Only W3 passed TG, so A12 reduced integration to two steps. M5 was a panel pre-check: a context family would continue only if it added ≥ 0.005 on top of Human1 + B-loo. M6 was a 2 × 2 factorial, {H1, H1 + B-loo} × {M0, M2(f = 1)}, scored on training and then in one validation look. Its candidate order was corrected to DESIGN §9's "medium < tasks < combination" before scoring (D10; notes/reviews/M6_harness_and_spec_review.md).

**3.9 W6, flux information.** Ridge models on matched-capacity PCA blocks, scored by LOLO, tested whether pFBA flux or KO profiles of the G\* and P1 arms add to RNA or protein. The three analyses were subtypes, a priori metabolic drugs and metabolomics. The decision rule was a Holm-corrected Wilcoxon test, with permuted and synthetic controls.

**3.10 W7, auto-research.** Gene-prior ML hybrids with GEM features were screened by nested LOLO against the current best. Every idea was logged before it ran (autoresearch/ideas.tsv). The loop had a 4 h time box (A12) and stopped after 10 consecutive failures. P2 is the best validated idea, or the simpler of two within one paired SE (D13).

**3.11 W8, resource.** P1 has no fitted parameters, so P1 models differ between lines only by medium. One SBML model was built per (parent, medium signature), with six gates (results/w8/README.md).

**3.12 Leakage control and verification** (DESIGN §11–12).
- **Label access.** Labels pass only through `labels.py`. A PreToolUse hook blocks direct reads and is probed at every session start, and a runtime tripwire and a label-copy scanner back it up.
- **Looks.** FROZEN.tsv SHA pins are checked before each look; the ledgers are hash-chained; validation looks reveal scalars only.
- **Exposure.** Every agent that saw label aggregates is listed in an exposure register (notes/label_exposure.md).
- **Review.** Every result relied on was recomputed by a non-author agent (notes/reviews/).

The GPR transform of Machado et al. (2016) was assessed and not adopted (A8). Its gene-level KOs equal the Boolean KOs exactly (notes/methods/GPRTransform_assessment.md).

## 4. Results

### 4.1 Generic models: Human1 beats Human2, and both are far below a gene prior

The table gives lineage-equal tAP_line on L_prob, all U genes, expression lines; "minus CE_fit" removes the fold's common essentials.
- **Training:** LOLO on 686 lines in 12 lineages (results/w1/train_grid.tsv).
- **Validation:** the single W1 look, 198 lines in 4 lineages (results/w1/val_look_W1_M2_baseline.tsv).

| Method | Training | Training, minus CE_fit | Validation | Validation, minus CE_fit |
|---|---|---|---|---|
| Generic Human1 (G\*) | <!--n:w1.tr.h1-->0.300 | <!--n:w1.tr.h1.noce-->0.145 | <!--n:w1.val.h1-->0.302 | <!--n:w1.val.h1.noce-->0.149 |
| Generic Human2 | <!--n:w1.tr.h2-->0.292 | <!--n:w1.tr.h2.noce-->0.142 | <!--n:w1.val.h2-->0.295 | <!--n:w1.val.h2.noce-->0.147 |
| Human1 + P57 task calls | <!--n:w1.tr.h1task-->0.285 | <!--n:w1.tr.h1task.noce-->0.155 | <!--n:w1.val.h1task-->0.287 | <!--n:w1.val.h1task.noce-->0.158 |
| B0 constant (null) | <!--n:w1.tr.b0-->0.124 | <!--n:w1.tr.b0.noce-->0.071 | <!--n:w1.val.b0-->0.126 | <!--n:w1.val.b0.noce-->0.073 |
| B3 expression | <!--n:w1.tr.b3-->0.326 | <!--n:w1.tr.b3.noce-->0.204 | <!--n:w1.val.b3-->0.331 | <!--n:w1.val.b3.noce-->0.208 |
| B8+ (G\*, ties by higher expression) | <!--n:w1.tr.b8p-->0.424 | <!--n:w1.tr.b8p.noce-->0.214 | <!--n:w1.val.b8p-->0.430 | <!--n:w1.val.b8p.noce-->0.222 |
| B8− (ties by lower expression) | <!--n:w1.tr.b8m-->0.224 | <!--n:w1.tr.b8m.noce-->0.108 | <!--n:w1.val.b8m-->0.224 | <!--n:w1.val.b8m.noce-->0.109 |
| B8′ (unexpressed genes zeroed) | <!--n:w1.tr.b8q-->0.307 | <!--n:w1.tr.b8q.noce-->0.150 | <!--n:w1.val.b8q-->0.310 | <!--n:w1.val.b8q.noce-->0.155 |
| B2 gene prior | <!--n:w1.tr.b2-->0.900 | <!--n:w1.tr.b2.noce-->0.686 | <!--n:w1.val.b2-->0.892 | <!--n:w1.val.b2.noce-->0.671 |
| B2m medium-class prior | <!--n:w1.tr.b2m-->0.898 | <!--n:w1.tr.b2m.noce-->0.686 | <!--n:w1.val.b2m-->0.892 | <!--n:w1.val.b2m.noce-->0.676 |
| B6 gradient boosting (M2 defaults) | <!--n:w1.tr.b6-->0.888 | <!--n:w1.tr.b6.noce-->0.673 | not in the W1 look | |

- **H1.2 is not supported.** Human2 − Human1 is <!--n:w1.tr.h2h1-->−0.0079 [−0.0083, −0.0075] on training (results/w1/train_grid.tsv), negative in all 12 lineages. On validation it is <!--n:w1.val.h2h1-->−0.0071 [−0.0079, −0.0063] (results/w1/val_look_W1_M2_baseline_gstar.json), which fixed G\* = <!--n:w1.val.gstar-->Human1. H1.1 holds.
- **H1.3 holds overwhelmingly.** B2 is about three times the generic GEM. Expression alone (B3) beats both generic models, and the B8+ tie-break gains more than any GEM lever in this study.
- **Paper protocols** (results/w1/train_binary.tsv, results/w1/val_look_W1_M2_baseline.tsv).
  - P57 × L-ge06 MCC, Human1 / Human2: <!--n:w1.tr.mcc.h1.p57-->0.329 / <!--n:w1.tr.mcc.h2.p57-->0.331 on training; <!--n:w1.val.mcc.h1.p57-->0.330 / <!--n:w1.val.mcc.h2.p57-->0.333 on validation.
  - PGR (biomass only): <!--n:w1.tr.mcc.h1.pgr-->0.334 / <!--n:w1.tr.mcc.h2.pgr-->0.331 on training; <!--n:w1.val.mcc.h1.pgr-->0.337 / <!--n:w1.val.mcc.h2.pgr-->0.336 on validation.
  - These lie within the published 0.25–0.44 range (LITERATURE §2.1). Under P57 Human2 edges ahead, by a negligible margin.
- **Recall ceiling.** Only <!--n:w1.addr.h2.n-->1,098 U genes are addressable in Human2, meaning their KO closes a reaction. <!--n:w1.addr.h2.out-->42% of positive training gene × line pairs lie outside that set (<!--n:w1.addr.h1.out-->45% for Human1) (results/w1/train_addressable.tsv).

### 4.2 Context-specific models (W4, Stage A): every family negative

Nested panel-LOLO Δ̄_lin vs G\* on 30 training lines in 10 lineages (results/w4/stage_a_analysis/train_selection.json). <!--n:w4.nbuilds-->1,380 builds took <!--n:w4.buildcpu-->108.7 CPU-h with no failure (results/w4/stage_a/MANIFEST.tsv).

| Family | Δ̄_lin vs Human1 [95% CI] | Folds > 0 | vs Human2 (build parent) | Flags |
|---|---|---|---|---|
| F0 thresholding + repair | <!--n:w4.f0-->−0.0080 [−0.0103, −0.0053] | <!--n:w4.f0.folds-->0/10 | <!--n:w4.f0.h2-->−0.0014 | gate not passed |
| F2 GIMME | <!--n:w4.f2-->−0.0067 [−0.0090, −0.0042] | <!--n:w4.f2.folds-->1/10 | <!--n:w4.f2.h2-->−0.0001 | gate not passed; carried to M6 |
| F3 FASTCORE | <!--n:w4.f3-->−0.0102 [−0.0133, −0.0072] | <!--n:w4.f3.folds-->1/10 | <!--n:w4.f3.h2-->−0.0036 | gate not passed |
| F4 E-Flux | <!--n:w4.f4-->−0.0029 [−0.0041, −0.0015] | <!--n:w4.f4.folds-->1/10 | <!--n:w4.f4.h2-->+0.0038 | gate not passed; method-dominated; viability gate failed (units artefact) |
| F5 label-calibrated presence | <!--n:w4.f5-->−0.0114 [−0.0131, −0.0095] | <!--n:w4.f5.folds-->0/10 | <!--n:w4.f5.h2-->−0.0048 | gate not passed; label-calibrated |

CORDA, R1, 10 lines, descriptive: Δ <!--n:w4.corda-->−0.0127 vs Human1 (results/w4/stage_a_analysis/train_corda_a6prime.tsv).

- **No family passes TG**, against Human1 or against its own Human2 parent (results/w4/stage_a_analysis/README.md):
  - less pruning is always better;
  - the genomic constraint changes tAP by at most 0.0002;
  - RNA + protein is worse than RNA alone.
  Expression-driven pruning removes correct generic calls faster than it adds line-specific ones.
- **H4.1 holds at every level: method dominates tissue** (results/w4/stage_a_analysis/h4_tests.json). R² for method vs lineage:
  - reaction content: <!--n:w4.h41.rxn.method-->0.588 vs <!--n:w4.h41.rxn.lineage-->0.052;
  - task profiles: <!--n:w4.h41.task.method-->0.254 vs <!--n:w4.h41.task.lineage-->0.073;
  - binary KO profiles: <!--n:w4.h41.ko.method-->0.199 vs <!--n:w4.h41.ko.lineage-->0.085.
- **H4.2.** F0, F2 and F3 carry tissue signal beyond permuted expression. F4 is method-dominated: its binary KO profile equals the generic one.
- **Why F2 R1 was carried.** F4's viability failure is a units artefact: E-Flux scales growth to about 1, and the GR task needs ≥ 1. The PI kept the frozen gate, so F2 R1, flagged "fails TG", was the configuration carried forward.

### 4.3 GPR curation (W5): negative on both parents

Training LOLO Δ̄_lin vs the generic parent (results/w5/train_summary.json).
- **W5a (literal rules):** <!--n:w5.h1-->−0.0285 [−0.0289, −0.0282] on Human1 and <!--n:w5.h2-->−0.0293 [−0.0295, −0.0291] on Human2, 0/12 folds. This is worse than every matched random-edit control (Human1: <!--n:w5.neg.h1-->−0.0145 to −0.0057).
- **W5b** equals W5a: no co-essential pair reached r ≥ 0.4.
- **Cause.** The AND → OR paralog rule merged distant RNA polymerase I/II/III subunits inside one lumped GPR, which the edit reviewer rejected on biology. With paralogs restricted to ≥ 30% identity, W5a gives <!--n:w5.pid30.h1-->−0.0021 [−0.0024, −0.0018] (Human1) and <!--n:w5.pid30.h2-->−0.0031 [−0.0032, −0.0030] (Human2).
- **Thymidylate synthase fix** (Human1): <!--n:w5.tyms-->+0.0007 [0.0004, 0.0010], <!--n:w5.tyms.folds-->9/12 folds. Consistent, but below the gate.

The curation lever was not selected, so test contrasts (iii) and (v) read "lever not selected".

### 4.4 Medium (W2): a realistic medium helps, but not line-specifically

<!--n:w2.nattr-->466 training lines with a documented medium sheet; Human1; Δ̄_lin vs Ham's M0 (results/w2/train/shortlist.json).
- **Documented medium plus serum, M2(f = 1):** <!--n:w2.m2-->+0.0159 [0.0146, 0.0172], <!--n:w2.m2.folds-->11/12 folds (one fold has too few lines to count), <!--n:w2.m2.noce-->+0.0182 minus CE_fit. It passed TG. Base medium alone (M1) gives <!--n:w2.m1-->+0.0089 (results/w2/train/step_scores.tsv).
- **W2-SHUF.** Documented minus shuffled medium is <!--n:w2.shuf-->+0.0005 [−0.0000, 0.0010]; a shuffled medium keeps <!--n:w2.shuf.share-->97% of the gain (results/w2/train/shuffled_control.tsv). The gain means "a realistic culture medium instead of Ham's", not line-specific medium: Ham's hypoxanthine and thymidine mask purine and pyrimidine synthesis, while serum carnitine, cholesterol and glutathione rescue artefact calls. This direction agrees with the medium effects reported in DepMap [10]. mΔAP_gene ≈ 0.
- **W2-MOPEN.** The open-medium negative control gives <!--n:w2.mopen-->+0.0195 on all genes but <!--n:w2.mopen.noce-->−0.0168 minus CE_fit (results/w2/train/step_scores.tsv).
- **Physiological media** score below the documented medium: HPLM <!--n:w2.hplm-->+0.0116, Plasmax <!--n:w2.plasmax-->+0.0100.
- **Error attributable to medium** (results/w2/train/attribution.json). The documented medium reaches <!--n:w2.oracle-->0.979 of a per-line training oracle's gain, but net-corrects only <!--n:w2.netflip-->2.7% of the generic model's binary errors. A small attributable share was pre-registered as a plausible outcome.
- **W2-CONST** (post hoc, descriptive). One constant medium, RPMI + 10% FBS, scored slightly above the documented media (results/w2/README.md).

### 4.5 Tasks, objectives and biomass (W3): one structural fix works

Training LOLO Δ̄_lin vs the generic parent, Human1, M0 (results/w3/train/shortlist.json, results/w3/train/mode_scores.tsv).
- **B-loo:cofactor_pool_biomass** removes the cofactor/vitamin pool from biomass MAR13082.
  - Δ̄_lin <!--n:w3.tr.bloo-->+0.0343 [0.0336, 0.0349], <!--n:w3.tr.bloo.folds-->12/12 folds; <!--n:w3.tr.bloo.noce-->+0.0213 minus CE_fit; Human2 <!--n:w3.tr.bloo.h2-->+0.0369.
  - It removes biomass-artefact false positives (retinoid and vitamin-D catabolism, carnitine, BH4, CoQ).
  - It also loses some genuine requirements (heme, PDXK, TPK1, RFK).
- **W3-D1.** No pre-registered negative control applies. The variant loses to a B2-informed symmetric removal comparator (<!--n:w3.tr.sym-->−0.020).
- **W3-MEV.** Under serum, HMGCR, MVK, PMVK and MVD become non-lethal, because ubiquinol leaves with the pool.
- **Metabolite-pool leave-out:** <!--n:w3.tr.metpool-->+0.0056, withdrawn before the look. The winner's-curse term over 134 modes, <!--n:w3.wc-->0.049, exceeds every training effect.
- **Task modes lower all-gene tAP but raise it minus CE_fit:**
  - T3: <!--n:w3.tr.t3-->−0.0149 (<!--n:w3.tr.t3.noce-->+0.0094 minus CE_fit);
  - T5: <!--n:w3.tr.t5-->−0.0288 (<!--n:w3.tr.t5.noce-->+0.0172).
  All task modes lose to the matched-call-rate control.
- **Objectives do nothing:** MOMA <!--n:w3.tr.moma-->−0.0022; task capacity <!--n:w3.tr.ocap-->+0.0013; biosynthetic demands <!--n:w3.tr.dbs-->−0.0012.
- **Validation look 2: gate passed.** B-loo on Human1 gives <!--n:w3.val-->+0.0322 [0.0308, 0.0334] vs G\* (tAP <!--n:w3.val.tap-->0.334) and <!--n:w3.val.noce-->+0.0191 [0.0175, 0.0205] minus CE_fit (results/w3/val/val_look_W3_M4_shortlist_decision.json, results/w3/val/val_look_W3_M4_shortlist.tsv). The training estimate replicated: the first lever to beat G\* on held-out lineages.

### 4.6 Context on top of the biomass fix (M5): small, not selected

On the 30-line panel, on top of Human1 + B-loo (results/m5/train/train_decision.json):
- GIMME adds <!--n:m5.f2-->+0.0031 [0.0017, 0.0045] (<!--n:m5.f2.folds-->8/10 folds).
- E-Flux adds <!--n:m5.f4-->+0.0044 [0.0035, 0.0052] (<!--n:m5.f4.folds-->10/10 folds; viability gate failed; method-dominated).

Both are below the +0.005 bar, so Stage B and C were skipped. The families turn slightly positive on top of the biomass fix, which suggests that extraction partly removes the same artefacts.

### 4.7 Integration factorial (M6) and P1

| Cell | Training Δ̄_lin vs G\* | Training tAP | Validation Δ̄_lin vs G\* | Validation tAP |
|---|---|---|---|---|
| G\* = Human1 / M0 | — | <!--n:m6.tr.tap.gstar-->0.300 | — | <!--n:m6.val.tap.gstar-->0.302 |
| (1) Human1 + B-loo / M0 | <!--n:m6.tr.c1-->+0.0343 [0.0336, 0.0349] | <!--n:m6.tr.tap.c1-->0.334 | <!--n:m6.val.c1-->+0.0322 [0.0308, 0.0334] | <!--n:m6.val.tap.c1-->0.334 |
| (2) Human1 / M2(f = 1) | <!--n:m6.tr.c2-->+0.0160 [0.0152, 0.0168] | <!--n:m6.tr.tap.c2-->0.316 | <!--n:m6.val.c2-->+0.0171 [0.0156, 0.0186] | <!--n:m6.val.tap.c2-->0.319 |
| (3) Human1 + B-loo / M2(f = 1) = **P1** | <!--n:m6.tr.c3-->+0.0534 [0.0523, 0.0545] | <!--n:m6.tr.tap.c3-->0.353 | <!--n:m6.val.c3-->+0.0517 [0.0499, 0.0538] | <!--n:m6.val.tap.c3-->0.353 |
| Human2 + B-loo / M2(f = 1), descriptive | | <!--n:m6.tr.tap.h2c3-->0.348 | <!--n:m6.val.h2c3-->+0.0460 [0.0441, 0.0480] | <!--n:m6.val.tap.h2c3-->0.348 |

Sources: results/m6/train/summary.json and results/m6/train/cells.tsv (training, 12/12 folds positive for every Human1 cell); results/m6/val/val_look_M6_factorial_decision.json and results/m6/val/val_look_M6_factorial.tsv (validation).

- **Selection (look 3, D11).** All three cells passed the validation gate, and only (3) was within one paired SE of the best, so P1 = <!--n:m6.val.selected-->M6:H1+Bloo/M2(f=1) (results/m6/val/val_look_M6_factorial_decision.json).
  - Winner's-curse term <!--n:m6.val.wc-->0.0015; adjusted estimate <!--n:m6.val.adj-->+0.0503.
  - On validation P1 also gains <!--n:m6.val.p1.noce-->+0.0395 [0.0367, 0.0424] minus CE_fit and <!--n:m6.val.p1.expr-->+0.0478 [0.0459, 0.0499] on expressed pairs. Its per-gene ranking across lines is unchanged (mΔAP_gene <!--n:m6.val.p1.mdap-->−0.0003) (results/m6/val/val_look_M6_factorial.tsv).
- **Interaction.** The two levers are slightly super-additive on training: <!--n:m6.tr.inter-->+0.0031 [0.0027, 0.0034] above the additive prediction of <!--n:m6.tr.additive-->+0.0503 (results/m6/train/interaction.tsv).
- **Flags.** P1 carries W3-D1, W3-MEV (on every serum medium, i.e. nearly every line), W2-SHUF and W2-MOPEN. Documentation slip: the look's rule string lists the pre-correction candidate order. Only one candidate was within one SE, so selection is unaffected (D11).
- **Against the baselines** (validation; descriptive, because the numbers come from two looks on the same 198 lines). P1 (<!--n:m6.val.tap.c3-->0.353, results/m6/val/val_look_M6_factorial.tsv) is below B8+ (<!--n:w1.val.b8p-->0.430) and far below B2 (<!--n:w1.val.b2-->0.892) (results/w1/val_look_W1_M2_baseline.tsv).
- **Paper protocols, training** (results/m6/train/binary.tsv, results/m6/train/binary_paired.tsv).
  - P57 × L-ge06 MCC: <!--n:m6.tr.mcc.p1.p57-->0.351 (Δ <!--n:m6.tr.mcc.d.p57-->+0.0227 [0.0223, 0.0230] vs G\*).
  - Fixed cutoff ratio_r < 0.01: <!--n:m6.tr.mcc.p1.fix-->0.392 vs <!--n:m6.tr.mcc.gstar.fix-->0.334 (Δ <!--n:m6.tr.mcc.d.fix-->+0.0577 [0.0563, 0.0590]).
  - The literal P57 protocol cannot see the medium, because each task uses its own inputs.

### 4.8 Confirmatory chain on the test tier

| Contrast | Estimate [95% CI] | One-sided p | Claim (§6.4) |
|---|---|---|---|
| C1, P1 − G\*, all genes | [T1: Δ̄_lin, CI — filled by scripts/report_build.py after T1] | [T1: p] | [T1: tier; "conditional on the 7 test lineages"] |
| C1b, P1 − G\*, minus CE_fit | [T1: Δ̄_lin, CI] | [T1: p] | [T1: tier vs SESOI_C1b] |
| C3, P1 per-gene ranking vs expression | [T1: mΔAP_gene difference, CI, permutation percentile] | [T1: pass/fail] | [T1: "detectable" or not; expressed-pairs sign] |
| C4, P2 − G\* | [T2: Δ̄_lin, CI] | [T2: p] | [T2: tier; confirmatory only if C3 passed] |
| C5, P2 − B6 (pre-registered M2 configuration) | [T2: Δ̄_lin, CI] | [T2: p] | [T2: tier] |

The T1 look will also report the following:
- [T1: transfer — lineage t-interval (df 6), exact 2^7 sign-flip p, test lineages with Δ > 0, transfer wording].
- [T1: estimates — P1 vs B2, B2m and B2c; P1 vs B8+, B8− and B8′; within-parent contrasts; Δ̄_line; headroom closed, if its denominator's CI excludes 0; failed-build counts].
- [T1: attribution contrasts (i) P1 vs G\*, (ii) medium alone, (iv) constraints alone; (iii) and (v) "lever not selected"].

**Reading C4 and C5** (D12; notes/reviews/M7_W7_review.md).
- C4 compares gene-prior ML with generic FBA. Because B2 − G\* is about +0.6 on training, C4 is expected to pass, and it is no evidence for GEMs.
- C5 measures only the monotone constraint that separates P2 from the pre-registered B6.

### 4.9 Binary endpoints vs the Human1/Human2 papers

Training and validation values are in §4.1 and §4.7. [T1: MCC, F1, precision and recall on test for G\*, generic Human2, P1(H1) and P1(H2), under P57 × L-ge06, PGR × L-ge06 and the fixed cutoff, with paired differences vs G\* — filled by scripts/report_build.py after T1].

### 4.10 Breakdowns per gene, lineage and task

[T2 unseal: per-lineage tAP, Δ vs G\* and lift for the 7 test lineages; per-gene tAP for P1, G\*, B2 and B3; the 57-task table; library and organoid strata. These are sealed at T1 and revealed after P2 is frozen]. Lineages are never ranked by raw AP (DESIGN §6.2).

### 4.11 Flux states carry no information beyond omics and medium (W6)

The verdict is "not supported" for both hypotheses, both arms and all three analyses (results/w6/provisional/analysis/decisions.tsv). The run is labelled provisional because it predates the P1 freeze, but its Human1 + B-loo arm is exactly P1 (results/w6/README.md).

| Analysis (units) | RNA + flux − RNA, median Δ, G\* / P1 | Holm p, G\* / P1 | Flux − mapped RNA, median Δ, G\* / P1 |
|---|---|---|---|
| Subtypes (<!--n:w6.n.subtypes-->150) | <!--n:w6.gstar.flux.subtypes.med-->0.000 / <!--n:w6.p1.flux.subtypes.med-->0.000 | <!--n:w6.gstar.flux.subtypes.p-->1.00 / <!--n:w6.p1.flux.subtypes.p-->1.00 | <!--n:w6.gstar.sim.subtypes.med-->−0.044 / <!--n:w6.p1.sim.subtypes.med-->−0.044 |
| Drugs (<!--n:w6.n.drugs-->106) | <!--n:w6.gstar.flux.drugs.med-->0.000 / <!--n:w6.p1.flux.drugs.med-->0.000 | <!--n:w6.gstar.flux.drugs.p-->0.38 / <!--n:w6.p1.flux.drugs.p-->1.00 | <!--n:w6.gstar.sim.drugs.med-->−0.239 / <!--n:w6.p1.sim.drugs.med-->−0.240 |
| Metabolomics (<!--n:w6.n.metabolomics-->225) | <!--n:w6.gstar.flux.metabolomics.med-->−0.007 / <!--n:w6.p1.flux.metabolomics.med-->−0.008 | <!--n:w6.gstar.flux.metabolomics.p-->1.00 / <!--n:w6.p1.flux.metabolomics.p-->1.00 | <!--n:w6.gstar.sim.metabolomics.med-->−0.399 / <!--n:w6.p1.sim.metabolomics.med-->−0.398 |

- **Positive control.** The synthetic positive control is recovered: median Δρ <!--n:w6.gstar.synth-->+0.42 for G\* and <!--n:w6.p1.synth-->+0.45 for P1 (results/w6/provisional/analysis/secondary.tsv).
- **Cause of the null.** Generic-parent flux depends only on the assay medium: one state for all of PRISM, two in GDSC, at most 20 for culture media. <!--n:w6.drugs.zero-->77/106 drug units have Δ exactly 0 (results/w6/provisional/analysis/decisions.tsv).
- **Interpretation.** The largest per-unit gains encode medium content (notes/reviews/M7_W6_review.md). The null is a property of line-invariant models, not of the test.

### 4.12 Auto-research (W7): GEM features add nothing to gene-prior ML

Training uses nested 12-fold LOLO on 686 lines (autoresearch/results.tsv). <!--n:w7.n.logged-->17 ideas were logged (autoresearch/ideas.tsv), <!--n:w7.n.scored-->15 were scored and <!--n:w7.n.promoted-->1 was promoted. The loop stopped after 10 consecutive failures. The current best at the start was B6, at tAP <!--n:w7.tr.b6-->0.889.
- **Every GEM-containing idea failed.** Gains over the current best:
  - B6 + GEM features: <!--n:w7.tr.02.d-->−0.0030;
  - the promoted model + GEM features: <!--n:w7.tr.12.d-->+0.0001;
  - sign-constrained GEM features: <!--n:w7.tr.14.d-->+0.0001;
  - stacking P1 with B6: <!--n:w7.tr.04.d-->+0.0000;
  - the artefact-only biomass variant as a feature: <!--n:w7.tr.08.d-->−0.0165.
- **The artefact-only biomass variant itself** (W7-03; post hoc, exploratory) keeps the mevalonate core lethal under serum and scores <!--n:w7.tr.03.tap-->0.330 (<!--n:w7.tr.03.vsp1-->−0.0230 vs P1).
- **Promoted: W7-05,** which is B6 with a monotone constraint on the gene priors.
  - Training: <!--n:w7.tr.05-->0.903, which is <!--n:w7.tr.05.vsb6-->+0.0144 [0.0133, 0.0154] vs B6 and <!--n:w7.tr.05.vsb2-->+0.0028 vs B2.
  - Gene-grouped CV: <!--n:w7.gcv.05-->0.185 (B6 <!--n:w7.gcv.b6-->0.169; results/w7/screen_genecv.tsv). Its accuracy is almost all gene identity.
  - Within-fold expression-neighbour priors gave the largest line-specific additions (<!--n:w7.tr.15.d-->+0.0019, <!--n:w7.tr.17.d-->+0.0023; exploratory), still below the gate.
- **Validation (looks 4–5; D13, D14).**
  - W7-05 passes the gate vs G\* with <!--n:w7.val.05.vsg-->+0.594 [0.590, 0.598].
  - Against the runner-up B6 it scores <!--n:w7.val.05-->0.896 vs <!--n:w7.val.b6-->0.895, a paired difference of <!--n:w7.val.diff-->+0.0007 (SE <!--n:w7.val.se-->0.0002). So P2 = <!--n:w7.val.selected-->W7:W7-05 (results/w7/val/val_look_W7_P2_select_decision.json, results/w7/val/val_look_W7_P2_select.tsv).
  - Its lead over B6 shrank about 20-fold from training to validation.
  - Beside the required baselines: B2 <!--n:w1.val.b2-->0.892, B2m <!--n:w1.val.b2m-->0.892, B8+ <!--n:w1.val.b8p-->0.430 and B8− <!--n:w1.val.b8m-->0.224 (results/w1/val_look_W1_M2_baseline.tsv).
- **D12 wording.** P2 contains no GEM output: it is a tuned non-GEM gene-prior model, and W7's GEM result is negative. The W7 runs used `allow_dirty` against §11.4. The reviewer's bitwise reproduction mitigates this.

## 5. Resource (W8)

| | Value |
|---|---|
| Models (SBML L3V1 FBC v2; Human2 built first) | <!--n:w8.models-->40 = 2 parents × <!--n:w8.signatures-->20 medium signatures |
| Lines mapped (line × parent rows) | <!--n:w8.lines-->1,361 (<!--n:w8.rows-->2,722) |
| Models passing all six gates | <!--n:w8.pass-->39 (Human2 <!--n:w8.h2pass-->19/20, Human1 <!--n:w8.h1pass-->20/20) |
| Models flagged W3-MEV (serum media) | <!--n:w8.mev.models-->38 |
| Size; build cost | <!--n:w8.mb-->78 MB; <!--n:w8.cpu-->2.6 CPU-h, <!--n:w8.wall-->27 min wall |
| MEMOTE total, 4-model sample | Human2 <!--n:w8.memote.h2-->0.64, Human1 <!--n:w8.memote.h1-->0.63 |

Sources: results/w8/p1/summary.json, resource/models.tsv, results/w8/memote/memote_summary.json.

- **What each model records.** The SBML notes give the pipeline, Human-GEM release, medium signature and tier, W3 variant, git commit, input SHA-256s and flags. Closed reactions keep bounds (0, 0), and no reaction was repaired. `scripts/w8_line_sbml.py` writes per-line files on demand (resource/README.md).
- **MEMOTE.** The failing MEMOTE tests are properties of the parent models, not of the P1 edits (results/w8/README.md).
- **One gate failure.** One Human2 model, for a single-line RPMI + testosterone medium, fails the random-KO gate: a lethal KO gives a ratio of about −8e-9 against 0 (results/w8/p1/gates.jsonl). This is probably LP noise just above the 1e-9 absolute floor, but it was not adjudicated. The resource commit message (6a71f7f) says "all gates pass", which the committed summary contradicts.
- **Deviations from DESIGN W8.**
  - Scope: the screened lines with a medium assignment, not the design's 2,030 lines.
  - Packaging: per-signature files instead of per-line files.
  - Timing: built before T1. This is legitimate because P1 has no fitted parameters.
  - FROG was not run.

## 6. Failures and negative results

- **Human2 is not better than Human1** (H1.2).
- **Context-specific extraction** loses to the generic model in all five families and in CORDA. Method dominates tissue.
- **GPR curation:** the literal paralog rule is wrong for lumped multi-subunit GPRs.
- **Medium:** it fails the shuffled-medium control, and the open-medium control passes the gate.
- **Tasks and objectives** do not help on all genes.
- **Flux information** is null, and GEM features add nothing to gene-prior ML.
- **Positive controls that failed because of model structure:**
  - serum cholesterol cannot rescue the post-squalene sterol genes, because biomass needs vitamin D3 products reachable only via lathosterol;
  - uridine cannot rescue CAD, because dihydroorotate sits in the biomass pool;
  - the ETF capability control fails in Human1;
  - P1 loses mevalonate-core essentiality under serum (W3-MEV).
- **Infeasible or incomplete.** M3 calibration was infeasible as specified, and the W2 single-removal flip analysis did not complete (results/w2/README.md).
- **Failure account on test** (DESIGN §10.3). [T1: P1's test errors by pre-registered class, and pFBA traces of the 20 worst lineage × gene cells. Sealed at T1, reported after the T2 unseal].

## 7. Limitations

- **Addressable-set ceiling.** <!--n:w1.addr.h2.out-->42% of positive pairs lie in genes whose KO changes nothing in the parent (results/w1/train_addressable.tsv). Mechanisms outside FBA, such as MTAP/PRMT5, ferroptosis and non-metabolic roles, are invisible to every GEM arm.
- **P1 is line-specific only through the medium,** and W2-SHUF shows that this information is not line-specific either. Its cross-line per-gene ranking is unchanged, which is why the resource holds 40 models.
- **The gene prior dwarfs the GEMs.** B2 scores <!--n:w1.val.b2-->0.892 against <!--n:m6.val.tap.c3-->0.353 for P1 (results/w1/val_look_W1_M2_baseline.tsv, results/m6/val/val_look_M6_factorial.tsv). Most label signal is gene identity shared across lineages, so a GEM gain is a small correction to a ranking the prior already gets mostly right.
- **Near-binary scores.** Results depend on large tied blocks, and the B8+ expression tie-break beats every GEM lever.
- **Chronos vs CERES.** The Human1 paper applied gene effect < −0.6 to CERES (19Q3); we apply it to Chronos (26Q3), whose scale differs. The Human2 paper states no threshold. Our MCCs match the papers' protocol, not their label scale, and neither paper's per-line extraction was reproduced (D2).
- **Selection bias.** P1's <!--n:m6.val.p1.pt3-->+0.052 was chosen on the validation lineages (results/m6/val/val_look_M6_factorial_decision.json). W3's training-tier winner's-curse term (<!--n:w3.wc-->0.049) exceeded its own effects (results/w3/train/shortlist.json). T1 is the unbiased estimate.
- **Post hoc elements.** The W3 control reading, the random-removal control, the constant medium and the artefact-only biomass variant were all decided after seeing results.
- **Labels and bounds.** L_prob and DepMap's common-essential lists are fitted on all lines, and Human2's curation may have been DepMap-informed. Uptake bounds encode availability, not rates.
- **Power.** The test CI is conditional on 7 lineages. Transfer is weakly powered, and the method-vs-tissue panel had 30 lines.
- **Retired seed.** A retired seed's draw was published in round 1, so lineages it named may appear in the test tier. Only aggregate statistics on them were seen (notes/label_exposure.md).

## 8. Runtime and compute

**What was summed.**
- **Included:** `cpu_s` and `wall_s` of every row of every committed results/*/**/MANIFEST.tsv and PROVENANCE.tsv at HEAD, grouped by top directory. Each row is one script run, including reruns, pilots and looks.
- **Excluded:** runs killed before writing a row, design-phase work, review probes, the baseline KO sweep and the GPRTransform assessment (these have no manifest), and agent time.
- **Wall time** is summed run duration, not elapsed time.

| Workstream | Recorded runs | CPU-h | Summed run wall-h |
|---|---|---|---|
| W1 generic benchmark | <!--n:rt.w1.runs-->12 | <!--n:rt.w1.cpu-->7.0 | <!--n:rt.w1.wall-->1.4 |
| W2 medium | <!--n:rt.w2.runs-->22 | <!--n:rt.w2.cpu-->19.8 | <!--n:rt.w2.wall-->5.3 |
| W3 tasks and biomass | <!--n:rt.w3.runs-->25 | <!--n:rt.w3.cpu-->20.9 | <!--n:rt.w3.wall-->6.9 |
| W4 context models (Stage A, smoke tests) | <!--n:rt.w4.runs-->6 | <!--n:rt.w4.cpu-->130.9 | <!--n:rt.w4.wall-->8.0 |
| W5 GPR curation | <!--n:rt.w5.runs-->17 | <!--n:rt.w5.cpu-->0.7 | <!--n:rt.w5.wall-->0.5 |
| M5 panel pre-check | <!--n:rt.m5.runs-->2 | <!--n:rt.m5.cpu-->3.0 | <!--n:rt.m5.wall-->0.4 |
| M6 factorial, T1 preparation, look | <!--n:rt.m6.runs-->18 | <!--n:rt.m6.cpu-->1.5 | <!--n:rt.m6.wall-->0.3 |
| W6 flux information | <!--n:rt.w6.runs-->9 | <!--n:rt.w6.cpu-->2.6 | <!--n:rt.w6.wall-->1.3 |
| W7 auto-research | <!--n:rt.w7.runs-->23 | <!--n:rt.w7.cpu-->8.4 | <!--n:rt.w7.wall-->2.2 |
| W8 resource, MEMOTE | <!--n:rt.w8.runs-->4 | <!--n:rt.w8.cpu-->3.7 | <!--n:rt.w8.wall-->1.2 |
| **Total** | <!--n:rt.total.runs-->138 | <!--n:rt.total.cpu-->198.6 | <!--n:rt.total.wall-->27.5 |

Source: results/*/**/MANIFEST.tsv and PROVENANCE.tsv at HEAD, summed by `runtime_table` in `scripts/report_build.py`.

- **Elapsed time.** <!--n:rt.elapsed-->40.4 h from the first commit to the last commit that changed results/ (git log), against a target of 1–1.5 days (D2). The process review blamed review and fix rounds, not compute (notes/milestone_reviews/M3M4_process_review.md).
- **Stage A dominates.** Its builds took <!--n:w4.buildcpu-->108.7 CPU-h (results/w4/stage_a/MANIFEST.tsv); the process total also counts overhead.
- **W7's ledger** records <!--n:rt.budget.w7-->8.8 CPU-h, including aborted runs (results/budget.tsv).
- **Unrecorded losses, quoted from READMEs and not checked:** about 8.5 CPU-h of killed W3 runs, 1.3 of a stalled W2 run, and 7 of a BLAS-thread incident in W6 (results/w3/README.md, results/w2/README.md, RESEARCH_TREE.md item 8.1).
- **Cost per unit of accuracy** favours the gene prior: B2 costs seconds, while any GEM arm costs hours of LP solves.

## 9. Reproducibility

- **Environment.** `.venv` (Python 3.11, requirements.txt pinned; cobra 0.32.1), Gurobi 13 through a slot pool or GLPK, declared per family. MEMOTE runs in a separate `.venv_memote`.
- **Pins and ledgers.** FROZEN.tsv pins frozen code, inputs and values by SHA-256 and is checked at commit and before every look. The ledgers are hash-chained (results/validation_ledger.tsv, results/test_ledger.tsv).
- **Provenance.** Every production script writes a MANIFEST row: commit, dirty flag, input and output SHA-256, versions, seeds, solver, workers, time and RSS. Long jobs ran from pinned clean worktrees. W3 and W7 ran from a dirty tree; reviewers reproduced their outputs bitwise.
- **Seeds.** SEEDS = [11, 23, 37, 41, 53, 67, 79, 83, 97, 101] plus the split seed (results/splits/SEED_SOURCE.md). A test forbids unseeded randomness.
- **Tests.** `.venv/bin/python -m pytest -q tests/` runs known-answer metric tests, leakage-guard exploit tests and toy tests of each look harness.
- **This report.** `scripts/report_build.py --check` regenerates every marked number from HEAD, and `--update` rewrites them. `--require-filled` fails while test placeholders remain. `tests/test_report_build.py` tests the marker logic.
- **Data.** data/ is not committed. DATA.md lists the source, version, licence and checksum of every input.

## 10. Path taken

**Decisions** (RESEARCH_TREE.md).
- **D0–D2.** Human1 := v1.19.0 (last 1.x release) and Human2 := v2.1.1. ftINIT dropped; paper-protocol MCC/F1 added; compute target about 1–1.5 days.
- **D3–D5.** Design frozen after three review rounds; round 1 alone produced 34 PI decisions, among them tAP, B8 and C1b. Split seed 38 by commit-reveal. The a priori lists were re-derived blind, because the PI had relaxed an authorship rule post hoc (A2).
- **D6–D8.** G\* = Human1. Lighter guard bookkeeping (A6). CORDA kept as a 10-line comparison.
- **D9–D11.** Scope cuts after M4 (A12). The W2 shortlist went forward with its flags. The factorial cells and their order were fixed before scoring. P1 selected.
- **D12–D14.** A non-GEM P2 is admissible, with fixed wording; runner-up comparison; P2 = W7-05.

**Process lessons** (LESSONS.md; notes/milestone_reviews/).
- **Independent recomputation is the core guarantee.** One review found a float-tie AP instability of up to 0.018 that 36 passing checks had missed. The W4 code review caught three defects, any of which would have invalidated Stage A.
- **Fix control readings and eligible sets before looking.** W3-D1 exists because this was not done, and W2's eligible set could not be verified as fixed in advance (notes/reviews/M4_W2_review.md).
- **Read every result against B2 and B8,** not only G\*, and look for structural false positives before adding data.
- **Every long job inherits time limits and pinned worktrees.** W3 lost about 8.5 CPU-h without them. About 40% of wall time went to review and fix rounds (notes/milestone_reviews/M3M4_process_review.md).

**Process deviations recorded in reviews.**
- The user stated seed 38 before the freeze commit, and the PI briefly wrote it into an unpushed commit (A1; notes/reviews/M1_infra_review.md).
- A round-1 reviewer left an all-line label copy in shared scratch (notes/milestone_reviews/M1_process_review.md).
- The W2 agent deleted shared scratch files it had not created (LESSONS.md).
- W2's scoring code was committed only after its first label-reading run (notes/reviews/M4_W2_review.md).
- Stage A sliced training labels outside `cv.lolo`; a fold-scoped rerun reproduced every choice (notes/reviews/M3_W4_stageA_review.md).
- W3 and W7 ran from a dirty tree (notes/reviews/M7_W7_review.md).
- The M6 rule string is stale (D11).
- A frozen tripwire test contradicted the test harness and was corrected under A6 (notes/guard_changelog.md).
- The resource commit message overstates its gate results (§5).

## References

Numbered as in LITERATURE.md; only works cited here are listed.

1. Robinson JL et al. An atlas of human metabolism. *Sci Signal* 2020;13:eaaz1482. doi:10.1126/scisignal.aaz1482.
2. Luo J et al. Reconstruction of human metabolic models with large language models. *PNAS* 2026;123:e2516511123. doi:10.1073/pnas.2516511123.
7. Opdam S et al. A systematic evaluation of methods for tailoring genome-scale metabolic models. *Cell Syst* 2017;4:318–329. doi:10.1016/j.cels.2017.01.010.
10. Lagziel S, Lee WD, Shlomi T. Inferring cancer dependencies on metabolic genes from large-scale genetic screens. *BMC Biol* 2019;17:37. doi:10.1186/s12915-019-0654-4.
22. Dempster JM et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biol* 2021. doi:10.1186/s13059-021-02540-7.
34. Becker SA, Palsson BO. Context-specific metabolic networks are consistent with experiments (GIMME). *PLoS Comput Biol* 2008;4:e1000082. doi:10.1371/journal.pcbi.1000082.
40. Vlassis N, Pacheco MP, Sauter T. Fast reconstruction of compact context-specific metabolic network models (FASTCORE). *PLoS Comput Biol* 2014;10:e1003424. doi:10.1371/journal.pcbi.1003424.
43. Colijn C et al. Interpreting expression data with metabolic flux models (E-Flux). *PLoS Comput Biol* 2009;5:e1000489. doi:10.1371/journal.pcbi.1000489.
50. Richelle A et al. Model-based assessment of mammalian cell metabolic functionalities using omics data (CellFie). *Cell Rep Methods* 2021;1:100040. doi:10.1016/j.crmeth.2021.100040.
55. Cantor JR et al. Physiologic medium rewires cellular metabolism and reveals uric acid as an endogenous inhibitor of UMP synthase (HPLM). *Cell* 2017. doi:10.1016/j.cell.2017.03.023.
56. Vande Voorde J et al. Improving the metabolic fidelity of cancer models with a physiological cell culture medium (Plasmax). *Sci Adv* 2019;5:eaau7314. doi:10.1126/sciadv.aau7314.
64. Segrè D, Vitkup D, Church GM. Analysis of optimality in natural and perturbed metabolic networks (MOMA). *PNAS* 2002;99:15112. doi:10.1073/pnas.232349399.
92. McSherry F, Najork M. Computing information retrieval performance measures efficiently in the presence of tied scores. *ECIR 2008*, LNCS 4956.

Not in LITERATURE.md, cited as given in DESIGN §17 A8: Machado et al. 2016. *PLoS Comput Biol*. doi:10.1371/journal.pcbi.1005140.
