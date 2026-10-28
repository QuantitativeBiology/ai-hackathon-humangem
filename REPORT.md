# Context-specific genome-scale metabolic models and cancer gene essentiality: a pre-registered benchmark on DepMap 26Q3

> **Read first.** This article was written by the AI pipeline. The human-written disclaimer at the top of [README.md](README.md) applies to it: the project is an exercise in agentic AI for bioinformatics (Beltrao-lab hackathon; INESC-ID & IST), no result has yet been verified by a human, and the text is not a traditional publication.

> **Status, 2026-10-28: complete.** Training, validation and both test looks (T1, T2) are done, and the sealed T1 outputs have been unsealed. Every number taken from a committed table is regenerated from the file cited next to it by `scripts/report_build.py --check --require-filled`, which fails if the text and the file disagree or if a test placeholder remains. A few counts quoted from README prose are cited but not machine-checked.

## Abstract

**Background.** Published tests of genome-scale metabolic models (GEMs) against CRISPR screens report binary metrics after choosing among many settings, and none holds out whole tissues (LITERATURE.md §2.1). We asked whether medium-aware, curated, context-specific versions of Human1 (Human-GEM v1.19.0) and Human2 (v2.1.1) predict metabolic gene essentiality in DepMap 26Q3 better than the generic models on lineages they have never seen.

**Design.** The benchmark was pre-registered (DESIGN.md v3, frozen at 61b0beb; amendments A1–A15). <!--n:split.n.total-->1,361 screened lines were split by lineage into 7 test, 4 validation and 12 training lineages (results/splits/split_meta.json). The primary metric is tie-aware average precision (tAP) per line, averaged with equal weight per lineage. Four levers were screened on training: medium, biomass and tasks, context-specific extraction, and GPR curation. At most 40 scalar-only validation looks were allowed; <!--n:looks.used-->7 were used, including the T1 rehearsal (results/validation_ledger.tsv).

**Selection.** Generic Human1 beat Human2 on validation (tAP <!--n:w1.val.h1-->0.302 vs <!--n:w1.val.h2-->0.295), so G\* = Human1 (results/w1/val_look_W1_M2_baseline.tsv, results/w1/val_look_W1_M2_baseline_gstar.json). No context-specific family beat G\* on training. Evidence-based GPR curation lowered tAP. A documented culture medium helped only as much as a shuffled one. Removing the cofactor/vitamin pool from the biomass reaction was the one structural fix that helped. The selected pipeline P1 (Human1, biomass pool removed, documented medium) carries four flags: W3-D1, W3-MEV, W2-SHUF and W2-MOPEN.

**Test (7 locked lineages).**
- **Headline (C1).** P1 improves lineage-averaged tAP over G\* by <!--n:t1.c1-->+0.0497 [0.0482, 0.0513] (sign-flip p <!--n:t1.c1.p-->< 0.001), i.e. it **<!--n:t1.c1.claim-->meaningfully improves**, conditional on the 7 test lineages (results/m6/test/test_look_T1_gates.json; flags W3-D1, W3-MEV, W2-SHUF, W2-MOPEN). P1 is generic Human1 with a biomass edit and the documented medium: no extraction and no curation.
- **C1b.** Without common essentials the gain is <!--n:t1.c1b-->+0.0374 [0.0351, 0.0397], which also **<!--n:t1.c1b.claim-->meaningfully improves** against its lift-scaled SESOI.
- **Transfer.** The lineage t-interval <!--n:t1.tr.t-->[0.0474, 0.0521] excludes 0 and Δ > 0 in <!--n:t1.tr.k-->7/7 lineages, so P1 "generalises to new tissues", conditional on these 7 lineages (results/m6/test/test_look_T1_gates.json).
- **C3 failed.** P1 does not rank genes across lines better than expression, so it shows no detectable line-specific signal beyond expression.
- **C4 and C5 are estimates only.** P2 is non-GEM, so these say nothing about GEMs (D12). It exceeds G\* by an estimated <!--n:t2.c4-->+0.602 [0.598, 0.605] and the pre-registered B6 by <!--n:t2.c5-->+0.0138 [0.0116, 0.0167] (results/m6/test/test_look_T2_gates.json).

**In context.** On test P1 scores tAP <!--n:t1.tap.p1-->0.347 against <!--n:t1.tap.gstar-->0.297 for G\* (results/m6/test/test_look_T1.tsv). That is below the generic model with an expression tie-break (B8+, <!--n:un.w1.b8p.all-->0.418) and far below a label-only gene prior (B2, <!--n:un.w1.b2.all-->0.895) (results/m6/test/unsealed/T1/w1_benchmark.tsv).

**Other findings.**
- Simulated flux states carried no information beyond the omics and the medium.
- GEM features added nothing to gene-prior machine learning, so P2 is a non-GEM model.
- We release <!--n:res.models-->48 SBML models of P1, one per parent and medium signature (P1 has no line-specific parameters), mapped to <!--n:res.lines-->2,114 DepMap lines (resource/models.tsv, results/w8/unscreened/build_summary.json).

**Conclusion.** On unseen lineages, the confirmed GEM gain (flags W3-D1, W3-MEV, W2-SHUF, W2-MOPEN) comes from removing structural false positives of the generic model and from a realistic medium, not from line-specific modelling. The strongest predictor of CRISPR essentiality here is a gene prior learned from the labels, not a GEM.

## 1. Question and pre-registration

**Question.** PROMPT.md asks for optimised Human1/Human2 models that predict DepMap 26Q3 metabolic gene essentiality better than the generic model, judged by AUPRC per gene, lineage and overall. It also asks for an SBML resource, a flux-information test and an account of failure. DESIGN §1 maps these to eight workstreams: W1 baseline, W2 medium, W3 tasks and biomass, W4 context-specific models, W5 GPR curation, W6 flux information, W7 auto-research, W8 resource.

**Pre-registered hypotheses** (DESIGN §9).
- **W1.** H1.1: the generic models beat a random ranking. H1.2: Human2 beats Human1, as its paper reports [2]. H1.3: the gene prior B2 beats generic FBA.
- **W2–W5.** A documented medium beats Ham's medium and a shuffled medium. A priori tasks or biomass edits beat their negative controls. Curated GPRs beat the generic ones.
- **W4 (H4.1–H4.3).** Extraction method explains more variance in model content and KO profiles than lineage. Within each method, lineage explains more than under permuted expression; methods that fail are "method-dominated".
- **W6.** Flux states predict subtypes, drug response and metabolomics beyond the omics they were built from.

**Design.** DESIGN.md v3 was frozen at 61b0beb after three review rounds by seven reviewer agents. Every later change is a dated row in §17 (A1–A14), and every user or PI decision is logged in RESEARCH_TREE.md (D0–D16). The user excluded ftINIT, because it underperformed the generic model in the user's own analysis, and set a compute target of about 1–1.5 days (D2).

**Endpoints and claim language** (DESIGN §6.3–6.4).
- The primary endpoint, C1, is Δ̄_lin(P1 − G\*): the lineage-equal mean of per-line ΔtAP on the 7 locked test lineages.
- Gatekeeping runs C1 → C1b → C3 → C4 → C5, each one-sided at α = 0.025. C1b removes CE_fit (genes positive in ≥ 90% of fitting lines). C3 tests P1's per-gene ranking against expression. C4 and C5 compare P2 with G\* and with B6.
- SESOI is 0.01 for C1, C4 and C5. SESOI_C1b is <!--n:sesoi.c1b-->0.0057, fixed on training before any model was scored (results/w1/sesoi_c1b.json).
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

**Labels** (DESIGN §5; results/labels/build_report.json). The primary label, L_prob, is DepMap gene dependency > 0.5. Chronos [22] gene effect < −0.5 and < −1 are sensitivities, and < −0.6 (L-ge06) serves the Human1-paper protocol. L_prob stayed primary because its agreement with L-ge05 on training was κ = <!--n:labels.kappa-->0.934, above the pre-registered 0.5 switch (results/validation_ledger.tsv, row 1). Labels are read only through `src/nextgem/labels.py`, fold-scoped.

**Gene universe.** U is the <!--n:genes.U-->2,866 Human1 or Human2 genes that have a CRISPR column. U∩, the <!--n:genes.Ucap-->2,766 genes in both parents, is a sensitivity set (results/design/design_numbers.json).

## 3. Methods as run

**3.1 Scores and statistics.** The knockout (KO) score is e = 1 − ratio_r, with ratio_r = clip(round(KO/WT biomass flux, 6), 0, 1).
- **Metric.** Generic scores are near-binary, so we use tAP, the exact expected AP over random orderings within tied blocks (McSherry & Najork [92]).
- **Uncertainty.** A lineage-stratified line bootstrap, with methods paired on identical gene × line sets. The test looks add a sign-flip permutation (10,000 flips) and, for C3, a within-lineage label-permutation null (DESIGN §6.3, §6.5).

**3.2 Baselines** (DESIGN §8).
- B0: constant.
- B2: gene prior (the fraction of fitting lines in which the gene is positive). B2m: the same within medium class.
- B3: expression.
- B6: gradient boosting on gene priors, omics, paralog and neighbour features, and medium class.
- B8+ / B8−: the generic score with ties broken by higher or lower expression. B8′: the generic score with unexpressed genes zeroed.

B8 measures what any pipeline gains by ordering the tied block with an expression-correlated quantity. For test looks, every fitted component is refit on training ∪ validation.

**3.3 W1.** Single-gene KOs of both parents on Ham's medium (M0). The binary endpoints follow the Human1 [1] and Human2 [2] papers' call rules: P57 (the KO fails ≥ 1 of 57 metabolic tasks) and PGR (biomass only), against L-ge06 (DESIGN §6.8). Our P57 calls equal Human-GEM's `taskEssentialGenes` exactly (results/w1/README.md).

**3.4 W2, medium.** Each line received its documented screen medium, parsed into components. The medium library was built by two blind agents and frozen (A9). The variant ladder:
- M-open: all uptakes open (negative control).
- M0: Ham's medium.
- M1: base medium plus supplements and a non-limiting micronutrient floor.
- M2(f): M1 plus a serum model at f × requirement.
- HPLM and Plasmax: physiological counterfactuals.

The negative control shuffles media within lineage × serum status. The calibrated regime M3 was excluded as uncalibrated, because growth tracked the O2 bound (A9).

**3.5 W3, tasks, objectives, biomass.** The candidate task sets and controls were fixed a priori (notes/apriori_lists_final.md, A2). The modes:
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

CORDA (the `corda` 0.5.1 package; LITERATURE.md gives no reference) cost 1.5–1.7 CPU-h per build, so it ran as a 10-line descriptive comparison (A7, D8). A fixed repair rule restores growth and the viability tasks. A4 fixed the implementation details before any build was scored, and A10 fixed time limits. Method vs tissue was tested by PERMANOVA against permuted-expression builds, the failure mode reported by [7, 8].

**3.7 W5, GPR curation.** For the 200 GPR groups where single edits flip most calls, evidence edits were scripted:
- OR → AND for co-complex members;
- AND → OR for paralogs that share no complex;
- removal of branches absent pan-cancer.

W5b added training-fold co-essentiality. Human1's reversible thymidylate synthase reaction was made irreversible as a separate item.

**3.8 M5 and M6, integration.** Only W3 passed TG, so A12 reduced integration to two steps. M5 was a panel pre-check: a context family would continue only if it added ≥ 0.005 on top of Human1 + B-loo. M6 was a 2 × 2 factorial, {H1, H1 + B-loo} × {M0, M2(f = 1)}, scored on training and then in one validation look. Its candidate order was corrected to DESIGN §9's "medium < tasks < combination" before scoring (D10; notes/reviews/M6_harness_and_spec_review.md).

**3.9 W6, flux information.** Ridge models on matched-capacity PCA blocks, scored by LOLO, tested whether pFBA flux or KO profiles of the G\* and P1 arms add to RNA or protein. The three analyses were subtypes, a priori metabolic drugs and metabolomics. The decision rule was a Holm-corrected Wilcoxon test, with permuted and synthetic controls.

**3.10 W7, auto-research.** Gene-prior ML hybrids with GEM features were screened by nested LOLO against the current best. Every idea was logged before it ran (autoresearch/ideas.tsv). The loop had a 4 h time box (A12) and stopped after 10 consecutive failures. P2 is the best validated idea, or the simpler of two within one paired SE (D13).

**3.11 W8, resource.** P1 has no fitted parameters, so P1 models differ between lines only by medium. One SBML model was built per (parent, medium signature), with six gates. Coverage was then extended to the unscreened lines in DESIGN's scope, using W2's frozen medium rule (results/w8/README.md).

**3.12 Test looks** (DESIGN §10.2). A rehearsal ran the frozen test harness on validation, and it counted as a validation look. T1 revealed the C1–C3 scalars and sealed everything else. After the W7 stop and the P2 freeze, T2 revealed C4/C5. The T1 bundle was then unsealed, and the 11 sealed tables were re-derived with matching SHA-256 (D16).

**3.13 Leakage control and verification** (DESIGN §11–12).
- **Label access.** Labels pass only through `labels.py`. A PreToolUse hook blocks direct reads and is probed at every session start; a runtime tripwire and a label-copy scanner back it up.
- **Looks.** FROZEN.tsv SHA pins are checked before each look. The ledgers are hash-chained. Validation looks reveal scalars only.
- **Exposure.** Every agent that saw label aggregates is listed in an exposure register (notes/label_exposure.md).
- **Review.** Training and validation results were recomputed by non-author agents with their own code (notes/reviews/). The test tier could not be recomputed the same way, because the tripwire allows test-label reads only from the frozen `evaluate_test.py`. Four checks stand in for it:
  - the test harness was reviewed and re-checked on synthetic labels (notes/reviews/M6_evaluate_test_review.md);
  - its rehearsal on validation reproduced the independently reviewed M6 look exactly: C1 <!--n:rh.c1-->+0.0517 in the rehearsal (results/m6/rehearsal/val_look_REHEARSAL_P1_M6_gates.json) against <!--n:m6.val.c3.pt-->+0.0517 in the M6 look (results/m6/val/val_look_M6_factorial_decision.json);
  - the PI recomputed the lineage-equal C1 from the unsealed per-lineage values (notes/reviews/M9_T1_partial_recompute.md);
  - W8 had no separate review; its six gates and the reload checks against the M6 sweeps are automated.

The GPR transform of Machado et al. (2016) was assessed and not adopted (A8). Its gene-level KOs equal the Boolean KOs exactly (notes/methods/GPRTransform_assessment.md).

## 4. Results

### 4.1 Generic models: Human1 beats Human2, and both are far below a gene prior

The table gives lineage-equal tAP_line on L_prob, all U genes, expression lines; "minus CE_fit" removes the fold's common essentials.
- **Training:** LOLO on 686 lines in 12 lineages (results/w1/train_grid.tsv).
- **Validation:** the single W1 look, 198 lines in 4 lineages (results/w1/val_look_W1_M2_baseline.tsv).

Test values are in §4.10.

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
  - These lie within the published range of about 0.25–0.45 (LITERATURE.md §0, item 2).
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
- **Cause.** The AND → OR paralog rule merged distant RNA polymerase I/II/III subunits inside one lumped GPR, which the edit reviewer rejected on biology. With paralogs restricted to ≥ 30% identity (a sensitivity declared in code before scoring, not in DESIGN), W5a gives <!--n:w5.pid30.h1-->−0.0021 [−0.0024, −0.0018] (Human1) and <!--n:w5.pid30.h2-->−0.0031 [−0.0032, −0.0030] (Human2).
- **Thymidylate synthase fix** (Human1): <!--n:w5.tyms-->+0.0007 [0.0004, 0.0010], <!--n:w5.tyms.folds-->9/12 folds. Consistent, but below the gate.

The curation lever was not selected, so test contrasts (iii) and (v) read "lever not selected".

### 4.4 Medium (W2): a realistic medium helps, but not line-specifically

<!--n:w2.nattr-->466 training lines with a documented medium sheet; Human1; Δ̄_lin vs Ham's M0 (results/w2/train/shortlist.json).
- **Documented medium plus serum, M2(f = 1):** <!--n:w2.m2-->+0.0159 [0.0146, 0.0172], <!--n:w2.m2.folds-->11/12 folds (one fold has too few lines to count), <!--n:w2.m2.noce-->+0.0182 minus CE_fit. It passed TG. Base medium alone (M1) gives <!--n:w2.m1-->+0.0089 (results/w2/train/step_scores.tsv).
- **W2-SHUF.** Documented minus shuffled medium is <!--n:w2.shuf.precise-->+0.0005 [−0.00003, 0.0010], and a shuffled medium keeps <!--n:w2.shuf.share-->97% of the gain (results/w2/train/shuffled_control.tsv). The gain means "a realistic culture medium instead of Ham's", not line-specific medium. Ham's hypoxanthine and thymidine mask purine and pyrimidine synthesis, while serum carnitine, cholesterol and glutathione rescue artefact calls. Medium-dependent metabolic essentiality is reported in DepMap [10] and in physiological-medium screens, including purine synthesis rescued by hypoxanthine [11]. mΔAP_gene ≈ 0.
- **W2-MOPEN.** The open-medium negative control gives <!--n:w2.mopen-->+0.0195 on all genes but <!--n:w2.mopen.noce-->−0.0168 minus CE_fit (results/w2/train/step_scores.tsv).
- **Physiological media** score below the documented medium: HPLM <!--n:w2.hplm-->+0.0116, Plasmax <!--n:w2.plasmax-->+0.0100.
- **Error attributable to medium** (results/w2/train/attribution.json). The documented medium reaches <!--n:w2.oracle-->0.979 of a per-line training oracle's gain, but net-corrects only <!--n:w2.netflip-->2.7% of the generic model's binary errors.
- **W2-CONST** (post hoc, descriptive). One constant medium, RPMI + 10% FBS, scored slightly above the documented media (results/w2/README.md).

### 4.5 Tasks, objectives and biomass (W3): one structural fix works

Training LOLO Δ̄_lin vs the generic parent, Human1, M0 (results/w3/train/shortlist.json, results/w3/train/mode_scores.tsv).
- **B-loo:cofactor_pool_biomass** removes the cofactor/vitamin pool from biomass MAR13082.
  - Δ̄_lin <!--n:w3.tr.bloo-->+0.0343 [0.0336, 0.0349], <!--n:w3.tr.bloo.folds-->12/12 folds; <!--n:w3.tr.bloo.noce-->+0.0213 minus CE_fit; Human2 <!--n:w3.tr.bloo.h2-->+0.0369.
  - It removes biomass-artefact false positives (retinoid and vitamin-D catabolism, carnitine, BH4, CoQ).
  - It also loses some genuine requirements (heme, PDXK, TPK1, RFK).
- **W3-D1.** No pre-registered negative control applies. Under the symmetric control used in the first scoring run, no W3 mode, this one included, qualified (results/w3/README.md, D-1). The variant also loses to a B2-informed symmetric removal comparator (<!--n:w3.tr.sym-->−0.020).
- **W3-MEV.** Under serum, HMGCR, MVK, PMVK and MVD become non-lethal, because ubiquinol leaves with the pool.
- **Metabolite-pool leave-out:** <!--n:w3.tr.metpool-->+0.0056, withdrawn before the look. The winner's-curse term over 134 modes, <!--n:w3.wc-->0.049, exceeds every training effect.
- **Task modes lower all-gene tAP but raise it minus CE_fit:** T3 <!--n:w3.tr.t3-->−0.0149 (<!--n:w3.tr.t3.noce-->+0.0094 minus CE_fit); T5 <!--n:w3.tr.t5-->−0.0288 (<!--n:w3.tr.t5.noce-->+0.0172). All task modes lose to the matched-call-rate control.
- **Objectives do nothing:** MOMA <!--n:w3.tr.moma-->−0.0022; task capacity <!--n:w3.tr.ocap-->+0.0013; biosynthetic demands <!--n:w3.tr.dbs-->−0.0012.
- **Which contexts (exploratory, training only).** <!--n:w3.ctx.claims-->1,416 of <!--n:w3.ctx.tests-->7,209 context × mode tests reached BH q < 0.05 on training (results/w3/train/context_claims.tsv). Per the W3 README, almost all of them restate the mode's global sign; the 3 reversals are Human2 contexts. None was checked on validation, so no context-specific task claim is made.
- **Validation look 2: gate passed.** B-loo on Human1 gives <!--n:w3.val-->+0.0322 [0.0308, 0.0334] vs G\* (tAP <!--n:w3.val.tap-->0.334) and <!--n:w3.val.noce-->+0.0191 [0.0175, 0.0205] minus CE_fit (results/w3/val/val_look_W3_M4_shortlist_decision.json, results/w3/val/val_look_W3_M4_shortlist.tsv). This was the first lever to beat G\* on held-out lineages.

### 4.6 Context on top of the biomass fix (M5): small, not selected

On the 30-line panel, on top of Human1 + B-loo (results/m5/train/train_decision.json):
- GIMME adds <!--n:m5.f2-->+0.0031 [0.0017, 0.0045] (<!--n:m5.f2.folds-->8/10 folds).
- E-Flux adds <!--n:m5.f4-->+0.0044 [0.0035, 0.0052] (<!--n:m5.f4.folds-->10/10 folds; viability gate failed; method-dominated).

Both are below the +0.005 bar, so Stage B and C were skipped. As a result, no lineage, genomic-background or lineage × background group models were built or evaluated, although PROMPT.md asks for them. Genomic data entered only as the Stage A loss constraint, which changed tAP by at most 0.0002 (§4.2).

### 4.7 Integration factorial (M6) and P1

| Cell | Training Δ̄_lin vs G\* | Training tAP | Validation Δ̄_lin vs G\* | Validation tAP |
|---|---|---|---|---|
| G\* = Human1 / M0 | — | <!--n:m6.tr.tap.gstar-->0.300 | — | <!--n:m6.val.tap.gstar-->0.302 |
| (1) Human1 + B-loo / M0 | <!--n:m6.tr.c1-->+0.0343 [0.0336, 0.0349] | <!--n:m6.tr.tap.c1-->0.334 | <!--n:m6.val.c1-->+0.0322 [0.0308, 0.0334] | <!--n:m6.val.tap.c1-->0.334 |
| (2) Human1 / M2(f = 1) | <!--n:m6.tr.c2-->+0.0160 [0.0152, 0.0168] | <!--n:m6.tr.tap.c2-->0.316 | <!--n:m6.val.c2-->+0.0171 [0.0156, 0.0186] | <!--n:m6.val.tap.c2-->0.319 |
| (3) Human1 + B-loo / M2(f = 1) = **P1** | <!--n:m6.tr.c3-->+0.0534 [0.0523, 0.0545] | <!--n:m6.tr.tap.c3-->0.353 | <!--n:m6.val.c3-->+0.0517 [0.0499, 0.0538] | <!--n:m6.val.tap.c3-->0.353 |
| Human2 + B-loo / M2(f = 1), descriptive | | <!--n:m6.tr.tap.h2c3-->0.348 | <!--n:m6.val.h2c3-->+0.0460 [0.0441, 0.0480] | <!--n:m6.val.tap.h2c3-->0.348 |

Sources: results/m6/train/summary.json and results/m6/train/cells.tsv (training; 12/12 folds positive for every Human1 cell); results/m6/val/val_look_M6_factorial_decision.json and results/m6/val/val_look_M6_factorial.tsv (validation).

- **Selection (look 3, D11).** All three cells passed the validation gate, and only (3) was within one paired SE of the best, so P1 = <!--n:m6.val.selected-->M6:H1+Bloo/M2(f=1) (results/m6/val/val_look_M6_factorial_decision.json). The winner's-curse term is <!--n:m6.val.wc-->0.0015, and the adjusted estimate <!--n:m6.val.adj-->+0.0503.
- **Interaction.** The two levers are slightly super-additive on training: <!--n:m6.tr.inter-->+0.0031 [0.0027, 0.0034] above the additive prediction of <!--n:m6.tr.additive-->+0.0503 (results/m6/train/interaction.tsv).
- **Flags.** P1 carries:
  - W3-D1;
  - W3-MEV, on every serum medium, which means nearly every line;
  - W2-SHUF;
  - W2-MOPEN.
- **Documentation slip.** The look's rule string lists the pre-correction candidate order. Only one candidate was within one SE, so selection is unaffected (D11).

### 4.8 Confirmatory chain on the test tier (T1, T2)

The test tier has <!--n:t1.n.lines-->380 expression lines in 7 lineages, plus <!--n:t1.n.noexpr-->14 lines without expression that are scored for generic models only. Failed builds: <!--n:t1.failed-->0 in every arm (results/m6/test/test_look_T1.tsv).

| Contrast | Estimate [95% CI] | One-sided p | Status and claim (§6.4) |
|---|---|---|---|
| **C1**, P1 − G\*, all genes (headline) | <!--n:t1.c1-->+0.0497 [0.0482, 0.0513] | <!--n:t1.c1.p-->< 0.001 | confirmed; **<!--n:t1.c1.claim-->meaningfully improves**, conditional on the 7 test lineages; flags W3-D1, W3-MEV, W2-SHUF, W2-MOPEN |
| **C1b**, P1 − G\*, minus CE_fit | <!--n:t1.c1b-->+0.0374 [0.0351, 0.0397] | <!--n:t1.c1b.p-->< 0.001 | <!--n:t1.c1b.status-->confirmed; **<!--n:t1.c1b.claim-->meaningfully improves** (SESOI_C1b); flags W3-D1, W3-MEV, W2-SHUF, W2-MOPEN |
| **C3**, P1 per-gene ranking − B3 | <!--n:t1.c3-->−0.0114 [−0.0183, −0.0115] | <!--n:t1.c3.p-->0.64 | not confirmed; <!--n:t1.c3.claim-->no detectable improvement |
| C4, P2 − G\* | <!--n:t2.c4-->+0.602 [0.598, 0.605] | — (not tested) | <!--n:t2.c4.status-->estimate only (C3 failed) |
| C5, P2 − B6 (pre-registered M2 configuration) | <!--n:t2.c5-->+0.0138 [0.0116, 0.0167] | — (not tested) | <!--n:t2.c5.status-->estimate only (C3 failed) |

Sources: results/m6/test/test_look_T1_gates.json, results/m6/test/test_look_T2_gates.json.

- **Transfer (C1).** Δ > 0 in <!--n:t1.tr.k-->7/7 test lineages. The lineage t-interval (df 6) is <!--n:t1.tr.t-->[0.0474, 0.0521], and the exact sign-flip p over 2^7 configurations is <!--n:t1.tr.p-->0.0078 (results/m6/test/test_look_T1_gates.json). So "generalises to new tissues" <!--n:t1.generalises-->holds, conditional on these 7 lineages. Per-lineage gains are in §4.10.
- **Shrinkage.** P1's gain was <!--n:m6.tr.c3-->+0.0534 [0.0523, 0.0545] on training (results/m6/train/summary.json) and <!--n:m6.val.c3-->+0.0517 [0.0499, 0.0538] on validation (results/m6/val/val_look_M6_factorial_decision.json). On test it is <!--n:t1.c1.pt3-->+0.050 (results/m6/test/test_look_T1_gates.json), so the validation selection was not inflated.
- **C3 (results/m6/test/test_look_T1.tsv).**
  - P1's mean per-gene ΔAP over the null is <!--n:t1.c3.mdap-->−0.0002 on <!--n:t1.c3.ngenes-->963 estimable genes, below the permutation threshold of <!--n:t1.c3.null-->0.0001.
  - Expression (B3) ranks genes across lines better, so the contrast is negative. The expressed-pairs version agrees: <!--n:t1.c3.expr-->−0.0117 [−0.0189, −0.0115].
  - The percentile interval excludes the point estimate. The cause (gene dropping or line duplication in the bootstrap) is not established. Both C3 criteria fail independently (CI and permutation), so the verdict does not depend on it.
  - P1 therefore shows **no detectable line-specific information beyond expression**, which is consistent with W2-SHUF.
- **Estimates** (results/m6/test/test_look_T1.tsv).

  | Contrast | Δ̄_lin [95% CI] |
  |---|---|
  | P1 − B2 | <!--n:t1.est.b2-->−0.5477 [−0.5513, −0.5440] |
  | P1 − B2m | <!--n:t1.est.b2m-->−0.5493 [−0.5530, −0.5456] |
  | P1 − B2c | <!--n:t1.est.b2c-->−0.5418 [−0.5453, −0.5381] |
  | P1 − B8+ | <!--n:t1.est.b8p-->−0.0708 [−0.0731, −0.0687] |
  | P1 − B8− | <!--n:t1.est.b8m-->+0.1240 [0.1220, 0.1259] |
  | P1 − B8′ | <!--n:t1.est.b8q-->+0.0423 [0.0408, 0.0438] |
  | P1(H2) − generic Human2, within parent | <!--n:t1.est.wp.h2-->+0.0524 [0.0509, 0.0538] |
  | P1(H2) − G\* | <!--n:t1.est.p1h2-->+0.0440 [0.0424, 0.0454] |

  - "Adds information beyond generic + expression" is **<!--n:t1.adds-->not claimed** (results/m6/test/test_look_T1_gates.json): it requires all three B8 contrasts to be positive, and P1 loses to B8+.
  - Headroom closed, Δ(P1 − G\*)/Δ(B2 − G\*), is <!--n:t1.headroom-->8.3% [8.1%, 8.6%]. Its denominator is <!--n:t1.headroom.den-->+0.5974 [0.5938, 0.6009] (results/m6/test/test_look_T1.tsv).
- **Attribution on test** (results/m6/test/unsealed/T1/attribution.tsv).
  - (i) P1 vs G\* equals C1.
  - (ii) Medium alone, M0 → M2(f = 1) on generic Human1: <!--n:un.attr.ii-->+0.0127 [0.0115, 0.0140]. It carries W2-SHUF: realistic medium, not line-specific medium.
  - (iii) and (v): lever not selected (curation).
  - (iv) Constraints alone equals (i), because P1 has generic GPRs.
- **Sensitivities of C1** (results/m6/test/unsealed/T1/sensitivities.tsv). L-ge05: <!--n:un.sens.ge05-->+0.0503 [0.0485, 0.0519]. L-ge10: <!--n:un.sens.ge10-->+0.0473 [0.0460, 0.0488]. U∩: <!--n:un.sens.ucap-->+0.0539 [0.0522, 0.0554].

**Reading C4 and C5** (D12; notes/reviews/M7_W7_review.md).
- C4 compares gene-prior ML with generic FBA and says nothing about GEMs.
- C5 measures only the monotone constraint that separates P2 from the pre-registered B6.
- Both are estimates, because C3 failed. Each is positive in <!--n:t2.c4.k-->7/7 and <!--n:t2.c5.k-->7/7 lineages respectively (results/m6/test/test_look_T2_gates.json).
- On test, P2 scores tAP <!--n:t2.tap.p2-->0.899 against <!--n:t2.tap.b6-->0.885 for B6 (results/m6/test/test_look_T2_grid.tsv). B2 scores <!--n:un.w1.b2.all-->0.895, B2m <!--n:un.w1.b2m.all-->0.896, B8+ <!--n:un.w1.b8p.all-->0.418 and B8− <!--n:un.w1.b8m.all-->0.223 (results/m6/test/unsealed/T1/w1_benchmark.tsv). These come from two looks on the same lines, so the comparison is descriptive.
- C5's path (training +0.014, validation +0.0007, test +0.014) is discussed in D16 (RESEARCH_TREE.md). Validation compared W7-05 with W7's B6 (W7-01); test compared it with the pre-registered B6 refit on training ∪ validation, so the three numbers are not the same contrast.

### 4.9 Binary endpoints vs the Human1/Human2 papers (test)

MCC on the test expression lines, lineage-equal mean, parent gene set (results/m6/test/test_look_T1.tsv). P1 carries the flags W3-D1, W3-MEV, W2-SHUF, W2-MOPEN. Under P57 the Human2-paper rows are identical to the Human1-paper rows.

| Call rule × label | G\* (Human1) | Human2 | P1(H1) | P1(H2) | P1(H1) − G\* [95% CI] |
|---|---|---|---|---|---|
| P57 × L-ge06 | <!--n:t1.mcc.p57.h1-->0.322 | <!--n:t1.mcc.p57.h2-->0.325 | <!--n:t1.mcc.p57.p1-->0.345 | <!--n:t1.mcc.p57.p1h2-->0.346 | <!--n:t1.mcc.p57.d-->+0.0222 [0.0217, 0.0226] |
| PGR × L-ge06 | <!--n:t1.mcc.pgr.h1-->0.334 | <!--n:t1.mcc.pgr.h2-->0.330 | <!--n:t1.mcc.pgr.p1-->0.379 | <!--n:t1.mcc.pgr.p1h2-->0.381 | <!--n:t1.mcc.pgr.d-->+0.0457 [0.0449, 0.0465] |
| ratio_r < 0.01 × L-ge06 | <!--n:t1.mcc.fix06.h1-->0.334 | <!--n:t1.mcc.fix06.h2-->0.330 | <!--n:t1.mcc.fix06.p1-->0.387 | <!--n:t1.mcc.fix06.p1h2-->0.394 | |
| ratio_r < 0.01 × L_prob | <!--n:t1.mcc.fixp.h1-->0.320 | <!--n:t1.mcc.fixp.h2-->0.316 | <!--n:t1.mcc.fixp.p1-->0.372 | <!--n:t1.mcc.fixp.p1h2-->0.377 | <!--n:t1.mcc.fixp.d-->+0.0522 [0.0503, 0.0539] |

Under P57 × L-ge06, P1's precision rises from <!--n:t1.p57.precision.h1-->0.423 to <!--n:t1.p57.precision.p1-->0.473, while recall falls from <!--n:t1.p57.recall.h1-->0.343 to <!--n:t1.p57.recall.p1-->0.333. F1 goes from <!--n:t1.p57.f1.h1-->0.378 to <!--n:t1.p57.f1.p1-->0.390 (results/m6/test/test_look_T1.tsv). Under PGR, F1 goes from <!--n:t1.pgr.f1.h1-->0.363 to <!--n:t1.pgr.f1.p1-->0.379 (Human2 <!--n:t1.pgr.f1.h2-->0.364, P1(H2) <!--n:t1.pgr.f1.p1h2-->0.384). Under P57, F1 for Human2 and P1(H2) is <!--n:t1.p57.f1.h2-->0.380 and <!--n:t1.p57.f1.p1h2-->0.393. Accuracy is <!--n:t1.p57.accuracy.h1-->0.892 → <!--n:t1.p57.accuracy.p1-->0.901 and specificity <!--n:t1.p57.specificity.h1-->0.951 → <!--n:t1.p57.specificity.p1-->0.961 (G\* → P1, P57). With the papers' line-weighted mean, P57 MCC is <!--n:t1.p57.linemcc.h1-->0.326 for G\* and <!--n:t1.p57.linemcc.p1-->0.349 for P1 (Human2 <!--n:t1.p57.linemcc.h2-->0.329, P1(H2) <!--n:t1.p57.linemcc.p1h2-->0.350). The 21-threshold L-ge06 grid of the Human1-paper protocol is in results/m6/test/unsealed/T1/binary_threshold_grid.tsv. The biomass fix removes false positives; it adds no new true calls. All MCCs remain in the published range of about 0.25–0.45 (LITERATURE.md §0, item 2). These are estimates, and comparability with the papers is in protocol, not in label scale (§7).

### 4.10 Breakdowns per lineage, gene and task (unsealed T1 outputs)

**Per lineage** (results/m6/test/unsealed/T1/per_lineage.tsv). The table gives mean tAP_line with Δ vs G\*. Lineages are listed alphabetically, not ranked by raw AP (DESIGN §6.2).

| Lineage (expression lines) | G\* | P1 | Δ P1(H1) | Δ P1(H2) | Δ Human2 | B2 |
|---|---|---|---|---|---|---|
| Breast (<!--n:un.lin.breast.n-->65) | <!--n:un.lin.breast.gstar-->0.309 | <!--n:un.lin.breast.p1-->0.364 | <!--n:un.lin.breast.d.p1-->+0.0543 | <!--n:un.lin.breast.d.p1h2-->+0.0477 | <!--n:un.lin.breast.d.h2-->−0.0058 | <!--n:un.lin.breast.b2-->0.886 |
| CNS/Brain (<!--n:un.lin.cns.n-->126) | <!--n:un.lin.cns.gstar-->0.297 | <!--n:un.lin.cns.p1-->0.348 | <!--n:un.lin.cns.d.p1-->+0.0505 | <!--n:un.lin.cns.d.p1h2-->+0.0439 | <!--n:un.lin.cns.d.h2-->−0.0095 | <!--n:un.lin.cns.b2-->0.901 |
| Eye (<!--n:un.lin.eye.n-->15) | <!--n:un.lin.eye.gstar-->0.284 | <!--n:un.lin.eye.p1-->0.333 | <!--n:un.lin.eye.d.p1-->+0.0486 | <!--n:un.lin.eye.d.p1h2-->+0.0448 | <!--n:un.lin.eye.d.h2-->−0.0089 | <!--n:un.lin.eye.b2-->0.895 |
| Liver (<!--n:un.lin.liver.n-->25) | <!--n:un.lin.liver.gstar-->0.299 | <!--n:un.lin.liver.p1-->0.349 | <!--n:un.lin.liver.d.p1-->+0.0502 | <!--n:un.lin.liver.d.p1h2-->+0.0439 | <!--n:un.lin.liver.d.h2-->−0.0098 | <!--n:un.lin.liver.b2-->0.905 |
| Ovary/Fallopian Tube (<!--n:un.lin.ovary.n-->59) | <!--n:un.lin.ovary.gstar-->0.303 | <!--n:un.lin.ovary.p1-->0.353 | <!--n:un.lin.ovary.d.p1-->+0.0500 | <!--n:un.lin.ovary.d.p1h2-->+0.0446 | <!--n:un.lin.ovary.d.h2-->−0.0085 | <!--n:un.lin.ovary.b2-->0.904 |
| Pancreas (<!--n:un.lin.pancreas.n-->69) | <!--n:un.lin.pancreas.gstar-->0.297 | <!--n:un.lin.pancreas.p1-->0.346 | <!--n:un.lin.pancreas.d.p1-->+0.0486 | <!--n:un.lin.pancreas.d.p1h2-->+0.0437 | <!--n:un.lin.pancreas.d.h2-->−0.0070 | <!--n:un.lin.pancreas.b2-->0.878 |
| Pleura (<!--n:un.lin.pleura.n-->21) | <!--n:un.lin.pleura.gstar-->0.291 | <!--n:un.lin.pleura.p1-->0.337 | <!--n:un.lin.pleura.d.p1-->+0.0459 | <!--n:un.lin.pleura.d.p1h2-->+0.0391 | <!--n:un.lin.pleura.d.h2-->−0.0096 | <!--n:un.lin.pleura.b2-->0.894 |

The P1 gain is homogeneous across lineages (<!--n:un.lin.range-->+0.0459 to +0.0543, results/m6/test/unsealed/T1/per_lineage.tsv). Human2 is below Human1 in every lineage, as it was on training.

**Full W1 benchmark on test** (lineage-equal tAP_line, L_prob, U; results/m6/test/unsealed/T1/w1_benchmark.tsv; P2 and B6 from results/m6/test/test_look_T2_grid.tsv).

| Method | All genes | Minus CE_fit | Pooled over all pairs ("overall") |
|---|---|---|---|
| G\* (Human1) | <!--n:un.w1.h1.all-->0.297 | <!--n:un.w1.h1.noce-->0.146 | <!--n:un.w1.h1.pooled-->0.295 |
| Human2 | <!--n:un.w1.h2.all-->0.289 | <!--n:un.w1.h2.noce-->0.143 | <!--n:un.w1.h2.pooled-->0.287 |
| Human1 + P57 task calls | <!--n:un.w1.h1task.all-->0.280 | <!--n:un.w1.h1task.noce-->0.153 | <!--n:un.w1.h1task.pooled-->0.278 |
| Human1 / M2(f = 1), medium alone | <!--n:un.w1.m2.all-->0.310 | <!--n:un.w1.m2.noce-->0.161 | <!--n:un.w1.m2.pooled-->0.309 |
| **P1(H1)** | <!--n:un.w1.p1.all-->0.347 | <!--n:un.w1.p1.noce-->0.184 | <!--n:un.w1.p1.pooled-->0.346 |
| P1(H2) | <!--n:un.w1.p1h2.all-->0.341 | <!--n:un.w1.p1h2.noce-->0.180 | <!--n:un.w1.p1h2.pooled-->0.340 |
| B0 null | <!--n:un.w1.b0.all-->0.124 | <!--n:un.w1.b0.noce-->0.071 | <!--n:un.w1.b0.pooled-->0.121 |
| B3 expression | <!--n:un.w1.b3.all-->0.317 | <!--n:un.w1.b3.noce-->0.198 | <!--n:un.w1.b3.pooled-->0.314 |
| B8+ | <!--n:un.w1.b8p.all-->0.418 | <!--n:un.w1.b8p.noce-->0.212 | <!--n:un.w1.b8p.pooled-->0.416 |
| B8− | <!--n:un.w1.b8m.all-->0.223 | <!--n:un.w1.b8m.noce-->0.109 | <!--n:un.w1.b8m.pooled-->0.222 |
| B8′ | <!--n:un.w1.b8q.all-->0.305 | <!--n:un.w1.b8q.noce-->0.151 | <!--n:un.w1.b8q.pooled-->0.303 |
| B2 gene prior | <!--n:un.w1.b2.all-->0.895 | <!--n:un.w1.b2.noce-->0.681 | <!--n:un.w1.b2.pooled-->0.890 |
| B2m | <!--n:un.w1.b2m.all-->0.896 | <!--n:un.w1.b2m.noce-->0.687 | <!--n:un.w1.b2m.pooled-->0.893 |
| B6 (refit on training ∪ validation) | <!--n:t2.tap.b6-->0.885 | <!--n:t2.tap.b6.noce-->0.673 | <!--n:t2.pooled.b6-->0.879 |
| **P2** (W7-05, non-GEM) | <!--n:t2.tap.p2-->0.899 | <!--n:t2.tap.p2.noce-->0.693 | <!--n:t2.pooled.p2-->0.896 |

- **Human2 vs Human1.** Human2 − Human1 on test is <!--n:un.w1.h2h1-->−0.0084 [−0.0089, −0.0079] (results/m6/test/unsealed/T1/w1_benchmark.tsv), so H1.2 is not supported on test either.
- **All 394 lines.** Including the lines without expression, generic Human1 scores <!--n:un.all.h1-->0.297 and B2 <!--n:un.all.b2-->0.894 on <!--n:un.all.n-->394 lines (results/m6/test/unsealed/T1/generic_all_lines.tsv).

**Per gene** (results/m6/test/unsealed/T1/per_gene.tsv). The table covers <!--n:un.gene.n-->963 estimable genes, each with ≥ 5 positive and ≥ 5 negative test lines.

| Method | Mean ΔAP_gene over the null | Genes changed or above the null |
|---|---|---|
| P1 | <!--n:un.gene.p1.mean-->−0.0002 | <!--n:un.gene.p1.nz-->71 genes changed from line-constant; <!--n:un.gene.p1.pos-->25 above the null |
| B3 expression | <!--n:un.gene.b3.mean-->+0.0112 | <!--n:un.gene.b3.pos-->489 genes above the null |
| B8+ | <!--n:un.gene.b8p.mean-->+0.0122 | — |
| B2m | <!--n:un.gene.b2m.mean-->+0.0331 | — |

P1's across-line signal comes only from medium differences. P2 reaches <!--n:t2.p2.mdap-->+0.0446 (results/m6/test/test_look_T2_grid.tsv), mostly through its medium-class prior and omics features.

**Strata of C1** (results/m6/test/unsealed/T1/strata.tsv; Δ̄_lin P1 − G\*).

| Stratum (lines) | Δ̄_lin [95% CI] |
|---|---|
| Expressed pairs (<!--n:un.str.expr.n-->380) | <!--n:un.str.expr-->+0.0466 [0.0450, 0.0480] |
| Addressable in Human1 (<!--n:un.str.addr1.n-->380) | <!--n:un.str.addr1-->+0.0902 [0.0876, 0.0925] |
| Library: Avana (<!--n:un.str.avana.n-->273) | <!--n:un.str.avana-->+0.0517 [0.0500, 0.0532] |
| Library: Humagne (<!--n:un.str.humagne.n-->85) | <!--n:un.str.humagne-->+0.0412 [0.0386, 0.0441] |
| Library: KY (<!--n:un.str.ky.n-->22) | <!--n:un.str.ky-->+0.0487 [0.0460, 0.0515] |
| Organoids (<!--n:un.str.organoid.n-->40) | <!--n:un.str.organoid-->+0.0364 [0.0324, 0.0413] |
| CN ≥ 4 pairs only (<!--n:un.str.cn4.n-->15) | <!--n:un.str.cn4-->+0.0054 [−0.0129, 0.0235] |
| Multi-copy loci only (<!--n:un.str.multi.n-->283) | <!--n:un.str.multi-->+0.0000 [0.0000, 0.0000] |
| Housekeeping-adjacent only (<!--n:un.str.hk.n-->380) | <!--n:un.str.hk-->+0.0167 [0.0161, 0.0173] |
| Moonlighting genes only (<!--n:un.str.moon.n-->330) | <!--n:un.str.moon-->−0.0119 [−0.0173, −0.0074] |

- P1 helps in every library and model type.
- It does not change multi-copy loci.
- It **loses** on moonlighting genes, whose essentiality is non-metabolic.

**Per task** (results/m6/test/unsealed/T1/per_task.tsv; P1 + one task's failure calls vs P1 alone).
- <!--n:un.task.adding-->31 of <!--n:un.task.n-->57 tasks add calls.
- <!--n:un.task.ci_pos-->0 tasks improve tAP, and <!--n:un.task.ci_neg-->30 lower it with CI < 0. The worst costs <!--n:un.task.worst-->−0.0247.

This replicates the training result that task calls are less precise than the biomass-based ranking.

### 4.11 Flux states carry no information beyond omics and medium (W6)

The verdict is "not supported" for both hypotheses, both arms and all three analyses (results/w6/provisional/analysis/decisions.tsv). The run is labelled provisional because it predates the P1 freeze, but its Human1 + B-loo arm is exactly P1 (results/w6/README.md).

| Analysis (units) | RNA + flux − RNA, median Δ, G\* / P1 | Holm p, G\* / P1 | Flux − mapped RNA, median Δ, G\* / P1 |
|---|---|---|---|
| Subtypes (<!--n:w6.n.subtypes-->150) | <!--n:w6.gstar.flux.subtypes.med-->0.000 / <!--n:w6.p1.flux.subtypes.med-->0.000 | <!--n:w6.gstar.flux.subtypes.p-->1.00 / <!--n:w6.p1.flux.subtypes.p-->1.00 | <!--n:w6.gstar.sim.subtypes.med-->−0.044 / <!--n:w6.p1.sim.subtypes.med-->−0.044 |
| Drugs (<!--n:w6.n.drugs-->106) | <!--n:w6.gstar.flux.drugs.med-->0.000 / <!--n:w6.p1.flux.drugs.med-->0.000 | <!--n:w6.gstar.flux.drugs.p-->0.38 / <!--n:w6.p1.flux.drugs.p-->1.00 | <!--n:w6.gstar.sim.drugs.med-->−0.239 / <!--n:w6.p1.sim.drugs.med-->−0.240 |
| Metabolomics (<!--n:w6.n.metabolomics-->225) | <!--n:w6.gstar.flux.metabolomics.med-->−0.007 / <!--n:w6.p1.flux.metabolomics.med-->−0.008 | <!--n:w6.gstar.flux.metabolomics.p-->1.00 / <!--n:w6.p1.flux.metabolomics.p-->1.00 | <!--n:w6.gstar.sim.metabolomics.med-->−0.399 / <!--n:w6.p1.sim.metabolomics.med-->−0.398 |

- **Positive control.** The synthetic positive control is recovered: median Δρ <!--n:w6.gstar.synth-->+0.42 for G\* and <!--n:w6.p1.synth-->+0.45 for P1 (results/w6/provisional/analysis/secondary.tsv).
- **Cause of the null.** Generic-parent flux depends only on the assay medium: one state for all of PRISM, two in GDSC, at most 20 for culture media. <!--n:w6.drugs.zero-->77/106 drug units have Δ exactly 0 (results/w6/provisional/analysis/decisions.tsv).
- **Interpretation.** The null is a property of line-invariant models, not of the test. The C3 failure on test points the same way.

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
- **Test (T2)** is in §4.8.
- **D12 wording.** P2 contains no GEM output: it is a tuned non-GEM gene-prior model, and W7's GEM result is negative. The W7 runs used `allow_dirty` against §11.4. The reviewer's bitwise reproduction mitigates this.

## 5. Resource (W8)

| | Value |
|---|---|
| Models (SBML L3V1 FBC v2; Human2 built first) | <!--n:res.models-->48 = 2 parents × <!--n:res.perparent-->24 medium signatures |
| Lines covered | <!--n:res.lines-->2,114, of which <!--n:res.uns.lines-->753 are unscreened (<!--n:res.rows-->4,228 line × parent rows) |
| Models passing all six gates | <!--n:res.pass-->47; the 8 models added for unscreened lines pass <!--n:res.uns.new-->8/8 |
| Models flagged W3-MEV (serum media) | <!--n:res.mev-->42 |
| Size | <!--n:res.mb-->93.3 MB |
| MEMOTE total score (<!--n:res.memote.n-->10 models) | Human2 <!--n:res.memote.h2-->0.638, Human1 <!--n:res.memote.h1-->0.630; identical across media |
| FROG reports | <!--n:res.frog.n-->2 models (Human1 and Human2, classic RPMI; both all optimal, see results/w8/README.md). Human2: <!--n:res.frog.fva-->12,748 FVA rows, <!--n:res.frog.gd-->2,870 gene and <!--n:res.frog.rd-->12,748 reaction deletions, all optimal; <!--n:res.frog.cpu-->0.66 CPU-h |

Sources: resource/models.tsv, results/w8/unscreened/build_summary.json, results/w8/memote/memote_summary.json, results/w8/memote/memote_summary_release6.json, results/w8/frog/frog_summary.json.

- **Contents of each model.** The SBML notes record the pipeline, Human-GEM release, medium signature and tier, W3 variant, git commit, input SHA-256s and flags. Closed reactions keep bounds (0, 0), and no reaction was repaired. `scripts/w8_line_sbml.py` writes per-line files on demand (resource/README.md).
- **The one gate failure is LP noise.** One Human2 model, for a single-line RPMI + testosterone medium, misses the random-KO gate. Its reloaded KO growth is <!--n:res.noise.default-->−8.0e−09 under GLPK default tolerances, and <!--n:res.noise.exact-->0.0 with exact rational arithmetic (results/w8/ko_noise/human2_20cc933a4b852b6a_summary.json). The gate rule was left unchanged and the miss is documented. The original resource commit message claimed that all gates passed; RESEARCH_TREE.md corrects it.
- **The unscreened extension.** It adds DESIGN's unscreened scope, using W2's frozen medium rule; non-documented lines get the flagged classic RPMI + 10% FBS. Its 8 new files were built from a dirty tree (results/w8/README.md).
- **Deviations.** No extraction-based context-specific models are released, because every extraction family failed TG and the context lever was not selected (§4.2, §4.6); the released models are P1. The resource ships per-signature files rather than per-line files, so it does not suggest more line specificity than P1 has. The Human1 FROG report has not been run. MEMOTE failures are properties of the parent models.

## 6. Failures and negative results

**Failure account on test** (DESIGN §10.3; results/m6/test/unsealed/T1/failure_account.tsv). P1 makes <!--n:un.fail.all_errors.fn-->100,383 false-negative and <!--n:un.fail.all_errors.fp-->10,934 false-positive gene × line calls, using fixed-cutoff calls (ratio_r < 0.01) against L_prob. A call can fall in several classes, so the shares overlap.

| Error class | False negatives (share) | False positives (share) |
|---|---|---|
| Complex subunit | <!--n:un.fail.complex_subunit.fn-->48,946 (<!--n:un.fail.complex_subunit.ffn-->48.8%) | <!--n:un.fail.complex_subunit.fp-->5,038 (<!--n:un.fail.complex_subunit.ffp-->46.1%) |
| Common essential | <!--n:un.fail.common_essential.fn-->39,252 (<!--n:un.fail.common_essential.ffn-->39.1%) | <!--n:un.fail.common_essential.fp-->526 (<!--n:un.fail.common_essential.ffp-->4.8%) |
| Paralog-buffered (same compartment) | <!--n:un.fail.paralog_buffered.fn-->37,154 (<!--n:un.fail.paralog_buffered.ffn-->37.0%) | <!--n:un.fail.paralog_buffered.fp-->0 (<!--n:un.fail.paralog_buffered.ffp-->0.0%) |
| Housekeeping-adjacent | <!--n:un.fail.housekeeping_adjacent.fn-->11,643 (<!--n:un.fail.housekeeping_adjacent.ffn-->11.6%) | <!--n:un.fail.housekeeping_adjacent.fp-->4,932 (<!--n:un.fail.housekeeping_adjacent.ffp-->45.1%) |
| Pre-registered FBA failure | <!--n:un.fail.expected_fba_failure.fn-->4,790 (<!--n:un.fail.expected_fba_failure.ffn-->4.8%) | <!--n:un.fail.expected_fba_failure.fp-->0 (<!--n:un.fail.expected_fba_failure.ffp-->0.0%) |
| Medium-sensitive | <!--n:un.fail.medium_sensitive.fn-->3,951 (<!--n:un.fail.medium_sensitive.ffn-->3.9%) | <!--n:un.fail.medium_sensitive.fp-->3,109 (<!--n:un.fail.medium_sensitive.ffp-->28.4%) |
| Missing from model | <!--n:un.fail.missing_from_model.fn-->3,099 (<!--n:un.fail.missing_from_model.ffn-->3.1%) | <!--n:un.fail.missing_from_model.fp-->0 (<!--n:un.fail.missing_from_model.ffp-->0.0%) |
| Unexpressed | <!--n:un.fail.unexpressed.fn-->2,670 (<!--n:un.fail.unexpressed.ffn-->2.7%) | <!--n:un.fail.unexpressed.fp-->80 (<!--n:un.fail.unexpressed.ffp-->0.7%) |
| Biomass-artefact lethal | <!--n:un.fail.biomass_artefact_lethal.fn-->490 (<!--n:un.fail.biomass_artefact_lethal.ffn-->0.5%) | <!--n:un.fail.biomass_artefact_lethal.fp-->1,500 (<!--n:un.fail.biomass_artefact_lethal.ffp-->13.7%) |
| CN-amplified | <!--n:un.fail.cn_amplified.fn-->197 (<!--n:un.fail.cn_amplified.ffn-->0.2%) | <!--n:un.fail.cn_amplified.fp-->4 (<!--n:un.fail.cn_amplified.ffp-->0.0%) |

- **Which errors dominate** (results/m6/test/unsealed/T1/failure_account.tsv). Misses are mostly complex subunits, common essentials and paralog-buffered genes. False positives are mostly complex, housekeeping-adjacent and medium-sensitive calls. Remaining biomass artefacts still account for <!--n:un.fail.biomass_artefact_lethal.ffp-->13.7% of false positives.
- **Classes not annotated.** "Gene-less bypass" and "model-contradicted loss" could not be annotated, because P1 has no context or loss step.
- **Upper bound** (results/m6/test/unsealed/T1/failure_account.tsv). Only <!--n:un.fail.upper-->54.7% of test positives lie in addressable genes, which bounds any GEM's recall.
- **Worst cells** (results/m6/test/unsealed/T1/worst_cells.tsv, results/baseline/human1_gene_mapping.tsv, results/baseline/human2_gene_mapping.tsv). All <!--n:un.worst.n-->20 worst lineage × gene cells are in <!--n:un.worst.lineages-->CNS/Brain, which ranks first only because it is the largest lineage. They are missed essentials, with <!--n:un.worst.fp-->0 false positives among them, in every line: <!--n:un.worst.genes-->ALDOA, ATP2A2, ATP6V0C, COPB1, DAD1, DHPS, DPAGT1, DUT, EEF2, GAPDH, HSPA9, POLD1, POLE, POLE2, RFC5, SOD1, UBA1, UBE2I, OGT, RNGTT. These are glycolytic, replication, translation, glycosylation and ubiquitin-machinery genes. Many are essential for non-metabolic roles or sit in lumped multi-subunit reactions. The pre-registered pFBA traces of these cells are not in the unsealed bundle.

**Negative results** (training and validation; §4).
- Human2 is not better than Human1.
- Context-specific extraction failed in all five families and in CORDA, and method dominates tissue.
- GPR curation failed: the literal paralog rule is wrong for lumped GPRs.
- Medium failed the shuffled-medium control.
- Tasks and objectives did not help.
- Flux information was null.
- GEM features added nothing to gene-prior ML.
- On test, C3 failed, and P1 does not beat B8+.

**Positive controls that failed because of model structure.**
- Serum cholesterol cannot rescue the post-squalene sterol genes.
- Uridine cannot rescue CAD.
- The ETF capability control fails in Human1.
- P1 loses mevalonate-core essentiality under serum (W3-MEV).

M3 calibration was infeasible as specified, and the W2 single-removal flip analysis did not complete (results/w2/README.md).

## 7. Limitations

- **Addressable-set ceiling.** For Human1, the parent of P1 and G\*, <!--n:w1.addr.h1.out-->45% of positive training pairs lie in genes whose KO changes nothing (Human2: <!--n:w1.addr.h2.out-->42%; results/w1/train_addressable.tsv). On test the Human1 share is <!--n:un.fail.outside-->45.3% (results/m6/test/unsealed/T1/failure_account.tsv). Mechanisms outside FBA are invisible to every GEM arm, as the worst cells show.
- **P1 is line-specific only through the medium.** C3 failed and W2-SHUF holds, so even that shows no detectable line-specific information: the confirmed C1 gain is a better generic ranking, not personalisation. This is why the resource holds <!--n:res.models-->48 models, not one per line (resource/models.tsv).
- **The gene prior dwarfs the GEMs.** On test, P1 closes only <!--n:t1.headroom-->8.3% [8.1%, 8.6%] of the headroom between G\* and B2 (results/m6/test/test_look_T1.tsv). Most label signal is gene identity shared across lineages.
- **Near-binary scores.** Results depend on large tied blocks, and the B8+ expression tie-break beats P1 on test.
- **Chronos vs CERES.** The Human1 paper applied gene effect < −0.6 to CERES (19Q3); we apply it to Chronos (26Q3), whose scale differs. The Human2 paper states no threshold. Our MCCs match the papers' protocol, not their label scale, and neither paper's per-line extraction was reproduced (D2).
- **Flags.** W3-D1 (control reading fixed after training results), W3-MEV (a positive-control failure built into P1), W2-SHUF and W2-MOPEN all apply to the confirmed C1 claim.
- **Post hoc elements.** The random-removal control, the constant medium and the artefact-only biomass variant were added after results were seen.
- **Labels and bounds.** L_prob and DepMap's common-essential lists are fitted on all lines, and Human2's curation may have been DepMap-informed. Uptake bounds encode availability, not rates.
- **Power and scope.** The C1 CI and the transfer claim are conditional on the 7 test lineages. The method-vs-tissue panel had only 30 lines.
- **Retired seed.** A retired seed's draw was published in round 1. Only aggregate statistics on its lineages were seen (notes/label_exposure.md).

## 8. Runtime and compute

**What was summed.**
- **Included:** `cpu_s` and `wall_s` of every row of every committed results/*/**/MANIFEST.tsv and PROVENANCE.tsv at HEAD, grouped by top directory. Each row is one script run, including reruns, pilots and looks; the T1/T2 harness runs fall under M6.
- **Excluded:** runs killed before writing a row, design-phase work, review probes, the baseline KO sweep, the GPRTransform assessment and the MEMOTE/FROG runs in `.venv_memote`, which have no manifest rows. Agent time is also excluded.
- **Wall time** is summed run duration, not elapsed time.

| Workstream | Recorded runs | CPU-h | Summed run wall-h |
|---|---|---|---|
| W1 generic benchmark | <!--n:rt.w1.runs-->12 | <!--n:rt.w1.cpu-->7.0 | <!--n:rt.w1.wall-->1.4 |
| W2 medium | <!--n:rt.w2.runs-->22 | <!--n:rt.w2.cpu-->19.8 | <!--n:rt.w2.wall-->5.3 |
| W3 tasks and biomass | <!--n:rt.w3.runs-->25 | <!--n:rt.w3.cpu-->20.9 | <!--n:rt.w3.wall-->6.9 |
| W4 context models (Stage A, smoke tests) | <!--n:rt.w4.runs-->6 | <!--n:rt.w4.cpu-->130.9 | <!--n:rt.w4.wall-->8.0 |
| W5 GPR curation | <!--n:rt.w5.runs-->17 | <!--n:rt.w5.cpu-->0.7 | <!--n:rt.w5.wall-->0.5 |
| M5 panel pre-check | <!--n:rt.m5.runs-->2 | <!--n:rt.m5.cpu-->3.0 | <!--n:rt.m5.wall-->0.4 |
| M6 factorial, T1/T2 looks, rehearsal | <!--n:rt.m6.runs-->23 | <!--n:rt.m6.cpu-->2.7 | <!--n:rt.m6.wall-->0.8 |
| W6 flux information | <!--n:rt.w6.runs-->9 | <!--n:rt.w6.cpu-->2.6 | <!--n:rt.w6.wall-->1.3 |
| W7 auto-research | <!--n:rt.w7.runs-->23 | <!--n:rt.w7.cpu-->8.4 | <!--n:rt.w7.wall-->2.2 |
| W8 resource | <!--n:rt.w8.runs-->13 | <!--n:rt.w8.cpu-->5.2 | <!--n:rt.w8.wall-->1.6 |
| **Total** | <!--n:rt.total.runs-->174 | <!--n:rt.total.cpu-->210.4 | <!--n:rt.total.wall-->30.0 |

Source: results/*/**/MANIFEST.tsv and PROVENANCE.tsv at HEAD, summed by `runtime_table` in `scripts/report_build.py`.

- **Elapsed time.** <!--n:rt.elapsed-->46.6 h from the first commit to the last commit that changed results/ (git log), against a target of 1–1.5 days (D2). The process review blamed review and fix rounds, not compute (notes/milestone_reviews/M3M4_process_review.md).
- **Stage A dominates.** Its builds took <!--n:w4.buildcpu-->108.7 CPU-h (results/w4/stage_a/MANIFEST.tsv). The context families it built were all rejected.
- **W7's ledger** records <!--n:rt.budget.w7-->8.8 CPU-h, including aborted runs (results/budget.tsv).
- **Not in the manifests.** The FROG report took <!--n:res.frog.cpu-->0.66 CPU-h (results/w8/frog/frog_summary.json). Unrecorded losses quoted from READMEs, not checked: about 8.5 CPU-h of killed W3 runs, 1.3 of a stalled W2 run, and 7 of a BLAS-thread incident in W6 (results/w3/README.md, results/w2/README.md, RESEARCH_TREE.md item 8.1).
- **Runtime alongside accuracy** (test tAP in §4.10).
  - G\*: one generic KO sweep took about 53 s wall on 15 workers (results/baseline/README.md; no manifest row).
  - P1: the M6 sweeps of the biomass variant on 19 media × 2 parents took <!--n:rt.p1.sweeps-->0.63 CPU-h (results/m6/MANIFEST.tsv).
  - B6 / P2: one full nested-LOLO training screen of W7-05 took <!--n:rt.w7.05.screen-->0.36 CPU-h (results/budget.tsv).
  - Test looks: T1, which scores every method on the test lines, took <!--n:rt.t1.cpu-->2,638 CPU-s. T2, including the refit of B6 and P2 on training ∪ validation, took <!--n:rt.t2.cpu-->154 CPU-s (results/m6/test/MANIFEST.tsv).
  - B2 is a column mean computed in seconds, and it is the most accurate per unit of compute by a wide margin.

## 9. Reproducibility

- **Environment.** `.venv` (Python 3.11, requirements.txt pinned; cobra 0.32.1), Gurobi 13 through a slot pool or GLPK, declared per family. MEMOTE and FROG run in a separate `.venv_memote`.
- **Pins and ledgers.** FROZEN.tsv pins frozen code, inputs and values by SHA-256 and is checked at commit and before every look. The ledgers are hash-chained (results/validation_ledger.tsv, results/test_ledger.tsv). The unsealed T1 tables were re-derived and matched their sealed SHA-256 (D16).
- **Provenance.** Every production script writes a MANIFEST row: commit, dirty flag, input and output SHA-256, versions, seeds, solver, workers, time and RSS. Long jobs ran from pinned clean worktrees. W3, W7 and the W8 unscreened extension ran from a dirty tree; reviewers reproduced W3 and W7 bitwise.
- **Seeds.** SEEDS = [11, 23, 37, 41, 53, 67, 79, 83, 97, 101] plus the split seed (results/splits/SEED_SOURCE.md).
- **Tests.** `.venv/bin/python -m pytest -q tests/` runs known-answer metric tests, leakage-guard exploit tests and toy tests of each look harness.
- **This report.** `scripts/report_build.py --check --require-filled` regenerates every marked number from HEAD and fails on any disagreement or open placeholder. `--update` rewrites the marked numbers. `tests/test_report_build.py` tests the marker logic.
- **Data.** data/ is not committed. DATA.md lists the source, version, licence and checksum of every input.

## 10. Path taken

**Decisions** (RESEARCH_TREE.md).
- **D0–D2.** Human1 := v1.19.0, the last 1.x release; Human2 := v2.1.1. ftINIT dropped. Paper-protocol MCC/F1 added. Compute target about 1–1.5 days.
- **D3–D5.** Design frozen after three review rounds; round 1 alone produced 34 PI decisions, among them tAP, B8 and C1b. Split seed 38 by commit-reveal. The a priori lists were re-derived blind, because the PI had relaxed an authorship rule post hoc (A2).
- **D6–D8.** G\* = Human1. Lighter guard bookkeeping (A6). CORDA kept as a 10-line comparison.
- **D9–D11.** Scope cuts after M4 (A12). The W2 shortlist went forward with its flags. Factorial cells and their order were fixed before scoring. P1 selected.
- **D12–D14.** A non-GEM P2 is admissible, with fixed wording. Runner-up comparison. P2 = W7-05.
- **D15 (T1).** C1 headline; C1b confirmed; C3 not confirmed, so C4/C5 became estimates; transfer holds.
- **D16 (T2, unseal).** C4/C5 estimated; all 11 sealed T1 tables re-derived with matching SHA-256 and published.

**Process lessons** (LESSONS.md; notes/milestone_reviews/).
- **Independent recomputation is the core guarantee.** One review found a float-tie AP instability of up to 0.018 that 36 passing checks had missed. The W4 code review caught three defects, any of which would have invalidated Stage A.
- **Fix control readings and eligible sets before looking.** W3-D1 exists because this was not done, and W2's eligible set could not be verified as fixed in advance (notes/reviews/M4_W2_review.md).
- **Read every result against B2 and B8,** not only G\*. On test the GEM gain was confirmed, yet it is still below B8+.
- **Look for structural false positives before adding data.** The biomass edit alone was not scored on test. Descriptively, P1 minus medium alone is <!--n:un.w1.p1.all-->0.347 − <!--n:un.w1.m2.all-->0.310 in test tAP (results/m6/test/unsealed/T1/w1_benchmark.tsv), consistent with the training and validation effect, whereas every data-driven extraction lost.
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
- The first resource commit message overstated its gate results (§5).
- The first rehearsal attempt aborted on an inherited PYTHONPATH and was rolled back without using a look (RESEARCH_TREE.md).

## 11. Follow-ups and remaining ideas

**Pipeline improvements.**
- **The artefact-only biomass variant** (W7-03; now evaluated as F2, see §12) removes only the retinoid, vitamin-D, tocopherol and carnitine groups. It keeps the mevalonate core lethal under serum, so it fixes W3-MEV, at the cost of part of P1's gain. It is the first candidate for a P1 successor and needs its own validation look.
- **Line-specific constraints beyond medium** (genomic loss and E-Flux now tested as F4, see §12; the rest remain open). C3 shows that P1 carries no detectable line-specific signal. Copy number and LOH (here only a CN proxy), proteomics-based capacity bounds, uptake rates and enzyme constraints are untested as P1-level levers.
- **Group models** (Stage C, skipped under A12): lineage, genomic-background (e.g. MTAP_del, ASS1_silenced) and lineage × background models, with their pre-registered within-lineage contrasts.
- **Curated GPRs revisited** (the ≥ 30%-identity rule on P1 is now tested as F1, see §12; the UniProt path remains open) with the ≥ 30%-identity paralog rule, the UniProt removal path (file still requested from the user) and the edit-review rejections applied, re-scored as a fresh pre-registered arm.

**Analyses still owed.**
- **The pFBA traces of the 20 worst test cells** (DESIGN §10.3): now done as F3, see §12.
- **W3 context claims** checked for the same sign on validation, as DESIGN requires before any claim.

**Stated deviations.** These pre-registered outputs were not produced:
- the §6.4 "top-k" unit, i.e. the change in true positives among each line's top-k genes;
- the optimistic 16-lineage LOLO estimate of the frozen P1 protocol (§6.4).


## 12. Post-hoc follow-ups (DESIGN §17 A15)

These were run after the test unseal at the user's request. They are labelled post hoc and exploratory, were evaluated on training and validation only (no new test evaluation), and do not change the T1/T2 results above. Agents for F1, F2 and F4 were kept blind to the test outputs. Each rule was fixed in RESEARCH_TREE.md (D17–D19) before its first label-reading run.

- **F1: stricter paralog curation on top of P1** (results/f1/README.md). AND→OR merges were limited to Ensembl paralogs with ≥ 30% identity in the same compartment, excluding lumped, large-complex and transport reactions, plus the TYMS fix. That gives 20 edits; all 12 RNA-polymerase merges from W5 are excluded. Training Δ̄ vs P1 = −0.00063 [−0.00071, −0.00055], 0/12 folds, so it **fails** the training gate. The loss comes from the RRM2/RRM2B merge. Without common essentials: +0.0021, below SESOI. An independent recompute agreed (results/f2/train/contrasts.tsv).
- **F2: artefact-only biomass on the documented medium** (results/f2/README.md). Validation look 7/40 (results/m10/val/val_look_M10_F2_decision.json): +0.0310 [0.0293, 0.0326] vs G*, against P1's +0.0517 in the same look, so the validated cost of this biologically safer variant is ≈ 0.021 (training −0.023). It keeps HMGCR/MVK/PMVK/MVD lethal on every medium. The training decomposition by gene group (results/f2/decomposition.tsv) gives:
  - the mevalonate artefact costs P1 −0.0014, i.e. it does not help;
  - un-calling the 35 cofactor-pool genes gives about 62% of P1's gain;
  - F2 keeps 20 of them lethal (heme, CoQ, BH4, PLP, TPK1, RFK …), which explains nearly all of its gap.
- **F3: pFBA traces of the 20 worst T1 cells** (results/f3/README.md; descriptive). All 20 are CNS/Brain false negatives with KO ratio 1.0 on every medium. They break down as:
  - 13 functions Human1 does not represent;
  - 4 DNA-polymerase subunits hidden by one lumped DNA-synthesis reaction with a 17-way OR rule;
  - 2 alternative pathways (GAPDH, ALDOA);
  - 1 gene absent from Human1.
  None is fixable through the medium. The 20 are a tie-broken subset of many 126/126-wrong cells.
- **F4: line-specific constraints on top of P1** (results/f4/README.md).
  - **Genomic loss:** genes lost in 665/686 training lines; only 22 KO calls change. Δ̄ vs P1 = −0.00013, 3/12 folds, so it **fails**.
  - **E-Flux expression bounds:** the 30-line panel gives −0.0007, so it **fails** the pre-check.
  Beyond the medium, no tested line-specific constraint improves P1.

## References

Numbered as in LITERATURE.md; only works cited here are listed.

1. Robinson JL et al. An atlas of human metabolism. *Sci Signal* 2020;13:eaaz1482. doi:10.1126/scisignal.aaz1482.
2. Luo J et al. Reconstruction of human metabolic models with large language models. *PNAS* 2026;123:e2516511123. doi:10.1073/pnas.2516511123.
7. Opdam S et al. A systematic evaluation of methods for tailoring genome-scale metabolic models. *Cell Syst* 2017;4:318–329. doi:10.1016/j.cels.2017.01.010.
8. Richelle A et al. Increasing consensus of context-specific metabolic models by integrating data-inferred cell functions. *PLoS Comput Biol* 2019;15:e1006867. doi:10.1371/journal.pcbi.1006867.
10. Lagziel S, Lee WD, Shlomi T. Inferring cancer dependencies on metabolic genes from large-scale genetic screens. *BMC Biol* 2019;17:37. doi:10.1186/s12915-019-0654-4.
11. Rossiter NJ et al. CRISPR screens in physiologic medium reveal conditionally essential genes in human cells. *Cell Metab* 2021;33:1248–1263. doi:10.1016/j.cmet.2021.02.005.
22. Dempster JM et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biol* 2021. doi:10.1186/s13059-021-02540-7.
34. Becker SA, Palsson BO. Context-specific metabolic networks are consistent with experiments (GIMME). *PLoS Comput Biol* 2008;4:e1000082. doi:10.1371/journal.pcbi.1000082.
40. Vlassis N, Pacheco MP, Sauter T. Fast reconstruction of compact context-specific metabolic network models (FASTCORE). *PLoS Comput Biol* 2014;10:e1003424. doi:10.1371/journal.pcbi.1003424.
43. Colijn C et al. Interpreting expression data with metabolic flux models (E-Flux). *PLoS Comput Biol* 2009;5:e1000489. doi:10.1371/journal.pcbi.1000489.
50. Richelle A et al. Model-based assessment of mammalian cell metabolic functionalities using omics data (CellFie). *Cell Rep Methods* 2021;1:100040. doi:10.1016/j.crmeth.2021.100040.
55. Cantor JR et al. Physiologic medium rewires cellular metabolism and reveals uric acid as an endogenous inhibitor of UMP synthase (HPLM). *Cell* 2017. doi:10.1016/j.cell.2017.03.023.
56. Vande Voorde J et al. Improving the metabolic fidelity of cancer models with a physiological cell culture medium (Plasmax). *Sci Adv* 2019;5:eaau7314. doi:10.1126/sciadv.aau7314.
64. Segrè D, Vitkup D, Church GM. Analysis of optimality in natural and perturbed metabolic networks (MOMA). *PNAS* 2002;99:15112. doi:10.1073/pnas.232349399.
92. McSherry F, Najork M. Computing information retrieval performance measures efficiently in the presence of tied scores. *ECIR 2008*, LNCS 4956. (Title and venue as given in LITERATURE.md, which flags them for verification against the Springer record; not yet verified.)

Not in LITERATURE.md, cited as given in DESIGN §17 A8: Machado et al. 2016. *PLoS Comput Biol*. doi:10.1371/journal.pcbi.1005140.