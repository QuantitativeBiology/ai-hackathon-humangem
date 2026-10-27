# Context-specific genome-scale metabolic models and cancer gene essentiality: a pre-registered benchmark on DepMap 26Q3

> **Status: skeleton (2026-10-27).** This is a draft written by the AI pipeline; see the disclaimer in README.md. Each section names the files or scripts that will fill it, and `scripts/report_build.py` assembles the final numbers. The only results quoted here are the training and validation aggregates already public in README.md. Test-tier results stay sealed until the design's look T1.

## Abstract
TODO (last). Filled from §4 by `scripts/report_build.py`, using the §6.4 wording rules (claim tier, SESOI 0.01, transfer wording) and the flags "gate not passed", "method-dominated", "label-calibrated" and "exploratory".

## 1. Question and pre-registration
- Can context-specific, medium-aware, curated versions of Human1 (Human-GEM v1.19.0) and Human2 (v2.1.1) predict metabolic gene essentiality in DepMap 26Q3 CRISPR screens better than the generic model, on held-out lineages?
- The pre-registered design is DESIGN.md v3, frozen at 61b0beb after three review rounds. Amendments A1–A10 are listed in §17, with user decisions D1–D8 in RESEARCH_TREE.md.
- Primary endpoint: the lineage-equal-weighted per-line ΔtAP (tie-aware AUPRC) of the selected pipeline P1 vs G\* on 7 locked test lineages. The confirmatory chain is C1→C1b→C3→C4→C5.
- TODO: summarise the hypotheses H1.x–H4.x from DESIGN §2–3.

## 2. Data and splits
- DepMap Public 26Q3: CRISPR gene effect and dependency, expression, copy number and models. Also Human-GEM, ProCan and Nusinow proteomics, CCLE metabolomics, PRISM, GDSC, CORUM, Complex Portal and Ensembl paralogs (DATA.md).
- Split, seed 38 (commit-reveal): 7 test lineages with 394 lines, 4 validation lineages with 200 lines, and 767 training lines (results/splits/).
- Labels: L_prob (primary), L-ge06 (paper protocol). Label building and sealing are described in results/labels/build_report.json.
- TODO: the cohort table from `results/splits/split_summary.tsv`.

## 3. Methods
- **3.1 Generic baseline (W1).** Single-gene KO on Ham's medium; ratio_r; tAP; baselines B0–B6 and B8±; §6.8 binary endpoints under the Human1/Human2 paper protocols (P57 task essentiality, PGR). Source: results/w1/README.md.
- **3.2 Medium (W2).** Frozen media library (two-agent extraction), micronutrient floor, serum model, M0–M3 ladder. Source: results/w2/README.md, A9.
- **3.3 Objectives, tasks, biomass (W3).** TODO: results/w3/.
- **3.4 Context-specific models (W4).** F0 thresholding, F2 GIMME, F3 FASTCORE, F4 E-Flux, F5 label-calibrated presence, and F3b CORDA (10-line comparison, A7). Repair rule, genomic constraint and Stage A on a 30-line panel. Source: results/w4/README.md, A4.
- **3.5 GPR curation (W5).** TODO: results/w5/.
- **3.6 Flux-state information (W6).** TODO.
- **3.7 Auto-research loop (W7).** TODO.
- **3.8 Resource (W8).** SBML export and validation gates. Source: src/nextgem/resource.py.
- **3.9 Leakage controls and statistics.** labels.py-only access, a PreToolUse hook, a tripwire, the FROZEN.tsv manifest, hash-chained ledgers, scalar-only validation looks, and the exposure register. Lineage-stratified bootstrap and the sign-flip permutation.
- **Method assessed and not adopted:** GPRTransform (Machado et al. 2016); see notes/methods/GPRTransform_assessment.md and A8.

## 4. Results
- **4.1 Generic models.** Tie-aware AUPRC per line, lineage-weighted:
  - training: Human1 0.300, Human2 0.292;
  - validation: Human1 0.302, Human2 0.295;
  - Δ(H2−H1) on validation: −0.0071 [−0.0079, −0.0063];
  - so **G\* = Human1**.

  MCC (P57 × gene effect < −0.6, training): 0.329 / 0.331. Gene-prior baseline B2: 0.90. Source: results/w1/.
- **4.2 Stage A, context-specific families (training panel).** TODO, after the independent review: results/w4/stage_a_analysis/.
- **4.3 Medium, task and GPR levers.** TODO: W2, W3 and W5 shortlists and validation looks.
- **4.4 Confirmatory chain C1→C1b→C3→C4→C5.** TODO: results/t1/, results/t2/ (sealed).
- **4.5 §6.8 binary endpoints vs the Human1/Human2 papers.** TODO.
- **4.6 Breakdowns per gene, lineage and task.** TODO: sealed until T2.

## 5. Resource
TODO: resource/manifest.tsv (2,030 models, Human2 first), the validation-gate pass rates, and the SBML conventions (closed reactions are kept with bounds (0,0); the notes record the medium tier, pipeline, commit, input SHAs and repaired reactions).

## 6. Failures and negative results
- Human2 does not improve on Human1 (H1.2 not supported).
- M3 growth calibration is infeasible as specified.
- The sterol-branch and CAD/uridine positive controls fail because of model structure.
- CORDA is too costly for the full panel.
- Gene-level E-Flux adds no information.
- TODO: Stage A, W3, W5, W7.

## 7. Limitations
TODO: the addressable-set recall ceiling (~42% of positive pairs are outside it); near-binary KO scores and ties; uptake bounds are availability, not rates; media annotations; lineage-level power.

## 8. Runtime and compute
TODO: aggregate the `results/*/PROVENANCE.tsv` and MANIFEST files. So far W1 took ≈ 3.5 CPU-h, the W2 sweeps ≈ 17 CPU-h and W4 Stage A 108.7 CPU-h.

## 9. Reproducibility
TODO: `scripts/reproduce.py`, the pinned environment (requirements.txt), FROZEN.tsv, the ledgers and the git commits per result.

## References
TODO: from LITERATURE.md.
