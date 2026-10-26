# Research tree

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-26, M0: project start)

Project started from PROMPT.md. Infrastructure exists: Python 3.11 env (`.venv`, cobra 0.32.1, Gurobi 13 WLS academic), Human2 = Human-GEM v2.1.1 (Luo et al. PNAS 2026) in `data/Human-GEM-main/`, DepMap Public 26Q3 (CRISPR, expression, CN, mutations, model annotations) in `data/depmap/`. FBA on Human2 with Ham's medium works (growth 74.76; CAD/UMPS KOs lethal). Human1 (last v1.x release) not yet downloaded; proteomics, metabolomics, drug screens and culture-media annotation not yet sourced.

Current phase: literature and dataset survey (→ LITERATURE.md), with the generic-model baseline pipeline (FBA + single-gene KO on Human1/Human2) built in parallel, since the baseline is fixed by the goal and uses no fitted parameters. DESIGN.md and its review rounds follow the survey; no DepMap labels are used to fit anything before the held-out splits are fixed in DESIGN.md.

## Decisions
- **D0 (2026-10-26)**: Human2 := Human-GEM v2.1.1 (local); Human1 := latest Human-GEM v1.x release (to fetch). Baseline infrastructure started before DESIGN.md is final because it uses only generic models and no DepMap labels for fitting. Git: commit locally; push only when the user asks. Time budget: multi-day effort (user, 2026-10-26); report runtime with every result.
- **D1 (2026-10-26)**: results/ is committed to git (PROMPT.md requirement); files >50 MB go under results/**/tmp/ or *.large.* (ignored) and must be regenerable by a script.

## Tree
- [~] 0. Setup & record keeping (PROMPT.md, RESEARCH_TREE.md, LESSONS.md, DATA.md)
- [~] 1. Literature & datasets survey → LITERATURE.md
  - [x] 1.1 Literature survey (Opus, 16 min, ~8.6k words, 88 refs) → LITERATURE.md. PI spot-checks passed (v1.19.0 = last v1.x tag; 419/2,798 screened Human2 genes are DepMap common essentials).
    - [!] No leakage-free baseline exists: published GEM-vs-CRISPR is MCC ~0.25–0.44, best-of-settings, no lineage hold-out; pooled metrics inflated by gene priors (ROC 0.94 pooled vs 0.61 per gene) → per-gene AUPRC + gene-prior baselines must be primary.
    - [ ] Method-vs-tissue clustering never tested on Human1/2 + ftINIT at DepMap scale → variance partitioning + permuted-expression control.
    - [ ] Per-screen medium available locally (ScreenSequenceMap → ModelCondition.GrowthMedia); fetch CRISPRGeneDependency.csv.
    - [ ] Paywalled papers to request from user (list in LITERATURE.md, batched with data requests).
  - [x] 1.2 Dataset survey/download (Sonnet, 49 min): fetched Nusinow proteomics (375 lines), CCLE metabolomics (928 lines), PRISM 24Q2 + 20Q2 secondary, GDSC1/2 r8.5, CORUM 5.1 + Complex Portal, Ensembl 116 paralogs, Cell Model Passports list (~0.9 GB). Blocked items → notes/USER_DOWNLOADS.md (asked user).
    - [!] Medium annotation: 2,097/2,554 models (1,282/1,361 CRISPR lines); base media RPMI 907, DMEM 291, MEM/EMEM 157, DMEM/F12 104, IMDM 72, F12 56, McCoy's 31, L-15 23; ~300 organoid/proprietary (unmappable); 239 serum-free. Model.csv OnboardedMedia/FormulationID columns are swapped (PI-verified).
- [ ] 2. DESIGN.md draft
  - [ ] 2.1 Review round 1 (external PI; prior work; statistics; data; algorithms/compute; AI-agent failure modes; domain experts)
  - [ ] 2.2 Review round 2 (if needed)
  - [ ] 2.3 Review round 3 (if needed)
- [~] 3. Baseline: generic Human1 & Human2 (pipeline agent, Opus; predictions only, no DepMap scoring until DESIGN.md fixes labels/splits)
  - [x] 3.1 Human1 = v1.19.0 (commit 022ed5c) fetched; Entrez mapping 2,874/2,887 (H1), 2,857/2,870 (H2) one-to-one; 2 unmapped (GK3P, SLC37A4 alt ENSG); join on Entrez (symbols renamed for 7–10 genes)
  - [~] 3.2 Single-gene KO, Ham's medium, Gurobi (16 WLS sessions OK): H1 149 / H2 154 lethal (ratio<0.01; none in 0.01–0.5); ~50 s per model on 15 workers. 36/37 checks pass. Independent review in progress (notes/review_baseline_ko.md).
    - [!] Human1 UMPS/DHODH non-lethal on Ham's: thymidylate synthase MAR04644 reversible in H1 (irreversible in H2) → UMP via thymidine bypass. Kept as released; a model difference to report.
    - Predicted ratios are near-binary → essentiality scores for AUPRC need a continuous score (e.g. ratio + tie-breaking / flux-based scores), to be fixed in DESIGN.md.
  - [ ] 3.3 Benchmark vs DepMap CRISPR (AUPRC overall / per gene / per lineage)
- [ ] 4. Medium (exchange bounds) — per-line culture media; error attributable to medium
- [ ] 5. Objectives / metabolic tasks (a-priori task set)
- [ ] 6. Context-specific models (omics integration; method-vs-tissue clustering test)
- [ ] 7. GPR curation (proteomics, paralogs); curation vs constraint contribution
- [ ] 8. Flux-state information content (subtypes, drug screens, metabolomics)
- [ ] 9. Resource: SBML models per DepMap 26Q3 line
- [ ] 10. REPORT.md
