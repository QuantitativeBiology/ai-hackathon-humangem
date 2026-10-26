# Research tree

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-26, M0: project start)

Project started from PROMPT.md. Infrastructure exists: Python 3.11 env (`.venv`, cobra 0.32.1, Gurobi 13 WLS academic), Human2 = Human-GEM v2.1.1 (Luo et al. PNAS 2026) in `data/Human-GEM-main/`, DepMap Public 26Q3 (CRISPR, expression, CN, mutations, model annotations) in `data/depmap/`. FBA on Human2 with Ham's medium works (growth 74.76; CAD/UMPS KOs lethal). Human1 (last v1.x release) not yet downloaded; proteomics, metabolomics, drug screens and culture-media annotation not yet sourced.

Current phase: literature and dataset survey (→ LITERATURE.md), with the generic-model baseline pipeline (FBA + single-gene KO on Human1/Human2) built in parallel, since the baseline is fixed by the goal and uses no fitted parameters. DESIGN.md and its review rounds follow the survey; no DepMap labels are used to fit anything before the held-out splits are fixed in DESIGN.md.

## Decisions
- **D0 (2026-10-26)**: Human2 := Human-GEM v2.1.1 (local); Human1 := latest Human-GEM v1.x release (to fetch). Baseline infrastructure started before DESIGN.md is final because it uses only generic models and no DepMap labels for fitting. Git: commit locally; push only when the user asks. Time budget: multi-day effort (user, 2026-10-26); report runtime with every result.

## Tree
- [~] 0. Setup & record keeping (PROMPT.md, RESEARCH_TREE.md, LESSONS.md, DATA.md)
- [~] 1. Literature & datasets survey → LITERATURE.md
  - [x] 1.1 Literature survey (Opus, 16 min, ~8.6k words, 88 refs) → LITERATURE.md. PI spot-checks passed (v1.19.0 = last v1.x tag; 419/2,798 screened Human2 genes are DepMap common essentials).
    - [!] No leakage-free baseline exists: published GEM-vs-CRISPR is MCC ~0.25–0.44, best-of-settings, no lineage hold-out; pooled metrics inflated by gene priors (ROC 0.94 pooled vs 0.61 per gene) → per-gene AUPRC + gene-prior baselines must be primary.
    - [ ] Method-vs-tissue clustering never tested on Human1/2 + ftINIT at DepMap scale → variance partitioning + permuted-expression control.
    - [ ] Per-screen medium available locally (ScreenSequenceMap → ModelCondition.GrowthMedia); fetch CRISPRGeneDependency.csv.
    - [ ] Paywalled papers to request from user (list in LITERATURE.md, batched with data requests).
  - [~] 1.2 Dataset survey/download agent (Sonnet) → data/, DATA.md, scripts/download_*.py
- [ ] 2. DESIGN.md draft
  - [ ] 2.1 Review round 1 (external PI; prior work; statistics; data; algorithms/compute; AI-agent failure modes; domain experts)
  - [ ] 2.2 Review round 2 (if needed)
  - [ ] 2.3 Review round 3 (if needed)
- [~] 3. Baseline: generic Human1 & Human2 (pipeline agent, Opus; predictions only, no DepMap scoring until DESIGN.md fixes labels/splits)
  - [ ] 3.1 Fetch Human1; harmonise gene IDs to DepMap (Entrez)
  - [ ] 3.2 FBA + single-gene KO of all metabolic genes; toy-model verification
  - [ ] 3.3 Benchmark vs DepMap CRISPR (AUPRC overall / per gene / per lineage)
- [ ] 4. Medium (exchange bounds) — per-line culture media; error attributable to medium
- [ ] 5. Objectives / metabolic tasks (a-priori task set)
- [ ] 6. Context-specific models (omics integration; method-vs-tissue clustering test)
- [ ] 7. GPR curation (proteomics, paralogs); curation vs constraint contribution
- [ ] 8. Flux-state information content (subtypes, drug screens, metabolomics)
- [ ] 9. Resource: SBML models per DepMap 26Q3 line
- [ ] 10. REPORT.md
