(Human-written text starts)
This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, to which we (INESC-ID & IST - Universidade de Lisboa) contribute, where we used Claude to run bioinformatics projects (https://github.com/evocellnet/ai-hackathon/). All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication and bioinformatics tool.


[Project-specific text] The goal of this specific project was to improve the prediction of cancer essential genes and metabolic synthetic lethal interaction using human genome-scale metabolic models. It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents.


Currently none of the results have been verified by a human. We plan to verify the agent’s methodology and benchmarking, and we will update this section accordingly.


The rest of this repository is entirely AI generated. It includes code for implementation of the predictive model and a AI-written report in the style of a journal article, detailing the final method, the path followed to reach it, and benchmarking against other methods.
(Human-written text ends)

---

# nextgem: context-specific Human-GEM models for cancer gene essentiality

Goal: optimise the human genome-scale metabolic models **Human1** (Human-GEM v1.19.0) and **Human2** (Human-GEM v2.1.1) so they predict metabolic gene essentiality in **DepMap Public 26Q3** CRISPR screens better than the generic model. Performance is measured by AUPRC per gene, per lineage and overall. MCC/F1 are also reported under the protocols of the Human1/Human2 papers, so results are comparable with them.

The work has four strands:
- build context-specific models from transcriptomics, proteomics and genomics;
- fix three limitations of the generic models: medium/exchange bounds, biomass and metabolic tasks, and GPR curation;
- test how much information flux states carry beyond omics;
- deliver SBML models for the DepMap 26Q3 lines.

Whole lineages are held out for testing. The full brief is in [PROMPT.md](PROMPT.md).

*Status as of 2026-10-28, 12:30 UTC: complete. All pre-registered looks done (T1, T2, unseal); REPORT.md final (non-author reviewed); post-hoc follow-ups F1–F4 done. A plain-language explainer of the method and how it was found is in [docs/explainer/](docs/explainer/README.md), and a paper-length research briefing for a broad audience in [docs/briefing/BRIEFING.md](docs/briefing/BRIEFING.md).*

## Status at a glance

| Milestone | State |
|---|---|
| M0–M2 Setup, pre-registration ([DESIGN.md](DESIGN.md), amendments A1–A14), generic baseline | done: G\* = Human1 |
| M3–M4 Medium (W2), biomass/tasks (W3), GPR curation (W5), context families (W4 Stage A) | done: W2 and W3 pass the training gate; W4 and W5 do not |
| M5 Context models on top of the passing levers | done: not selected (< +0.005) |
| M6 Factorial → **P1** = Human1 + biomass cofactor-pool removal + documented medium | done (validation look 3) |
| M7 W6 flux information; W7 auto-research → **P2**; test looks **T1, T2**, unseal | done |
| M8 SBML resource ([resource/](resource/)) | done: 48 models, 2,114 lines, 47/48 pass all gates |
| M9 [REPORT.md](REPORT.md) | done: non-author reviewed, fixes applied |
| Post-hoc follow-ups F1–F4 (A15) | done: curation and line-specific constraints add nothing; a safer biomass variant costs ≈0.02 ([REPORT §12](REPORT.md)) |

Validation looks used: 7 of 40 (W7: 2 of 20; look 7 is the post-hoc F2 look). Test looks: T1, T2 and the unseal; no third look.

## Key results

**Headline (T1, first test number, 7 held-out lineages; 394 test lines, of which 380 have expression and are scored):** P1 improves tie-aware AUPRC over the generic Human1 model from 0.297 to 0.347: **Δ̄ = +0.050 [0.048, 0.051]**, sign-flip p < 0.001, "meaningfully improves" (SESOI 0.01). The CI is conditional on the 7 test lineages. Without common essentials (C1b): +0.037 [0.035, 0.040], confirmed. The gain is positive in 7/7 test lineages, so the pre-registered "generalises to new tissues" criterion holds.

| | Validation (4 lineages) | Test (7 lineages) |
|---|---|---|
| Generic Human1 (G\*) tAP | 0.302 | 0.297 |
| P1 tAP | 0.353 | 0.347 |
| P1 − G\* | +0.052 | +0.050 |
| Medium alone − G\* | +0.017 | +0.013 |
| Gene prior B2 tAP | 0.892 | 0.895 |

- **What P1 is.** Two structural fixes to Human1. (1) The biomass reaction's cofactor pool is removed, which deletes false "essential" calls created by the biomass reaction itself (vitamin-D/retinoid catabolism and carnitine artefacts, but also genuine cofactor genes such as heme, CoQ and BH4 synthesis, which DepMap rarely finds essential). (2) Each line gets its documented culture medium instead of Ham's.
- **What P1 is not.**
  - It is not line-specific beyond the medium. Shuffled media give about 96–97% of the medium gain (flag W2-SHUF).
  - It does not rank genes across lines better than expression: C3 failed on test (−0.011).
  - It does not beat all three generic + expression baselines.
  - Under serum, the mevalonate-core genes lose essentiality (flag W3-MEV); a safer artefact-only biomass variant fixes this; on validation it gains +0.031 instead of +0.052 (post hoc, REPORT §12).
- **Context-specific models** (GIMME, FASTCORE, E-Flux, thresholding, CORDA, label-calibrated presence) all lose to the generic model on training lineages. On top of P1's biomass edit they add < 0.005, so they were not selected.
- **GPR curation** as pre-registered loses 0.029. Its paralog AND→OR merges joined non-interchangeable RNA polymerase subunits.
- **Flux states (W6)** add no information beyond omics + medium for subtypes, drug response or metabolomics; the positive control is recovered.
- **Gene priors dominate.** GEM features add nothing to gene-prior ML (W7). P2, the best ML model (B6 with a monotone prior constraint, no GEM features), reaches test tAP of about 0.90. T2: P2 − G\* = +0.602, P2 − B6 = +0.014. Both are estimates only, because C3 failed. P1 closes 8.3% of the gap between G\* and B2.
- **Paper-protocol MCC** (Human1 paper, P57 × gene effect < −0.6, test): G\* 0.322 → P1 0.345. With growth-only calls: 0.334 → 0.379.
- **SBML resource:** 48 gzipped SBML models (Human1 and Human2 P1, one per medium signature), with a manifest mapping all 2,114 DepMap lines. 47/48 pass all six gates; the one miss is LP noise (−8e-9; exactly 0 in exact arithmetic). MEMOTE run on 10 models; FROG run on both classic-RPMI models (Human1 and Human2); every model rebuilt from clean pinned commits.

## Known pitfalls and open issues
- Most generic KO scores are near-binary (lethal or 1), so ties dominate AUPRC. Every metric is tie-aware, and tie-break baselines (B8±) are mandatory.
- Some a priori positive controls are non-informative as built (GOT1/2, SLC7A11 and CAD/uridine), because the models bypass them.
- Several design gaps were resolved by documented amendments A1–A14 (DESIGN.md §17). Small changes to the guards are logged in [notes/guard_changelog.md](notes/guard_changelog.md) (user decision D7).
- Compute budget: the user's target was ≈1–1.5 days; compute met it (≈30 h summed run time, ≈210 CPU-h); elapsed time from first to last commit was ≈46 h, of which ≈13 h were usage-limit stalls and waits (docs/explainer/journey.md); earlier "2.2–2.3 days" figures were agent estimates. Stage A is ≈110 CPU-h, above its ≈100 CPU-h row because of the CORDA comparison and one build per process.

## Leakage controls (read before touching data)
- DepMap labels (CRISPR gene effect and dependency, common-essential lists) are read **only** through `src/nextgem/labels.py`, fold-scoped. A PreToolUse hook blocks direct reads, and a runtime tripwire blocks them inside Python.
- Validation and test labels are sealed. Validation is scored only by the frozen `scripts/evaluate_val.py`. Each look is scalar-only, logged in a hash-chained ledger and capped at 40 in total (6 used).
- [FROZEN.tsv](FROZEN.tsv) pins code, inputs and values by SHA-256. It is checked at commit, at push and before every look.
- Every agent or person who has seen label summaries is recorded in [notes/label_exposure.md](notes/label_exposure.md).

## Repository layout
- `src/nextgem/`: package (model loading, knockouts, labels, metrics, splits, media, context-model families, provenance).
- `scripts/`: runnable steps (`w1_*` baseline, `w2_*` medium, `w4_*` context models, `evaluate_val.py`, guards).
- `tests/`: 567 tests, including known-answer tests for metrics and leakage guards.
- `results/`: committed outputs per workstream (`w1/`, `w2/`, `w4/`, `w5/`, `gprt/`, `splits/`, …). Each has a README and a provenance manifest.
- `notes/`: design reviews, milestone process reviews, independent result reviews and method assessments.
- `data/` (git-ignored): Human-GEM, DepMap 26Q3 and other inputs. Provenance is in `data/README.md` and `DATA.md`.

## Setup
- Python 3.11 env at `.venv` (`requirements.txt` pinned).
- Solvers: Gurobi 13 (WLS academic licence; GLPK fallback). Gurobi sessions are rationed through `src/nextgem/gurobi_pool.py`, 8 slots.
- Sanity check: `.venv/bin/python scripts/fba_sanity_check.py` should give growth ≈ 74.76 on Ham's medium.
- Tests: `.venv/bin/python -m pytest -q tests/`.

Process rules for agents are in [CLAUDE.md](CLAUDE.md) and [LESSONS.md](LESSONS.md).
