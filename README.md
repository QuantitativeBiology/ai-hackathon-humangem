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

*Status as of 2026-10-27, 19:00 UTC.*

## Status at a glance

| Milestone | State |
|---|---|
| M0 Setup, literature ([LITERATURE.md](LITERATURE.md)), data ([DATA.md](DATA.md)) | done |
| M1 Pre-registered design ([DESIGN.md](DESIGN.md), frozen v3 after 3 review rounds); leakage infrastructure; lineage split (seed 38: 7 test / 4 validation / 12 training lineages) | done |
| M2 Generic-model baseline (W1), training plus one validation look | done: **G\* = Human1** |
| M3 Pilots: medium (W2, library frozen), context-model families (W4), Stage A on the 30-line training panel | builds done; scoring next |
| M4–M9 Biomass/tasks (W3), GPR curation (W5), flux information (W6), auto-research loop (W7), SBML resource (W8), test looks T1/T2, REPORT.md | planned |

The live research record is [RESEARCH_TREE.md](RESEARCH_TREE.md). It holds the decisions D0–D8 and the findings tree; resume from there.

## Key results so far

**Generic models** (tie-aware AUPRC, mean per line, lineage-weighted):

| | Training (686 lines, 12 lineages) | Validation (198 lines, 4 lineages) | MCC, P57 × gene effect < −0.6 (paper protocol, training) |
|---|---|---|---|
| Human1 | 0.300 | 0.302 | 0.329 |
| Human2 | 0.292 | 0.295 | 0.331 |

- **Human1 vs Human2.** Human1 leads in all 12 training lineages and on validation: Δ(H2−H1) = −0.0071 [−0.0079, −0.0063]. So "Human2 improves on Human1" is not supported. MCC is tied and within the published range of ~0.25–0.44.
- **Strong gene-prior baseline.** A label-only gene-prior baseline (B2) reaches AUPRC 0.90, about 3× the generic GEM. Every context-specific gain is read against B2 and the expression tie-break baseline B8, not only against G\*.
- **Recall ceiling.** About 42% of essential gene × line pairs involve genes outside the model's GPR-addressable set.
- **Medium.** Medium choice moves few genes (mean |Δratio| ≈ 0.001 per line). The design's growth-rate calibration is infeasible as specified, so that variant is reported as uncalibrated.
- **Sterol branch.** Serum cholesterol cannot rescue the sterol-synthesis genes in either model. Biomass needs vitamin D3 products that can only be made from lathosterol, so this is a model-structure issue for W3/W5.
- **GPRTransform** (Machado et al. 2016). Gene-level KO equals the Boolean KO exactly on both models, so it is not adopted as a family. The original isozyme encoding inflates capacity (+25% WT growth). A gene-pFBA cost tie-break is kept as an exploratory score.
- **CORDA.** It costs 1.5–1.7 CPU-h per model, so it runs as a 10-line descriptive comparison only (user decision D8).

## Running now
- **Stage A finished** (2026-10-27, 16:31 UTC): 1,380 distinct context-model builds plus genomic re-sweeps on the 30-line training panel, 0 failed and 0 timeouts, 108.7 CPU-h, 6.2 h wall, from the reviewed commit `c280413`. Outputs are in `results/w4/stage_a/`. Next comes label-blind gating and nested panel LOLO selection (training labels only), then independent review.
- **W2 medium library frozen** after the two-agent extraction and re-check (A9).
- **Starting:** W3 (objectives, tasks, biomass), W5 (GPR curation) and the W8 SBML export plus REPORT.md skeleton.

## Known pitfalls and open issues
- Most generic KO scores are near-binary (lethal or 1), so ties dominate AUPRC. Every metric is tie-aware, and tie-break baselines (B8±) are mandatory.
- Some a priori positive controls are non-informative as built (GOT1/2, SLC7A11 and CAD/uridine), because the models bypass them.
- Several design gaps were resolved by documented amendments A1–A10 (DESIGN.md §17). Small changes to the guards are logged in [notes/guard_changelog.md](notes/guard_changelog.md) (user decision D7).
- Compute budget: ≈1–1.5 days in total. Stage A is ≈110 CPU-h, above its ≈100 CPU-h row because of the CORDA comparison and one build per process.
- Not yet started: SBML export, REPORT.md and the auto-research loop.

## Leakage controls (read before touching data)
- DepMap labels (CRISPR gene effect and dependency, common-essential lists) are read **only** through `src/nextgem/labels.py`, fold-scoped. A PreToolUse hook blocks direct reads, and a runtime tripwire blocks them inside Python.
- Validation and test labels are sealed. Validation is scored only by the frozen `scripts/evaluate_val.py`. Each look is scalar-only, logged in a hash-chained ledger and capped at 40 in total (1 used).
- [FROZEN.tsv](FROZEN.tsv) pins code, inputs and values by SHA-256. It is checked at commit, at push and before every look.
- Every agent or person who has seen label summaries is recorded in [notes/label_exposure.md](notes/label_exposure.md).

## Repository layout
- `src/nextgem/`: package (model loading, knockouts, labels, metrics, splits, media, context-model families, provenance).
- `scripts/`: runnable steps (`w1_*` baseline, `w2_*` medium, `w4_*` context models, `evaluate_val.py`, guards).
- `tests/`: 461 tests, including known-answer tests for metrics and leakage guards.
- `results/`: committed outputs per workstream (`w1/`, `w2/`, `w4/`, `gprt/`, `splits/`, …). Each has a README and a provenance manifest.
- `notes/`: design reviews, milestone process reviews, independent result reviews and method assessments.
- `data/` (git-ignored): Human-GEM, DepMap 26Q3 and other inputs. Provenance is in `data/README.md` and `DATA.md`.

## Setup
- Python 3.11 env at `.venv` (`requirements.txt` pinned).
- Solvers: Gurobi 13 (WLS academic licence; GLPK fallback). Gurobi sessions are rationed through `src/nextgem/gurobi_pool.py`, 8 slots.
- Sanity check: `.venv/bin/python scripts/fba_sanity_check.py` should give growth ≈ 74.76 on Ham's medium.
- Tests: `.venv/bin/python -m pytest -q tests/`.

Process rules for agents are in [CLAUDE.md](CLAUDE.md) and [LESSONS.md](LESSONS.md).
