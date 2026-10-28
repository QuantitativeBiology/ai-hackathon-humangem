# Two fixes to Human1: predicting cancer metabolic gene essentiality

*Written by the AI pipeline (Claude); not yet verified by a human. See the disclaimer at the top of the [repository README](../../README.md). An interactive version is in [index.html](index.html).*

**In one sentence.** Given a cancer cell line's culture medium, the method (called **P1** in the repository) scores each of 2,866 metabolic genes by how much simulated growth the edited Human1 metabolic model loses when the gene is knocked out, and uses that score to rank the genes the line most likely needs.

**In three sentences.** P1 is the published Human1 model with two changes: its growth ("biomass") reaction no longer demands a small pool of vitamin and cofactor derivatives, and the model may only take up what the line's culture medium and serum supply instead of a nucleoside-rich test medium. Knocking out each gene with flux balance analysis then gives a short, mostly tied list of predicted essential genes per line. On 7 cancer tissues held out from all development, this raised tie-aware AUPRC against DepMap CRISPR essentiality from 0.297 to 0.347 (+0.050, positive in 7/7 tissues), by removing errors shared by all lines; it adds no line-specific signal and stays far below a simple gene-frequency prior (0.895).

## How it differs, and why that matters

Standard approaches either run knockouts on the generic model with a default medium, or tailor a model to each line from expression data. P1 does neither: it keeps the network and its gene rules as released, uses no omics, and changes what the model must achieve and what it may eat. In this project every expression-tailored family *lowered* accuracy (−0.003 to −0.013), while the single biomass edit gained +0.032 on held-out tissues. The largest fixable error was in the objective, not in missing context.

## Overview

```
Human1 v1.19.0 ─┬─ 1 generic network, gene rules as released
                ├─ 2 remove the cofactor pool from the biomass reaction (all lines)
line medium ────┼─ 3 uptake = documented base medium + serum (classic RPMI + 10% FBS if unknown)
                └─ 4 knock out each gene: e = 1 − KO growth / normal growth   ──►  ranked genes per line
                   5 read as a ranking with ties, against gene-prior and expression baselines
```

| Step | What | Why (evidence) | Without it |
|---|---|---|---|
| [1](method.md#step-1) | Generic Human1, no pruning, no GPR edits | Human1 beat Human2 on every tier; expression pruning and GPR curation both lost accuracy | −0.003 to −0.029 |
| [2](method.md#step-2) | Remove the cofactor/vitamin pool from biomass | The pool makes 35 genes lethal that DepMap finds essential in only 14% of lines | gain falls from +0.052 to +0.017 (validation) |
| [3](method.md#step-3) | Line's medium + serum instead of Ham's | Ham's nucleosides hide 8 de novo nucleotide genes; serum rescues 13 rarely essential calls | gain falls to +0.032 (validation) |
| [4](method.md#step-4) | Knock out every gene, score growth lost (rounded) | Loss of biomass is the FBA proxy for fitness; rounding makes reruns identical | float noise moved AP by 0.018 |
| [5](method.md#step-5) | Read as a tied ranking vs strong baselines | 94% of scores tie; an expression tie-break alone is worth +0.12 | "gains" that are only tie-breaking |

## Worked example (real run)

NCI-H460 (lung, training tier, "RPMI + 10% FBS"), re-run with `scripts/explainer_worked_example.py`; all scores equal the committed P1 table exactly.

| | Reference (Ham's) | + pool removed | P1 (+ medium) |
|---|---|---|---|
| Genes called lethal | 147 | 112 | 106 |
| CYP26A1 / BBOX1 (vitamin catabolism, carnitine) | lethal | not lethal | not lethal |
| PPAT / DHODH / UMPS / TYMS (de novo nucleotides) | not lethal | not lethal | lethal |
| HMGCR (mevalonate; essential in 90% of lines) | lethal | lethal | **not lethal** (known cost) |
| NCI-H460 tAP | 0.284 | 0.329 | 0.336 |

Details: [method.md § worked example](method.md#worked-example).

## Results

| Lineage-equal tAP | Training (12 lineages) | Validation (4) | Test (7) |
|---|---|---|---|
| Reference: generic Human1 on Ham's | 0.300 | 0.302 | 0.297 |
| **P1** | **0.353** | **0.353** | **0.347** |
| P1 − reference [95% CI] | +0.053 [0.052, 0.055] | +0.052 [0.050, 0.054] | +0.050 [0.048, 0.051] |
| Reference, ties broken by expression | 0.424 | 0.430 | 0.418 |
| Gene prior | 0.900 | 0.892 | 0.895 |

Pre-registered design, one test look for the headline, positive in 7/7 test tissues, +0.037 without common essentials. The line-specificity test failed. Post-hoc follow-ups (stricter GPR curation, a safer biomass variant, genomic-loss and E-Flux constraints) all leave P1 unchanged ([follow-ups](method.md#follow-ups)).

![Accuracy by pipeline and tier](figures/fig1_tap_tiers.svg)

## Limitations in brief

- Not line-specific: a shuffled or constant medium does as well; expression ranks lines better for a given gene.
- Known error: the mevalonate core (HMGCR, MVK, PMVK, MVD) becomes non-lethal under serum. A safer variant fixes it but costs 0.021.
- About 45% of essential gene × line pairs are invisible to any FBA knockout in Human1; the worst errors are housekeeping genes outside the model.
- Far below a gene-frequency prior; use P1 to study what a metabolic model gets right, not as the best essentiality predictor.

## The documents

- [method.md](method.md): Deliverable A in full: steps, worked example, pseudocode, results, limitations, credits.
- [journey.md](journey.md): how the agents arrived at the method: eight turning points, an epilogue, the research process, dead ends.
- [insights.md](insights.md): scientific, methodological and engineering insights, ranked.
- [discrepancies.md](discrepancies.md): discrepancies between documents and results, open questions, unsure points.
- [figures/](figures/): five SVG figures (theme-aware), generated by `scripts/make_figures.py` from committed files.

**Reproduce.** Figures: `/data/tmpA/nextgem/.venv/bin/python scripts/make_figures.py --out md/figures`. Gene-group table: `scripts/gene_group_summary.py --out data`. Worked example: `scripts/explainer_worked_example.py --line ACH-000463 --out data` (in the main repository). None reads a label file.
