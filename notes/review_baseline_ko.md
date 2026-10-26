# Review: baseline single-gene knockout pipeline (Human1 v1.19.0, Human2 v2.1.1)

Reviewer: independent agent (did not write the pipeline). Date: 2026-10-26.
Reviewed: `src/nextgem/{model,knockouts,genes,__init__}.py`, `scripts/{download_human1,run_baseline_ko,verify_knockouts}.py`, `results/baseline/*`.
Scratch code (not part of the project): `/tmp/claude-10105/-data-tmpA-nextgem/c5c0f45b-d54c-4e57-bdef-3d0657510b12/scratchpad/review_baseline/` (`indep_ko.py`, `cobra_sgd_sample.py`, `umps.py`, `rerun.py`).
No DepMap CRISPR values were read; only the header of `CRISPRGeneEffect.csv`.

## Verdict: ACCEPT WITH FIXES

The knockout numbers are correct: an independent full recomputation matches every gene in both models, and so does a second path (cobra `single_gene_deletion`, GLPK) on a stratified sample. Gene mapping is correct, and the Human1 download is pinned and its checksum matches upstream. I found no blockers. One major issue (float noise in ties) has to be fixed before any AUPRC is computed. The rest are minor robustness and wording fixes.

## Independent recomputation

| Check | Human1 | Human2 |
|---|---|---|
| A. Own code path, **all genes**. Own GR-task parser; uptake closed on all 1-metabolite boundary reactions (not only `model.exchanges`); own `ast` GPR evaluator; manual bound closure and restore; serial Gurobi, 1 thread | 2887/2887 genes | 2870/2870 genes |
| WT growth (A vs parquet) | 65.57338319775573 (identical) | 74.75676660379064 (identical) |
| max abs ratio difference (A) | 4.6e-13 | 2.2e-12 |
| lethal (<0.01) / partial / neutral class agreement (A) | 2887/2887 | 2870/2870 |
| B. cobra `single_gene_deletion`, **GLPK**, stratified sample (seed 42: 25 lethal + 25 partial + 25 neutral + CAD/UMPS/DHODH/TYMS/UPP1/UPRT/TYMP/DTYMK) | 83 genes (27 lethal / 27 partial / 29 neutral) | 83 genes (29 / 26 / 28) |
| max abs ratio difference (B), class agreement | 5.5e-11, 83/83 | 5.0e-11, 83/83 |

- No state leaks between workers. The parquet comes from the 15-worker parallel run, and it matches my serial recomputation to 1e-12. In each worker, the medium and objective come from the pickled model, which is pickled after `configure()`.
- With the uptake of every boundary reaction closed (not only cobra's `exchanges`), WT is identical, so there are no open sinks or demands. No other boundary reactions exist.
- Genes whose KO closes no reaction: Human1 1801, Human2 1750. All of them have ratio 1.
- All solves are optimal. No infeasible KOs occur in the genome-scale models, so the "infeasible → 0" path is exercised only by the toy model (`toy_infeasible_ko`). I read the code and it is correct: `_solve` returns 0.0 for INFEASIBLE.
- Every model gene is in the outputs (no drops). No model gene has zero reactions.

## Point-by-point

1. **Knockout logic.** This is correct.
   - `Gene.knock_out()` in cobra 0.32.1 sets `functional = False`, which is reverted on context exit. It closes only reactions whose GPR evaluates False with the set of non-functional genes, so AND, OR and nested rules are handled. My evaluator agrees for all genes in both models. (`verify_knockouts.py` checks only Human2's GPRs; my run covers Human1.)
   - The ratio uses the WT from the same configured model.
   - The tolerances are Gurobi defaults. The lethal threshold of 0.01 is far from any value: no ratio lies in (0, 0.5), and the smallest non-lethal ratio is 0.505.
2. **Recomputation.** See the table above.
3. **Gene mapping.**
   - For every gene, `entrez_all` and `symbol` equal the model's own `genes.tsv`, for both releases.
   - The renamed genes checked against the DepMap headers are present under the same Entrez ID: TAZ→`TAFAZZIN (6901)`, GBA→`GBA1 (2629)`, NDUFA4→`COXFA4 (4697)`, SLC22A18→`SLC67A1 (5002)`, ZADH2→`PTGR3`, NAT8L→`ASPNAT`, TRMT12→`TYW2`, FUT10/11→`POFUT3/4`, B3GNTL1→`QTGAL`.
   - In the one-to-many cases the primary (first) Entrez ID is the protein-coding gene. KLHL23 and USP17L2 are not in DepMap under any ID.
   - 62 Human1 and 63 Human2 genes are absent from DepMap. I checked: none of their symbols appears in the DepMap header under a different Entrez ID, so no symbol rescue is missed. They include POLR2A, UQCRH, FABP5, the 13 MT- genes, AMY*, USP17L* and CDY*, all genuinely absent from CRISPRGeneEffect.
   - There are no symbol collisions.
4. **Reproducibility.**
   - I re-ran `run_baseline_ko.py` unchanged, with only OUT redirected to scratch; wall time 2 min 22 s.
   - The mapping TSVs and `model_differences.txt` are byte-identical, and `summary.tsv` is identical apart from runtimes.
   - The KO parquets are **not bitwise identical**. Ratios differ by up to 5e-13; statuses and lethal sets are identical. Only 312 of 2887 (Human1) and 285 of 2870 (Human2) ratios are bitwise equal. See issue M1.
   - Human1 is pinned: tag v1.19.0 → commit 022ed5c2 (checked via the GitHub API). Re-downloading the SBML from the commit and from the tag gives SHA-256 `5aa161fb…`, which equals both the local file and `data/checksums.tsv`. All Human1 files match `checksums.tsv`, and the DATA.md row is present.
5. **Human1 UMPS discrepancy: CONFIRMED.**
   - MAR04644 (TYMS) has bounds (-1000, 1000) in Human1 and (0, 1000) in Human2.
   - In the pFBA solution of the UMPS KO in Human1, the reactions run in this order:
     - thymidine uptake
     - TK (MAR04483)
     - reverse TYMS, dTMP → dUMP (MAR04644, −6.36)
     - dUMP → deoxyuridine (MAR04642)
     - UPP reversed, deoxyuridine → uracil (MAR04680)
     - UPRT, uracil → UMP (MAR04343, +6.36)
   - Positive control for this mechanism: removing thymidine from the medium makes UMPS lethal (ratio 0).
   - Setting MAR04644 lb = 0 changes the lethal set by exactly {UMPS, DHODH}. The README's "changes nothing else" is slightly overstated: WT drops 65.573 → 65.412, and 65 other ratios shift by ≤0.002 (OXPHOS/ATP synthase, TYMP); TYMS goes from 0.9975 to 1.
6. **Bias risks.** See issues M1, m1 and m2 below.

## Issues

**M1 (major, fix before any benchmark): float noise inside ties makes AUPRC non-reproducible.**
- 2662 (Human1) and 2623 (Human2) genes have ratio ≈ 1, but this is spread over about 60 distinct floats in [1 − 1e-13, 1 + 2e-13], and the values change from run to run.
- I simulated AP with 200 random label vectors at 10% prevalence (no DepMap data). Between the original run and my rerun, AP differs by up to 0.018 (mean 0.0035). After rounding the ratios to 1e-6, the difference is exactly 0.
- Fix:
  - Add a rounded, clipped score column, e.g. `ratio = clip(round(ratio, 6), 0, 1)`, and use it in every downstream metric.
  - Document tie handling. sklearn `average_precision_score` treats equal scores as a single threshold, which is the correct behaviour.
  - Report the precision/recall operating point of the binary lethal call alongside AUPRC, since the score is effectively 3-level (0, partial, 1).

**m1 (minor, interpretation): "partial" genes are mostly artefacts of unbounded uptake.**
- With every Ham's input at 1000, thymidine is used as a bulk carbon source: TYMP flux is 993, and WT falls from 65.6 to 46.7 without thymidine.
- The partial stratum (0.5–1) is OXPHOS, ATP synthase, TYMP and SUCLG1. Another 23 (Human1) and 31 (Human2) "partial" genes have ratios ≥ 0.9924, many > 0.99999, which are negligible effects.
- Fix: in DESIGN.md, pre-register that the baseline score is the rounded ratio (or a binary lethal call), and treat partial ratios as medium-dependent. The medium module should quantify this.

**m2 (minor): no NaN policy at the scoring stage.**
- `essential_genes()` and the summary silently treat NaN ratios (non-optimal solver status, or a WT that does not grow) as non-essential. None occur now, but context-specific models may produce them.
- Fix: `run_baseline_ko.py` should fail, or loudly flag, if any status is not optimal or infeasible. The benchmark should state its NaN handling explicitly (exclude and count them, or treat them as non-lethal), and not do either silently.

**m3 (minor, footgun): wrong default task file for Human1.**
- `configure()`, `hams_medium()` and `single_gene_knockouts()` default `tasks_file` to the **Human2** task file. Calling them on Human1 without `tasks_file` silently uses Human2's GR task. That is harmless today because the GR tasks are identical, but it is fragile.
- Fix: derive the task file from the model (for example a `gem` argument), or make the argument required.

**m4 (minor): the download script never verifies checksums.**
- `download_human1.py` writes files without checking them against an expected SHA-256, and skips existing files without verifying them either. The checksums are recorded only afterwards, by `checksum_data.py`.
- Fix: hardcode the expected SHA-256 of the 8 files (they are provenance, not data) and verify them both after download and for files that already exist.

**m5 (minor, README accuracy).**
- "Setting MAR04644 lb 0 ... changes nothing else" should say "changes no other lethal call (WT −0.25%, 65 ratios shift ≤0.002)".
- Also report Human1's count of genes with no closure (1801).

**m6 (minor, verification coverage).**
- `gpr_closure` runs only on Human2. Extend it to Human1.
- Add a genome-scale infeasible-KO test, for example a KO combined with a forced lower bound, so that the status → 0 path is checked outside the toy model.

## Observations (not pipeline bugs; for DESIGN.md)
- The lethal sets differ between models. Lethal only in Human1: AAAS, ABCB6, CRLS1, NPC1L1, PGS1, PTS, SCARB1, SHMT2. Lethal only in Human2: ACAA2, ACADL, ACADSB, AMACR, CYP4F2, DHODH, ECHS1, HADHA, HADHB, HSD17B7, SHMT1, UGP2, UMPS.
- Fatty-acid β-oxidation genes being lethal in Human2 on Ham's is surprising and worth a mechanistic look before biological interpretation.
- The benchmark denominators should be the genes in the model that also have an Entrez ID in DepMap (Human1 2825, Human2 2807). The 13 MT- genes and POLR2A are not scorable, and this should be stated.
