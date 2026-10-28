# W6: information in flux states beyond the inputs (provisional run = P1 arm)

Author: W6 engineer (Claude Opus), 2026-10-28. Design: DESIGN.md §9 W6, §6.5 (W6 tests), §7 (seeds), §11, §13 (W6 budget
40 CPU-h), amendments A2 (the a priori lists are `notes/apriori_lists_final.md`), A9 (frozen W2 media library), A11 (W3
B-loo flags), A12 (scope). Nothing is committed by the author. DESIGN.md, FROZEN.tsv and the ledgers were not edited.

**Status: provisional run complete** (pinned clean worktree `/data/tmpA/nextgem_runs/71c2a7b`, commit 71c2a7b, 2.06
CPU-h, manifest rows dirty=False). The outputs, labelled "provisional (pre-P1)" in their file headers, are in
`provisional/`; `provisional_pilot/` holds the earlier label-free pilot only.

**P1 equivalence (coordinator, tree D11).** P1 = Human1 + W3 B-loo:cofactor_pool_biomass on M2(f = 1), gate-free,
generic parent, no context family. The arm `H1+B-loo` uses exactly this bound scheme: W2 ladder step "M2" (not "M2g"),
`gate_open = true`, which means no Asp/Glu gate is applied. Its flux features are therefore P1's. The arm `H1` (generic
Human1, same bounds) is **G\***. In this README **P1 = arm `H1+B-loo`** and **G\* = arm `H1`**; the directory and
file names keep the arm names. Flag W3-MEV applies to P1 (statin and FDPS targets); W3-D1 is carried.

## 0. Results (both arms; `provisional/analysis/decisions.tsv`, `secondary.tsv`)

**Decision (one-sided Wilcoxon per analysis, Holm across the 3 analyses, α = 0.025, above the permuted control):
"not supported" for both hypotheses, both arms, all 3 analyses.**

| Arm | Hypothesis | Subtypes (150 tasks, ΔF1) | Drugs (106 units, Δρ) | Metabolomics (225, Δρ) |
|---|---|---|---|---|
| G\* | flux adds info: RNA+X_flux − RNA | median 0 (40 tasks exactly 0), p_Holm 1 | median 0 (77/106 exactly 0), p 0.13, p_Holm 0.38 | median −0.007, p_Holm 1 |
| P1 | flux adds info | median 0, p_Holm 1 | median 0 (77/106 exactly 0), p 0.34, p_Holm 1 | median −0.008, p_Holm 1 |
| G\* | simulation > mapping: X_flux − X_rxnRNA | −0.044 | −0.239 | −0.399 (all p = 1) |
| P1 | X_flux − X_rxnRNA | −0.044 | −0.240 | −0.398 (all p = 1) |
| G\* / P1 | X_KO − X_rxnRNA | −0.051 / −0.047 | −0.240 / −0.240 | −0.405 / −0.408 (all p = 1) |

**Permuted-flux control (SEEDS[0:5]).**
- H_flux permuted medians: 0 for subtypes and drugs; −0.013 to −0.023 for metabolomics.
- Metabolomics is the only "above permuted" case (real −0.007 > permuted −0.02). Its median is negative, so the
  verdict is not supported.

**Synthetic positive control** (flux-only tanh signal, 20 units): recovered. Δρ(RNA + flux − RNA) median +0.42 (G\*)
and +0.45 (P1), p = 9.5e-7.

**Synthetic negative control** (a block that is a linear function of RNA): no gain. Median Δ is −0.007 for drugs
(p 0.97), +0.0003 for metabolomics (p 0.14) and 0 for subtypes (p 0.25).

**Generic KO check.** The G\* KO e-vectors equal the W2 sweeps bitwise in all 20 bounds states.

**Why the result is null.**
- Generic flux depends only on the medium: one state for every PRISM line, 2 in GDSC, at most 20 for the culture
  media. With one state per assay, Δ = 0 exactly.
- 77 of 106 drug units have Δ exactly 0: PRISM, plus GDSC and CTRP units where the medium state is constant within the
  unit's PCA. The non-zero units are GDSC (RPMI vs DMEM:F12) and CTRP and metabolomics (culture media).
- Comparators. Median held-out ρ: RNA 0.11 for drugs (CTRP 0.37, GDSC2 0.22, PRISM 0.03); for metabolomics, RNA 0.28,
  X_rxnRNA 0.23, X_flux −0.18.
- **Note on H_sim magnitudes.** The pooled LOLO ρ of a predictor without line information is negative, not 0. Its
  per-fold intercept anti-correlates with the held-out lineage's mean. So the H_sim Δ values partly measure "no
  information" against "RNA-derived information".

**Secondary** (unadjusted, exploratory):
- RNA + X_KO − RNA: drugs median 0 (p 0.006 G\*, 0.002 P1); metabolomics +0.0005 G\* (p 0.002), +0.003 P1
  (p 9e-6).
- This is still medium-only information, and it is smaller than RNA + X_rxnRNA on metabolomics (+0.008, p 5e-5).
- PROT + X_flux shows no gain. PCA-50 matches PCA-20. Off-model drugs show no gain.

**Fit-free test and positive-control pairs.**
- Target e is constant within every lineage for 100/104 (G\*) and 97/104 (P1) target × unit rows.
- For the positive-control pairs: 46 (G\*) and 43 (P1) are non_informative and 10 have no eligible unit.
- Evaluated pairs:
  - TYMS–5-FU, GDSC2: ρ_e +0.23 vs ρ_expr −0.11. This is a medium effect: DMEM:F12 contains thymidine.
  - TYMS–5-FU, CTRP: ρ_e −0.04.
  - ATP5F1A–oligomycin, CTRP: ρ_e −0.13 vs ρ_expr +0.17.
  - P1 only, HMGCR–statins, CTRP: ρ_e +0.01 to +0.10 vs ρ_expr −0.11 to −0.21. Here e varies only between serum and
    serum-free media; this is W3-MEV.
- Target expression across all fit-free rows: median ρ_expr −0.015.

**Interpretation.** With generic-parent P1, flux states carry no information beyond the omics inputs and the medium.
The pipeline can detect such information when it is present (the synthetic positive is recovered). The null is a
property of line-invariant generic models, not of the test.

**Cosmetic issue.** `w6_simulate.py` read the compound TSV back with NaN `targets`, so `features/checks.json` lists
47 targets including a spurious "nan". It maps to no gene and is skipped; the features are unaffected. The fix is
`target_genes` skipping NaN (not changed, to keep the run's commit).

**Label hygiene.** The session-start probe `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **blocked** by the
PreToolUse hook. W6 reads no CRISPR label, no label-derived source (§11.2 denylist) and nothing in `data/locked/`.
`tests/test_w6.py::test_w6_code_has_no_unseeded_randomness` asserts that the W6 code never names a label file or the
labels API. W6 uses every cohort line regardless of tier, because its targets are drug screens, metabolomics and
subtypes, not labels (§11.2: "Drug screens are W6 targets only").

## 1. How to run (one command)

```
# pinned clean worktree (§11.4), data/ symlinked as for c280413:
cd /data/tmpA/nextgem_runs/<commit>
PYTHONPATH=$PWD/src /data/tmpA/nextgem/.venv/bin/python scripts/w6_run_all.py --config provisional   # pre-P1
PYTHONPATH=$PWD/src /data/tmpA/nextgem/.venv/bin/python scripts/w6_run_all.py --config p1            # after M6
```

- `w6_run_all.py` runs three steps: `w6_inputs.py` (only if its caches are missing, or with `--refresh-inputs`), then
  `w6_simulate.py`, then `w6_analyze.py`. Each step writes one provenance row to `MANIFEST.tsv`.
- **Guards.**
  - The run refuses a dirty `src/`, `scripts/` or `tests/`. `--allow-dirty` exists for pilots only and is recorded in
    the manifest row.
  - It refuses to run if the imported `nextgem` is not the package next to the script.
  - Each process uses 1 BLAS thread, with at most 6 worker processes.
  - Every job declares a CPU budget and stops at 2× that budget.
- **Solver.** GLPK only; no Gurobi session is opened. Every LP has a 300 s limit and one retry (`w3.solve_status`). A
  non-optimal status after the retry aborts the job.
- **P1 re-run.** Copy `config/p1_TEMPLATE.json` to `config/p1.json` and fill it from the frozen P1 spec (see §6).

## 2. What is built

| Item | Where |
|---|---|
| Library: compound table, screen matching, assay media, state simulation (WT pFBA, KO e-vector, KO pFBA of drug targets), comparator blocks, state-space PCA, LOLO ridge, subtype kNN, controls, decision statistics | `src/nextgem/w6.py` |
| Inputs (label-free): compounds, matches, responses, assay media, omics and reaction-expression blocks, subtype task list (written **before** any fitting), coverage | `scripts/w6_inputs.py` → `inputs/` (+ caches in `tmp/inputs/`) |
| Simulation, ≤ 6 GLPK workers | `scripts/w6_simulate.py` → `<tag>/features/` (+ `tmp/<tag>/koflux_long.parquet`) |
| Analyses, controls, decisions | `scripts/w6_analyze.py` → `<tag>/analysis/` |
| One command | `scripts/w6_run_all.py`; configs in `config/` |
| Tests: 20 toy and known-answer tests, all passing | `tests/test_w6.py` |

**Features.** A *state* is one arm × one set of medium bounds. Media are built with the frozen W2 library through
`media.uptake_bounds` and deduplicated by bounds signature.

- **X_flux:** asinh(v_pFBA / μ). This is the WT pFBA flux normalised by biomass flux; asinh is a fixed,
  fold-free variance stabiliser.
- **X_KO:** e = 1 − ratio_r, from the same KO engine as W3 (`w3.ko_sweep`).
- **X_KOflux:** the KO pFBA flux of the 46 curated drug-target genes (compute SF6), normalised by the **WT** biomass
  flux so that lethal KOs stay finite, then reduced by PCA. One LP is solved per distinct closure set (all complex-I
  subunits share one LP, represented by MT-ND1). An exact shortcut reuses the WT pFBA solution when the KO closes no
  reaction that carries flux. A test asserts that the shortcut is exact.
- **X_rxn:** reaction presence. For generic models every reaction is present, so the block is constant.

**Comparators.**

- **X_RNA:** log2(TPM+1) of 19,218 protein-coding genes for the 1,292 cohort lines with expression.
- **X_prot:** ProCan-DepMapSanger proteins detected in ≥ 80% of the 598 cohort lines (3,631 proteins). A missing value
  takes the line's own 1st percentile, so imputation uses no cross-line statistic.
- **X_rxnRNA / X_rxnProt:** GPR-mapped E-Flux capacity (AND = min, OR = sum, as `context.capacity`), log2(cap + 1).
  This gives 7,927 reactions for RNA and 1,652 for protein. An OR with an unmeasured child is unknown, so protein
  coverage is sparse.

**Assay media** (DESIGN W6; counts in `inputs/assay_media.tsv`).

| Analysis | Medium | Source |
|---|---|---|
| PRISM | RPMI 1640 + 10% FBS, W2 key `RPMI` | PRISM SOP on disk (`data/prism/PRISM_SOP_Supplemental_Materials.pdf`): RPMI 1640 without phenol red (Gibco 11835055) + 10% FBS |
| GDSC2 | `RPMI` or `DMEM:F12` | GDSC "Screen Medium": R → RPMI, D/F12 → DMEM:F12 |
| Metabolomics | the line's culture medium | `results/w2/media_library/line_media.tsv` (22 keys) |
| CTRP, subtypes | the line's culture medium | DESIGN is silent. CTRP screened each line in its own growth medium (interpretation) |

**Generic models are line-invariant except through the medium.** Every flux-type block of a generic arm (X_flux, X_KO,
X_KOflux, X_rxn) is a function of the assay medium only.

- Under PRISM, every line has the same state, so Δ(RNA + flux − RNA) = 0 exactly. A test asserts this.
- Under GDSC there are 2 states. Under the culture medium there are at most 20 distinct states (22 keys).
- The flux blocks therefore carry at most medium identity. No line-specific constraint outside the spec was added.
  This is part of the answer, not a defect. Line-specific flux needs a line-specific P1, such as a context family,
  which A12 has not selected.

## 3. Data coverage (`inputs/coverage.json`, `inputs/compound_matches.tsv`, `inputs/subtype_tasks.tsv`)

| Analysis | Coverage |
|---|---|
| Cohort | 1,361 lines with a W2 medium. 1,292 have RNA and 598 have ProCan proteomics. LOLO folds use the split's lineage column |
| 1. Subtypes | **150 tasks** in 23 lineages: 17 OncotreeSubtype and 133 `OmicsInferredMolecularSubtypes` (PAM50 merged into one categorical; other columns binary). A task needs ≥ 2 classes of ≥ 5 lines with RNA. Tasks have 11–123 lines (median 64) |
| 2. Drugs | 105 compounds in §4, verbatim, 8 of them off-model controls. **121 eligible units** (≥ 100 lines with features): PRISM 85, GDSC2 17, CTRP 19. That is 106 primary units from 77 compounds, plus 15 off-model-control units. Lines per unit: PRISM 231–717, GDSC2 294–530, CTRP 287–588 (medians 443 / 443 / 574) |
| 3. Metabolomics | 225 metabolites; 652 cohort lines, 651 with RNA |
| Drug targets (X_KOflux) | 46 genes, all resolved in Human1 v1.19.0 |

**PRISM source per compound.** The 20Q2 secondary AUC is used if it has ≥ 100 lines with features (47 units).
Otherwise the 24Q2 primary log-fold change is used (38 units). GDSC2 and CTRP use AUC.

**Matching.** Compounds are matched by normalised name, then aliases, then the '/'-separated parts of the list token.
When several screen IDs match, the one with the most lines is used.

## 4. Analyses, controls and decision statistics (as pre-registered)

**Matched capacity.**
- Every block is reduced to PCA-20, with PCA-50 as a sensitivity analysis. PCA is fitted on training lines only, in
  state space from the block Gram matrix. This is exact with repeated states; a test checks it against sklearn.
- Drugs and metabolites use ridge regression with the lineage one-hot as a covariate. Alpha is chosen by exact
  leave-one-out over 17 values; a test checks it against sklearn `RidgeCV`. Folds are LOLO.
- Subtypes use stratified 5-fold CV within lineage, kNN with k = 5, and SEEDS[0:3]. The metrics are macro-F1 (the
  decision metric) and silhouette.
- LOLO is deterministic, so seeds have no effect there.

**Metrics.**
- **Drugs:** held-out Spearman ρ, pooled across held-out lineages (decision metric), and the mean within-lineage ρ
  (secondary).
- **Fit-free test:** the within-lineage ρ between the target's e (from X_KO under the assay medium) and sensitivity
  (−AUC or −LFC), next to the same ρ for the target's expression.
- **Positive-control pairs:** transcribed from §4 (43 gene × compound rows).

**Controls.**
- **Permuted flux:** the line → state assignment is shuffled within lineage (SEEDS[0:5]).
- **Synthetic positive:** a tanh-injected signal carried only by the flux block (20 units, seeds spawned from
  SEEDS[0]).
- **Synthetic negative:** a 200-column random linear map of model-gene expression ("flux as a linear function of
  expression").
- All three have known-answer tests in `tests/test_w6.py`:
  - the injected signal is recovered, and permutation removes it;
  - the linear-of-expression block adds nothing;
  - a line-constant block adds exactly nothing.

**Decision** (`<tag>/analysis/decisions.tsv`, per arm). For each analysis there is one one-sided Wilcoxon signed-rank
test of the per-unit Δ, with Holm correction across the 3 analyses at α = 0.025 one-sided (the §6.3 convention).
"Supported" requires all of: Holm p < α, median Δ > 0, and the real median above every permuted-seed median.

| Hypothesis | Δ |
|---|---|
| Flux adds information | metric(RNA + X_flux) − metric(RNA) |
| Simulation adds beyond network mapping | metric(X_flux) − metric(X_rxnRNA); also X_KO − X_rxnRNA |

Secondary contrasts are listed in `secondary.tsv`, unadjusted:
- RNA + X_KO, RNA + X_KOflux, RNA + X_rxnRNA, RNA + X_rxn, each vs RNA;
- PROT + X_flux and PROT + X_rxnProt, each vs PROT (units with ≥ 100 ProCan lines; subtype tasks re-checked for
  eligibility);
- the PCA-50 versions;
- the off-model drugs;
- silhouette;
- the synthetic negative and positive controls.

## 5. Cost (MANIFEST.tsv; CPU-h)

| Run | CPU-h | Note |
|---|---|---|
| `w6_inputs.py` (pilot, dirty) | 0.14 | |
| `w6_simulate.py` pilots: 2 runs of 4 states, 6 targets, plus 1 crashed run with no row (~0.04) | 0.25 | about 60 s per state with 2 KO-flux LPs. The generic arm's KO equals the W2 sweeps bitwise (`provisional_pilot/features/checks.json`) |
| `w6_analyze.py` pilots: 3 runs | 0.18 | 6 drug, 10 metabolite and 10 subtype units and 4 synthetic units per arm. 130 CPU-s after the numpy LOO ridge replaced RidgeCV (366 s before) |
| **Killed analysis pilot (no MANIFEST row)** | **≈ 7** (≥ 6.75 measured by `ps`) | **Incident.** OpenBLAS defaulted to ~40 threads on this 48-core host. The process ran at ~4,000% CPU, mostly thread spin, and was killed. Fixed: every W6 script pins 1 BLAS thread through environment variables and `threadpoolctl` |
| Scratch timing probe and tests | ≈ 0.1 | |
| **Provisional run, pinned worktree 71c2a7b** (inputs 0.02, simulate 1.30, analyze 0.74) | **2.06** | 40 states at about 105 s each (about 21 KO-flux LPs per state); analysis 2,650 s |
| Scanner (`--roots results/w6` + scratchpad): 0 hits | — | |
| **Spent so far** | **≈ 9.8** | of the 40 CPU-h W6 budget. No separate P1 re-run is needed (arm H1+B-loo = P1) |

**Projection for the full provisional run: ≈ 1.9 CPU-h**, so it needs the pinned worktree.
- Inputs: 0.14.
- Simulation: 40 states (20 bounds signatures × 2 arms) at about 110 s each, including about 15 KO-flux LPs, plus
  worker start-up: ≈ 1.25.
- Analysis: (121 + 225) units × 2 arms × ≈ 2.5 s, plus 150 subtype tasks × 2 arms × ≈ 1 s, plus the omics PCA cache:
  ≈ 0.5.
- A P1 re-run with one arm costs about 1 CPU-h, or about 2 if P1 has two arms.
- Total W6 projected: ≈ 11 CPU-h, within the 40 CPU-h budget.

## 6. Interpretations and deviations (for PI ratification)

1. **PRISM medium.** Taken from the SOP PDF on disk. Phenol red is not a nutrient, so the medium maps to W2 classic
   `RPMI`. Serum % is not used for scaling (A9f).
2. **Bounds scheme for the provisional arms: W2 M2 with f = 1, gate open.** All assays use FBS media. P1's W2 choice
   is open, so this is a config field.
   - **W3-MEV applies:** under serum media, B-loo removes mevalonate-core essentiality. HMGCR (the statins) and FDPS
     targets are affected in the H1 + B-loo arm.
   - W3-D1 is carried with the arm.
3. **CTRP and subtypes use the line's culture medium.** DESIGN names media only for PRISM, GDSC and metabolomics.
4. **Compound list, read verbatim from §4.**
   - A parser of the §4 table is checked against a transcription (105 tokens). The off-model-control class is
     excluded from the primary drug test and reported separately.
   - Where the table gives no target, the target is taken from the class or from literature, and `target_source` says
     which. Literature targets: raltitrexed → TYMS; thiopurines → PPAT for "de novo purine"; azaserine → PPAT and PFAS;
     lonidamine → HK2; 3-bromopyruvate → HK2 and GAPDH; vorasidenib → IDH1 and IDH2.
   - Fixed resolutions: "complex I" → MT-ND1; "ATP synthase" → ATP5F1A, as the positive-control pair names it.
   - "trifluridine+tipiracil" is screened as trifluridine. "mycophenolic acid / mofetil" uses the first alias found.
5. **"Above the permuted control"** is read as: the real median Δ is above the median of every permuted seed. It is a
   check, not a second test, so there is still one Wilcoxon test per analysis.
6. **Units.**
   - Drugs: compound × dataset.
   - Molecular-subtype tasks: one task per binary inferred-subtype column per lineage, with PAM50 as one categorical
     task.
   - Held-out ρ is pooled across held-out lineages.
7. **H_sim is read as standalone blocks**, X_flux (or X_KO) vs X_rxnRNA. The combination version (RNA + X_flux vs
   RNA + X_rxnRNA) can be read from `units_*.tsv`.
8. **X_KOflux is normalised by WT, not KO, biomass**, so that lethal KOs stay finite.
9. **PCA for LOLO is fitted on every cohort line that has the block, minus the held-out lineage.** It is not fitted on
   the unit's response lines only. This is fold-safe and lets the PCA be shared across units.

**Pilot observations** (2 media states only; plumbing, **not results**):
- With 2 states, X_flux, X_KO and X_KOflux all reduce to one PCA component, the medium indicator. As a result the H1
  and H1 + B-loo arms give identical numbers, and so do the FLUX and KO contrasts.
- The synthetic positive is recovered: median Δρ +0.14 to +0.16 over 4 units, p = 0.0625, which is the smallest p
  possible with n = 4.
- In GDSC2, the generic TYMS e varies within lineage only through the screen medium. DMEM:F12 contains thymidine. Any
  fit-free ρ for antifolates or TYMS drugs is therefore a medium effect and must be read as such.

## 7. What the P1 re-run needs

- **The frozen P1 spec, put into `config/p1.json`:** the parent, the W3 variant (name in `variants_<parent>.json`), the
  W2 bound scheme (M1, or M2 with f) and the gate.
- **If P1 stays generic** (H1 + B-loo ± a W2 bound scheme), the pipeline runs as is: one command, ≈ 1 CPU-h.
- **If P1 includes a context family** (A12(b) panel pre-check → Stage B), states become line-specific.
  - That needs a per-line state builder: the W4 build per line under each assay medium, then the same
    `simulate_state`. `RunConfig` refuses `context_family` until this exists.
  - X_rxn then becomes informative.
  - Cost: about 1,300 lines × 3 assay maps is too much. Restrict it to analysis lines and the medium keys used. Pilot
    ≤ 10 lines first (§11.5). Budget about 25 CPU-h within the remaining ≈ 30.
- **The pinned worktree:** `data/` symlinks and `PYTHONPATH=<worktree>/src`.
