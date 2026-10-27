# W4 Stage-A analysis (training tier only; DESIGN §9 W4, §13; amendments A4, A7, A8, A10)

**Author:** W4 engineer agent. **Date:** 2026-10-27.
**Inputs:** the Stage-A production run from the pinned worktree at c280413 (outputs committed in 7b22d8b): 1,380 builds plus 600 A4 re-sweeps, 0 failed, 0 timeouts.
**Status:** not committed; no validation look was run (A6).

**Label hygiene**
- The probe `head -c 20 data/depmap/CRISPRGeneEffect.csv` was **blocked** at session start.
- Labels were read only through `labels.train()` (panel lines) and `cv.lolo` / `labels.fold()` (CE_fit per held-out lineage).
- `logs/label_access.log` has 42 W4 entries, all `allow train ... fold-scoped training API`; there is no validation or test access.
- The label-copy scanner reports **0 hits** over this directory.

## Scripts and outputs
- **`scripts/w4_stage_a_labelfree.py`** (no labels): writes `build_gates.tsv`, `config_summary.tsv`, `a7_stability.tsv`, `a6prime_comparison.tsv`, `permuted_control.tsv`, `positive_controls.tsv`, `variance_partition.tsv` and `h4_tests.json`.
- **`scripts/w4_stage_a_train.py`** (training labels):
  - per-line, per-config, nested-selection, modality, CORDA, A7, permuted-control and verification tables: `train_per_line.tsv`, `train_config_scores.tsv`, `train_nested_selection.tsv`, `train_modality.tsv`, `train_corda_a6prime.tsv`, `train_a7_stability.tsv`, `train_permuted_control.tsv`, `train_verification.tsv`;
  - the selection and shortlist, `train_selection.json`.
- **`scripts/w4_val_proposal.py`** (no labels): `val_look_proposal.json` and `val_look_dryrun.json`.
- **Library:** `src/nextgem/w4_analysis.py`, with known-answer tests in `tests/test_w4_analysis.py` (`test_w4a_*`).
- **Provenance:** `MANIFEST.tsv` (dirty = True: the code is uncommitted). Runtimes: label-free part 66 s wall (135 CPU-s, 3.0 GB); training part 35 s (48 CPU-s, 2.1 GB). No solver and no Gurobi session were used.

## 1. Label-free summary

### Gates (§9 W4: growth, the 10 viability tasks, repair ≤ 5% of reactions)
A configuration passes if ≥ 95% of its builds pass (*ours*, after the §11.4 5% failed-line flag). Results (`config_summary.tsv`):
- **F0, F2, F3 (every threshold, modality and genomic variant), F5 (all 24 grid points) and F3b: all builds pass.**
  - Every build grows.
  - Maximum repair: 3.1% of reactions (F0 R3), 1.1% (F3), 0.1% (F2) and 4.3% (F5).
- **F4 (E-Flux) fails the viability gate.**
  - Builds passing all 10 viability tasks: 33–47% on RNA, 0–40% on RNA+protein, 100% on protein-only.
  - **The only failing task is GR,** which requires a biomass output ≥ 1 unit. F4's E-Flux-scaled WT growth is about 1 unit (ratified scaling, README §5 item 5).
  - GR fails exactly when WT growth < 1: 144 of 144 builds below 1 fail, 96 of 96 at or above 1 pass. This is a **units artefact of the ratified scaling, not a biological failure**.
  - Under the frozen gate, F4 is excluded from selection and scored descriptively only. **PI decision needed** if F4 should be gated on scale-free task feasibility.

### Variance partition (method vs tissue; `variance_partition.tsv`, `h4_tests.json`)
- **Method:** marginal PERMANOVA R² on √Jaccard distances, 999 permutations from SEEDS[0].
- **Builds:** F0, F2, F3, F4 × {R1, R3} × {real, permuted}, 30 panel lines.
- **Expression reference:** R²(lineage) is 0.48 on the 30 lines' expression and 0.31 on permuted expression. 0.31 is the chance level for 10 groups of 3 (9/29).

| Level | R²(method) | R²(threshold) | R²(lineage), real builds | R²(lineage), permuted builds | H4.1 diff [percentile CI; basic CI] |
|---|---|---|---|---|---|
| reactions | 0.588 | 0.025 | 0.052 (p 0.001) | 0.036 (p 0.75) | 0.537 [0.501, 0.535; 0.538, 0.572] |
| 57-task profile | 0.254 | 0.077 | 0.073 (p 0.001) | 0.036 (p 0.97) | 0.180 [0.128, 0.189; 0.171, 0.233] |
| binary KO profile | 0.199 | 0.068 | 0.085 (p 0.001) | 0.069 (p 0.68) | 0.113 [0.047, 0.097; 0.130, 0.180] |

- **H4.1** (R²(method) > R²(lineage), line bootstrap B = 2,000): **holds at all three levels.** Building method dominates tissue.
  - The percentile interval sits below the point estimate: resampled duplicate lines are identical points inside one lineage, which inflates R²(lineage) in every replicate. The basic (pivot) interval corrects that bias. Both are reported, and the verdict is the same.
- **H4.2** (within method, R²(lineage | real) > the 95th percentile of R²(lineage) on the permuted-expression builds under 999 lineage relabellings):
  - **F0, F2, F3: tissue-informed** at all three levels;
  - **F4: method-dominated** at all three levels. E-Flux floors almost all bounds, and its binary KO profile equals the generic model's (README §3.3).
- **H4.3** (H4.1 and H4.2 on KO profiles): H4.1 holds; H4.2 holds for F0, F2 and F3.
- **CORDA sensitivity (A7; the 10 A6′ lines, one per lineage, so lineage = line):** R²(method) is 0.62 without CORDA and 0.58 with it (reactions); 0.23 and 0.20 (KO). Adding CORDA changes the method-vs-line picture little.

### Controls
- **Permuted expression** (`permuted_control.tsv`):
  - real vs permuted reaction-set Jaccard: 0.93 (F0 R1), 0.93 (F2 R1), 0.91 (F3 R1); 0.77–0.80 at R3; 1.00 for F4;
  - lethal-set Jaccard: 0.91–0.95 at R1, 0.76–0.79 at R3;
  - every permuted build passes the gates.
- **Positive controls** (`positive_controls.tsv`, R1; mean fraction of the gene's reactions open per line):
  - TYR (11 reactions): 0.33 in the 3 Skin lines against 0.003 in the other 27; DCT: 0.67 against 0.04 (F0, F2);
  - the urea cycle is closed: OTC in 0/30 lines and ARG1 in 1/30 for F0 and F2;
  - F3 keeps OTC in 14/30, because FASTCORE's support includes it;
  - F4 keeps everything (no pruning).
- **A7 stability** (`a7_stability.tsv`; 20 lines, 2 per panel lineage, noise drawn per line):
  - GIMME objective perturbation: reaction Jaccard between seeds 0.9999, lethal-set Jaccard vs unperturbed 0.999;
  - FASTCORE expression noise (N(0, 0.25)): 0.974 and 0.986.
- **A6′ CORDA vs the others** (`a6prime_comparison.tsv`, descriptive, not eligible for P1):
  - reaction-set Jaccard at R1: 0.58 vs F0/F2, 0.82 vs F3, 0.54 vs F4;
  - CORDA costs 86–107 CPU-min per build, against 2–3 CPU-min for the LP families.

## 2. Training-tier scoring (panel: 30 training lines, 10 lineages; L_prob, all genes, U)

**References:** G\* = Human1 (W1 look); also the Stage-A parent, generic Human2 on M0.
- tAP_line, lineage-equal mean: **H1 0.2982**, **H2 0.2916**.
- Verification: my generic per-line tAP equals W1's `train_per_line.tsv` to ≤ 4.9e-7 (W1 stores 6 decimals).

### Per configuration (`train_config_scores.tsv`; paired line bootstrap B = 2,000, lineage-stratified, SEEDS[0])

| Configuration | tAP_line | Δ vs H1 [CI] | Δ vs H2 | mΔAP_gene | MCC (ratio < 0.01 × L_prob) |
|---|---|---|---|---|---|
| generic H1 (G\*) | 0.2982 | — | +0.0066 | 0.000 | 0.316 |
| generic H2 (parent) | 0.2916 | −0.0066 | — | 0.000 | 0.314 |
| F0 R1 | 0.2903 | −0.0080 [−0.0103, …] | −0.0014 | +0.0025 | 0.310 |
| **F2 R1** | **0.2915** | **−0.0067** [−0.0090, …] | **−0.0001** | **+0.0069** | 0.309 |
| F3 R1 | 0.2880 | −0.0102 | −0.0036 | +0.0045 | 0.306 |
| F4 R1 (gate-failed, descriptive) | 0.2954 | −0.0028 | **+0.0038 [+0.0032, …], 10/10 folds** | −0.0001 | 0.314 |
| F5, best grid point (w_p 0, w_c 0, τ −1.5; label-calibrated) | 0.2869 | −0.0114 | −0.0047 | −0.0032 | 0.310 |
| F3b CORDA R1 (10 lines; descriptive, not eligible for P1) | 0.2909 | −0.0127 (same 10 lines) | −0.0033 | +0.0039 | 0.305 |
| B3 (expression; W1) | 0.3252 | +0.0270 | +0.0336 | −0.0029 (gene unit) | — |
| B8′ (H1; W1) | 0.3068 | +0.0086 | | | — |
| B8+ (H1; W1) | 0.4212 | +0.1230 | | | — |
| B2 (gene prior; W1, LOLO) | 0.9041 | +0.6059 | | | — |

- **Threshold:** every family gets worse with more pruning (R1 > R2 > R3; F5's most permissive τ is best).
- **Genomic constraint (A4):** changes tAP by ≤ 0.0002.
- **Minus CE_fit and expressed pairs:** the same picture. F2 R1 vs H1 is −0.0007 without CE_fit and −0.0096 on expressed pairs.
- **mΔAP_gene:** the context models carry a small line-specific per-gene signal, above B3 by +0.010 for F2 R1 (−0.0029 for B3, 0 for the generic models). This is descriptive, not a C3 test.
- **No context model approaches B8′, B8+ or B2.**
- **§6.8 binary endpoints:**
  - fixed cutoff × L_prob: F2 R1 0.309 against G\* 0.316 (Δ −0.0065);
  - PGR (Human1-paper biomass-only; PGR = the FBA lethal set) × L_ge06: F2 R1 0.320 against 0.326;
  - P57 calls were not computed for line models (§6.8 requires them for P1 only).

### Selection (frozen rule): label-blind gates, then nested panel LOLO, ≥ 8-of-10 bar, vs G\* (`train_selection.json`, `train_nested_selection.tsv`)
- **Candidates:** for each family, its gate-passing RNA configurations {R1, R2, R3} × {no genomic, genomic}, or the 24 grid points for F5.
- **Nested choice:** for each held-out panel lineage, the configuration is chosen on the other 9 lineages (ties go to the simplicity order). The panel TG is computed on the held-out values.

| Family | Gate-passing configs | Nested Δ̄_lin vs H1 [CI] | Folds positive | TG vs H1 | Nested vs H2 [CI], folds | Configuration chosen (every fold) |
|---|---|---|---|---|---|---|
| F0 | 6 | −0.0080 [−0.0103, −0.0053] | 0/10 | fail | −0.0014 [−0.0030, 0.0004], 3/10 | R1, no genomic |
| F2 | 6 | −0.0067 [−0.0090, −0.0042] | 1/10 | fail | −0.0001 [−0.0016, 0.0015], 5/10 | R1, no genomic |
| F3 | 6 | −0.0102 [−0.0133, −0.0072] | 1/10 | fail | −0.0036 [−0.0062, −0.0012], 3/10 | R1, no genomic |
| F4 (descriptive: gate failed) | 0 | −0.0029 [−0.0041, −0.0015] | 1/10 | fail | +0.0038 [+0.0032, +0.0043], 10/10 (below the 0.005 bar) | R1 / R2 |
| F5 (label-calibrated) | 24 | −0.0114 [−0.0131, −0.0095] | 0/10 | fail | −0.0048 [−0.0062, −0.0033], 0/10 | w_p 0, τ −1.5 (w_c 0 or 1) |

- **Thresholds kept (≤ 2 per family, full-panel rank):** R1 and R2 (no genomic) for F0, F2, F3 and F4.
- **Modality** (`train_modality.tsv`; 10 ProCan panel lines, 6 lineages; descriptive):
  - **protein-only beats RNA for F0** (R1 +0.0055 [+0.0035, +0.0074]; R2 +0.0024) but stays below H1 (−0.0034);
  - for F2, F3 and F4 the protein arms are within noise;
  - **RNA+protein is worse than RNA in every family** (≈ −0.007 to −0.011);
  - the modality kept is **RNA**: the only modality available for every line. F0 protein-only is noted as descriptive.
- **CORDA** (`train_corda_a6prime.tsv`; 10 lines):
  - R1: −0.0127 vs H1, −0.0006 vs F0, −0.0021 vs F2;
  - R3: −0.0165 vs H1;
  - descriptive only (A7).
- **A7 tAP stability** (`train_a7_stability.tsv`): GIMME's tAP is identical across the 3 seeds; FASTCORE's per-line SD is 0.0018 (maximum range 0.008).
- **Permuted control (tAP)** (`train_permuted_control.tsv`):
  - real beats permuted for F2 (R1 +0.0021 [+0.0008, +0.0034]; R3 +0.0041), F0 R3 (+0.0037) and F3 R3 (+0.0034);
  - F4 ≈ 0;
  - so the expression input carries signal, but not enough to beat the generic G\*.

## 3. Shortlist for Stage B (frozen §9 W4 rule)
- **No family passes the panel TG, against G\* or against the H2 parent.** So no configuration qualifies for validation or test builds: these are made only for TG-passing configurations.
- **Carried to M6, flagged "fails TG"** (§9 W4: "the best configuration is carried to M6 even if it fails TG"): **F2 GIMME, R1, RNA, no genomic constraint**. It is the best gate-passing configuration: nested −0.0067 vs H1, −0.0001 vs H2.
- **Proposed Stage B (for the PI):** build F2 R1 (optionally with F2 R2, the second kept threshold) on all training lines on M0. The definitive 12-fold training TG (≥ 9 of 12) is then computed there. No validation or test builds are made unless it passes.

## 4. Proposed validation look (A6: config entry + label-free dry run; NOT run)
- **Proposal and config SHA:** `val_look_proposal.json`, look `W4_M5_stageA_carry`, config SHA `5a1e9872…`.
- **Dry run:** `val_look_dryrun.json`. It reads no labels.
  - **Not ready.** The Stage-B validation builds do not exist.
  - The W4 code and inputs are not frozen or committed (10 pre-flight problems).
  - `evaluate_val.py` has only the W1 look type. A W4 scoring path must be added, which changes a frozen harness file and so needs a full independent code review (A6 (2)).
- Under the frozen rule there is no TG-passing configuration to look at. **Recommendation: do not spend a validation look on W4 now.**

## 5. Deviations and interpretations (for the PI)
1. **Config gate threshold:** ≥ 95% of builds (*ours*).
2. **F4 gate failure** is a units artefact (GR requires biomass ≥ 1; E-Flux growth ≈ 1). F4 is scored descriptively. It is the family closest to G\* and the only one positive vs the H2 parent (10/10 folds), but its +0.0038 is below the 0.005 TG bar anyway.
3. **TG reference:** G\* = Human1 for the decision, as the best generic pipeline. The W2, W3 and W5 shortlisted generic pipelines are not available yet, so the reference is generic H1 on Ham's M0. Δ vs the H2 parent is reported throughout.
4. **"≤ 2 thresholds per family plus the modality":** read as the two best thresholds on the full panel, plus a modality comparison on the 10 ProCan lines (6 lineages, so the ≥ 8-of-10 bar cannot apply). RNA is kept.
5. **H4.2 null:** R²(lineage) on the permuted-expression builds under 999 lineage relabellings. Its 95th percentile (or the observed permuted R², if larger) is the bar.
6. **H4.1 interval:** the percentile interval as frozen, plus a basic interval because of the duplicate-resampling bias.
7. **§6.8:** P57 not computed for line models; PGR = the FBA lethal set on Ham's (equal on the generic parents, W1).
8. **CORDA:** its partition is a sensitivity analysis on the 10 A6′ lines, where lineage = line (A7).

## Independent review and PI notes (2026-10-27)
- **Review:** notes/reviews/M3_W4_stageA_review.md (GO). The reviewer's own code reproduces per-line tAP (to 2e-16 on 44 configurations), the nested Δ̄ for F2/F3/F4/F5, every per-fold choice, the variance-partition R² and the permutation nulls. F0 and CORDA were not recomputed.
- **cv.lolo:** this analysis slices a `labels.train()` frame rather than calling `cv.lolo`. The reviewer's harness version, run through `cv.lolo` / `labels.fold(k)`, reproduces every choice, so the selection is fold-scoped.
- **Gate rule:** the "≥ 95% of builds" config rule changes nothing. Every non-F4 config passes 100% of its builds, and every F4 RNA config passes ≤ 47%.
- **F4 (E-Flux):** it fails only the GR viability task, a scaling artefact: GR fails in exactly the 196/196 builds with growth < 1. The PI kept the frozen gate (no post-hoc change after seeing scores). F4 would not be selected without the gate either: vs H1 −0.0029, CI < 0, 1/10 folds; vs H2 +0.0038 < 0.005. Because the gate is kept, the configuration carried to M6 as "fails TG" is **F2 R1** (−0.0067). Without the gate it would have been F4 R1.
- **H4.1 intervals:** the frozen percentile interval is primary. A duplicate-free half-sampling interval (2 of 3 lines per lineage) is a sensitivity. The "basic" interval above is not a valid correction. The shift comes mainly from the small-sample upward bias of R²(lineage), which makes the comparison conservative. H4.1 holds under every interval (KO-level lower bounds 0.048 / 0.073).
- **H4.2:** two calls are borderline, F3 at the KO level (0.182 vs bar 0.180) and F4 at the task level. F4 is method-dominated by construction.
- **Provenance:** the outputs were regenerated from a clean pinned worktree at 37a39f3 (MANIFEST rows dirty = False). They are byte-identical to the author's run, apart from the runtime field and the derived config SHAs of the look proposal.
