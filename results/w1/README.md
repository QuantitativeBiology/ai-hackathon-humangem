# W1 — generic Human1 / Human2 benchmark, baselines and M1 leftovers (milestone M2, training tier)

Author: W1 analyst (Claude Opus), 2026-10-26. Design: DESIGN.md v3 (frozen) + §17 amendments A1–A3; a priori lists per A2 (`notes/apriori_lists_final.md`). Nothing here is committed (the PI commits).

**Status.** Everything on the **training tier** is done. The single **validation look is pending**:
- `scripts/evaluate_val.py` is written and was reviewed (`notes/reviews/M2_W1_review.md`, ACCEPT WITH FIXES); the review fixes are applied.
- It has not been frozen or run on real labels (coordinator's instruction).
- Its pre-flight refuses to start until its code and inputs are in the committed `FROZEN.tsv` and clean in git.

**G\* therefore stays pending**, as do all validation numbers.

Label hygiene:
- The raw-label probe `head -c 20 data/depmap/CRISPRGeneEffect.csv` was **blocked** by the hook at session start.
- Labels were read only through `nextgem.labels`: `train()`, and `fold(k)` via `nextgem.cv.lolo`.
- No validation or test labels were read.
- The label-copy scanner reports **0 hits** over 727 files (project and `/tmp/claude-*`).

## Headline (training tier, 12 eligible lineages, 686 lines with expression; label L_prob)

| | Human1 v1.19.0 | Human2 v2.1.1 | H2 − H1 [95% CI] |
|---|---|---|---|
| tAP_line, lineage-equal mean (all U genes) | 0.300 [0.298, 0.301] | 0.292 [0.290, 0.293] | −0.0079 [−0.0083, −0.0075] |
| lift over E[AP_random] = 0.124 (line); pooled lift | 2.42; 2.44 | 2.36; 2.37 | |
| tAP_line without CE_fit (LOLO), null 0.071 | 0.145 (lift 2.05) | 0.142 (lift 2.01) | −0.0031 [−0.0035, −0.0026] |
| MCC, Human1-paper protocol, P57 × L-ge06 | 0.329 [0.326, 0.331] | 0.331 [0.329, 0.334] | +0.0028 [0.0025, 0.0032] |
| MCC, Human1-paper protocol, PGR (biomass) × L-ge06 | 0.334 [0.332, 0.336] | 0.331 [0.329, 0.333] | −0.0027 [−0.0032, −0.0022] |
| MCC, Human2-paper protocol, P57 × L-ge06 | 0.329 | 0.331 | +0.0028 (identical computation to the Human1-paper P57 row) |
| MCC, our fixed cutoff ratio_r < 0.01 × L_prob | 0.322 [0.320, 0.323] | 0.318 [0.317, 0.320] | −0.0033 [−0.0038, −0.0028] |

- **M2 gate (W1 decision): PASS.**
  - H2 lift on training is 2.36 (lineage-equal mean) and 2.37 (pooled), both > 1.
  - No label-free audit is required.
- **H1.1** (lift > 1 for both parents) holds.
- **H1.2** (H2 > H1) is **not supported on training**.
  - H1 is ahead by 0.008 tAP, consistently in all 12 lineages, in all three labels and in every non-NA stratum except the housekeeping-adjacent genes.
  - Under the P57 paper protocol, H2 is ahead by a negligible +0.003 MCC.
- **H1.3** (gene prior B2 > generic FBA) holds overwhelmingly: B2 tAP is 0.900, against 0.30 for the generic models.
- **SESOI_C1b (frozen from training): 0.005723.**
  - It equals 0.01 × 0.07069 / 0.12353, with L_prob, all 767 training lines and CE_fit from all training lines (167 genes).
  - Source: `sesoi_c1b.json`.
- **G\***: pending the validation look. G\* is decided on validation only: Human1 iff Δ̄_lin(H2 − H1) < 0, and a tie goes to Human2.
  - On training, H1 leads.

## What was done (DESIGN §14 M2 row and the M1 leftovers)

| Item | Script | Output |
|---|---|---|
| Label-free inputs: expression, CN, mutation (default entries); gene features (addressable sets, isozyme counts, H2 enzyme neighbours, Ensembl paralogs); line features (medium class, serum flag) | `scripts/w1_inputs.py` | `gene_features.tsv`, `line_features.tsv`, `inputs_summary.json`, `tmp/*.parquet` |
| M1 leftover: gene panels generated from `notes/apriori_lists_final.md` §1 | `scripts/w1_gene_panels.py` | **`data/controls/gene_panels.tsv`** (DESIGN path; SHA-256 `e812724f…a4e78318`), `gene_panels_summary.tsv`, `gene_panels_unresolved.tsv` |
| M1 leftover: evaluable background groups per tier (label-free) | `scripts/w1_backgrounds.py` | `background_groups.tsv`, `background_members.tsv` |
| P57 / PGR task-based KO calls on both generic parents (§6.8) | `scripts/w1_p57.py` | `p57_human{1,2}.parquet`, `.json` |
| B6 (HistGradientBoosting, LOLO; plus a full-training model for the validation look) | `scripts/w1_b6.py` | `b6_run.json`, `tmp/b6_oof.parquet`, `tmp/b6_full.pkl` |
| tAP grid: all methods, labels, gene sets, universes and strata; per line / gene / lineage / pooled; SESOI_C1b; addressable set; gate; BCa check; per-gene error table | `scripts/w1_score_train.py` (library `src/nextgem/w1.py`) | `train_grid.tsv`, `train_per_line.tsv`, `train_per_lineage.tsv`, `train_per_gene.tsv`, `train_addressable.tsv`, `sesoi_c1b.json`, `train_summary.json` |
| §6.8 binary endpoints (both paper protocols, 21-threshold grid, fixed cutoff, controls, context baselines) | `scripts/w1_binary.py` | `train_binary.tsv`, `train_binary_grid.tsv`, `train_binary_paired.tsv` |
| §12 controls and second code path | `scripts/w1_verify.py` [`--task-equality`] | `verification.tsv` (34/34 pass), `verification_tasks.tsv` (2/2 pass) |
| Validation look (written, reviewed, fixes applied, toy-tested, **not run**) | `scripts/evaluate_val.py`, `tests/test_evaluate_val.py` | — |
| Provenance wrapper (DESIGN §11.4), used by every W1 script and by the look | `src/nextgem/run.py`, `tests/test_run.py` | `MANIFEST.tsv` |
| Driver, timing and determinism check; README tables | `scripts/w1_run_all.py`, `scripts/w1_tables.py` | `runtime.tsv`, `rerun_sha_check.tsv`, `MANIFEST.tsv` |

### How it was computed (choices the design leaves to the implementer are marked *ours*)

**Evaluation set**
- Training lines of the 12 eligible lineages.
- Primary: the 686 lines with a default expression profile (§4.3).
- Secondary: all 737 lines, for the generic models and label-only baselines.
- The 30 rare-lineage training lines are used for fitting only; they are never a fold and never evaluated (*ours*). This follows from "LOLO over the 12 eligible lineages".

**Gene universe**
- U = 2,866 genes; U∩ = 2,766.
- A gene absent from a parent scores e = 0 (§3.3).
- Generic score: e = 1 − ratio_r. Task-augmented score (§6.1): e = 1 if P57, else 1 − ratio_r.

**Fitting (LOLO through `nextgem.cv.lolo`)**
- B2, B2m, B2c, CE_fit, the selective-gene set and B6 are all fitted LOLO: a line in lineage k is scored with components fitted on training minus k, including the rare lines.
- *Interpretation:* the brief asked for common essentials "derived from train ∪ val". DESIGN §5.3 assigns train ∪ val to **test** evaluations only. For the training tier the fitting lines are training minus fold, and for validation looks they are training, so that is what is implemented. The train ∪ val list belongs to `evaluate_test.py` at T1.

**Baselines (§8)**
- B0: constant.
- B1: uniform random, SEEDS[0:5], seeded per (seed, line); reported per seed and as the mean.
- B2: gene prior.
- B2m: prior within `base_medium_class`, falling back to B2 for classes with fewer than 10 fitting lines.
- B2c: −mean Chronos.
- B3: log2(TPM+1).
- B8+ / B8−: exact lexicographic rank of (e, ±expression), for **both** parents, because G\* is pending.
- B8′: e with genes at TPM < 1 set to 0.
- B6 features:
  - B2 and B2m, stacked leave-own-lineage-out inside the fitting set (*ours*, so that the learner never sees its own label inside a prior);
  - expression, CN, damaging and hotspot mutation;
  - mean expression of H2 enzyme neighbours, i.e. genes sharing a non-currency metabolite (*ours*: a fixed currency list plus metabolites in ≥ 100 reactions);
  - maximum Ensembl-116 paralog expression, paralogs with %id ≥ 30 (*ours*);
  - isozyme-gene count;
  - medium class and serum flag.
- B6 hyperparameters are fixed at M2 (the grid is W7's): 200 iterations, learning rate 0.1, 31 leaves, no early stopping, random_state 11.
- B6 is trained and evaluated on L_prob only.
- B6's medium class is the split's `base_medium_class`, and its serum flag is parsed from the formulation text. W2's tiered annotation does not exist yet, so B6 must be refitted once it does (§8 requires the same annotation as W2).

**Metrics**
- Frozen `nextgem.metrics`: tAP, E[AP_random], the §5.2 NA rules and the lineage-stratified bootstrap (B = 2,000, SEEDS[0]).
- Deltas are paired per line.
- AUROC and sklearn step AP are secondary or sensitivity rows, computed for the primary cells only.

**P57 calls**
- Same algorithm as Human-GEM `taskEssentialGenes.py`, with strict solver statuses and every touched task tested.
- In both releases, all 57 tasks pass on the wild type.
- The calls equal `find_task_essential_categories` exactly: 226/226 genes (H1) and 229/229 (H2).
- PGR equals the FBA lethal set (149 / 154).
- The task file equals Robinson 2020 data file S3 on all 57 (id, description) pairs. Inputs and outputs are equal up to two metabolite renamings.

**Binary gene set**
- The parent's genes with a CRISPR column (U ∩ parent) and a non-NaN label; U-wide is a sensitivity row.
- NaN-denominator lines are excluded and counted; zero-filled is the sensitivity row.

**Background groups**
- Definitions are those of the final list, §2 (A2).
- Within a tier, a lineage is eligible if it has ≥ 3 members and ≥ 3 non-members. A group is evaluable if ≥ 2 lineages are eligible and they hold ≥ 10 members.
- "TPM < x" is read as TPM (*ours*). NEG-3 burden is the damaging-allele row sum; its 90th percentile is 98, whereas the list file says ≤ 88 (definition difference).

## Runtime

Final clean re-run (`scripts/w1_run_all.py`, 2026-10-26 20:54–21:35):

| Step | Wall (s) | CPU (s) | Solver / workers |
|---|---|---|---|
| `w1_inputs.py` | 54 | 52 | — |
| `w1_gene_panels.py` | 37 | 41 | — (GLPK structure only) |
| `w1_backgrounds.py` | 34 | 33 | — |
| `w1_p57.py --model human1` | 380 | 1,251 | Gurobi, 7 workers + parent = 8 sessions |
| `w1_p57.py --model human2` | 421 | 1,357 | Gurobi, 8 sessions |
| `w1_b6.py` (13 HGB fits, about 11 s each) | 167 | 939 | 8 threads |
| `w1_score_train.py` (345 cells, 7,335 rows) | 498 | 6,215 | 14 processes |
| `w1_binary.py` | 51 | 75 | — |
| `w1_verify.py` | 67 | 69 | — |
| `w1_verify.py --task-equality` | 759 | 2,564 | Gurobi, 8 sessions |
| **Total** | **≈ 41 min** | **≈ 12,600 s ≈ 3.5 CPU-h** | ≤ 14 cores in the clean re-run (22 briefly during development); ≤ 8 Gurobi sessions |

- Peak RSS was ≈ 2.3 GB per process.
- No flock slot pool exists (DESIGN §11.4 gap): the cap of 8 Gurobi sessions is enforced by the process count in `w1_p57.py`. Solver failures raise `SolverInfraError` and abort.
- P57 for the generic parents costs ≈ 0.35 CPU-h per model, far below the "few hours" bound, so it was run in full.

**Determinism.** A full re-run reproduced every output bitwise, including the P57 parquets, the B6 out-of-fold scores and every training table. Two files changed, both explained:
- `train_grid.tsv`, because a column (sklearn AP) was added between the runs;
- `verification_tasks.tsv`, whose detail field records elapsed seconds.

See `rerun_sha_check.tsv`.

## Results (training tier). All tables are printed by `scripts/w1_tables.py` from the output files

### Primary cell: L_prob, all genes, U, 686 expression lines

| method | tAP_line [95% CI] | E[AP_random] | lift | Δ vs H1 [CI] | Δ vs H2 [CI] | mΔAP_gene | pooled tAP (lift) | AUROC | sklearn AP |
|---|---|---|---|---|---|---|---|---|---|
| H1 | 0.300 [0.298, 0.301] | 0.124 | 2.42 | — | 0.0079 [0.0075, 0.0083] | 0.0000 | 0.296 (2.44) | 0.636 | 0.258 |
| H2 | 0.292 [0.290, 0.293] | 0.124 | 2.36 | −0.0079 [−0.0083, −0.0075] | — | 0.0000 | 0.288 (2.37) | 0.633 | 0.252 |
| H1_task | 0.285 [0.282, 0.287] | 0.124 | 2.30 | −0.0149 [−0.0159, −0.0140] | −0.0070 | 0.0000 | 0.281 (2.32) | 0.638 | 0.244 |
| H2_task | 0.280 [0.278, 0.283] | 0.124 | 2.27 | −0.0191 | −0.0112 [−0.0120, −0.0104] | 0.0000 | 0.277 (2.28) | 0.636 | 0.242 |
| B0 | 0.124 [0.123, 0.125] | 0.124 | 1.00 | −0.1760 | −0.1680 | 0.0000 | 0.121 (1.00) | 0.500 | 0.121 |
| B1 (mean of 5 seeds) | 0.123 [0.122, 0.125] | 0.124 | 1.00 | −0.1761 | −0.1682 | −0.0002 | 0.121 (1.00) | 0.500 | 0.123 |
| B2 | 0.900 [0.897, 0.903] | 0.124 | 7.28 | 0.6006 | 0.6085 [0.6062, 0.6107] | −0.0327 | 0.894 (7.37) | 0.978 | 0.900 |
| B2m | 0.898 [0.895, 0.900] | 0.124 | 7.26 | 0.5980 | 0.6059 | 0.0090 | 0.893 (7.36) | 0.974 | 0.893 |
| B2c | 0.895 [0.893, 0.898] | 0.124 | 7.24 | 0.5957 | 0.6036 | −0.0270 | 0.890 (7.33) | 0.975 | 0.895 |
| B3 (expression) | 0.326 [0.322, 0.329] | 0.124 | 2.63 | 0.0261 | 0.0340 [0.0313, 0.0368] | 0.0081 | 0.320 (2.64) | 0.798 | 0.326 |
| B6 | 0.888 [0.885, 0.891] | 0.124 | 7.18 | 0.5884 | 0.5963 [0.5938, 0.5988] | 0.0219 | 0.878 (7.24) | 0.978 | 0.889 |
| B8+ (H1) | 0.424 [0.422, 0.427] | 0.124 | 3.43 | 0.1247 | 0.1326 | 0.0078 | 0.419 (3.46) | 0.828 | 0.424 |
| B8− (H1) | 0.224 [0.223, 0.226] | 0.124 | 1.81 | −0.0755 | −0.0675 | −0.0021 | 0.222 (1.83) | 0.443 | 0.223 |
| B8′ (H1) | 0.307 [0.306, 0.309] | 0.124 | 2.49 | 0.0078 | 0.0158 | −0.0000 | 0.304 (2.50) | 0.637 | 0.264 |
| B8+ (H2) | 0.413 [0.411, 0.416] | 0.124 | 3.35 | 0.1139 | 0.1219 [0.1204, 0.1233] | 0.0078 | 0.409 (3.37) | 0.823 | 0.413 |
| B8− (H2) | 0.217 [0.215, 0.218] | 0.124 | 1.75 | −0.0829 | −0.0750 [−0.0755, −0.0745] | −0.0021 | 0.215 (1.77) | 0.444 | 0.216 |
| B8′ (H2) | 0.301 [0.299, 0.302] | 0.124 | 2.43 | 0.0009 | 0.0089 [0.0087, 0.0090] | −0.0000 | 0.297 (2.45) | 0.636 | 0.259 |

**Without CE_fit (LOLO; L_prob; null 0.071).** Lift is shown in parentheses.

| Group | Method | tAP (lift) |
|---|---|---|
| Generic | H1 | 0.145 (2.05) |
| | H2 | 0.142 (2.01) |
| | H1_task | 0.155 (2.19) |
| | H2_task | 0.153 (2.16) |
| Label-free | B3 | 0.204 (2.89) |
| | B8+ (H2) | 0.207 (2.93) |
| | B8− (H2) | 0.105 |
| | B8′ (H2) | 0.148 |
| Label-fitted | B2 | 0.686 (9.70) |
| | B2m | 0.686 |
| | B2c | 0.675 |
| | B6 | 0.673 (9.53) |

**Without the DepMap common-essential list** (reporting only; null 0.022):

| Method | tAP (lift) |
|---|---|
| H1 | 0.028 (1.24) |
| H2 | 0.027 (1.22) |
| B3 | 0.046 (2.06) |
| B8+ (H2) | 0.042 (1.87) |
| B2 | 0.275 (12.3) |
| B6 | 0.261 (11.6) |

### Labels × gene sets (U, expression lines; tAP_line with lift)

| label | gene set | H1 | H2 | H2 − H1 [CI] | B2 | B8+ (H2) |
|---|---|---|---|---|---|---|
| L_prob | all | 0.300 (2.42) | 0.292 (2.36) | −0.0079 [−0.0083, −0.0075] | 0.900 (7.28) | 0.413 |
| L_prob | noCE | 0.145 (2.05) | 0.142 (2.01) | −0.0031 [−0.0035, −0.0026] | 0.686 (9.70) | 0.207 |
| L_prob | noDepMapCE | 0.028 (1.24) | 0.027 (1.22) | −0.0006 [−0.0008, −0.0004] | 0.275 (12.26) | 0.042 |
| L_ge05 | all | 0.294 (2.50) | 0.285 (2.43) | −0.0082 [−0.0086, −0.0078] | 0.900 (7.67) | 0.404 |
| L_ge05 | noCE | 0.136 (2.10) | 0.133 (2.05) | −0.0032 [−0.0036, −0.0028] | 0.672 (10.36) | 0.192 |
| L_ge05 | noDepMapCE | 0.025 (1.24) | 0.024 (1.21) | −0.0006 [−0.0008, −0.0004] | 0.256 (12.64) | 0.036 |
| L_ge10 | all | 0.203 (3.46) | 0.196 (3.34) | −0.0070 [−0.0073, −0.0067] | 0.902 (15.36) | 0.297 |
| L_ge10 | noCE | 0.089 (2.65) | 0.085 (2.53) | −0.0041 [−0.0044, −0.0037] | 0.718 (21.37) | 0.125 |
| L_ge10 | noDepMapCE (**NA** by §5.2) | 0.008 (1.25) | 0.007 (1.19) | −0.0004 | 0.208 (33.90) | 0.012 |

Other cuts:
- **U∩:** H1 0.308, H2 0.299, H2 − H1 −0.0091 [−0.0095, −0.0087].
- **All 737 eligible-lineage lines**, including those without expression: H1 0.299, H2 0.291, B2 0.898.

### Strata (§5.4; L_prob, all genes, U)

| stratum | lines estimable | NA | H1 | H2 | H2 − H1 [CI] | B2 | E[AP_random] |
|---|---|---|---|---|---|---|---|
| none | 686/686 | no | 0.300 | 0.292 | −0.0079 [−0.0083, −0.0075] | 0.900 | 0.124 |
| expressed pairs | 686/686 | no | 0.346 | 0.338 | −0.0080 [−0.0084, −0.0076] | 0.908 | 0.172 |
| selective genes | 686/686 | no | 0.444 | 0.436 | −0.0078 [−0.0085, −0.0072] | 0.789 | 0.379 |
| addressable (H2) | 686/686 | no | 0.433 | 0.415 | −0.0173 [−0.0179, −0.0167] | 0.921 | 0.186 |
| addressable (H1) | 686/686 | no | 0.437 | 0.424 | −0.0132 [−0.0138, −0.0126] | 0.918 | 0.183 |
| library Avana | 573/573 | no | 0.300 | 0.292 | −0.0082 | 0.906 | 0.124 |
| library Humagne | 84/84 | no | 0.305 | 0.297 | −0.0079 | 0.882 | 0.125 |
| library KY | 28/28 | no | 0.286 | 0.280 | −0.0063 | 0.857 | 0.119 |
| organoid | 55/55 | no | 0.311 | 0.306 | −0.0052 | 0.855 | 0.127 |
| CN ≥ 4 pairs only | 27/686 | **NA** | 0.729 | 0.719 | — | 0.861 | 0.682 |
| CN ≥ 4 pairs excluded | 686/686 | no | 0.300 | 0.292 | −0.0080 | 0.901 | 0.123 |
| multi-copy loci only (§1.6b) | 480/686 | no | 0.128 | 0.127 | −0.0012 [−0.0043, 0.0023] | 0.808 | 0.135 |
| multi-copy excluded | 686/686 | no | 0.305 | 0.298 | −0.0072 | 0.902 | 0.125 |
| housekeeping-adjacent only (§1.6g) | 686/686 | no | 0.767 | 0.774 | **+0.0064** [0.0061, 0.0067] | 0.984 | 0.595 |
| housekeeping excluded | 686/686 | no | 0.155 | 0.150 | −0.0057 | 0.854 | 0.096 |
| moonlighting only (§1.6f) | 604/686 | **NA** (< 30 genes) | 0.547 | 0.547 | 0 | 0.963 | 0.547 |
| pseudo-reaction GPR genes excluded (§1.4d) | 686/686 | no | 0.203 | 0.195 | −0.0075 | 0.879 | 0.109 |

Per lineage (`train_per_lineage.tsv`):
- H2 − H1 is negative in all 12 training lineages, ranging from −0.0056 (Bowel) to −0.0096 (Bone).
- Mean tAP of H2 lies between 0.278 and 0.301 and of B2 between 0.877 and 0.919. Lineages are not ranked by raw AP (§6.2).

### Addressable set and recall ceiling (§5.4)

The addressable set holds the genes whose KO closes ≥ 1 reaction: 1,064 U genes in H1 and 1,098 in H2.

| Label | Positive training pairs | Ceiling, H1 | Ceiling, H2 | Ceiling without CE_fit, H1 | Ceiling without CE_fit, H2 |
|---|---|---|---|---|---|
| L_prob | 235,187 | 0.549 | 0.577 | 0.536 | 0.562 |
| L_ge05 | — | 0.550 | 0.579 | 0.539 | 0.565 |
| L_ge10 | — | 0.563 | 0.600 | 0.478 | 0.510 |

- The ceiling is the share of positives inside the addressable set, so about 42–45% of training positive pairs lie in genes whose KO changes nothing in the parent.
- Addressable-set mAP: H2 0.415, H1 0.437 (strata table).

### §6.8 binary endpoints (training expression lines; parent gene set; lineage-equal mean [95% CI]; line mean)

| protocol | model | call | label | MCC | F1 | precision | recall | specificity | accuracy | called |
|---|---|---|---|---|---|---|---|---|---|---|
| Human1-paper | H1 | P57 | L-ge06 | 0.329 [0.326, 0.331]; 0.330 | 0.384 | 0.431 | 0.347 | 0.951 | 0.893 | 216 |
| Human1-paper | H1 | PGR | L-ge06 | 0.334 [0.332, 0.336]; 0.334 | 0.364 | 0.513 | 0.283 | 0.972 | 0.905 | 147 |
| Human2-paper | H1 | P57 | L-ge06 | 0.329 [0.326, 0.331]; 0.330 | 0.384 | 0.431 | 0.347 | — | 0.893 | 216 |
| fixed cutoff | H1 | FBA < 0.01 | L_prob | 0.322 [0.320, 0.323]; 0.323 | 0.344 | 0.563 | 0.248 | 0.974 | 0.886 | 147 |
| Human1-paper | H2 | P57 | L-ge06 | 0.331 [0.329, 0.334]; 0.333 | 0.387 | 0.423 | 0.358 | 0.950 | 0.894 | 219 |
| Human1-paper | H2 | PGR | L-ge06 | 0.331 [0.329, 0.333]; 0.332 | 0.365 | 0.493 | 0.291 | 0.969 | 0.906 | 152 |
| Human2-paper | H2 | P57 | L-ge06 | 0.331 [0.329, 0.334]; 0.333 | 0.387 | 0.423 | 0.358 | — | 0.894 | 219 |
| fixed cutoff | H2 | FBA < 0.01 | L_prob | 0.318 [0.317, 0.320]; 0.319 | 0.346 | 0.543 | 0.254 | 0.972 | 0.887 | 152 |
| context: not a paper protocol | B2 (fitting-rate threshold) | — | L-ge06 | 0.787 [0.784, 0.789] | | | | | | |
| context: not a paper protocol | B3 (fitting-rate threshold) | — | L-ge06 | 0.233 [0.230, 0.237] | | | | | | |

Notes on this table:
- No line had an undefined MCC for H1 or H2 (n_nan = 0).
- The "called" column counts parent genes that have a CRISPR column. That is why it is lower than the model-wide P57 counts of 226 and 229.
- The 21-threshold grid (−1.0 … 0.0) is in `train_binary_grid.tsv`.
- U-wide sensitivity rows and the all-lines rows are in `train_binary.tsv`.

**Paired H2 − H1 MCC** (CIs in `train_binary_paired.tsv`):

| Call | MCC difference |
|---|---|
| P57 | +0.0028 |
| PGR | −0.0027 |
| Fixed cutoff | −0.0033 |

For P57, H2's higher recall (+0.011) outweighs its lower precision (−0.008).

### Per-gene error table (ext-PI SF1; `train_per_gene.tsv`)

The table gives, for each U gene:
- the generic calls;
- the fraction of training lines positive, PPV if called and miss rate if uncalled;
- whether the gene is in CE_fit;
- the matching panels;
- tAP_gene for H1, H2, B2, B3 and B6.

H2 has 72 FP-majority genes: called lethal but positive in < 50% of lines. 49 of them are in a §1.4 false-positive panel.

H2 has 225 FN-majority genes: uncalled but positive in ≥ 50% of lines.
- 108 of them are in CE_fit.
- 151 are in a §1.2, §1.5 or §1.6 panel.

The majority split is descriptive and *ours*.

### Background groups (M1 leftover; screened lines; counts are members)

| group | train | validation | test | all screened |
|---|---|---|---|---|
| MTAP_loss | 99, evaluable (10 lineages) | 43, evaluable (3) | 117, evaluable (6) | 259 |
| ASS1_silenced | 81, evaluable (7) | 46, evaluable (3) | 34, evaluable (4) | 161 |
| CDS1_null | 209, evaluable (11) | 111, evaluable (3) | 121, evaluable (7) | 441 |
| NRF2_active | 38, evaluable (4) | 5, no | 6, no | 49 |
| STK11_loss | 49, evaluable (3) | 2, no | 8, no | 59 |
| PHGDH_amp | 21, evaluable (4) | 3, no | 6, no | 30 |
| LDHB_null | 21, evaluable (4) | 2, no | 16, evaluable (2) | 39 |
| NAPRT_silenced | 17, no (2 lineages, < 10 members) | 2, no | 11, evaluable (2) | 30 |
| ME2_loss | 10, no | 2, no | 12, evaluable (2) | 24 |
| SLC7A11_null | 17, no | 11, evaluable (2) | 12, no | 40 |
| ASNS_low / SQLE_null / HK1low_HK2high / VHL_loss | not evaluable in any tier | | | 10 / 11 / 33 / 26 |
| ENO1_del, IDH_mut | in-silico only (0 and 8 lines) | | | |
| NEG-1 CDKN2A_del_MTAP_intact | 74, evaluable | 6, no | 43, evaluable | 123 |
| NEG-3 TTN_damaging_nonhyper | 62, evaluable | 15, evaluable | 13, evaluable | 90 |
| NEG-4 SMARCA4_loss | 41, evaluable | 7, no | 16, evaluable | 64 |
| NEG-5 KRAS_hotspot (secondary) | 120, evaluable | 7, no | 82, evaluable | 209 |

NEG-2 (size-matched random within-lineage sets) is available wherever its target group is evaluable.

**Evaluable on test (the M1 list):** MTAP_loss, ASS1_silenced, CDS1_null, LDHB_null, NAPRT_silenced, ME2_loss, plus the controls NEG-1, NEG-3, NEG-4 and NEG-5.

### Gene panels (M1 leftover)

`data/controls/gene_panels.tsv` covers §1.1–§1.6 of the final list: 2,859 rows, 2,279 genes and 97 subpanels. Each row carries:
- section, panel and subpanel;
- symbol and Entrez;
- membership in H1, H2 and U;
- role (trigger/target, modifier, seed/rule) and tag (strong, moderate, model-masked, weaker, uncertain);
- how the row was resolved (token, range, wildcard, `gpr_of:MAR…`, computed rule);
- the source SHA.

Eleven descriptive phrases name no genes and are listed, unresolved, in `gene_panels_unresolved.tsv`:
- "the core complex I and complex III subunits";
- "ETC subunits";
- "OXPHOS partials";
- the pseudogene half of the §1.6b ≥ 95% rule, because no local pseudogene identity table exists.

## Verification (§12, W1 rows) — `verification.tsv` 34/34 and `verification_tasks.tsv` 2/2 pass

**Second code path** (H1 and H2, training, primary cell):
- Per-line tAP was recomputed from scratch for all 686 lines: pyarrow read of the KO parquet, its own Entrez join, its own numpy label pivot, and its own tAP implementation.
- Lineage-equal and pooled tAP were recomputed the same way.
- All agree with the outputs to the files' 6-significant-digit precision (max |diff| ≤ 5e-7), with an identical NA pattern.
- On 5 lines per model, a Monte-Carlo check (2,000 random tie-breaks with sklearn AP) agrees within 4 SE.
- Lineage-equal MCC for H2 P57 × L-ge06 via `sklearn.matthews_corrcoef` equals the output.

**Positive controls:**
- Oracle scorer: tAP = 1 on all 686 lines. Oracle caller: MCC = 1.
- CAD and UMPS are lethal in H2. CAD is lethal in H1, and UMPS is not (known reversible-TYMS issue).
- P57 equals Human-GEM `taskEssentialGenes` on both parents, with identical category sets.

**Negative controls:**
- B0 tAP equals E[AP_random] exactly.
- B1 lift is between 0.996 and 1.003 for all 5 seeds.
- With labels shuffled across genes within each line, H2 lift is 1.002 (real: 2.36).
- Constant caller: MCC NaN under the primary rule and 0 in the sensitivity row.
- A fold-scoped memorising scorer gets exactly E[AP_random] on every held-out fold (max deviation 3e-17), so no held-out label reaches fitting through the harness.
- LAP3 (H2) and TGDS (H1) have ratio_r = 1.
- Orientation: mean prevalence of H2-lethal genes is 0.543, against 0.097 for the rest.

**Toy tests of `evaluate_val.py`** (`tests/test_evaluate_val.py`, 21 pass) and of the provenance wrapper (`tests/test_run.py`, 3 pass):
- label-free preparation, including a failing B6 scorer, fails before any ledger row is written;
- the config JSON must exist before `started`;
- the `started` row is committed before validation labels are read;
- outputs are scalar-only. The hardened guard rejects:
  - duplicate scalar families;
  - more rows than expected;
  - line IDs, lineage names and genes (Entrez or symbol) in cells or column names;
  - any numeric column that is not a whitelisted aggregate;
- a look id is never repeated;
- the cap of 40 is enforced before any label is read;
- a crashed look resumes with a ledger `output` row `{"event": "resumed", "attempt": n}`, written before labels are re-read, and the cap counts it once;
- G\* follows the sign rule, with a tie going to Human2;
- both B6 branches are exercised: excluded by default, and included with a stub scorer;
- pre-flight refuses each of: a missing FROZEN row, a SHA mismatch, an untracked file and a modified file. Git-ignored inputs are checked by SHA;
- provenance: full SHA-256 of inputs and outputs, versions, no row on failure, and a dirty tree refused unless allowed.

The full suite passes: **430 tests**, run on 2026-10-27 without appending to `results/test_runs.tsv`.

## The pending validation look (for the PI)

What the look does, in order (`scripts/evaluate_val.py`):

1. **Pre-flight**, before anything is written. The look is refused unless all of the following hold:
   - every required file has a `file` row in the **committed** `HEAD:FROZEN.tsv`, with a matching SHA-256;
   - every required file is clean in git. Git-ignored inputs under `data/` and `results/**/tmp/` are checked by SHA only;
   - `scripts/check_frozen.py` passes **with** the test-run check.

   The required files are:
   - `scripts/evaluate_val.py` and `src/nextgem/{labels,w1,metrics,run}.py`;
   - both KO parquets and both P57 parquets;
   - `line_features.tsv`, `gene_features.tsv`, `tmp/expr.parquet` and `tmp/cn.parquet`;
   - `data/controls/gene_panels.tsv` and `sesoi_c1b.json`.
2. `nextgem.run.provenance`, which refuses dirty `src/`, `scripts/` or `tests/`.
3. Everything that needs no validation label:
   - the label-free context;
   - training-only fitting;
   - **all validation score matrices**;
   - `val_look_W1_M2_baseline_config.json`, with the input and code SHA-256 and the config SHA.
4. Ledger `started` (committed). A resumed look appends the `resumed` marker row here.
5. Only then:
   - `labels.validation()` and scoring;
   - `val_look_W1_M2_baseline.tsv` (scalars only) and `…_gstar.json`;
   - ledger `output` and `completed`, and a MANIFEST row.

**B6 is not in the W1 look** (`INCLUDE_B6 = False`). It is untuned and provisional, and W7 owns its validation looks.
- DESIGN lists B0–B6 among the M2 deliverables (§14) and says the W1 grid is computed on training and validation (W1 Method), but it never names B6 in the W1 validation look explicitly.
- The flag is part of the config, so it is covered by the config SHA, if the PI decides otherwise.

**G\* rule.** G\* = Human1 iff Δ̄_lin(H2 − H1) < 0, otherwise Human2.
- The contrast is computed on L_prob, all genes, U and the validation expression lines, following DESIGN §6.3 ("the better … on validation Δ̄_lin").
- **A tie goes to Human2**, the default parent (§3.1).
- The rule is documented in `decide_gstar()`.

**To run:**
1. Print the candidate rows: `.venv/bin/python scripts/evaluate_val.py --look W1_M2_baseline --print-frozen-candidates`.
2. Add them to FROZEN.tsv, and commit the code, the inputs and FROZEN.tsv.
3. Check with `--preflight-only`.
4. Run the look once: `.venv/bin/python scripts/evaluate_val.py --look W1_M2_baseline`

## Caveats and things that surprised me

- **H1 beats H2 on training tAP** (−0.008, tight CI, all 12 lineages, all labels). This is the opposite of H1.2 and of the Luo 2026 claim. Under the paper P57 protocol the order flips, but only by +0.003 MCC. The difference is small in absolute terms, and the two models are near-equivalent; the binary MCC is ≈ 0.33 for both.
- **Task-augmented scores lower tAP on all genes but raise it without common essentials.** H2_task is −0.011 vs H2 on all genes and +0.011 without CE_fit. The 75–77 extra P57 genes are, on average, less often essential than the FBA-lethal block they join.
- **The B6 positive control fails at M2 defaults.** The §12 W7 control requires B6 ≥ B2 on training; B6 scores 0.888 against 0.900 for B2. B6 does carry line-specific information (mΔAP_gene +0.022, against −0.033 for B2), but it loses on per-line ranking. B6 is untuned (the grid is W7's), so this is flagged rather than fixed; tuning it here would be selection on training labels outside W7.
- **B2's negative mΔAP_gene (−0.033) is a LOLO artefact.** A lineage with high prevalence gets a lower out-of-fold prior, which anti-correlates the prior with held-out labels across folds. It is not an error, but it matters for C3, which compares P1 against B3 and not B2.
- **The gene prior dominates.** B2 lifts are 7.3 with all genes and 9.7 without CE_fit, against about 2.4 and 2.0 for the generic models. With the DepMap common-essential list removed, generic lift drops to 1.2. Expression-only B3 (0.326) beats both generic models, and B8+ (expression tie-break) adds +0.12 over G-candidates, which is the "free" gain B8 was designed to expose.
- **About 42% of positive pairs are outside the addressable set** (H2), which bounds the recall of any pipeline that keeps the parent's GPR closure.
- **BCa check (§6.5).** The skew of the H2 − H1 line deltas is 0.15, so BCa is not triggered. B2 − H2 has skew −1.68. The decisive check, on P1 − G\*, has to be repeated when P1 exists.
- **Moonlighting stratum.** It is NA because fewer than 30 genes are estimable. H1 and H2 score its genes identically (no difference), because all of them sit in the same tied block in both models.
- **Not done in this M2 pass** (W1 items not in the brief, or cut-list items):
  - the failure-mode traces (H1 TYMS diagnostic, the pFBA trace of H2 β-oxidation, CYP4F2/BBOX1/TMLHE);
  - the v1.3.0 sensitivity (cut-list item 1).
  The training strata by CN, expression, library and CE are in the grid.
- **Process gaps for the PI:**
  - `nextgem.run.provenance()` now exists, in a minimal form (`src/nextgem/run.py`).
    - `MANIFEST.tsv` was regenerated through it. Each script run has one row with: full SHA-256 of inputs and outputs; Python, cobra, optlang and gurobipy versions; solver version; seeds; wall and CPU time; peak RSS.
    - Rows record `dirty=True` because the W1 code is uncommitted; the runs used `NEXTGEM_ALLOW_DIRTY=1`, which is approved for M2.
    - The wrapper does not hash label files. For label-reading scripts, the label files' FROZEN.tsv hashes are recorded in the `extra` field instead.
    - The superseded first MANIFEST is kept in `tmp/MANIFEST_v1_superseded.tsv`.
  - No flock Gurobi slot pool exists.
  - `metrics.py`, `evaluate_val.py`, `w1.py`, `run.py` and the look's inputs are not frozen yet; `--print-frozen-candidates` lists the candidate rows. `metrics.py` was never in FROZEN.tsv, contrary to the earlier version of this README.
  - `ledger.py` is frozen and has no `resumed` event, so a resumed look is logged as an `output` row with a `resumed` detail. A dedicated event would need an amendment.
  - B6's medium inputs must switch to the W2 annotation when it exists.
- **Label exposure of this agent:**
  - training-tier aggregate label statistics (per-line and per-lineage tAP; per-gene training prevalence in `train_per_gene.tsv`);
  - no validation or test labels, and no raw files.
  - Under the PI rule this is not a disqualifying exposure, because it is training-only. The register entry is for the PI to add.

## Review response (notes/reviews/M2_W1_review.md, 2026-10-27)

| Item | Status |
|---|---|
| 1. Label-free steps before `started` | Done. `real_context`, the training fit, every validation score matrix and the config JSON now run before `ledger.look()`. Only `labels.validation()` and the scoring come after it. |
| 2. Pre-flight (committed FROZEN.tsv and clean git) | Done. `preflight()` covers the required code and inputs; `--print-frozen-candidates` prints the missing rows and `--preflight-only` runs the check. Today it correctly refuses the look, because the rows are not yet in FROZEN.tsv. |
| 3. `check_frozen.py` with the test-run check | Done. `--no-test-runs` is removed, and the check currently passes. |
| 4. B6 in the W1 look | Removed: `INCLUDE_B6 = False`, recorded in the config. DESIGN does not name B6 in the W1 validation look explicitly (§14 M2 lists B0–B6 as M2 deliverables; W1 says "full grid on training and validation"). Toy tests cover both branches. |
| 5. Config JSON before `started` | Done: `val_look_<id>_config.json`, which includes `config_sha`. `run_look` refuses to start if the file is missing or stale. |
| 6. Resumed looks logged | Done, via an `output` row with detail `{"event": "resumed", "attempt": n}`, written before the label re-read. `ledger.py` is frozen and has no `resumed` event. The cap counts the look once (tested at 39 + 1 looks). |
| 7. Hardened scalar guard | Done: one row per scalar family; row caps per section; no line ID, lineage or gene in any cell or column name; numeric columns must be on the whitelist. Each check has a test. |
| 8. G\* on the sign of Δ̄_lin(H2 − H1) | Done; a tie goes to Human2 (`decide_gstar`). |
| 9. `nextgem.run.provenance()` | Done (`src/nextgem/run.py`, with 3 tests). `MANIFEST.tsv` was regenerated by a full re-run (`NEXTGEM_ALLOW_DIRTY=1 scripts/w1_run_all.py`, 41 min, Gurobi ≤ 8). |

**No training number changed.** In the re-run, every tracked training output is bitwise identical (`rerun_sha_check.tsv`). The one exception is `verification_tasks.tsv`, whose detail text records elapsed seconds; both checks still pass.

**Current look config SHA:** `331a2c8d1b25c509fb12ac6e1771fb792b9d48f13426f698ebd092b0801a60e5`. It changes if any frozen input or code file changes.
