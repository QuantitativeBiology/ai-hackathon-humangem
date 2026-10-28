# M5 — A12(b) panel pre-check: context families on H1 + W3 B-loo (training tier)

**Author:** M5 engineer agent. **Date:** 2026-10-28. **Design:** DESIGN.md §17 amendment A12(b), §9 W4 (Stage A), A4, A10, A11.
**Code:** commit 4729547 (`scripts/m5_build.py`, `scripts/m5_score_train.py`, `tests/test_m5_precheck.py`, 10 tests pass).
Both production steps ran from the pinned clean worktree `/data/tmpA/nextgem_runs/4729547`; both MANIFEST rows have dirty = False.

## Verdict (pre-registered rule A12 b)
A family enters Stage B only if its combined Δ̄_lin vs H1 + B-loo is ≥ +0.005 on the 30-line panel.

| Family (on H1 + B-loo) | Δ̄_lin vs H1 + B-loo [95% CI] | Folds positive | Build gates | Enters Stage B |
|---|---|---|---|---|
| F2 GIMME R1 | **+0.0031** [+0.0017, +0.0045] | 8/10 | 30/30 pass | **no** (< 0.005) |
| F4 E-Flux R1 | **+0.0044** [+0.0035, +0.0052] | 10/10 | 8/30 pass (GR units artefact) | **no** (< 0.005) |

**No family enters Stage B.** Under A12(b), the context lever is not selected: the M6 factorial drops its context axis, and Stage B and Stage C are skipped.

The F4 point estimate is below the bar, although its upper CI bound (+0.0052) is above it. F4 also fails the Stage-A viability gate again, for the same GR-units reason as in Stage A, so it would stay descriptive even if it were above the bar.

## What was built (`builds/`, label-free; `scripts/m5_build.py`)
- **Parent:** Human1 with the W3 variant B-loo:cofactor_pool_biomass. The variant (MAM10012c removed from MAR13082) comes from `results/w3/structure/variants_human1.json` and is applied permanently with `nextgem.w3.Variant`.
  - The medium is M0 (Ham's, from Human1's GR task).
  - Parent wild-type growth is 65.5774.
- **Lines and configurations:** the 30 Stage-A panel lines (`results/w4/stage_a_panel.tsv`, 10 training lineages × 3 lines). The configurations are exactly Stage A's (`train_selection.json`):
  - F2 `R1|rna|real|none`, the configuration carried to M6;
  - F4 `R1|rna|real|none`, F4's full-panel best;
  - genomic constraint: none, as in both carried configurations;
  - the A4 KO gene-state convention.
- **Execution:**
  - Same W4 builders (`context.make_parent`, `context_run.worker_run`, `context_run.Runner`), with only the parent swapped.
  - GLPK, the declared solver for F2 and F4 (A4 c). No Gurobi session was opened.
  - 8 workers, one fresh process per build.
  - A10 limits: 300 s per LP with one retry, CPU kill at 1,800 s per job.
- **Difference from Stage A:** Stage A's parent was Human2. A12(b) specifies the H1 + B-loo parent. R1 is a global threshold (TPM 1), so no Human2-fitted statistic enters.
- **Result:** 61 builds (1 generic + 30 F2 + 30 F4): 0 failed, 0 timeouts, 0 non-optimal KO solves.

**Generic check:** the H1 + B-loo KO sweep built here equals W3's validated sweep row `M0|B-loo:cofactor_pool_biomass` for all 2,887 genes (max |Δ| 0; 114 lethal genes in both).

**Repair rule with the modified biomass (F2):**
- It works. The growth LP was optimal in 30/30 builds and all 300/300 viability-task LPs (GR + 9 ER) solved.
- **0 reactions were re-opened in all 30 builds.** GIMME keeps biomass at 0.9 × max, as in Stage A, where F2 repair was ≤ 0.1%.
- Wild-type growth is 0.9985–1.000 of the parent. The KO sweep calls 118–128 lethal genes (median 120.5).
- Tasks passed: 38–57 of 57 (generic parent: 57). All gates pass (growth, viability, repair ≤ 5%).

**F4 (not repaired, by design):**
- Wild-type growth is 0.010–0.021 of the parent (0.66–1.39 in scaled units).
- The GR task (biomass ≥ 1) fails in exactly the 22/30 builds with growth < 1; the other 8 pass all 10 viability tasks. This is the Stage-A units artefact.
- The lethal set equals the generic one in every build (114 genes), so F4's signal is in its partial ratios only.

## Scoring (`train/`; `scripts/m5_score_train.py`)
- **Machinery:** Stage A's scoring (`w4_stage_a_train.load`, `Scorer`, `delta_stats`; `w4_analysis.nested_lolo_select`, `tg_panel`).
  - Scores are e = 1 − ratio_r over U, with Human1 genes mapped by the W1 generic-table Entrez convention.
  - Primary cell: L_prob, all genes.
  - Paired line bootstrap, lineage-stratified, B = 2,000, SEEDS[0].
- **Nested panel LOLO:** one fixed configuration per family, so the inner choice is that configuration in every fold, and the outer Δ equals the plain Δ (asserted).
- **Labels:** `labels.train()` (panel lines) and `cv.lolo` (CE_fit per held-out panel lineage) only.

| Contrast | Δ̄_lin [CI] | Folds + | Minus CE_fit | Expressed pairs | ΔMCC (L_prob) |
|---|---|---|---|---|---|
| (2) F2 + B-loo vs (1) H1 + B-loo | +0.0031 [+0.0017, +0.0045] | 8/10 | +0.0011 | +0.0013 | −0.0039 |
| (3) F4 + B-loo vs (1) H1 + B-loo | +0.0044 [+0.0035, +0.0052] | 10/10 | +0.0008 | +0.0037 | +0.0000 |
| (1) H1 + B-loo vs generic H1 (sanity) | **+0.0332** [+0.0310, +0.0354] | 10/10 | +0.0203 | +0.0286 | +0.0390 |
| F2 + B-loo vs generic H1 | +0.0363 [+0.0334, +0.0391] | 10/10 | +0.0213 | +0.0299 | +0.0351 |
| F4 + B-loo vs generic H1 | +0.0376 [+0.0353, +0.0401] | 10/10 | +0.0211 | +0.0323 | +0.0390 |

- **Sanity check:** the B-loo gain on the panel (+0.0332) matches the W3 training gain (+0.0343, 12 folds) and its validation look (+0.032).
- **Per-fold means:** `train_folds.tsv`. F2 is negative in Cervix (−0.0024) and Peripheral Nervous System (−0.0021). F4 is positive in all 10 lineages (+0.0004 to +0.0074).
- **Per-configuration lineage-equal tAP_line** (`train_config_scores.tsv`):
  - generic H1: 0.2982
  - H1 + B-loo: 0.3314
  - F2 + B-loo: 0.3345
  - F4 + B-loo: 0.3358
- **mΔAP_gene:** F2 +0.0019, F4 −0.0015.
- **Context:** the context families add about +0.003 to +0.004 on top of B-loo. This is similar to their Stage-A effect vs their own parent (F2 −0.0001, F4 +0.0038 vs H2) and below the pre-registered bar.
- **Verification** (`train_verification.tsv`; all pass):
  - Human1 Entrez map reproduces W1's e_H1 exactly;
  - generic H1 tAP_line equals W1's `train_per_line` (max |Δ| 4.9e-7);
  - H1 + B-loo tAP_line equals W3's per-line values for this variant (max |Δ| 5.0e-7);
  - the KO sweep equals W3's (0).

## Runtime (MANIFEST.tsv)
| Step | Wall | CPU |
|---|---|---|
| 2-line cost pilot (dirty tree, `tmp/pilot`, not committed) | 174 s | 0.21 CPU-h |
| `m5_build.py` (8 GLPK workers) | 1,469 s | 2.96 CPU-h in total (2.35 CPU-h inside builds; F2 112 s and F4 167 s mean per build) |
| `m5_score_train.py` | 35 s | 0.01 CPU-h |

Peak RSS is 1.42 GB per worker.

## Label hygiene and exposure
- The session-start probe `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **blocked**.
- Labels were read only by `m5_score_train.py` through `w4_stage_a_train.load`: `labels.train()` plus `cv.lolo` CE_fit, training tier and panel lines only. There was no validation, test or `data/locked/` access.
- **Training-tier aggregates viewed by the author:**
  - the contrast table above (Δ̄_lin, CIs, folds positive, minus-CE_fit / expressed-pairs / MCC deltas);
  - the per-configuration lineage-equal tAP_line, mΔAP_gene and MCC;
  - the per-fold (lineage) mean Δ for the 5 contrasts;
  - the verification max |Δ|.
- One earlier accidental view of a single row of `results/w3/train/per_line.tsv` (recorded by the coordinator). `train/train_per_line.tsv` was written by the script but not opened.
- The label-copy scanner reports 0 hits over `results/m5`.
