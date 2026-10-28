# F2: artefact-only biomass leave-out on documented medium (post hoc / exploratory (A15))

**Post hoc / exploratory (A15).** This follow-up comes after the T2 unseal (DESIGN §17 A15; RESEARCH_TREE D17, D18). It never replaces T1/T2 or the headline, and it gets no test evaluation. The F2 agent did not open any T1/T2 output, `results/m6/test/`, `results/f3/` or `data/locked/`.

Author: F2 engineer agent (Claude Opus), 2026-10-28; a non-author of M6, W7 and F1. Status: **complete; validation look run as M10_F2 (look 7/40): F2 +0.0310 [0.0293, 0.0326] vs G*; see results/m10/val/.** Code committed at 17c04d2 before any label read. The training recompute and the decomposition ran on real training labels from the clean pinned worktree of 17c04d2 (MANIFEST dirty=False); the spec was then rewritten with `--pin-train`.

## Label hygiene

- Session-start probe `head -c 50` of the raw CRISPR gene-effect file: **BLOCKED** by the label hook.
- Labels were read only by `f2_score_train.py` and `f2_decompose.py` (committed code, clean tree), through `labels.train()` and `labels.fold(k)`: training tier only. Development dry runs used `--synthetic-labels` (output in `tmp/`, git-ignored).
- `evaluate_val.py --dry-run` (the reviewed harness) read training labels for synthetic positive rates, as designed (A6); it printed nothing label-derived.
- Exposure carried in: D17, D18, the F1 result note, the F1 README (training aggregates), the M6 README (training aggregates, W3-MEV), the W3 README (training aggregates), W7's results.tsv row for W7-03 and the W1 G* file (validation Δ̄(H2 − H1) and tAP_lin of H1/H2, already in RESEARCH_TREE D6).

## F2 definition (D18)

F2 = Human1 + artefact-only B-loo (W7-03: `B-loo:retinoid_catabolites`, `B-loo:vitaminD_catabolites`, `B-loo:tocopherol_catabolites`, `B-loo:carnitine` removed from the cofactor pool MAR10065; ubiquinol, heme, lipoate, cobamide and the real cofactors kept), on the documented medium M2(f=1) of each line (P1's medium; `results/m6/scores/line_assignment.tsv`, asserted equal to `results/w2/controls/line_sweeps.tsv`). Scores: e = 1 − ratio_r of the W7-03 sweep (`results/w7/sweeps/ratio_r_human1.parquet`) for the line's M2(f=1) signature, rounded to 6 decimals; a U gene absent from Human1 scores 0.

## 1. Label-free scores and checks (`scripts/f2_scores.py`; `scores_checks.json`, `mev_check.tsv`)

- **Gene map.** The own ratio → e map reproduces M6's sweep-keyed P1 e table (`results/m6/scores/sweep_e.parquet`) bitwise on all 19 M6 B-loo sweeps.
- **W3-MEV (task b).** HMGCR, MVK, PMVK and MVD are lethal under F2 (ratio_r = 0, max 0.0) on **19/19 serum M2(f=1) signatures** and on the serum-free one: 1,361/1,361 assigned lines. Under P1 the core is lethal on 0/19 serum signatures (min core ratio 0.548–0.989; only the 2 serum-free lines keep it). The re-derived ratios equal `results/w7/sweeps/mev.tsv` exactly.
- **Validation scores (task c).** `results/m10/val/f2_val_scores.parquet`: 198 validation expression lines × 2,866 U genes (567,468 rows), columns `H1`, `H2`, `M6:H1+Bloo/M2(f=1)` (P1) and `F2:H1+artBloo/M2(f=1)` (F2). H1, H2 and P1 equal `results/m6/val/m6_val_scores.parquet` bitwise. F2 differs from P1 on 15,098 pairs; fixed-cutoff calls (e ≥ 0.99): F2 25,836, P1 21,092. SHA-256 `933f7e92836d8fe96edc52772eb9033fa7b9c626fbf6044d5a6c11b5f53321da`.

## 2. Validation look spec (task d, e; `scripts/f2_val_spec.py`)

`results/m10/val/look_spec_M10_F2.json`: kind shortlist, workstream = milestone = "M10" (A13), reference H1 (= G*), candidates [F2], P1 descriptive, L_prob, genesets all/noCE, U, strata none/expressed_pairs, B 2,000, seed 11. Decision rule: the §8 gate vs G* (m = 1). F2 − P1 is read as Δ̄_lin(F2 − H1) − Δ̄_lin(P1 − H1) on the same lines (exact point estimate; the harness gives no CI for it). Provenance pins the label-free inputs, the code and (written with `--pin-train`, which refuses unless every recomputed contrast agrees) `train/summary.json` and `train/contrasts.tsv`.

- `evaluate_val.py --look M10_F2 --print-config` (main tree): exit 0, **config_sha `c665f0c5f0899e92769384319d0bfa5b3aaa149b78e131d0d1a3ca89f5a7acd0`** (the pre-recompute label-free spec had `269b58aa…b28f`).
- `--dry-run` (synthetic validation labels): exit 0, same config_sha.
- The real look was **not** run.

## 3. Training recompute (tasks 1 and 2a; `scripts/f2_score_train.py`; `train/`)

Own code (`scripts/f2_common.py`): tie-aware AP (block form of DESIGN §6.1, unit-tested against brute force and `nextgem.metrics.tap`), lineage-equal Δ̄_lin, lineage-stratified bootstrap (B = 2,000, SEEDS[0]), TG fold rule. 686 non-rare training lines with expression, 12 LOLO lineages, L_prob, all U genes; noCE uses CE_fit from `labels.fold(k)` of the line's lineage. Contrasts fixed in the code: F1 − P1, F2 − G*, F2 − P1, P1 − G* (control). Verification built in: P1's and G*'s training scores equal M6's bitwise (passes label-free), and with real labels their per-line tAP must match `results/m6/train/per_line.tsv` (≤ 1e-6). Agreement rule: estimate equal at the reported precision, same folds count, CI ends within 0.0002.

Verification passed: P1 and G* scores equal M6 bitwise; their per-line tAP equals `results/m6/train/per_line.tsv` (max 5.0e-7, the file's rounding), all and noCE. **All four contrasts agree.**

| Contrast | Reported | Recomputed Δ̄_lin [95% CI], folds | noCE |
|---|---|---|---|
| F1 (curated P1) − P1 | −0.00063 [−0.00071, −0.00055], 0/12 | −0.00063 [−0.00071, −0.00055], 0/12 | +0.00206, 12/12 |
| F2 − G* | +0.0304 [0.0294, 0.0313], 12/12 | +0.0304 [0.0294, 0.0312], 12/12 | +0.0312, 12/12 |
| F2 − P1 | −0.0230 [−0.0236, −0.0224], 0/12 | −0.0230 [−0.0237, −0.0224], 0/12 | −0.0098, 0/12 |
| P1 − G* (control) | +0.0534 [0.0523, 0.0545], 12/12 | +0.0534 [0.0522, 0.0545], 12/12 | +0.0410, 12/12 |

## 4. Gene-group decomposition (`scripts/f2_decompose.py`; `decomposition.tsv`)

Groups (label-free, disjoint, in this priority): (a) mevalonate core (4); (b) W3 sterol panel, Human1 (12: CYP51A1, DHCR24, DHCR7, EBP, FDFT1, HSD17B7, LSS, MSMO1, NSDHL, SC5D, SQLE, TM7SF2); (c) the 35 W3 cofactor-pool B-loo calls (35); (d) W2 flip lists, Human1 (`flips/flip_genes.tsv` ∪ `flips/cached_doc_flips.tsv`), minus (a)–(c) (13); (e) the other 2,802 U genes. For P1 and F2, each group's contribution is Δ̄_lin(hybrid − G*), where the hybrid has G*'s scores except on that group. The sum of the five is reported next to the total, with the residual. Output: `decomposition.tsv` (aggregates only), `decomposition.json` (group lists).

Training, L_prob, all U genes, 686 lines, 12 lineages; contribution = Δ̄_lin(hybrid − G*) [95% CI], folds > 0.

| Group | P1 | F2 | F2 − P1 |
|---|---|---|---|
| (a) mevalonate core (4) | **−0.0014 [−0.0016, −0.0013], 1/12** | 0 (unchanged vs G*) | +0.0014 |
| (b) sterol panel (12) | +0.0072 [0.0071, 0.0073], 12/12 | +0.0072 (identical) | 0 |
| (c) cofactor-pool calls (35) | +0.0331 [0.0325, 0.0337], 12/12 | +0.0127 [0.0125, 0.0129], 12/12 | −0.0204 |
| (d) medium flips (13) | +0.0061 [0.0054, 0.0068], 11/12 | +0.0061 (identical) | 0 |
| (e) other (2,802) | +0.0029, 12/12 | +0.0026, 12/12 | −0.0003 |
| sum of groups | +0.0479 | +0.0286 | −0.0193 |
| total Δ̄_lin vs G* | +0.0534 | +0.0304 | −0.0230 |
| residual (total − sum) | +0.0055 (10%) | +0.0018 (6%) | −0.0037 |

Interpretation: P1's W3-MEV artefact does **not** contribute to its gain. Making the mevalonate core non-lethal costs P1 −0.0014 on its own (1/12 folds positive), so P1's training gain would be slightly larger without it, and F2 recovers this. Most of P1's gain (+0.033, 62%) comes from un-calling the 35 cofactor-pool genes. F2 keeps only +0.013 of that, because it keeps 20 of the 35 calls lethal on the main RPMI-family medium (label-free check). Those 20 are the real-cofactor genes: heme (ALAD, CPOX, HMBS, UROD, UROS, ABCB6), CoQ (COQ2/3/6/7), BH4 (GCH1, PTS, SPR), PLP (PDXK, PNPO), TPK1, RFK, MMAA, SLC5A6 and SLC25A26. This group alone accounts for −0.020 of the −0.023 F2 − P1 gap. The sterol and medium-flip gains are identical in P1 and F2. The positive residual (super-additivity, as in M6's interaction) is smaller for F2.

## Compute

Label-free score build 8 s; synthetic dry runs 6–9 s each; real training recompute and decomposition about 1 min CPU together; harness dry runs 8 s each. No LP was solved (all sweeps were reused). Total under 0.05 CPU-h.

## Files

| File | Content |
|---|---|
| `scores_checks.json` | gene-map check, MEV summary, validation-score checks and SHA |
| `mev_check.tsv` | per M2(f=1) signature: lines, F2 and P1 mevalonate-core status |
| `train/summary.json`, `train/contrasts.tsv` | training recompute (aggregates) |
| `decomposition.tsv`, `decomposition.json` | gene-group decomposition (aggregates) and group lists |
| `MANIFEST.tsv` | provenance rows (development runs dirty; the real runs and the spec clean at 17c04d2) |
| `../m10/val/f2_val_scores.parquet` | validation scores for the look |
| `../m10/val/look_spec_M10_F2.json` | the look spec |
