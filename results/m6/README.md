# M6: integration factorial toward P1 (training; validation look M6_factorial run as look 3/40, see val/val_look_M6_factorial_decision.json; P1 selected, RESEARCH_TREE D11)

**Author:** M6 engineer agent (Claude Opus). **Date:** 2026-10-28.
**Design:** DESIGN.md §10.1, §10.2, §6.3–§6.8 and §17 (A10, A11, A12c). Cells were fixed in RESEARCH_TREE D10 before any factorial scoring.
**Status:**
- Nothing here is committed by the author.
- DESIGN.md, FROZEN.tsv, the ledgers, `scripts/evaluate_val.py` and `scripts/evaluate_test.py` were not edited.
- No validation look was run.

**Label hygiene:**
- The session-start probe (`head -c 50` of the raw CRISPR gene-effect file) was **blocked** by the hook.
- Labels were read only by `scripts/m6_score_train.py`, through `labels.train()` and `labels.fold(k)` via `cv.lolo`. The tripwire log shows 50 `train` opens and no `val`, `test` or `data/locked/` access.
- Test-line scores are model outputs, not labels (see "Scores" below).

## Cells (D10)
| Cell | Parent | Biomass | Medium | Role |
|---|---|---|---|---|
| G* | H1 | as released | M0 (Ham's) | reference (= W1 generic H1, asserted) |
| (1) | H1 | W3 B-loo:cofactor_pool_biomass | M0 | candidate (flag W3-D1) |
| (2) | H1 | as released | M2(f=1) | candidate (flags W2-SHUF, W2-MOPEN) |
| (3) | H1 | B-loo | M2(f=1) | candidate (W3-D1, W2-SHUF, W2-MOPEN, **W3-MEV**) |
| H2, H2 (1)–(3) | H2 | | | descriptive |

**Medium M2(f=1)** means:
- the W2 documented base medium, plus supplements, plus the serum model at f = 1, gate-free;
- assigned per line from `results/w2/controls/line_sweeps.tsv`;
- all 686 training lines are scored. The 220 lines outside the attribution set get W2's integrated-pipeline medium: imputed, fallback and unmappable lines get classic RPMI + 10% FBS, and onboarded lines get their onboarded medium.
- There is no context axis, because the context lever was not selected at M5 (A12b).

## 1. Label-free KO sweeps (`sweeps/`; `scripts/m6_sweeps.py`)
- **Distinct media.** M2(f=1) has 20 distinct signatures per parent: 19 serum media and 1 serum-free (DMEM:F12 serum-free, 2 lines).
- **Reused sweeps.**
  - Plain parents on M0 and on all 20 M2(f=1) signatures: the W2 cache.
  - B-loo on M0, and on the classic-RPMI signature (= W3's RPMI proxy): the W3 cache.
- **New sweeps:** 38 (19 per parent): B-loo on each remaining M2(f=1) signature.
  - Bounds were rebuilt with `w2_engine.build_sweeps` and matched to the W2 sweep IDs by signature; every ID was found.
  - The variant was applied with `nextgem.w3.Variant` inside `w3.ko_sweep`.
- **Solver and A10 limits.** GLPK, 8 workers. Each LP had a 300 s limit plus one retry, and each sweep job a 1,800 s CPU limit (`ITIMER_PROF`).
  - **42/42 jobs were ok: 0 failed builds and 0 timeouts.**
- **Spot checks (`spot_checks.tsv`).** Four reused sweeps were recomputed and all are **bitwise equal** (max |Δ| = 0):
  - H1 and H2 plain on DMEM M2(f=1), against the W2 cache (solved on Gurobi there);
  - H1 B-loo on classic RPMI, against W3;
  - H2 B-loo on M0, against W3.
- **Cost:** 2,265 CPU-s (0.63 CPU-h), 6.1 min wall. With the 2-sweep pilot (150 CPU-s), the total is 0.67 CPU-h. This is under 1 CPU-h, so the sweeps ran from the working tree (`dirty=True` in MANIFEST).

## 2. Scores (`scores/`, `val/`; `scripts/m6_scores.py`; label-free)
- **`sweep_e.parquet`:** e = 1 − ratio_r for all 84 cell sweeps, keyed by sweep (not by line). It uses the W1 Entrez convention: an absent gene gets 0, rounded to 6 decimals.
- **`line_assignment.tsv`:** ModelID → (parent, step) → W2 sweep ID, medium key and stratum, for all 1,361 lines. It has no tier column.
- **`scores_{train,validation,test}.large.parquet`:** line_gene layout (ModelID, entrez, 8 method columns) for all 767 / 200 / 394 lines of each tier.
  - These files are git-ignored (`*.large.*`) and regenerable by `m6_scores.py`. Their SHA-256s are in **`SHA256SUMS.tsv`**, which is tracked.
  - The test file is kept untracked because the pre-commit test-keyed rule would refuse a test-ModelID × gene table. Together, `sweep_e.parquet` and `line_assignment.tsv` determine every test score.
- **`val/m6_val_scores.parquet`:** the look's file, covering the 198 validation expression lines. SHA-256 `888d6ca268717d65cd0b629f2ff20ce58568d0fb23e34ab3ce1083b3adacdc65`.
- **Asserted:**
  - T0|M0 equals `results/baseline` bitwise;
  - the H1 and H2 columns equal the W1 generic scores in every tier;
  - every split line has a sweep for every cell.

## 3. Training scoring (`train/`; `scripts/m6_score_train.py`)
**Set-up.** 686 lines in 12 lineages; L_prob over all U genes; paired lineage-stratified line bootstrap, B = 2,000, SEEDS[0]; CE_fit is fold-scoped (LOLO).

**Machinery reused.** `w1.method_summary`, `w1.delta_summary`, `w1.ce_mask`, `w1.apply_stratum` (expressed pairs) and `w2_score.folds_positive` (the TG fold rule).

**Verification (`verification.tsv`; all pass).** The G* tAP_line equals W1's `train_per_line` for all genes and for noCE (≤ 5e-7, the TSV's 6-digit rounding).

| Cell | tAP_lin [95% CI] | noCE | expressed pairs | Δ̄_lin vs G* [CI] | folds + | TG |
|---|---|---|---|---|---|---|
| G* H1/M0 | 0.2996 [0.2978, 0.3013] | 0.1451 | 0.3464 | — | — | — |
| (1) H1+B-loo/M0 | 0.3338 | 0.1665 | 0.3769 | **+0.0343** [0.0336, 0.0349] | 12/12 | pass |
| (2) H1/M2(f=1) | 0.3156 | 0.1634 | 0.3618 | **+0.0160** [0.0152, 0.0168] | 12/12 | pass |
| (3) H1+B-loo/M2(f=1) | 0.3529 | 0.1861 | 0.3963 | **+0.0534** [0.0523, 0.0545] | 12/12 | pass |
| H2/M0 | 0.2916 | 0.1421 | 0.3385 | −0.0079 | 0/12 | — |
| H2+B-loo/M0 | 0.3285 | 0.1649 | 0.3711 | +0.0289 | 12/12 | (descriptive) |
| H2/M2(f=1) | 0.3047 | 0.1567 | 0.3520 | +0.0051 | 11/12 | (descriptive) |
| H2+B-loo/M2(f=1) | 0.3476 | 0.1829 | 0.3914 | +0.0480 | 12/12 | (descriptive) |

**Δ̄_lin vs G* without CE_fit and on expressed pairs:**
- (1) +0.0214 / +0.0305;
- (2) +0.0183 / +0.0154;
- (3) +0.0410 / +0.0499.

**mΔAP_gene:** 0 for the M0 cells (line-constant); −0.0004 to −0.0005 for the M2 cells. The medium adds no per-gene ranking across lines, as W2 found.

**Pairwise contrasts (H1, all genes; `contrasts.tsv`):**
- (2) − (1) = −0.0182 [−0.0192, −0.0171], 0/12 folds;
- (3) − (1) = +0.0191 [0.0181, 0.0201], 12/12 folds;
- (3) − (2) = +0.0373 [0.0366, 0.0380], 12/12 folds;
- noCE: (3) − (1) +0.0196; (3) − (2) +0.0227.

**Training one-SE rule (descriptive; selection happens on validation).** All three candidates pass TG. The best is (3), and neither (1) nor (2) is within one paired SE of it (differences 0.019 and 0.037; SE 0.0005 and 0.0004). The training pick would therefore be (3).

### 2×2 interaction, B-loo × medium (`interaction.tsv`)
The interaction is defined as I = [(3) − (2)] − [(1) − G*] per line, then averaged with lineage-equal weights.

| Parent | Gene set | I [95% CI] | Folds I > 0 | Δ B-loo / Δ M2 / Δ joint (additive prediction) |
|---|---|---|---|---|
| H1 | all | **+0.0031** [0.0027, 0.0034] | 12/12 | 0.0343 / 0.0160 / 0.0534 (0.0503) |
| H1 | noCE | +0.0014 [0.0008, 0.0019] | 9/12 (3 negative) | 0.0214 / 0.0183 / 0.0410 (0.0396) |
| H1 | expressed pairs | +0.0040 [0.0037, 0.0044] | 12/12 | |
| H2 | all | +0.0060 [0.0056, 0.0063] | 12/12 | 0.0369 / 0.0131 / 0.0560 (0.0500) |

The two levers are **slightly super-additive**, not sub-additive as the W2 README expected. The W2 overlap (serum carnitine rescuing some cofactor-pool artefacts) is outweighed. One likely cause is that W3-MEV and the B-loo sterol relief make further sterol/isoprenoid genes non-lethal under serum. This is not decomposed here.

### §6.8 binary endpoints (`binary.tsv`, `binary_paired.tsv`)
**Set-up.** L-ge06; the parent's genes with a CRISPR column; 686 lines; lineage-equal mean; 0 NaN lines.

**Call rules:**
- **P57 / PGR (paper protocol).** Each task uses its own inputs, so GR = growth on Ham's under the cell's biomass.
  - The 56 non-GR tasks close every boundary and never use the biomass. Their failure sets are the W3 capacity-matrix sets for every cell.
  - Asserted: P57(T0) and PGR(T0) equal W1's calls, and GR from the T0|M0 sweep (KO growth < 1, the GR task's LB) equals the W3 GR task column. The G* MCC equals W1's `train_binary` (≤ 4e-7).
  - Consequence: under these rules the medium does not enter, so (2) ≡ G* and (3) ≡ (1) by construction. P57 per cell is therefore cheap: no new task runs are needed.
- **Cell-medium sensitivity (ours).** The GR task is evaluated on the line's own medium.
- **Fixed cutoff.** ratio_r < 0.01.

| Cell | MCC P57 | MCC PGR | MCC P57 cell-medium | MCC FBA<0.01 (L-ge06) | F1 P57 | F1 FBA<0.01 |
|---|---|---|---|---|---|---|
| G* H1 | 0.3286 | 0.3339 | 0.3286 | 0.3339 | 0.3836 | 0.3636 |
| (1) | 0.3513 (+0.0227 [0.0223, 0.0230]) | 0.3788 (+0.0449) | 0.3513 | 0.3788 (+0.0449) | 0.3962 (+0.0126) | 0.3786 (+0.0150) |
| (2) | 0.3286 (0 by construction) | 0.3339 (0) | 0.3333 (+0.0047) | 0.3463 (+0.0124) | 0.3836 | 0.3719 (+0.0083) |
| (3) | 0.3513 (+0.0227) | 0.3788 (+0.0449) | 0.3543 (+0.0257) | **0.3916 (+0.0577 [0.0563, 0.0590])** | 0.3962 | 0.3842 (+0.0206) |
| generic H2 | 0.3315 | 0.3312 | 0.3315 | 0.3312 | 0.3868 | 0.3650 |
| H2 (1) / (2) / (3) | 0.3536 / 0.3315 / 0.3536 | 0.3816 / 0.3312 / 0.3816 | 0.3536 / 0.3360 / 0.3567 | 0.3816 / 0.3456 / 0.3995 | 0.3998 / 0.3868 / 0.3998 | 0.3842 / 0.3749 / 0.3925 |

Precision, recall, accuracy, CIs and the L_prob rows are in the TSVs.

## 4. W3-MEV diagnostic (`scores/w3_mev.tsv`; label-free)
**Cell (3) loses mevalonate-core essentiality in every serum medium, in both parents.** In all 19 of 19 serum M2(f=1) signatures, HMGCR, MVK, PMVK and MVD are non-lethal:
- KO ratio 0.548–0.989 in H1;
- 0.505–0.989 in H2.

This covers 765/767 training, 200/200 validation and 394/394 test lines.

The core stays lethal:
- in cell (3) on the one serum-free medium (2 training lines);
- in cell (2) on every medium;
- in cell (1) and G* on M0.

So **flag W3-MEV applies to cell (3) on effectively every line.** It is a positive-control failure: the mevalonate core is a known essential pathway. If (3) is selected, P1 carries it.

## 5. Validation look (`val/look_spec_M6_factorial.json`; `scripts/m6_val_spec.py`): run as look 3/40 (D11). Candidate order in the run follows DESIGN §9 (medium first), as corrected before scoring; text below that lists the biomass cell first predates the correction
**Spec contents:**
- kind "shortlist", layout line_gene;
- methods: H1, H2, `M6:H1+Bloo/M0`, `M6:H1/M2(f=1)`, `M6:H1+Bloo/M2(f=1)` and the three H2 cells (descriptive);
- candidates (1), (2), (3) in D10 order; reference H1;
- L_prob; all and noCE gene sets; U; strata none and expressed_pairs; B = 2,000; seed 11;
- the §8 gate, then the one-SE rule;
- 32 provenance pins: the scores' inputs, the W2/W3/M5 inputs, the M6 code and the training outputs;
- the flags in `notes`.

**Harness fit: the look is refused.**
- `--print-config` and `--dry-run` both exit 2 with "workstream must be W<n> and match the spec's directory results/m6/val (DESIGN §17 A11g; look caps)". The harness requires `workstream` to match `W\d+` *and* equal the spec directory's name. "M6" in `results/m6/val` fails the first test.
- **The harness was not changed.** The options for the PI are:
  - (a) A harness amendment allowing `M<n>` with its own look cap. This is a harness change, so A6(2) requires a full review.
  - (b) File the spec as a W2 or W3 look (`results/w2/val`, workstream W2), since D10 makes this look W2's validation answer. That counts against that workstream's cap.
- **Provisional in-process simulation.** I ran the harness in-process with only the `W\d+` pattern widened to `[WM]\d+` (`results/m6/tmp/sim_dry_run.json`; the harness file is unchanged).
  - `--print-config` gives config_sha **`40e4d6b6ebb79bb3bbcf3df2d1a299bf1ad971d5aa2f1f76fc8cfd7e34d8be53`**.
  - `--dry-run` on synthetic labels passes: 198 lines, 8 methods, 32 scalar rows, scalar guard pass, no validation label read, 0 ledger rows.
  - The real config SHA will change with any harness change, because the config embeds the harness code SHA, and with any re-pin after the commit.
- Toy test `tests/test_m6_factorial.py::test_look_spec_structure_and_workstream_rule` shows the spec passes every static `check_spec` rule except the workstream rule.

## Deviations and notes for the PI
1. **Harness workstream rule** (above). The spec is written as requested but cannot be used as is.
2. **Test-tier scores are untracked** (`*.large.*`), with SHA-256s in `SHA256SUMS.tsv`. The pre-commit test-keyed rule would refuse a tracked test-ModelID × gene table. The sweep-keyed `sweep_e.parquet` plus `line_assignment.tsv` (both tracked) regenerate them.
3. **P57/PGR under the literal paper protocol are medium-independent.** Cell (2) equals G* and cell (3) equals cell (1) by construction. A cell-medium GR sensitivity is reported alongside, and is ours.
4. **A10 per-job CPU limit is cooperative between LPs** (`ITIMER_PROF`): an LP in progress finishes first, bounded by 2 × 300 s. No job came near the limit; the maximum was 56 s.
5. **All 686 lines are scored, not only the 466-line attribution stratum** that W2 used for its decision. The cell (2) Δ of +0.0160 is consistent with W2's attribution-stratum +0.0159 and its all-line +0.0160.
6. **Code ran from the dirty tree** (`allow_dirty`; every MANIFEST row has `dirty=True`). Every run was under 1 CPU-h, and the largest was the sweeps at 0.63 CPU-h.
   - Before the look, the code must be committed and the spec and scores regenerated from the clean commit. That requires re-running `m6_scores.py`, `m6_score_train.py` and `m6_val_spec.py`, about 0.1 CPU-h. The sweeps can be re-run too for a clean MANIFEST row (0.6 CPU-h, under 1 CPU-h).
7. **One unrecorded run.** The first real `m6_score_train.py` run aborted on its own verification: the G* MCC differed from W1's TSV by 4e-7, which is TSV rounding, and the tolerance had been 1e-9. After setting the tolerance to 1e-6 (the same as the tAP check), it was re-run. The aborted run wrote no MANIFEST row (≈ 245 CPU-s); its outputs were overwritten by identical code paths.
8. **Label-copy scanner** (`--roots results/m6 <session scratchpad>`): **0 hits in results/m6.** There are 115 name-only hits, all toy pytest fixtures under `scratchpad/pt*/` in this session's shared scratchpad. They were written by another process (the harness/test agents) and were not created by M6, so they were not deleted.

## Runtime (MANIFEST.tsv)
| Step | Wall | CPU | Workers / solver |
|---|---|---|---|
| `m6_sweeps.py --pilot 2` | 92 s | 150 s | 2, GLPK |
| `m6_sweeps.py` (38 new + 4 spot checks) | 364 s | 2,265 s | 8, GLPK |
| `m6_scores.py` | 9 s | 8 s | — |
| `m6_score_train.py --synthetic-labels` (dry run) | 66 s | 202 s | — |
| `m6_score_train.py` | 88 s | 245 s | — |
| `m6_val_spec.py` | 2 s | 7 s | — |
| aborted scoring run (no row) | ≈ 85 s | ≈ 245 s | — |
| tests (`tests/test_m6_factorial.py`, 16; with the harness suites, 94 pass) | 106 s | — | — |

The total is ≈ 3,120 CPU-s, or **≈ 0.87 CPU-h**.
