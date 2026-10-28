# M4 W3 review: objectives, tasks, biomass (independent non-author review)

Reviewer: independent agent (Claude Opus), 2026-10-28. Scope: results/w3 (uncommitted), src/nextgem/w3.py, scripts/w3_*.py, tests/test_w3.py, against DESIGN §9 W3, the TG (§9: ≥ 9/12 folds), §17 A2/A6 and apriori_lists_final §1.4/§5. Nothing committed; no author file, DESIGN.md, FROZEN.tsv or ledger edited.

Outputs in results/w3/review/: `rv_ko.py`, `rv_score.py`, `rv_ratio_r_<parent>.parquet` (label-free), `rv_ko_meta.tsv`, `rv_scores.tsv`, `rv_random_removal.tsv`, `rv_changed_calls.tsv` (aggregates and gene lists only).

**Probe:** `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **BLOCKED** by the PreToolUse label guard.
**Compute:** ≈ 0.6 CPU-h, GLPK, ≤ 6 processes.

## 1. Independent recomputation (check 1)
**KO sweeps** (`rv_ko.py`). Own path:
- cobra `Gene.knock_out()` plus a GLPK LP for **every** gene. No flux-support shortcut, no w3/w2_engine import.
- M0 = `hams_medium()`, asserted equal to the author's M0. The variant is the `remove` spec in `variants_<parent>.json`.
- Result: ratio_r **bitwise identical** to `sweeps/ratio_r_<parent>.parquet` for T0 and both variants in both parents (max |Δ| = 0). Lethal counts: H1 149 / 114 / 139, H2 154 / 115 / 145 (T0 / cofactor / metabolite).
- The derived e vectors equal `val/w3_val_scores.parquet` exactly (all 6 columns).

**Scoring** (`rv_score.py`). labels.train() only; metrics.per_line / bootstrap_mean; no w1/w3 scoring code; 686 lines, 12 folds, B = 2,000, SEEDS[0]. Calls are counted as removed / added.

| parent | mode | tAP_lin | Δ̄_lin vs T0 [CI] | folds | removed / added |
|---|---|---|---|---|---|
| H1 | T0 | 0.2996 [0.2978, 0.3013] | — | — | — |
| H1 | B-loo:cofactor_pool | 0.3338 | +0.0343 [0.0336, 0.0349] | 12/12 | 35 / 0 |
| H1 | B-loo:metabolite_pool | 0.3051 | +0.0056 [0.0051, 0.0059] | 12/12 | 10 / 0 |
| H2 | T0 | 0.2916 [0.2898, 0.2934] | — | — | — |
| H2 | B-loo:cofactor_pool | 0.3285 | +0.0369 [0.0363, 0.0375] | 12/12 | 39 / 0 |
| H2 | B-loo:metabolite_pool | 0.2969 | +0.0053 [0.0050, 0.0056] | 12/12 | 9 / 0 |

**Agreement with mode_scores.tsv: exact to 4 dp** on every tAP, Δ, CI and fold count. Not recomputed: the noCE cells and T0 = W1 (that equality is asserted in the author's code). Each variant also shifts 52–56 partial ratios by ≤ 1e-3 through the WT-growth denominator.

## 2. POST-HOC label-free random-removal control (not a gate)
Per line and seed (SEEDS[0:5], rng [seed, line]): set k random T0 lethal calls (e ≥ 0.99), drawn from the parent's T0 calls, to 0.

| | k | control vs T0 | variant vs control (seed-mean) [CI] | seeds | folds |
|---|---|---|---|---|---|
| H1 cofactor | 35 | −0.0281 | **+0.0624** [0.0616, 0.0632] | +0.0621…+0.0630 | 12/12 |
| H1 metabolite | 10 | −0.0078 | **+0.0134** [0.0129, 0.0138] | +0.0128…+0.0135 | 12/12 |
| H2 cofactor | 39 | −0.0293 | +0.0662 [0.0654, 0.0670] | +0.0656…+0.0667 | 12/12 |
| H2 metabolite | 9 | −0.0066 | +0.0119 [0.0115, 0.0122] | +0.0116…+0.0120 | 12/12 |

**Random removal hurts:** T0 calls are far more precise than the base rate. So neither gain is a call-rate effect; it depends on *which* calls the edit removes.

**The B2 control is a different kind of comparator.** The author's symmetric control drops the lowest-B2 calls. B2 is a training-label prior (tAP ≈ 0.90), which makes that control near-oracle: it gains +0.054 over T0 at k = 35 (H1). Losing to it (−0.020 / −0.008) means only that a label-free edit chooses calls less well than the labels' own prior. It is not evidence of a call-rate artefact and not a failed null.

## 3. D-1 ruling
**What the texts say.** DESIGN: "matched call rate. The same number of extra calls is assigned to the genes with the highest B2". Final list §5 (A2): "a call-rate-matched comparator … guards against task-set overfitting". Both define the control only for modes that **add** calls. For a removal-only B-loo, the literal control is the identity (≡ T0), so it is **not applicable**, as pseudo-tasks are not.

**Ruling.**
- The literal reading is the correct pre-registered reading. The symmetric / B2-removal rule is the author's own addition.
- The B-loo candidates therefore pass "both negative controls" only **vacuously**: they qualify on TG alone.
- The reading was fixed after both control readings had been seen on training aggregates. This is a forking-paths decision and needs a PI ruling as a §17 row.

**Mitigation.** The cofactor-pool edit is a priori motivated: final list §1.4a names the MAR10065 dead-ends as FP drivers, and list T10 pre-registered a "de-artefacted biomass". L-T10 itself (cofactor + glycogen out, plus list-T3 constraints) fails TG: H1 +0.0022, 7/12 folds.

**Required flag** (for shortlist.json, the val config and the README):
> **W3-D1 (control-contingent):** passes TG; no pre-registered negative control applies to biomass leave-outs (the matched-call-rate control is defined for added calls only; pseudo-tasks n/a). The control reading was fixed after the training aggregates were seen. Loses to the B2-informed symmetric removal comparator (H1 −0.020 / −0.008); beats a post-hoc label-free random-removal control (+0.062 / +0.013). No claim may rest on "beats controls".

## 4. Other deviations
- **D-2** (pseudo controls on H1 only): acceptable. G* = H1, and pseudo-tasks do not apply to B-loo. Keep H2 `beats_pseudo = NA` explicit.
- **D-3** (T-set naming): acceptable and clearly mapped.
- **D-4** (T4_11 ≡ T0): correct, since the Growth group is GR only; the test checks 11/10 groups.
- **RPMI proxy:** acceptable as secondary, but it exposes a control failure (§6).
- **No pinned worktree; MANIFEST rows dirty = True; no code SHA in MANIFEST:** a real gap. The first scoring run (23:12, symmetric primary) cannot be reconstructed, because w3_score_train.py changed at 23:13. My bitwise recomputation covers the headline numbers only.

## 5. Enumeration and ordering (check 2)
**Enumeration.** MAR13082 has 9 reactants in both models; minus ATP and H2O that leaves 7 direct components (verified on the loaded models). The MAR10065 groups partition the pool exactly: H1 38 = 20 retinoid + 4 vitamin-D + 2 tocopherol + 5 singles + 7 unassigned; H2 36 = 20 + 2 + 2 + 5 + 7 (unassigned: biotin, FADH2, PLP, riboflavin, BH4, THF, TPP). The variant JSONs equal `bloo_components.tsv` species for species.

**Order.** Structure runs 19:25, 19:32, 21:08 (current SHAs = the 21:08 MANIFEST row) → M3_W5_review §4 (label-free CAD/dihydroorotate finding) ≈ 20:28 → sweeps 21:32 / 21:52 → first W3 label open 23:08:17 (w3_score_train.py).

So all components, including dihydroorotate, were fixed about 2 h before any W3 label access. Dihydroorotate is outside DESIGN's list (it is a MAR10064 species), but it is label-free and was added before scoring; it belongs in the §17 row. It scored −0.0004 / 0.0000 and is not shortlisted.

## 6. Biology (check 3)
**H1 cofactor-pool B-loo un-calls 35 genes.**
- **Artefacts (a priori §1.4a), about 20 genes:** CYP26A1/B1, CYP24A1, CYP1A2, BBOX1, TMLHE, ALDH9A1, SHMT2, MMAA, GCH1/SPR/PTS, SLC5A6, COQ2/3/6/7, NPC1L1, SCARB1, SLC25A26. Removing these calls is the intended correction.
- **Real requirements lost, about 15 genes:**
  - heme synthesis: ALAD, HMBS, UROS, UROD, CPOX, ABCB6;
  - cofactor activation: PDXK/PNPO → PLP, TPK1 → TPP, RFK → FMN/FAD;
  - signal peptidase: SPCS1/2/3, SEC11A/C. These are true essentials, lethal here through a pool-precursor route (mechanism not traced).
- PDXK/TPK1/RFK are *explained* (PLP, TPP and riboflavin sit in the same pool, so whole-pool removal drops them), not *justified*.

**New finding: a positive-control failure (REJECTS a README claim).** With serum (RPMI M2 f = 1), B-loo:cofactor_pool makes **HMGCR, MVK, PMVK and MVD non-lethal (ratio 0.985) in both parents** (sweeps row `RPMI_M2f1|B-loo:cofactor_pool_biomass`).
- Cause: ubiquinol, the only non-sterol isoprenoid demand, leaves with the pool.
- DESIGN §9 W2 says loss of core essentiality under serum "flags a modelling error". README §2 says the core "also stays lethal under the shortlisted cofactor- and metabolite-pool B-loo"; that is true on M0 only.
- It does not affect the M0 look, but it will appear in the M6 W2 × W3 factorial.

**Verdict on the edit.** Whole-pool removal also removes legitimate requirements (isoprenoid/CoQ, heme, PLP, TPP, FAD): blunt, though informative.

**Narrower variant (label-free; do NOT substitute into this look).** "Artefact groups only": retinoid + vitamin-D + tocopherol catabolites + carnitine, keeping ubiquinol, heme, lipoate, cobamide and the 7 unassigned vitamins. Per-group aggregates have been seen, so it is post hoc: it needs a §17 row and an exploratory label, and must be fixed before scoring.

**Metabolite-pool B-loo.**
- Un-calls PSAT1 and GCLC/GCLM/GSS: plausible FPs in serine/cystine media.
- Also un-calls COASY, PPCS, NADK, GNPNAT1, PGM3 and, in H1, CAD: real CoA, NADP, UDP-GlcNAc and pyrimidine demands.
- Δ 0.0056 sits just above the 0.005 floor; the noCE cell is −0.0003 (author's value).
- The winner's-curse optimism, 0.049, exceeds **both** candidates' effects. The README says only candidate 2's.
- Weak on biology and on statistics.

## 7. Leakage (check 4)
- **Scanner:** `scan_label_copies.py --roots results/w3` scanned 380 files, 0 hits.
- **Tripwire, since 19:00:**
  - W3 scripts opened only `train` (w3_score_train.py 50, w3_val_dryrun.py 2; 23:08–23:20).
  - No W3 `val`, `test` or `locked` opens. The val/test hash opens are check_frozen.py; the DENYs are pytest toy fixtures.
  - Unattributed `train` opens (script = None): 20:22:14 (pid 1187213), plus bursts at 19:30–19:37 that coincide with W5 runs. None is a W3 script; the coordinator should attribute them.
- **Val config:** scalars_only, B = 2,000, seed 11, L_prob, U, strata none/expressed. Code SHAs match the current files, and config_sha recomputes.
- **Dry run:** synthetic Bernoulli labels, scalar guard pass, 0 ledger rows.
- **`w3_val_scores.parquet`:** 2,866 × (entrez + 6 float model-score columns). No ModelID, no labels.
- **Harness:** the look needs a new `W3_M4_shortlist` branch in evaluate_val.py, which requires a full independent code review under A6(2).

## 8. FROZEN candidates (check 5)
- All 6 file-row SHAs match.
- 52 value rows (README says 51): 32 B-loo rows equal `bloo_components.tsv`; 20 B-add rows, 2 NaN = logged skips.
- Suggest adding a row for `sweeps/ratio_r_<parent>.parquet` (the scored vectors).

## 9. Tests (check 6)
`pytest tests/test_w3.py`: 8 passed.
- **Meaningful:** shortcut = brute force on a toy model for 3 variant kinds; variant revert and JSON round-trip; pseudo-task determinism; capacity / mode_e / CellFie rules; group counts.
- **Gaps:** no test of the scoring or control logic (`mcr_matrix`; `fold_stats` and its < 5-line rule; TG; one-SE shortlist); no test of `pool_remove`; one vacuous assert (`… not in kos or True`).

## 10. Verdicts
| Claim | Verdict |
|---|---|
| T0 0.2996 / 0.2916 (= W1) | CONFIRMED (recomputed; W1 equality is the author's assertion) |
| cofactor B-loo TG pass, +0.0343 [0.0336, 0.0349], 12/12; H2 +0.0369 | CONFIRMED numerically; QUALIFIED as a shortlist entry (W3-D1, W3-MEV, lost real demands) |
| 35 removed calls, incl. PDXK/TPK1/RFK | CONFIRMED (`rv_changed_calls.tsv`) |
| metabolite B-loo +0.0056 [0.0051, 0.0059], gain in CE | CONFIRMED numerically (noCE not recomputed); QUALIFIED: weak |
| shortlist = these two; other modes don't qualify | QUALIFIED: correct under the literal control; other modes read from mode_scores, not rescored |
| P57 equality 8/8 | CONFIRMED from p57_equality.tsv (not re-run) |
| sterol / vitamin-D rescue; mevalonate core lethal everywhere | CONFIRMED for the 4 sterol conditions; **REJECTED** for the extension to cofactor B-loo under serum (0.985) |
| ETF control fails H1, passes H2 | CONFIRMED (train_summary) |
| dihydroorotate added before scoring; enables CAD uridine rescue | CONFIRMED (order); rescue read from the pyrimidine tables |
| no val/test labels read; val artefacts label-free | CONFIRMED |

## 11. Required fixes before commit
1. Add the W3-D1 flag (§3) to shortlist.json, the val config (a `flags` field) and the README. The PI records the D-1 reading and the dihydroorotate addition as a §17 row.
2. Correct README §2 and train_summary: under serum, the cofactor B-loo loses mevalonate-core essentiality (0.985, both parents). Add this as flag **W3-MEV**.
3. README: the winner's-curse optimism (0.049) exceeds both effects; there are 52 value rows, not 51.
4. Val config: the gate candidates are the [H1] variants only; the [H2] columns are descriptive.
5. Commit the W3 code, re-run w3_score_train.py from the pinned commit (≈ 10 min) and confirm identical outputs.
6. Add tests for `mcr_matrix`, `fold_stats` (< 5-line rule) and `pool_remove`; drop the vacuous assert.
7. The evaluate_val.py `W3_M4_shortlist` branch gets a full independent code review before the look.

## 12. Recommendation on the validation look
A **single** look is DESIGN-compliant (≤ 2 TG passers).
- **Candidate 1:** forward it with the W3-D1 and W3-MEV flags.
- **Candidate 2:** I recommend the PI drop it. It sits at the TG floor, inside the winner's-curse band, its gain is CE-only, and it removes real CoA, NADP, UDP-GlcNAc and pyrimidine demands. If it is kept, it goes in the same single look, flagged.
- A validation result for candidate 1 does not endorse whole-pool removal. The narrower variant (§6) is a separate, pre-specified exploratory item.

## Label exposure (for notes/label_exposure.md)
Reviewer, 2026-10-28: training-tier aggregates only.
- **From `train/mode_scores.tsv`:** rows T0, B-loo:cofactor_pool_biomass, B-loo:metabolite_pool_biomass and L-T10, for H1 and H2. Columns: tap_lin + CI, Δ̄_lin + CI, folds_positive, delta_noCE, vs_mcr, mcr_delta_vs_T0, vs_mcr_sym + CI.
- **Also viewed:** `train/shortlist.json`, the README headline table, and my own recomputed tAP_lin, Δ̄_lin + CI and folds for these 6 cells, plus the random-removal Δs (§1–2).
- **Not viewed:** no per-line or gene × line label values printed or written; no validation or test labels; no raw file; nothing in data/locked/. 