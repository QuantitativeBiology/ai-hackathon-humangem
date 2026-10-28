---
reviewed: scripts/evaluate_test.py (sha256 989ba31c…cc95e), tests/test_evaluate_test.py (853c8f45…de11f), results/m6/p{1,2}_spec.template.json (142df6ba…, b90d0045…)
commit: 7969549
author: evaluate_test.py harness agent
reviewer: independent non-author reviewer (Claude Opus 5.5 subagent)
verdict: APPROVE WITH FIXES (F1 blocks the freeze; F2–F5 and F11 before the rehearsal; rest before T2 or report-only)
independent_code_path: own enumeration / own RNG recomputation of every statistic (scratch handcheck.py, 22/22 pass)
exposure: synthetic only (no label value, label summary, ledger detail or locked-dir content seen)
---

# M6 review — test-look harness `scripts/evaluate_test.py`

**Exposure.** Synthetic only. Hook probe `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **BLOCKED**. I read no label
file, no locked dir, no ledger detail column and no skew/score summary. Dry runs read training labels inside the harness
only; their stdout carries synthetic-label scalars, which I saw.

**What I ran.** `pytest tests/test_evaluate_test.py test_evaluate_val*.py test_metrics*.py test_labels.py test_frozen.py
test_ledger.py`: **226 passed** (BLAS 1 thread, NEXTGEM_NO_TEST_RECORD=1). Own-code toy recomputation: sign-flip p vs
exact 2^10 enumeration, unequal lineages (0.0249 vs 0.0234, < 1 SE; line-weighted would be 0.0195); bootstrap CI vs exact
resampling ([0.05, 0.85] both); C3 observed mΔAP_gene with own AP/E[AP_random] (to 1e-12), own within-lineage null (97.5th
pct 0.083 vs 0.079, MC); gatekeeping and claim tiers incl. SESOI_C1b; transfer t (df 6), exact 2^7 flip, k-of-7 with a
zero lineage mean; headroom and reportability: all agree. On a label-free toy spec in `results/m6/tmp/review_toy/`
(deleted): `--print-config` (T1, T2), `--preflight-only` (refuses correctly), `--print-frozen-candidates`, `--dry-run` T2
(0 ledger rows, 0 files, guard pass); the T1 dry run was stopped at 1 h 45 min CPU (see Runtime).

## Findings (by severity)

**F1 — HIGH (blocks the freeze). Frozen tripwire test contradicts point #4.** `tests/test_tripwire.py:203` asserts
`evaluate_test.py == "val:denied test:{test_ok}"`, i.e. a *frozen* evaluate_test.py must be denied validation reads
(A5: "trusted for its own tier only"; §11.1: validation labels only from evaluate_val.py). The harness reads
`labels.validation()` at `evaluate_test.py:1641` (T1/T2) and `:1631` (unseal); `tripwire.ENTRY_PERMS` allows it.
Scenario: PI adds the FROZEN row → next full run gives `val:ok` → test fails → `check_test_runs` fails → every look's
pre-flight (`check_frozen_full`) is blocked. Fix: a §17 row (A14) recording the train ∪ val refit read (§4.2/§5.3) and
the assertion changed to `val:{test_ok}`, in its own commit before the freeze.

**F2 — MEDIUM. Unseal relies on bitwise re-derivation that nothing pins.** `run_unseal` (`:1439-1473`) publishes a file
only if the re-computed %.6g TSV hash matches; there is no sealed reader. No BLAS thread limit, no library versions in
`cfg` (`:338-361`). Scenario: T1 with 8 OpenBLAS threads, unseal with 1 (or a numpy bump) → one 6th digit flips →
`per_gene.tsv` is never published and that sealed output is lost. The config is also not compared with the look's
`started` config_sha before test labels are re-read; unseal runs outside `provenance()`; the mismatch path is untested.
Fix: `threadpool_limits(1)` in real modes, versions + threads in `cfg`; refuse unseal unless the config_sha matches
(before any label read); wrap in provenance; add a mismatch → not-published test.

**F3 — MEDIUM. Look-state refusals come after the validation read.** `state_problems` (`:542-580`) does not refuse
T1/T2 when the id is already completed; the ledger refuses only at `started`, after `labels.validation()` (`:1641`) and
after `write_config` (`:1643`) overwrote the committed `results/m6/test/test_look_T1_config.json`. Scenario: an
accidental second `--look T1` with an edited spec reads validation labels (no ledger row) and clobbers a committed file;
same for a repeated rehearsal. Fix: refuse completed ids / differing config_sha in `state_problems`; write config after.

**F4 — MEDIUM. T3 is not tied to T1's P1.** T3 (`:573-576`) uses whatever `p1_spec.json` is committed; the
authorisation names only a test and a fix commit. Scenario: spec swapped to another pipeline + any committed test → a
second P1 is evaluated on test. Fix: authorisation names old/new spec SHAs; check spec_id/primary/methods/arms against
T1's spec (`git show <T1 head_commit>:results/m6/p1_spec.json`) and that the test is frozen; test acceptance/"used".

**F5 — MEDIUM-LOW. Rehearsal match is spec-only** (`:563`), not code. Scenario: harness fixed after the rehearsal → T1
runs an unrehearsed path. Fix: store the `cfg["code"]` hash in the started detail and require equality.

**F6 — LOW-MEDIUM. Per-task table ≠ the W3 M4 table.** `per_task_table` (`:1147-1171`) uses pooled pair precision and
tAP of the binary call; W3 (`scripts/w3_score_train.py:529-561`) uses lineage-equal per-line precision and tAP of "score
with task-failing genes set to 1" with Δ vs base. Training and test are not comparable. Mirror W3.

**F7 — LOW. Revealed set slightly exceeds §10.2:** tAP/E[AP_random]/lift of P1, G*, arms (`:946-953`; `e_random_lin` is
a test-prevalence summary), headroom (`:939`), non-primary arm vs G* (`:935`), C3 expressed pairs (`:918`, defensible),
non-G* generic binary rows and arm − parent MCC (`:852-870`). All scalar. Seal the `tap:*` rows or log an A6 note.

**F8 — LOW. Test quality.** `test_real_constants_are_enforced` (`:778`) is near-vacuous (`check_constants` compares `cfg`
with the constants it was built from; assert DESIGN literals). The C3 "perfect" branch reuses the same function.
Untested: unseal mismatch, T1 spec_sha branch, T3 acceptance, real `b6_fit_predict` (stubbed; my T2 dry run ran it on
real features), real sealed write through the tripwire (recorder only). All parametrised tests have explicit `ids=`.

**F9 — LOW.** A bare `evaluate_test.py` defaults to the real T1 (`:1560`). `bootstrap_method` is a free field: pin
`results/m6/t1/skew.json` and check it. Failed-build calls are all False (`:743-745`), so MCC is NaN and those lines are
excluded (`n_nan` and the "zero" rows show it). The selective-genes stratum is skipped in the sealed strata (`:1112`).

**F10 — INFO (T2).** W7 ideas (B6 + GEM features, stacking) are label-fitted, so `fitted_on_labels: false` (`:231`) blocks
T2 for them: refit plug-in + review needed. Pin `b6.hgb` to W7's values; refuse T2 unless the primary label is L_prob.

**F11 — MEDIUM (coordinator's C1 = C1b question). CE_fit is applied correctly but is empty in every test-shape dry
run.** C1b uses Y minus `prep["fitted"].ce[lab]` (`:879`, `:902`), with CE_fit from `fit_components(train ∪ validation)`
(`:718`): the real look is right. `dry_run` (`:1499-1501`) adds *synthetic* validation labels at a flat ~0.12 rate, so a
gene positive in 100% of ~700 training lines gets a pooled prior ≈ (700 + 200·0.12)/900 = 0.80 < 0.90 → CE_fit = ∅ and
C1b ≡ C1 (hence the same +0.00003 in `results/m6/tmp/t1/dry_run_T1.json`). Toy (own script, 20 true CEs): train only
20; train + real-like validation 20; train + the harness's synthetic validation **0**. The rehearsal's test-shape stage 1
(`:1614`) shares this; its validation stage (training-only fit) does not. Fix: draw synthetic validation labels per gene
at the training prior; assert a non-empty CE_fit and a smaller C1b gene set in the dry run.

**Checked and correct.** Pre-flight → provenance → refit → committed `started` → evaluation-tier read. No dry mode calls
`labels.validation`/`test`. Sealed bytes go only via `labels._write_protected` to the locked `sealed/<look>/` dir. The ledger
holds names and SHA-256 only. T2 needs the `W7_stop` row naming the P2 SHA; unseal needs a completed T2. B2/B2m/B2c/CE_fit
are refit on train ∪ val; B3 and B4/B5 equal the W1 generic scores to 1e-9; B8±/B8′ use G*; the B6 refit matches
`w1_b6.fit/predict`. §6.8 uses `ev.BINARY_PROTOCOLS` (as in W1).

**Runtime (author: T1 1–1.5 h) — underestimated.** Real shape (380 lines × 2,866 genes, 1 thread, load 11–19): bootstrap
0.49 s/rep, permutation 0.18 s/rep → **one C3 ≈ 88 min; C3 + expressed-pairs C3 ≈ 3 h; T1 ≈ 3–3.5 h**; rehearsal ≈ 4–5 h.
My toy T1 dry run was still running at 1 h 45 min CPU when stopped. Over 1 CPU-h, so §11.4 needs a pinned worktree (not
enforced). Parallelise the two C3 calls or the replicate batches. T2 dry run: 3.3 min, 2.3 GB.

## Rulings on the 16 points
1 Agree (PI-accepted), with F2. 2 Agree. 3 Agree, add the code match (F5). 4 Agree on substance (§4.2/§5.3 require it;
nothing else is fitted on validation; fit_sha recorded), but it needs the A14 row and test change (F1). 5 Agree: SEEDS[0]
spawned [boot, flip]/[boot, perm] per contrast (common streams across contrasts are fine), estimates B = 2 000; BCa field
acceptable if pinned to skew.json (F9). 6 Agree with F7. 7 Agree; pre-declare CI-lower-bound > 0 for all three as the claim
rule and the estimate version as descriptive. 8 Agree (§11.4 rules exact; binary caveat F9). 9 Agree (sealed, exploratory;
pFBA tracing after unseal). 10 Disagree as implemented (F6). 11 Agree (A12c); the real P1 uses M2(f=1), so (ii) medium
alone must be estimated, not "lever not selected" as in the template example. 12 Agree with F4. 13 Agree. 14 Agree
(cosmetic; the test ledger has no workstream caps). 15 Agree (safe refusal; see F10). 16 Agree (`label.primary` = L_prob,
`sesoi.C1b` value row; sesoi_c1b.json also pinned).

## FROZEN rows to add (after F1, and after F2–F5/F11 if code changes; recompute SHAs then)
- `scripts/evaluate_test.py` file 989ba31c790fce560f564fa2f4e56e3e52d83d97618820f20ae28bb86d8cc95e (test looks: code)
- `tests/test_evaluate_test.py` file 853c8f458dbcc238abc7a237eb395143d90547acdfd0eaa61b26d4dde02de11f (known-answer test)
- `tests/test_tripwire.py` refreshed under A14 (F1)
- before the rehearsal: `results/m6/p1_spec.json` and every file it names (`--print-frozen-candidates`)
- before T2: `scripts/w1_b6.py` fe84857e…41a8, `scripts/w1_common.py` e77fed31…24a6, `results/w1/tmp/{mut_dam,mut_hot}.parquet`,
  `results/w1/tmp/gene_neighbours.json`, `results/m6/p2_spec.json` and its files (not frozen today).

## Re-check (18dc2dc: evaluate_test.py 8f72687f…, tests 84948f0d…; F1 in 1ee9f6a) — final verdict **APPROVE** for T1
Exposure was synthetic only. `pytest test_evaluate_test.py test_tripwire.py test_evaluate_val*.py test_metrics*.py test_w7.py`: **253 passed**. `--look T1 --print-config` on `results/m6/p1_spec.json` (config_sha 326015bc…): threads 1 and versions recorded, (ii) = [M2[H1], H1], percentile bootstrap pinned to skew.json. `--look T1 --dry-run --workers 8`: rc 0, **14.0 min wall**, 3.5 GB, 0 ledger rows and 0 files, no validation or test read, CE_fit non-empty, C1b ≠ C1, guard pass.
- **F1–F11, #7 and #11 are verified fixed.** F1: the test now expects `val:{test_ok}`; logged under A6, since tripwire.py already granted the read. F2: config_sha re-checked before any label read; the unseal runs under provenance. F3: refusal happens before pre-flight, any read and any write. F4: T3 is tied to T1's spec (SHAs, identity fields, frozen regression test). F5: `code_sha` must match the rehearsal's. F6: the per-task table uses the W3 definitions. F8: the constants test asserts the DESIGN literals. F11: synthetic labels follow per-gene priors; `prepare` refuses an empty CE_fit. Advisable, not blocking: a one-line §17 note that §11.1 ("validation only from evaluate_val.py") excludes the train ∪ val refit.
- **Parallel C3 is bitwise equal to serial in my own toy** (7 unequal lineages, ties, 3% NaN, 11 batches, 113 permutations, workers 2 and 8). Seeds cannot be reordered: the counts are drawn once in the parent, permutation i always uses child i, and results are written by index.
- **Refit plug-in:** only W7-01 and W7-05 are accepted (checked in the spec and again at fit time), and only B6-feature classifiers (W7-02 and W7-17 are refused, tested). Fitting lines must belong to the split's training or validation tier; a test line raises `RefitLeakageError`, and predicted lines cannot be fitting lines. `w7.py` and `w1_b6.py` are pinned and frozen; the hyperparameters must equal `w7.HGB_M2`. Validation-derived quantities reach only `S[P2]`, alongside the pre-registered B2 and B6 refits.
- **Residual (low; decide before any T3):** `unseal_problems` refuses the whole unseal if any sealed look's config changed. A T3 fix changes code or spec, so T1's config would no longer match and the unseal would be blocked. Before authorising a T3, decide whether T3 supersedes T1's bundle or whether T1 is unsealed from a worktree pinned at its commit.

**FROZEN rows for T1** (path · sha256; every other required file is already frozen at its current SHA):
`scripts/evaluate_test.py` 8f72687f410a84bc6db54d1d796886834cbc752d5f17a08337e18646f8d75d97 · `tests/test_evaluate_test.py` 84948f0d44a063d82b508c1f75a3cc8d166abbd3faa20d3f96a0610d3f9187f8
`results/m6/p1_spec.json` 750054c3e78510483e57bb49c0348bcacfc79a1dd71d3b81a4032e2e03b12104 · `results/m6/t1/skew.json` 2a70c42bfd55a2fb449fab34488fb97a92e2e85850a8882a5c7a65074943b51f
`results/m6/t1/p1_scores_validation.parquet` 18665d5992b5b3fd2fea06f1fe1cb6e994e88943d12d1fa44f5acb5621dc30a0 · `results/m6/t1/p1_calls_validation.parquet` e7d39aa08ccd4ba48a7d0b14388ef804cdda05e3a1d6612f2333a39aeb64ef0f
`results/m6/tmp/t1/p1_scores_test.parquet` 658b6d9b3c9f8d7a8de416151014bb03bf3b56c55c2fa8388ca50bef66925fb2 · `results/m6/tmp/t1/p1_calls_test.parquet` 04ce66b3fc1fb72be546daa190b74c5398ca5dab8eb45cf769b1427aa250ada9
`results/m6/t1/p1_task_calls.parquet` b769f20e73bb9806a9d902aee5fbc67061a01c1e9b6735b33301b755cf2f36e0 · `results/m6/t1/p1_gene_annotations.parquet` 074275d6f5e627352b0d8e8a22823f46a41ea4e61864f0a317219a0ea573758d
`results/m6/t1/tasks_bloo_human1.parquet` 2cdab12418140b366de00fc97f58bffccd09c314056d2aff75220cd4149d25a8 · `results/m6/t1/tasks_bloo_human2.parquet` b1e25f128f22b00c9007cf8c5446507b70ea1bcfb12a9427d5f7f358cbad16ec
`results/m6/t1/tasks_bloo_human1_tasks.tsv` 62f94d17d43034de632c614a46164e657f77ab5b2ba472905dbe1e56956f33d0 · `results/m6/t1/tasks_bloo_human2_tasks.tsv` b88ae1873988808fa425f00990cb4016fe7463f9d2d61ffea5eba703e5403892
`results/m6/t1/build_summary.json` aadb2300a8f9e70ee3120a4306dd9944f950c04af027bdef5e5b44871d450af4 · `scripts/t1_prep_common.py` efd35abec57826b6559af3b16999fc06e6b8dd9ccbc75fd30593ba1f9f2fd57d
`scripts/t1_prep_tasks.py` 12b377eaa03be2d95b21668009c007bab7c7f7ac91a6ac6d93391e9e56a68541 · `scripts/t1_prep_build.py` 9b9133b4bfd07fddb6f0088c1d573dd36763534cb8ca789e026a931dc897db3d
`scripts/t1_prep_skew.py` 0d4b3d20f6cec0f30ac3dc333fbf51ebc54e092cab414d475e263d9c35ce1a18 · `scripts/t1_prep_spec.py` 8390867e7b59c016b42f6d58558f6abb609c638443f6146d05e5532155ce2e7c

**FROZEN rows for T2** (separately, after the W7 stop; `w7.py` and `w1_b6.py` are already frozen):
`scripts/w1_common.py` e77fed313e65d8245288c41124be4451dd262b35cc10a29359f839c6f8ce24a6 · `results/w1/tmp/mut_dam.parquet` 40120741609fe308eb956dff408c39fe6b07b13248c29a021d95fbc97cddcb22
`results/w1/tmp/mut_hot.parquet` f7221063c4969928ae8e339b3078fd32c78c0c3a8fff871854ea3319b548302f · `results/w1/tmp/gene_neighbours.json` b492702e51a20f64d58f8a4f3562d97e3d6442b3ac244d39f1ae977ca0fcc233 · `results/m6/p2_spec.json` (not yet written: fitted_on_labels true, refit {w7_idea, W7-05, method = primary}) and its score and failed-build files.
