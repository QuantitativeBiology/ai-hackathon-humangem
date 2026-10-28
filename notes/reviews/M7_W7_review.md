# M7 W7 review (auto-research loop; non-author), commit ea8823f

**Overall: PASS, with one process deviation (F1) and a required REPORT wording (F4).** The look W7_shortlist_1 may run once the 26 FROZEN rows below are committed under A6.
Hook probe: `head -c 50` of the raw CRISPR gene-effect CSV was **BLOCKED** by the PreToolUse label guard, as required.
Labels: training tier only, via `cv.lolo` (labels.fold) and `labels.train()`. I read no validation or test label, no raw label file and nothing in data/locked/. I ran no look; only --print-config, --print-frozen-candidates and --dry-run. CPU used ≈ 0.70 CPU-h, ≤ 5 processes.

## 1. Protocol: PASS (F1 deviation)
- **Logged before run.** All 15 scored ideas have `logged_at` earlier than the start of their MANIFEST run, and `w7_common.assert_logged` refuses unlogged ids in code. Two margins are tight: W7-15/16 were logged 03:40 and their run started ≈03:42; W7-17 was logged 03:49 and its run started ≈03:50. W7-08 was logged at 02:58, before W7-03 was scored at 03:14. W7-07 and W7-10 were logged but not run.
  - Limit: these timestamps are self-reported. Everything landed in the single commit ea8823f, so git cannot corroborate the order. Only the MANIFEST end times and the ideas.tsv mtime (03:49, the last row) do, and both are consistent.
- **F1 (deviation, should-fix).** All 21 W7 MANIFEST rows have `dirty=True`: w7_screen.py and w7_val_spec.py call `provenance(..., allow_dirty=True)`, against §11.4 ("refuses to run when src/scripts/tests are dirty"). The code that produced the screens was therefore not pinned at run time.
  - Mitigation: §2 and §3 below reproduce the promoted screen and the look scores bitwise from committed code.
  - Fix: set `allow_dirty=False` before any further W7 run; record the deviation in RESEARCH_TREE.
- **Whitelist respected.** No drug screens, external essentiality sets, DepMap CE lists or co-essentiality (grep of w7.py and scripts/w7_*.py; all inputs come from w1/w2/m6 outputs; paralogs come from Ensembl 116, not De Kegel/Ito).
  - Minor: the kNN prior (W7-15/16/17) is a within-fold, label-derived prior over expression-similar lines. It is a B2 generalisation, not gene×gene co-essentiality, but it is not literally on the whitelist. It was not promoted, so the point is moot.
  - No idea targets a lineage. W7-15's "haematopoietic" example is generic.
- **Stop rule correct.** The 10 consecutive failures after the W7-05 promotion are W7-06, 09, 11, 14, 13, 15, 16, 08, 17 and 12. Each was logged separately; counting W7-15 and W7-16 as two ideas is acceptable.
- **Time box and budget.** The pause is recorded (autoresearch/pauses.tsv, 63 min API limit). Total wall time 02:42–05:10 is under 4 h even with the pause included. W7 used 8.84 CPU-h in results/budget.tsv.
  - stop.json cites a 100 CPU-h cap (§13/§14) where §W7 says 300. 100 is the tighter cap, so this is harmless.

## 2. Independent recompute of W7-05: PASS (exact)
`results/w7/review/rv_w7_05.py` does not import w7.py or scripts/w7_*. It reuses the w1 B6 feature code and `w1_b6.stacked_priors`, which I read and confirmed are fit-set-only. It runs 12-fold LOLO × SEEDS[0:3] on 686 estimable training lines in 0.63 CPU-h (full 3 seeds, not reduced).

| | reviewer | author |
|---|---|---|
| tAP_lin W7-05 / B6 / B2 | 0.9029 / 0.8886 / 0.9001 | 0.9029 / 0.8886 / 0.9001 |
| Δ W7-05 − B6 | +0.0144 [0.0133, 0.0154], 12/12, TG pass | +0.0144 [0.0133, 0.0154], 12/12 |
| Δ W7-05 − B2 | +0.0028 [0.0022, 0.0034], 12/12 | +0.0028 [0.0022, 0.0034], 12/12 |
| per-seed W7-05 | 0.90290, 0.90291, 0.90300 | same |

The per-line tAP is identical to the author's perline_r1/r2b tables: max |Δ| = 0 over 2,058 rows for B6 and for W7-05, and 686 rows for B2.

## 3. Leakage in the validation scores: PASS
- **Refit check.** I refit B6 and W7-05 on `labels.train()` only (SEEDS[0:3], every fitting line asserted training-tier) and predicted the 198 validation lines from label-free features. The seed means are **bitwise equal** to the `W7:B6` and `W7:W7-05` columns of the look parquet (567,468 rows, 0 missing).
  - So these scores are a function of training labels alone.
- **Other columns.** H1, H2 and P1 are bitwise equal to results/m6/val/m6_val_scores.parquet. The parquet holds no label column.
- **Features of W7-05.** b2/b2m for validation lines come from `fit_components(train)` via `heldout_priors`. The other features are expression, CN, mutation, aggregates, isozyme count, medium class and serum. W7-05 has no kNN, deviation or within-validation prior, so no validation-derived feature can enter.
- **Tests.** tests/test_w7.py: 12/12 pass (run with NEXTGEM_NO_TEST_RECORD=1).
  - Coverage: the poisoning tests corrupt validation/test omics, aggregates and GEM features (and expression for the kNN ideas) and require bitwise-identical predictions. Fits containing non-training lines raise LeakageError.
  - Gaps (minor): the mono=True idea is not in a poisoning test (only a run test); line_features and the W2 medium features are not poisoned. The refit check above covers W7-05 directly.

## 4. A6 config review of look_spec_W7_shortlist_1.json: PASS
- **Workstream, candidate, reference, methods.** W7 under results/w7/val, so it counts against the W7 sub-cap of 20 (ledger: 0 W7 looks used, 3/40 in total); milestone M7 matches §14. Candidate W7:W7-05; reference H1 (G* from the W1 look); H1, H2, P1 and W7:B6 are descriptive.
- **Grid and seeds.** B 2000, seed 11, labels [L_prob], genesets [all, noCE], universes [U], strata [none, expressed_pairs]. These are identical to the W2_M4, W3_M4 and M6_factorial specs.
- **Provenance.** 27/27 pins equal the committed files; the ideas.tsv SHA 92fdbef… also matches. --print-config matches the spec.
- **--dry-run.** config_sha ed9b6876… equals stop.json. Output is identical to results/w7/tmp/dryrun_W7_shortlist_1.txt apart from wall_s. scalar_guard passes, no validation labels are read, 0 ledger rows are written and git status is unchanged.
- **FROZEN rows needed (26, from --print-frozen-candidates).**
  - Look inputs: look_spec_W7_shortlist_1.json (1b4c13b8…), w7_val_scores_W7_shortlist_1.parquet (0bcf12c8…). Code: src/nextgem/w7.py, scripts/w7_common.py, w7_screen.py, w7_summarise.py, w7_val_spec.py, w1_b6.py.
  - Results: screen_{g1,genecv,r1,r2a,r2b,r4..r13}.tsv (15 files), nested_W7-04_beta.json, sweeps/ratio_r_human1.parquet, sweeps/sweeps.json.
  - Already frozen and unchanged: w1.py, m6 SHA256SUMS.tsv, m6_val_scores.parquet, evaluate_val.py and metrics.py.
  - --preflight-only will fail until these rows are committed.

## 5. Is W7-05 a legitimate P2 candidate? Yes procedurally; F4 on wording
- **Why it is legitimate.** §W7 step 7 takes the best validated idea. A12a and the M3M4 process review put "B6 tuning" explicitly in W7 scope, and §8 places B6's grid in W7. The whitelist allows non-GEM features. W7-05 passed TG against the current best, which was B6, the stronger of P1 and B6 on training.
- **F4 (must, REPORT.md).** P2 = W7-05 contains **no GEM output**. It is B6 with its M2 hyperparameters plus a monotone-increasing constraint on b2/b2m: a tuned non-GEM gene-prior model. REPORT must state the following.
  - **(a) What C4 means.** C4 = P2 − G* compares non-GEM ML with generic FBA. Training B2 − G* is ≈ +0.60, so C4 is expected to pass and is no evidence for GEMs.
  - **(b) What C5 means.** C5 = P2 − B6 measures only the effect of the monotone constraint. "B6" in C5 must stay the pre-registered M2 configuration, named as such.
  - **(c) W7's GEM result is negative.**
      - H-a (W7-02 / W7-12 / W7-14): −0.0030, +0.0001, +0.0001 against the current best. H-d (W7-04): β = 0 in 8/12 folds, Δ +0.0000. The artefact-groups variant scored 0.3299. No GEM-containing idea was promoted.
  - **(d) P2 vs B2.** On training P2 beats B2 by only +0.0028, below the 0.01 SESOI. Report it next to B2, B2m and B8±.
  - **(e) Gene identity.** Gene-grouped CV gives tAP_lin 0.185 for W7-05 vs 0.169 for B6, so the gain rests on gene identity.
  - **(f) Counts and shrinkage.** 17 ideas logged, 15 scored, 1 promoted, 1 validation look. Shrinkage LOLO → validation → test, and σ·√(2 ln m) with m = 1.
  - The PI should confirm before T2 that a non-GEM P2 is acceptable under §6.3. Otherwise P2 = P1 (no validated GEM idea) and C4/C5 are reported as such.

## Training aggregates viewed (all training tier)
- autoresearch/results.tsv (every row: tAP_lin, Δ, CIs, folds); screen_r2b.tsv and screen_genecv.tsv (B2, Gstar, P1, W7-01, W7-05, genecv); my rv_w7_05_lolo.json (§2); per-line training tAP in memory for the comparison with perline_r1/r2b; the synthetic-label dry-run scalars. No validation or test aggregate.

Outputs (not committed): results/w7/review/rv_w7_05.py, rv_w7_05_lolo.json, rv_w7_05_valfit.json; results/w7/review/tmp/ (git-ignored per-line table and logs).
