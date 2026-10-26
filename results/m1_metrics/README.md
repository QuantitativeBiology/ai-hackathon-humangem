# M1 metrics library: `src/nextgem/metrics.py`

This library implements the metrics of DESIGN.md v3, which is frozen: §5.2 (the NA rules), §6.1–6.6 (tAP, units, contrasts, claim tiers, SESOIs, bootstrap), §6.8 (binary endpoints under the Human1 and Human2 paper protocols) and the §11.4 failed-build rule. It also covers the §12 verification row "tAP and metrics".

The functions are pure. They take arrays, Series or DataFrames, never open a file, and every random draw needs an explicit `seed`; `None` is refused. **Only toy and synthetic data were used. No label file was read.**

Files:
- `src/nextgem/metrics.py`: the library, one file, so FROZEN.tsv covers it as a single entry.
- `tests/test_metrics.py`: tAP, null, units, NA rules, pairing, failed lines and inference.
- `tests/test_metrics_binary.py`: §6.8 binary metrics.
- `scripts/verify_metrics.py`: the 10^5-jitter Monte-Carlo check, the bootstrap calibration simulation, the known-gap positive control and runtimes. It writes `verify_metrics.json` and `verify_metrics.log` (this folder).
- `pytest==9.1.1` (with `iniconfig`, `pluggy`) is installed in `.venv` and pinned in `requirements.txt`.

## API

Conventions: score and label frames are **lines × genes**. Labels are {0, 1}, with NaN for excluded pairs. Scores are finite on every non-NaN label pair, and a higher score means more essential. Ties are exact float equality, so scores are rounded upstream (§6.1).

**tAP and null (§6.1)**

| function | what it does |
|---|---|
| `tap(scores, labels)` | Primary tAP: the O(n_g) form of §6.1, vectorised. NaN if there is no positive. |
| `tap_hypergeometric(scores, labels)` | The defining hypergeometric order-statistic sum, with log-binomials. Cross-check only. |
| `tap_weighted(scores, labels, weights)` | Integer multiplicities (duplicated items), using a harmonic closed form of the O(n_g) sum. Used by the bootstrap and permutation code. |
| `tap_columns(S, Y)` | Vectorised tAP of every column. NaN labels are excluded. |
| `expected_ap_random(n, n_pos)` | E[AP_random] = [H_n + (p−1)/(n−1)·(n − H_n)]/n. This is exactly the tAP of a constant scorer. |
| `ap_step_sklearn` | Sensitivity row only. Returns NaN, not 0, when there is no positive. |
| `auroc` | Secondary metric. |
| `expected_tp_at_k` | The §6.4 "true positives in the top k" translation, tie-aware. |

**Units and summaries (§5.2, §6.2)**

| function | what it does |
|---|---|
| `per_line(S, L)` | One row per line: n, n_pos, prevalence, estimable, tap, e_random, lift, delta_null. A line with fewer than 5 positives is NA. |
| `per_gene(S, L)` | The same per gene, plus n_neg. A gene needs at least 5 positive and 5 negative lines. |
| `pooled(S, L, lines=None)` | Pooled tAP over pairs, optionally within one lineage. |
| `m_delta_ap_gene` | mΔAP_gene. |
| `cell_is_na` | The 50%-of-lines / 30-genes rule for a grid cell. |
| `paired_line_delta` | Per-line ΔtAP. |
| `per_lineage_summary` | Mean tAP_line, E[AP_random], π, lift and Δ vs a reference. |
| `lineage_equal_mean`, `line_weighted_mean`, `lineage_means`, `lineage_weights` | Δ̄_lin and Δ̄_line. |

**Pairing (§6.5) and failed builds (§11.4)**

| function | what it does |
|---|---|
| `assert_same_pairs(*frames)` | Raises `PairingError` on a set mismatch. Every unit function also raises it on a non-finite score at a scored pair, and on different n or n_pos between methods. |
| `no_information_scores(S, failed, unit)` | `"line"`: the failed line gets all genes tied, so tAP_line = E[AP_random]. `"gene"`: the failed line takes the per-gene median over the non-failed lines of the same tier. |

**Inference (§6.3–6.5)**

| function | what it does |
|---|---|
| `stratified_bootstrap_counts(lineage, B, seed)` | Resampling multiplicities, with lines drawn within their lineage. |
| `bootstrap_mean(values, lineage, B, seed, weighting="lineage"\|"line", method="percentile"\|"bca")` | Bootstrap CI of a mean. |
| `skew_triggers_bca` | The \|skew\| > 1 rule. |
| `sign_flip_test(delta, lineage, n_flip, seed)` | p = (1 + #{T\* ≥ T_obs})/(1 + n) on Δ̄_lin. |
| `lineage_t_interval` | df = K − 1, two-sided 95%. |
| `exact_sign_flip_lineage` | Exact p over the 2^K sign configurations. |
| `transfer_summary` | t-interval, exact p, k-of-K with k = ⌈0.7K⌉, and a `generalises` flag. |
| `claim_tier(passed, estimate, ci_low, sesoi, *, gate_status)` | The §6.4 wording tiers. `gate_status` (from `fixed_sequence`) is required: a contrast that is "estimate only" or "pending" gets no tier. `sesoi=None` (C3) allows only "detectable". |
| `sesoi_c1b` | The lift-scaled SESOI for C1b. |
| `fixed_sequence` | C1 → C1b → C3 → C4 → C5 gatekeeping, **enforced**: a sequence must be a prefix of `GATE_ORDER` in order (otherwise ValueError); a mapping is reordered into `GATE_ORDER`. It stops at the first failure. A contrast not yet run is `None`, shown as "pending". |
| `confirmatory_contrast(delta, lineage, *, contrast, sesoi, seed)` | For C1, C1b, C4 and C5: estimate, bootstrap CI (B = 10,000), sign-flip p (10,000 flips), test pass at p ≤ 0.025, Δ̄_line and the transfer summary. `contrast` and `sesoi` are **required, with no default**. C1, C4 and C5 must use the frozen 0.01 (`SESOI_BY_CONTRAST`). C1b must pass its own `sesoi_c1b(...)` value; 0.01 is refused. No claim is returned here. The seed is spawned into separate bootstrap and sign-flip streams. |
| `gatekept_claims({contrast: result or None})` | Runs `fixed_sequence` over the results and attaches each contrast's gate status and claim (via `claim_tier` with that status). This is the only path to claim wording. |
| `headroom_closed` | Shared-resample ratio CI, reportable only if the denominator's CI excludes 0. |
| `gene_bootstrap(S_a, S_b, L, lineage, B, seed)` | The C3 line bootstrap of mΔAP_gene(a) − mΔAP_gene(b). Eligibility is fixed on the observed data; in each replicate, genes with 0 positives or 0 negatives are dropped and the mean count dropped is reported. |
| `c3_permutation_null(S, L, lineage, n_perm, seed)` | Within-lineage label permutation (2,000 permutations, 97.5th percentile). |
| `permute_within_lineage` | The permutation used by the C3 null. |
| `c3_decision(S_P1, S_B3, L, lineage, *, seed)` | The combined C3 decision. It asserts that P1, B3 and the labels are paired, then passes iff the gene-bootstrap CI lower bound of P1 − B3 is > 0 **and** mΔAP_gene(P1) is above the 97.5th percentile of its permutation null. The seed is spawned into bootstrap and permutation streams. The result feeds `gatekept_claims`. |
| `c3_claim_condition(c3_all, c3_expressed)` | The §6.3 claim condition: C3 passes, and the expressed-pairs version has the same sign with CI lower bound > 0. |

**Binary endpoints (§6.1, §6.8)**

| function | what it does |
|---|---|
| `binary_metrics(tp, fp, tn, fn, zero_division=nan)` | MCC, F1, precision, recall, specificity and accuracy. A zero denominator gives NaN; `zero_division=0` gives the sensitivity row. |
| `confusion_counts` | Confusion counts for one line. |
| `binary_per_line(calls, truth)` | Per-line counts and metrics. |
| `labels_from_effect(effect, -0.6)` | L-ge06, with a strict `<`. |
| `binary_threshold_grid` | The 21-point grid, `THRESHOLD_GRID` = −1.0 … 0.0 in steps of 0.05. |
| `summarise_binary(per_line, lineage)` | Lineage-equal mean, line mean and median, under two NaN policies: exclude (primary, with a count) and zero (sensitivity). |
| `fba_calls(ratio_r)` | ratio_r < 0.01. |
| `rate_threshold`, `calls_at_threshold` | Non-FBA scorers dichotomised at the positive rate of the fitting lines. |

Frozen constants: `SESOI_BY_CONTRAST`, `GATE_STATUSES`, `MIN_POS_LINE`, `MIN_POS_GENE`, `MIN_NEG_GENE`, `CELL_MIN_*`, `ALPHA`, `SESOI`, `GATE_ORDER`, `B_CONFIRMATORY`, `B_OTHER`, `N_SIGN_FLIPS`, `N_C3_PERMUTATIONS`, `FBA_CUTOFF`, `L_GE06_THRESHOLD`, `THRESHOLD_GRID`.

## Test results

Run: `.venv/bin/python -m pytest tests/test_metrics.py tests/test_metrics_binary.py`

**61 passed, 0 failed** in 10 s wall, with 1.1 GB peak RSS (2026-10-26, after the review fixes below; it was 56 passed before them). The §12 checks map to tests as follows:

| §12 check | test(s) |
|---|---|
| Brute-force enumeration of tAP on toy rankings with ties | `test_hand_example` (17/24 by hand); `test_brute_force_enumeration_random_toys` (300 toys, n ≤ 7, every within-block ordering in exact Fractions, compared with all three tAP forms); `test_brute_force_single_big_block` |
| O(n_g) form = hypergeometric form | 200 random cases, plus a 2,866-gene near-binary case and a 2,650-item block with 250 positives, where C(n, p) overflows float. The vectorised and weighted forms equal `tap` and duplicated items. All agree to 1e-12. |
| Constant scorer → E[AP_random] exactly | Eight (n, p) pairs, from (1, 1) to (5000, 4999), at three constants. A line-constant scorer gives tAP_gene = E[AP_random], so mΔAP_gene = 0. Also matches the stats-review reference E(192, 5) = 0.0507. |
| Mean over 1,000 seeds of score + 1e-9 noise ≈ tAP | A toy case and a near-binary case, within 4 MC SE. The untied AP agrees with sklearn. |
| 10^5 jitters (the §6.1 wording) | In `verify_metrics.py`: tAP 0.240526 vs MC mean 0.240548 (SE 3.4e-5, z = 0.62). sklearn step AP is 0.2057, lower because of the tie penalty; E[AP_random] is 0.114. |
| Perfect and reversed rankings | tAP = 1 for a perfect ranking, including tied positives on top. The reversed ranking equals the closed form; untied tAP equals sklearn AP. |
| The ≥5-positive NA rule | Lines: 4 positives → NA, 5 → estimable, NaN pairs excluded. Genes: the 5/5 rule. Grid cells: the 50%-of-lines / 30-genes rule. |
| Set-mismatch assertion | `PairingError` on a different line or gene set, on a NaN score at a scored pair, on different per-line n between methods, on headroom and binary mismatches, and on a line without a lineage. Reordering rows or columns is aligned, not rejected. |
| Failed-line per-gene median rule | `test_failed_line_rules`: per line, tAP = E[AP_random]; per gene, the median of non-failed lines, entered as a tie. |
| MCC and F1 vs hand-computed confusion matrices | Four hand-computed matrices. Zero denominators give NaN; `zero_division=0` gives the sensitivity row. |
| MCC and F1 vs sklearn where defined | 280+ random cases with no zero denominator: all six metrics equal sklearn to 1e-12. |
| Constant-call negative control | MCC is NaN under the primary rule and 0 in the sensitivity row (sklearn's convention). A random caller gives mean MCC ≈ 0. |
| Inference | Stratified counts preserve lineage sizes and are deterministic. Bootstrap SE matches the analytic SE. BCa ≈ percentile on symmetric data. The sign-flip p is uniform under the null, and its floor is 1/(n + 1). The t-interval equals scipy's. The exact sign-flip gives 1/128 and 8/128. ⌈0.7K⌉ = 5 for K = 7, and is computed without float artefacts. Claim tiers, fixed-sequence gatekeeping, SESOI_C1b, headroom and the end-to-end C1 contrast are each tested. |
| C3 | The gene bootstrap equals brute-force duplication of resampled lines, for every one of 40 replicates. The permutation keeps the NaN pattern and the per-lineage counts. A perfect scorer passes (p = 1/201). Line-constant and lineage-only scorers do not pass: their null equals the observed value. |
| Review fixes (notes/reviews/M1_metrics_review.md) | Each of these tests fails on the pre-fix code. `test_fixed_sequence_enforces_gate_order` (F1): wrong order, skipped C1b and unknown contrasts are rejected. `test_confirmatory_contrast_requires_contrast_specific_sesoi` (F2): a missing SESOI, 0.01 for C1b and a non-frozen SESOI for C1/C4/C5 are rejected. `test_c3_decision_combines_bootstrap_and_null` (should-fix 1): both arms are required (identical P1 = B3 fails the bootstrap arm), pairing vs B3 is asserted, plus the claim condition. `test_claim_tier_respects_gatekeeping_state` and `test_gatekept_claims_downgrade_after_failure` (should-fix 2): a C4 that passes after a C3 failure is "estimate only". |
| Seeds (§7) | `metrics.py` contains no `default_rng()`, `random_state=None` or `np.random.seed`. |

## Bootstrap calibration (§12 negative control; stats r1 SF5)

`scripts/verify_metrics.py` ran 500 synthetic datasets per regime, seeded from SEEDS[0] = 11 with `SeedSequence.spawn`, on 8 workers.

Each dataset has:
- 7 lineages with sizes 126/75/50/35/28/20/15, 349 lines in all;
- two exchangeable near-binary GEM-like scorers: 150 genes at e = 1, 50 partial, and about 93% of genes tied at 0, with 60 per-line swaps (informative swap-ins);
- training-like prevalence: a median of 320 positives in 2,866 genes ("all genes"), or 37 in 2,450 ("no common essentials"), taken from the round-1 stats review aggregates;
- a lineage random effect u_k ~ N(0, τ²) added to every line Δ.

Because the two scorers are exchangeable, E[Δ | lineage] = u_k exactly. The conditional estimand is therefore mean_k u_k, and the population estimand is 0. The within-lineage SD of line Δ came out at s ≈ 0.018–0.021. The bootstrap uses B = 10,000 and the sign-flip uses 10,000 flips, the production settings.

Coverage, with Monte-Carlo SE about 0.01:

| regime | τ | gap | boot CI ∋ conditional | boot CI ∋ population | t-interval ∋ population | sign-flip reject rate (p ≤ 0.025) |
|---|---|---|---|---|---|---|
| all genes | 0 | 0 | **0.954** | 0.954 | **0.964** | 0.020 |
| all genes | 0.01 | 0 | **0.956** | 0.438 | **0.942** | 0.238 |
| all genes | 0.01 | 0.01 | **0.952** | 0.438 | **0.942** | 0.958 |
| no CE | 0 | 0 | **0.934** | 0.934 | **0.956** | 0.016 |
| no CE | 0.01 | 0 | **0.940** | 0.534 | **0.956** | 0.204 |
| no CE | 0.01 | 0.01 | **0.934** | 0.534 | **0.956** | 0.954 |

What the table shows:
- **The stratified CI covers the conditional estimand at 0.934–0.956.** That is inside or at the edge of the 0.93–0.97 band that SF5 expected. The no-CE regime sits about 1.5 MC SE below 0.95, consistent with the slight percentile-bootstrap shrinkage in small (15-line) lineages.
- **The lineage t-interval covers the population estimand at 0.942–0.964.**
- **The stratified CI does not cover the population estimand when τ > 0** (0.44–0.53), and the per-line sign-flip rejects 20–24% of the time when only the lineage effect is present. This is the documented behaviour: the C1 CI and p are conditional on the 7 test lineages, and only the t-interval with the ≥5/7 rule supports a "generalises" claim. That rule passed in 1.8–2.6% of the null runs.
- **Positive control (known gap δ = 0.01):** the mean estimate is 0.0098 (all genes) and 0.0099 (no CE); sign-flip power is 0.96 and 0.95; "generalises" is reached in about 52% of runs at τ = 0.01.

## Runtime (single thread unless stated)

| call | time |
|---|---|
| `per_line` / `per_gene`, 349 lines × 2,866 genes | 0.23 s each |
| `confirmatory_contrast` (B = 10,000, 10,000 flips) | 0.34 s |
| `gene_bootstrap`, C3, 2,866 eligible genes, B = 10,000 | about 20 min (extrapolated from B = 200: 23.5 s). This is the worst case; the real eligible count is smaller. |
| `c3_permutation_null`, 2,000 permutations | about 3 min (extrapolated from 50: 4.6 s) |
| Unit tests | 10 s |
| `verify_metrics.py`, 8 workers | 169 s in total (MC check 8 s, calibration 1,000 datasets × 3 scenarios about 130 s, benchmarks 29 s) |

## DESIGN ambiguities resolved in the implementation (DESIGN.md not changed)

1. **§12 says a constant-call method gives MCC 0, but §6.8 says a zero denominator gives NaN.** Implemented as NaN in the primary summaries (lines excluded and counted) and 0 in the sensitivity row (`zero_division=0`, which matches sklearn). The test checks both.
2. **Per-lineage and grid "lift"** is the ratio of means: mean tAP_line / mean E[AP_random] over the lineage's estimable lines.
3. **Claim tiers are nested and gated.** A claim is given only to a contrast that is confirmatory in the fixed sequence (`gatekept_claims`). "Meaningfully improves" also requires the test to pass and the estimate to be ≥ SESOI. A contrast passes when p ≤ 0.025. The C1 sign-flip flips each line independently and recomputes the lineage-equal Δ̄_lin. A 1e-12 tolerance is used in the ≥ comparisons and in the C3 "> 97.5th percentile" comparison, because line-constant scorers produce exact floating-point ties.
4. **The lineage t-interval is two-sided 95%,** so its lower bound corresponds to one-sided 0.025. The exact sign-flip p counts configurations ≥ the observed value, identity included, over 2^K.
5. **`confirmatory_contrast` spawns its seed into [bootstrap, sign-flip] streams.** The C3 null spawns one child per permutation.
6. **C3.** Eligibility (≥5/≥5) is fixed on the observed data. A replicate drops only genes with 0 positives or 0 negatives. The permutation null leaves NaN-label pairs in place, so per-gene counts and eligibility do not change.
7. **Failed-line per-gene median.** It is taken over the non-failed rows of the matrix passed in, which the caller must restrict to the same tier. With an even count the median can fall between two score levels; it is still entered as one tie value.
8. **tAP is NaN for zero positives** (sklearn returns 0). Pooled tAP has no NA threshold beyond ≥1 positive, because §5.2 gives none. NaN values, such as NA lines or NaN binary metrics, are dropped before the bootstrap resamples.
9. **Ties are exact float equality.** Rounding (ratio_r to 6 dp) is the caller's job, under §6.1.
10. **Non-FBA dichotomisation.** The threshold is the ⌈rate·n⌉-th largest fitting score, and a gene is called essential if its score is ≥ that threshold. Ties at the threshold can push the realised rate above the target.
11. **Calibration parameters are synthetic.** The SF5 inputs (training-like prevalence, lineage sizes) come from the aggregates quoted in the design reviews, not from any label read. τ = 0.01 equals the SESOI. s ≈ 0.02 is at the low end of the 0.02–0.05 range in §6.4.

## Review fixes applied (ACCEPT WITH FIXES, notes/reviews/M1_metrics_review.md)

- **F1.** `fixed_sequence` enforces `GATE_ORDER`.
- **F2.** `confirmatory_contrast` has no SESOI default. It needs `contrast` and checks the SESOI against the frozen table; C1b needs `sesoi_c1b`.
- **Should-fix 1.** `c3_decision` (plus `c3_claim_condition`) implements the C3 decision.
- **Should-fix 2.** `claim_tier` needs `gate_status`, and `gatekept_claims` is the single path to claim wording.

No number changed: the calibration and runtime figures above are unaffected, because the fixes touch only the decision and wording layer.

Notes from the review that needed no change:
- The negative-control report must say that a constant caller has MCC 0 only in the sensitivity row.
- `THRESHOLD_GRID` is rounded to 2 dp. It differs from MATLAB's unrounded linspace only for values that fall exactly on a threshold.
- Plan the compute for the C3 bootstrap (about 20 min worst case).
