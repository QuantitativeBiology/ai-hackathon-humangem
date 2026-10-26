---
artifact: src/nextgem/metrics.py (+ tests/test_metrics.py, tests/test_metrics_binary.py, scripts/verify_metrics.py, results/m1_metrics/)
design_ref: DESIGN.md FROZEN v3, section 6 and 6.8, section 12 row "tAP and metrics"
author: metrics builder
reviewer: independent reviewer (Claude agent; did not write the library)
date: 2026-10-26
verdict: ACCEPT WITH FIXES
independent_code_path: yes
label_access: none (synthetic data only; no label file opened)
files_edited: none (no commit)
---

# M1 metrics library review

Independent scripts (not committed) are in the session scratchpad `.../scratchpad/rev/r1.py` to `r4.py`. They import only `nextgem.metrics` as the system under test; every reference value comes from code I wrote from the DESIGN text and from first principles.

## Verdict

ACCEPT WITH FIXES. The statistics match the frozen text, and all numerical checks agree to machine precision or within Monte-Carlo error. There are no numerical errors. The fixes below are guard-rail gaps (a gatekeeping order that is not enforced, a SESOI that is easy to forget, no combined C3 decision), plus two documentation points. None changes a number already produced.

## (1) Tie-aware expected AP, independent implementation

Reference methods, none using the O(n_g) or hypergeometric formulas:
- A. Brute force over all within-block permutations (n <= 7): 250 cases.
- B. Exact enumeration of every placement of each block's p positives among its n slots, in exact Fractions (n = 8..16): 150 cases.
- C. Monte-Carlo over random tie-breaks, 4,000 draws (n = 30..200, 1-5 distinct score levels): 30 cases.
- D. E[AP_random] against `tap` of a constant scorer and against enumeration: 11 checks.

The case mix includes heavy ties, all-tied scores, no positives (both sides NaN) and all positives (AP = 1).

| check | cases | agreement | max abs error |
|---|---|---|---|
| `tap` vs permutation brute force | 250 | 250/250 | 2.5e-14 |
| `tap_hypergeometric` vs permutation brute force | 250 | 250/250 | 2.5e-14 |
| `tap` vs exact placement enumeration | 150 | 150/150 | 1.1e-16 |
| `tap` vs Monte-Carlo | 30 | 30/30 within 4.5 SE | 2.2e-3 (about 1-2 SE) |
| `expected_ap_random` vs constant scorer / enumeration | 11 | 11/11 | 5.6e-17 |

Total 691 cases, 0 mismatches. In addition, `tap_weighted` against explicit duplication of items agrees on 100 random cases (0 differences). That is the identity the bootstrap and permutation code rely on.

## (2) Binary metrics and the threshold grid (section 6.8)

- Scalar check: my own MCC, F1, precision, recall, specificity and accuracy against `binary_metrics`, on 2,000 random confusion matrices (12,000 comparisons, 305 of them zero-denominator NaN cases). There were 0 mismatches. `zero_division=0` also matches (undefined set to 0).
- Per-line and grid check: `binary_threshold_grid` against my own confusion counts on 40 synthetic lines with about 5% NaN effects, a constant-call line, a no-positive line and a fully-positive line. This covered 21 thresholds, 5,040 metric comparisons and all four counts: 0 mismatches. The grid is exactly 21 steps from -1.0 to 0.0, the label is a strict `<`, and a value of exactly -0.6 is not positive.
- `summarise_binary`: the lineage-equal mean and the line mean agree with my own under both NaN policies. The NaN count is correct.
- Zero denominator: NaN is excluded and counted, and the 0 row is a sensitivity row. This matches section 6.8 and the repository MATLAB script (0/0 = NaN). F1 = 2TP/(2TP+FP+FN) as in the script.

**Section 12 vs section 6.8.** Section 12 says "constant-call method gives MCC 0" as the negative control. Section 6.8 says a zero denominator is NaN, excluded and counted, with a sensitivity row at 0. DESIGN.md has no explicit precedence clause (grep for precedence, conflict and prevails finds none), so the builder had to choose. The choice is defensible and consistent with the frozen text:
- Section 6.8 is the normative specification of the metric. Section 12 is a verification plan.
- The resolution keeps the primary rule as written in 6.8 (NaN) and satisfies the 12 negative control through the explicit sensitivity row (MCC = 0, equal to sklearn's convention).
- It violates neither section.

Two caveats:
- A constant caller has MCC = 0 only in the sensitivity row, not in the primary summary. Any report of the negative control must say so.
- A **non-degenerate** random caller gives MCC about 0 under both policies (tested).

## (3) Inference against DESIGN

| item | finding |
|---|---|
| Lineage-equal weighting | w = 1/(K n_k), computed on non-NaN lines. Matches my independent mean of lineage means (diff 1e-17). |
| Stratified bootstrap | Multinomial resampling within lineage with n_k fixed. My loop-based resampler (B = 20,000) gives CI 0.02557-0.04455; the library gives 0.02534-0.04449 (B = 20,000) with the same estimate. Percentile intervals; BCa only through `skew_triggers_bca`. OK. |
| Sign-flip p | Per-line independent flips within lineage, recomputing the lineage-equal mean, p = (1+#{>=obs})/(1+n). Exact enumeration of all 2^14 flips gives p = 3.662e-4; the library Monte-Carlo (200,000 flips) gives 3.750e-4, SE 4.3e-5, so they agree. The floor is 1/(n+1). OK. |
| Seed handling | Seed spawned into bootstrap and sign-flip streams; `None` is refused. OK. |
| Lineage t-interval and exact sign-flip | df = K-1 interval equals my scipy construction (diff 0). The exact p over 2^7 equals my enumeration (0.046875). OK. |
| ceil(0.7K) rule | Equals `math.ceil(0.7K)` for K = 5, 7, 10, 15, 20 (4, 5, 7, 11, 14), without float artefacts. "Generalises" needs t lower bound > 0 and k >= ceil(0.7K). OK. |
| Claim tiers | Tiers are nested (passes and estimate >= SESOI for "improves"; CI lower bound >= SESOI for "meaningfully"). Boundary cases behave as in section 6.4. OK. |
| SESOI_C1b | Matches the formula (0.01 x mean E_random(U minus CE_fit) / mean E_random(U)). OK. |
| Gatekeeping | `fixed_sequence` stops at the first failure and marks later contrasts "estimate only". Pending handling is correct. See fix F1. |
| >=5-positive NA rule | A line with 4 positives is NA; one with 5 is estimable. Genes need >=5 positive and >=5 negative lines. The cell rule is NA below 50% of lines or below 30 genes. NaN pairs are excluded from the counts. OK. |
| Pairing | `assert_same_pairs`, `_align` and `paired_line_delta` raise `PairingError` on differing line sets, gene sets, n, n_pos or estimability, and on a non-finite score at a scored pair. I verified row, column and relabel mismatches. OK. |
| C3 | The permutation null is within lineage, per gene, with NaN pairs fixed, 2,000 permutations and the 97.5th percentile of P1 alone. The bootstrap drops per-replicate genes with 0 positives or 0 negatives and reports the mean number dropped. Matches section 6.2 and 6.3. See fix F3. |
| Headroom | Shared-resample ratio; reportable only if the denominator CI excludes 0. OK. |

## (4) Re-run

- Test suite: `pytest tests/test_metrics.py tests/test_metrics_binary.py` gives **56 passed** in 10 s, as reported.
- Short bootstrap calibration: I re-ran the builder's simulator (`one_sim`, imported, not `main`, so `verify_metrics.json` was not overwritten) with a different seed (SeedSequence(777)), 150 datasets per regime, 4 workers, production B = 10,000 and 10,000 flips. Runtime was 72 s.

| regime | tau | gap | boot CI covers conditional | boot CI covers population | t-interval covers population | sign-flip reject | "generalises" |
|---|---|---|---|---|---|---|---|
| all genes | 0 | 0 | 0.933 | 0.933 | 0.967 | 0.027 | 0.007 |
| all genes | 0.01 | 0 | 0.933 | 0.380 | 0.953 | 0.253 | 0.020 |
| all genes | 0.01 | 0.01 | 0.933 | 0.387 | 0.953 | 0.967 | 0.527 |
| no CE | 0 | 0 | 0.940 | 0.940 | 0.953 | 0.033 | 0.033 |
| no CE | 0.01 | 0 | 0.940 | 0.467 | 0.940 | 0.260 | 0.027 |
| no CE | 0.01 | 0.01 | 0.940 | 0.467 | 0.940 | 0.967 | 0.560 |

MC SE is about 0.02. Every figure is within about 1 SE of the builder's reported 500-dataset values: conditional coverage 0.934-0.956, t-interval coverage 0.942-0.964, sign-flip reject 0.016-0.020 at tau = 0 and 0.20-0.24 at tau = 0.01, "generalises" about 52% in the positive control. I also checked that the README and `verify_metrics.json` numbers match each other (they do).

Interpretation (the builder's reading is correct): the stratified CI and the sign-flip p are conditional on the 7 lineages. With a lineage random effect, the sign-flip rejects a null about 20-26% of the time, so those do not support a transfer claim. Only the t-interval with the 5/7 rule does (2-3% false-pass rate at tau = 0.01 and no gap). The documents say this. The synthetic design is limited: s is about 0.02, at the low end of section 6.4's range, and the two scorers are exchangeable. This was disclosed and is acceptable for a negative control.

## Must-fix and should-fix

**Must-fix (small; before the library is relied on for T1)**
- F1. `fixed_sequence` does not enforce `GATE_ORDER`. It accepts `[("C5", True), ("C1", True)]` and reports both as confirmed. Add an assertion that the contrast names are a prefix of `GATE_ORDER` in order (names and order), or at least an unknown-name check. As written, the frozen chain C1 -> C1b -> C3 -> C4 -> C5 is only a convention of the caller.
- F2. `confirmatory_contrast` defaults to `sesoi=0.01`. C1b needs `sesoi_c1b(...)`. Passing the wrong value silently produces an incorrect claim tier. Make `sesoi` a required argument, or add a `contrast` argument ("C1", "C1b", "C4", "C5") that rejects C1b without an explicit SESOI.

**Should-fix**
- F3. There is no function that returns the C3 decision. It is the conjunction of the gene_bootstrap CI lower bound > 0 (P1 - B3) and the P1 observed value > the 97.5th percentile from `c3_permutation_null`. The pieces are correct but the conjunction is left to the caller, as is the feed into `fixed_sequence`. A small `c3_decision()` would remove that room for error.
- F4. `claim_tier` and `confirmatory_contrast.passed` do not know about the gatekeeping state. A C4 that passes its test after a C3 failure would still be labelled "detectable improvement". The wording must be downgraded to "estimate only" by the caller, using `fixed_sequence`. Document this, or take the sequence status as an input.

**Notes (no change required)**
- The grid is `np.round(linspace(-1,0,21),2)`, whereas the MATLAB script uses the unrounded linspace. Only values that fall exactly on a threshold could differ (none in practice, since DepMap effects are continuous).
- The C3 line bootstrap is slow in the worst case (extrapolated about 20 min for B = 10,000 on 2,866 genes; the permutation null about 3 min). Plan the compute for it.
- Pooled tAP has no NA threshold (section 5.2 gives none). Per-lineage "lift" is a ratio of means; both are disclosed by the builder.
- `README.md` says calibration ran on 8 workers; fine as history (the CLAUDE.md-style limit applies to this review only).

## Result of the check on label safety

Nothing in `metrics.py` or the scripts reads from disk, and no label path appears in them. My scripts used only synthetic arrays. No label file was opened.
