# M3 / W4 Stage-A training-tier review (independent reviewer)

**Date:** 2026-10-27.
**Reviewed:** `results/w4/stage_a_analysis/` (README and tables) and the analysis code, which is uncommitted:
`scripts/w4_stage_a_labelfree.py`, `scripts/w4_stage_a_train.py`, `src/nextgem/w4_analysis.py` and
`tests/test_w4_analysis.py`. These run on the Stage-A outputs committed in 7b22d8b.
**Spec:** DESIGN §9 (TG), §9 W4 (Stage A, selection, method-vs-tissue test), §6.1–6.2, §11.1, §11.4, §13, and
amendments A4, A7 and A10.
**Reviewer's code and outputs:** `results/w4/review/`:
- `r13_stagea_train.py` → `r13_stagea_train.json`: training labels, aggregates only;
- `r14_variance_partition.py` → `r14_variance_partition.json`: label-free.

## 0. Probe, label exposure and compute
- **Probe.** `head -c 50` on the raw CRISPR gene-effect file was **BLOCKED** by the hook at session start.
- **Label exposure: training aggregates only.**
  - Labels were read only through `labels.train()` and `cv.lolo` / `labels.fold()`, from my own script
    `r13_stagea_train.py`. `logs/label_access.log` records 3 runs (PIDs 1173994, 1174464 and 1175510), all
    `allow train … fold-scoped training API`.
  - No `validation()` or `test()` call, nothing under `data/locked/`.
  - I wrote only aggregates (per configuration and per fold), no gene × line values.
  - The values I saw: training-panel tAP_line means, Δ̄_lin, fold means and CIs, the same aggregates as the author's
    README. I also saw the per-line tAP column of the author's `train_per_line.tsv`, as max-abs-difference
    comparisons only.
- **Other log entries.** The val/test "allow" entries in the log at 19:18 and 19:30 are another session's test
  fixtures ("deletion of a copy outside the data roots", toy roots), not this review.
- **Compute.** No solver and no Gurobi session. About 6 CPU-min for r13 (3 runs) and 5 CPU-min for r14.

## Verdicts

| # | Item | Verdict |
|---|---|---|
| 1 | Independent recomputation: F2, F3 and F5 (plus F4) nested-LOLO Δ̄ vs H1 | **PASS**: identical |
| 2 | Selection follows the frozen rules; no outer-fold leakage (incl. F5) | **PASS** (NOTE 2a, a process deviation) |
| 3a | Config gate "≥ 95% of builds pass" | **PASS as interpretation**: immaterial to every outcome |
| 3b | Gate / TG reference = generic H1 (G\*) | **PASS**: the conclusion is the same against H2 |
| 3c | "Corrected" basic interval for H4.1 | **NOTE**: the diagnosis is only partly right and the basic interval is not a valid fix; the H4.1 verdict is unaffected |
| 4 | F4: gate kept frozen and F4 reported descriptively | **PASS**, with one caveat (NOTE 4a) |
| 5 | Variance partition and permuted control | **PASS**: method sound, conclusion right (NOTEs 5a, 5b) |
| 6 | Label-copy scanner on `stage_a_analysis/` | **PASS**: 0 hits |

**GO for committing these results as the Stage-A record**, on condition C1 below (provenance), with the README
edits R1–R4.

## 1. Independent recomputation: PASS

**My own implementation** (no `metrics`, `w1` or `w4_analysis` code):
- tie-aware AP in the O(n_g) form of §6.1;
- score e = round(1 − ratio_r, 6) over U, with U genes without a parent gene scored 0;
- estimable if ≥ 5 positives; lineage-equal means; fold positive if ≥ 2 estimable lines and mean Δ > 0;
- a lineage-stratified line bootstrap (B = 2,000, SEEDS[0]).

**Nested LOLO run through the `cv.lolo` harness:**
- for each held-out panel lineage, `fit_predict` sees only `labels.fold(k)`, restricted to the other 9 lineages'
  panel lines, and picks the configuration with the best lineage-equal Δ vs H1, ties going to the simplicity
  order;
- the held-out labels are used only in the harness's `score` callback.

**Per-line tAP:** my values equal `train_per_line.tsv` for all 44 compared configurations to **2.2e-16**:
- the generic references H1 and H2;
- F2, F3 and F4 × {R1, R2, R3} × {no genomic, genomic};
- the 24 F5 grid points.

Generic tAP_line (lineage-equal): H1 0.29824, H2 0.29161.

| Family | Nested Δ̄_lin vs H1, own (author) | Own CI (author CI) | Folds positive, own (author) | Same choice in every fold | TG |
|---|---|---|---|---|---|
| F2 | −0.006709 (−0.006709), \|Δ\| 2e-18 | [−0.0090, −0.0044] ([−0.0090, −0.0042]) | 1 (1) | ✓ R1, no genomic | fail |
| F3 | −0.010215 (−0.010215) | [−0.0131, −0.0073] ([−0.0133, −0.0072]) | 1 (1) | ✓ R1, no genomic | fail |
| F5 | −0.011402 (−0.011402) | [−0.0131, −0.0097] ([−0.0131, −0.0095]) | 0 (0) | ✓ w_p 0, w_c 0, τ −1.5 (w_c 1 for Kidney) | fail |
| F4 (descriptive) | −0.002860 (−0.002860) | [−0.0041, −0.0016] | 1 (1) | ✓ R1 (R2 for Head and Neck) | fail |

- The CIs differ only through the bootstrap random stream.
- The F5 full-panel best grid point (w_p 0, w_c 0, τ −1.5) also matches.

## 2. Selection rules and leakage: PASS

**Order of steps** (`w4_stage_a_train.py`, `run`, nested block):
- **Gates first.** Candidates are the gate-passing RNA configurations: {R1, R2, R3} × {no genomic, genomic}, and
  the 24 F5 points. A family with no gate-passing configuration is scored descriptively, and its `tg_pass` is
  forced False.
- **Nested panel LOLO** with all three TG items: Δ̄_lin ≥ 0.005, CI lower bound > 0, and ≥ 8 of 10 folds positive.
  A fold with < 2 estimable lines counts as not positive.
- **≤ 2 thresholds per family.** These are kept by full-panel rank (R1, R2 everywhere). That is a training-tier
  choice for Stage B and does not enter the nested estimate.
- **Modality.** Compared on the 10 ProCan lines (6 lineages). The 8-of-10 bar cannot apply there, so it is
  descriptive and RNA is kept. This reading is reasonable.

**No outer-fold leakage.**
- `nested_lolo_select` uses only `D.loc[lin != k]` for fold k's choice.
- My harness-based implementation, where fold k's labels are structurally unavailable to the choice, reproduces
  every choice and every outer value exactly. So there is no leakage in outcome.
- For F5, the grid is fixed and label-free (A4 d), the z statistics and accepted losses are label-free, and the
  choice among the 24 points is nested exactly as for the LP families.
- CE_fit (secondary cells only) is fitted per held-out lineage through `cv.lolo`.

**NOTE 2a (process deviation, not an error).** §11.1 says "all CV goes through the `nextgem.cv.lolo` harness".
The author's nested selection instead slices a `labels.train()` frame. It is logically fold-scoped, and my
harness version gives identical results. For the record, either move the selection into `cv.lolo` or cite this
equivalence check in the README.

## 3. The author's interpretations

**(a) Config gate "≥ 95% of builds pass". PASS as an interpretation; immaterial.**
- DESIGN gives no config-level aggregation. It says only "label-blind gates come first (growth, viability tasks,
  repair size ≤ 5% of reactions)". The 95% figure mirrors §11.4's "> 5% failed lines" flag, which is a coherent
  reading.
- It changes nothing (`config_summary.tsv`):
  - every F0, F2, F3, F5 and F3b configuration has 100% of builds passing;
  - every F4 RNA and RNA+protein configuration has ≤ 47%; the F4 protein-only configurations have 100%, but they
    are not RNA candidates.
- Any threshold in (0.47, 1] gives identical gate outcomes. State this in the README (R4).

**(b) Gate / TG reference = generic H1 (G\*). PASS.**
- §9 W4 says "TG is computed against the best shortlisted generic pipeline". The W2, W3 and W5 shortlists do not
  exist yet, so G\* (the better generic parent at M2) is the natural reference.
- Because the Stage-A parent is Human2, Δ vs H1 mixes the parent difference (H2 − H1 = −0.0066 here) with the
  context effect. The author reports Δ vs H2 throughout, which is right.
- The conclusion does not depend on the reference. Against H2, the best gate-passing family is F2 at −0.0001
  (5/10 folds); the best of all, F4, is +0.0038 (10/10 folds), still below the 0.005 bar. No family passes TG
  against either reference.
- A later, stronger shortlisted generic reference can only make the TG harder, so the negative result stands
  a fortiori.

**(c) "Corrected" interval for duplicate lines in H4.1. NOTE: the diagnosis is partly right, the fix is not
valid, and the verdict is unaffected.**
- My recomputation (`r14`) reproduces the observed R² values exactly. Reactions: method 0.5884, lineage 0.0519,
  difference 0.5365. Binary KO: 0.1987, 0.0853, 0.1134.
- At the KO level, the line-bootstrap percentile interval [0.048, 0.097] excludes the point estimate. It is
  shifted down.
- A **duplicate-free** half-sampling (2 of 3 lines per lineage without replacement, 2,000 draws) is also shifted
  down: KO [0.073, 0.096], mean 0.084.
- So the shift is not only duplicate resampling. It is the small-sample upward bias of marginal R² for a 9-df
  factor (lineage) when there are fewer lines per lineage. R²(lineage) has a chance floor of roughly (G−1)/(n−1)
  at the line level, and this bias makes the frozen H4.1 comparison **conservative**: it favours lineage.
- The "basic" interval assumes a pivotal bias. At the KO level it lies entirely **above** the point estimate
  ([0.130, 0.180] against 0.113), which shows that it is not a valid correction here.
- **Recommendation (R1):**
  - report the frozen percentile interval as primary;
  - add the half-sampling interval as a duplicate-free sensitivity;
  - drop the basic interval, or label it as not interpretable;
  - state that marginal R² favours the 9-df lineage factor.
- **H4.1 holds at every level under every interval:**

| Level | Percentile lower bound | Half-sampling lower bound | Basic lower bound |
|---|---|---|---|
| Reactions | 0.501 | 0.511 | 0.538 |
| Tasks (author) | 0.128 | — | 0.171 |
| KO | 0.048 | 0.073 | 0.130 |

## 4. F4 and the frozen gate: PASS, with a caveat

**The failure is a scale artefact, confirmed label-free on all 360 F4 build variants:**
- GR fails exactly when WT growth < 1: 196 of 196 builds below 1 fail GR, and 164 of 164 at or above 1 pass;
- all 9 ER viability tasks pass in every F4 build;
- F4 WT growth is 0.75–1.63 on RNA and 2.1–16.5 on protein-only (where 100% pass).

(The README's 144/96 counts refer to a subset; the conclusion is the same.)

**F4 would not be selected even without the gate.**
- Nested vs H1: −0.00286, CI [−0.0041, −0.0016], 1/10 folds. All three TG items fail (own recomputation).
- Nested vs H2: +0.0038, 10/10 folds, CI > 0, but below the 0.005 minimum.
- So F4 cannot pass TG against either reference. Keeping the frozen gate is the right anti-forking decision.

**NOTE 4a (caveat for the PI).** §9 W4 carries "the best configuration … to M6 even if it fails TG".
- With the gate, that is **F2 R1** (−0.0067 vs H1).
- **Without the gate it would be F4 R1** (−0.0029 vs H1, the best nested value of any family).
- So the gate decides *which* failing configuration is carried, though not whether anything passes TG.
- Under the PI's decision (gate unchanged), F2 R1 is carried. The README should state this dependence explicitly
  (R2).

## 5. Variance partition and permuted control: PASS

**Method.**
- **R² definition.** Marginal PERMANOVA R² (Anderson 2001) on Jaccard dissimilarity, i.e. squared √Jaccard,
  which is Euclidean-embeddable. Lineage and line p-values come from permutations at the line level; each line
  keeps all its builds.
- **Reproduction.** I recomputed it with my own code (`r14`):
  - observed R² for method, threshold and lineage are identical;
  - the H4.2 values for F0, F2 and F3 are identical: reactions real 0.2598 / 0.2603 / 0.2420 against permuted
    0.1839 / 0.1834 / 0.1803;
  - my null 95th percentiles (lineage relabelled across the lines of the permuted-expression builds, 999 draws)
    agree with the author's to ≤ 0.001.
- **H4.2 null.** The author takes max(95th percentile of the null, observed permuted R²). With only one
  permuted-expression realisation built (SEEDS[0]), this is a reasonable and conservative reading of "the 95th
  percentile of R²(lineage | permuted)".
- **Expression reference.** Lineage R² is 0.48 on real expression and 0.31 on permuted expression, which equals
  the chance level 9/29. The permutation therefore destroys lineage structure as intended.

**Conclusion.** "Method explains more variance than tissue" (H4.1, H4.3) is right at all three levels. The
direction is robust, and it is conservative, since marginal R² favours lineage (§3c). F0, F2 and F3 are
tissue-informed (H4.2); F4 is method-dominated.

**Permuted control on tAP.** Real beats permuted for F2 (R1 and R3), F0 R3 and F3 R3, but no family reaches G\*.
Expression carries some line signal; it is not enough.

**NOTE 5a. Two H4.2 calls are borderline.**
- F3 at the KO level: real 0.1823 against the bar 0.1802 (author) or 0.1795 (mine).
- F4 at the task level: 0.3797 against 0.3868.
- Flag them as "near the bar".

**NOTE 5b. F4 is degenerate at two levels.**
- At the reaction level, all `none` builds open every reaction, so R² is NaN.
- At the KO level, the binary profile equals the generic model's.
- "Method-dominated" is right, but the README should say it holds by construction, not as a measured result.

## 6. Label-copy scanner: PASS
- `scan_label_copies.py --roots results/w4/stage_a_analysis --test-keyed`: **21 files, 0 hits**, exit 0.
- The same over `results/w4/review` (my outputs): 66 files, 0 hits.
- The author's 6 analysis tests pass (`tests/test_w4_analysis.py`, not recorded).

## Conditions and requested edits

**C1 (provenance, required before the commit).**
- `stage_a_analysis/MANIFEST.tsv` rows have dirty = True: the analysis code is uncommitted and ran with
  `--allow-dirty`.
- Commit `scripts/w4_stage_a_labelfree.py`, `scripts/w4_stage_a_train.py`, `src/nextgem/w4_analysis.py` and
  `tests/test_w4_analysis.py`.
- Re-run both scripts from a clean tree, or a pinned worktree at that commit: about 2 min, no solver.
- Confirm the outputs are byte-identical apart from the MANIFEST, then commit the outputs with the clean-tree
  MANIFEST row.

**Requested README edits (non-blocking):**
- **R1:** the H4.1 intervals, as in §3c.
- **R2:** the F4 carry dependence on the gate (NOTE 4a).
- **R3:** the cv.lolo equivalence (NOTE 2a), and the borderline or degenerate H4.2 calls (NOTEs 5a, 5b).
- **R4:** the 95% config rule is immaterial for any threshold in (0.47, 1] (§3a).

**Bottom line.** The headline is correct. I reproduced F2, F3, F4 and F5 independently; the F0 and CORDA values
are the author's and were not recomputed:
- no context family passes the Stage-A TG;
- every nested Δ̄ vs G\* is negative (F0 −0.0080, F2 −0.0067, F3 −0.0102, F4 −0.0029 descriptive, F5 −0.0114,
  CORDA −0.0127 on 10 lines);
- method dominates tissue.

**GO** to commit the Stage-A record once C1 is done.
