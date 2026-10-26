> **PI note (2026-10-26):** this review quotes all-line aggregate label statistics. Reading it is an aggregate exposure (see `notes/label_exposure.md`); agents authoring or approving a priori lists should avoid it.

# Design review, round 2: evaluation statistics

Reviewer: biostatistics agent (the round-1 statistics reviewer; not the author). Reviewed: DESIGN.md v1 (commit d84a593), together with notes/design_review_r1/PI_decisions.md, the §17 revision history and Appendix D.

**What I computed.** Only synthetic checks of the tAP formula and simple binomial and SE arithmetic. I read no CRISPR label values and scored no model. My scratch scripts are not project files.

**Verdict: REVISE, one new must-fix (NM1).** It is cheap to fix. Once it is done, I have no statistics objection to freezing at M1, and I do not need a round 3. The should-fix items below can be settled by the author without another review.

---

## 1. Round-1 must-fixes: verification

| R1 item | v1 location | Status |
|---|---|---|
| MF1 tie-aware AP | §6.1, §12 | **Adequate.** The formula is correct (see §1.1). The null is E[AP_random], and sklearn is kept as a sensitivity row. The known-answer tests include the jitter-invariance test. |
| MF2 C3 on noise | §6.2, §6.3 | **Adequate.** It is defined as ΔtAP vs E[AP_random], paired vs B3, with a within-lineage permutation null. There is one threshold inconsistency (SF-a). |
| MF3 common essentials | §5.2, §5.3, §6.3 | **Adequate.** C1b uses CE_fit from fitting lines and the same label. The ≥5-positive rule and the NA cell rule (<50% of lines or <30 genes estimable) are good. |
| MF4 multiplicity | §6.3 | **Adequate.** Fixed sequence, one-sided α = 0.025, matching the 95% CI lower bound, and C2 is estimate-only. The order itself is a PI choice; see SF-c and R2-1. |
| MF5 transfer | §6.2, §6.4 | **Adequate.** The lineage-equal-weighted primary and the stratified bootstrap are consistent with each other. The claim is explicitly conditional, there is a lineage t-interval (df = 6) with an exact 2^7 sign-flip, and the ≥5/7 rule. The minimum exact p is now 1/128 = 0.0078, so a lineage-level claim is attainable. |
| MF6 test-look sealing | §10.2, §4.2, W1 | **Adequate.** T1 reveals only scalars. Per-lineage, per-gene, attribution and failure-account outputs are sealed (hash in the ledger) until P2 is frozen. Validation looks are scalar-only for every workstream. The "k-of-7" count and the lineage t-interval disclose a little lineage-level information at T1. That is acceptable, because they are pre-declared transfer scalars. |

### 1.1 Spot-check of the tAP definition (§6.1)
The definition is **mathematically correct**:
- The j-th positive of a tied block (n_g items, p_g positives, preceded by N items and P positives) sits at in-block rank r with probability C(r−1, j−1)·C(n_g−r, p_g−j)/C(n_g, p_g).
- Its precision there is (P+j)/(N+r).
- Summing over blocks, j and r, then dividing by the total number of positives, gives E[AP] under uniformly random within-block order. This is the uninterpolated AP and equals sklearn's AP when there are no ties.

**Synthetic checks:**
- Against brute-force enumeration of all within-block orderings on 300 random toy cases (n ≤ 9, three score levels): max |difference| = 1.6e-15.
- On a 2,866-gene near-binary case (top block 150, 320 positives): the formula gives 0.29387, Monte Carlo over 2,000 jitters gives 0.29396, and sklearn step AP gives 0.24900.
- E[AP_random](192, 5) = 0.0507, matching the round-1 simulation.

**Implementation notes for `metrics.py`** (not design changes):
- **Overflow.** C(2650, 250)-sized binomials overflow float. Use log-gamma, or use this equivalent O(n_g) form, which I verified identical on all the cases above (it holds because precision is linear in the hypergeometric count):

  contribution of block = Σ_{r=1..n_g} (p_g/n_g) · [P + 1 + (r−1)(p_g−1)/(n_g−1)] / (N + r), with the bracketed fraction taken as 0 when n_g = 1.

- **Cross-check.** Keep both forms and assert that they agree. Blocks with p_g = 0 contribute nothing, but they still advance N.

---

## 2. Answers to the open questions in my area

**R2-4. κ rule for switching the primary label: acceptable.**
- κ(L-prob, L-ge05) is agreement between two labels. It involves no score and no model, so choosing the label on it cannot select for performance.
- **Conditions:**
  - (i) Compute it once at M1, inside `labels.py`, on training lines and pooled over all U pairs. Log it in the validation ledger as a non-look event before any metric is computed.
  - (ii) The rule is mechanical, as now: < 0.5 means L-ge05 is primary. No PI discretion.
  - (iii) Report κ also without CE_fit, for information only. It does not enter the rule.

**R2-10. SESOI = 0.01: accept 0.01 for C1, but fix its basis wording and define C1b's SESOI now (SF-b).**
- **Rationale (ii) misreads my simulation.** The minimum-detectable-effect numbers show what is *detectable*, not what is *meaningful*.
- **Basis I would write instead:** "0.01 absolute tAP ≈ 9% of the null E[AP_random] ≈ 0.115 on all genes, about one additional correctly top-ranked positive per ~100 positives. Fixed by convention before any scoring."
- **Updated precision for the new design.** For an illustrative 7-lineage tier of about 350 lines (sizes 15–126), the equal-weighted SE is ≈ 0.066·s, against ≈ 0.053·s line-weighted, where s is the within-lineage SD of line ΔtAP. The MDE (80% power, one-sided α = 0.025) is ≈ 0.19·s: 0.004 at s = 0.02 and 0.009 at s = 0.05. "Improves" is therefore attainable, and "meaningfully improves" needs Δ ≳ 0.01 + 0.13·s.
- **Do not switch to a relative SESOI fixed from training at M2.** It adds a label-informed degree of freedom for no gain. Report "headroom closed" as planned.

**R2-7. TG "≥ 9 of 12 folds": accept.**
- 9 = ⌈0.7·12⌉, which matches the test-tier ⌈0.7·k⌉ rule.
- **Null pass rate:** 0.073, against 0.21 for the old 9/14 rule.
- **Power:** if a true gain is positive in a fraction p of lineages, the rule passes with probability 0.49, 0.80 and 0.97 at p = 0.7, 0.8 and 0.9. Combined with the mean ≥ 0.005 and LB > 0 criteria, this is a sensible screening bar.
- **Fold definition.** Folds with < 5 estimable lines should count as "not positive", not be dropped, so that the denominator stays 12. State this.

**R2-2. No-information score for failed line models: agree. It is the conservative intention-to-treat choice.**
- **Why not the alternatives.** Excluding failed lines rewards failing on hard lines. Falling back to the generic parent makes failure cost nothing and hides it.
- **Three additions are needed (SF-d):**
  1. **The rule applies identically to every model-building method.** That includes B7 (ftINIT time-outs or gap > 1e-4) and the H1 arm. Otherwise C6 and the within-parent contrasts are biased.
  2. **Define the per-gene value.** "All genes tied" defines tAP_line but not the value used in tAP_gene across lines. Rule: a failed line takes, for each gene, the method's **median score for that gene over its non-failed lines in the same tier**. That carries no line-specific information, and it is placed as a tie, not at the bottom.
  3. **Report two sensitivities:** per-protocol (failed lines excluded for all methods) and generic fallback. Report the failure count per method and lineage at T1 as a revealed scalar.

**R2-8. B8 reverse tie-break: yes, report it.**
- Both directions are label-free and cheap. B8 is a guard, so the conservative rule is that the claim "adds information beyond generic + expression" requires P1 to exceed **both** B8+ and B8− (and B8′).
- Do not choose the direction on test.
- If one direction must be named as "the" B8 for tables, choose it on training at M2, and fix it in FROZEN.tsv.

**R2-1. C6 (P1 vs B7): put it in the sequence directly after C1b, provided one comparator is pre-specified.**
- §6.4 requires the abstract sentence to give C1, C1b and C6 together, so C6 should carry confirmatory status.
- B7g and B7t are two comparators. Fix one now: the better of the two on validation Δ̄_lin, chosen at M3 like G\*. Report the other as an estimate.
- Proposed order: **C1 → C1b → C6 → C3 → C4 → C5**. See SF-c on whether C3 should block C4 and C5.

**R2-11. Aggregate label exposures: not acceptable as stated, because the test tier is already effectively known (NM1).**
- Once NM1 is fixed by re-seeding, the exposures listed in §11.5 are acceptable. They are all-line or old-tier aggregates (prevalences, positive counts, eligible-gene counts), with no model performance.

---

## 3. New must-fix introduced by v1

### NM1. The split seed is compromised: the "unseen" v1 test tier can be read off a committed round-1 file (§4.1 "The real draw is unseen", §7, §11.5, R2-11)
- **What is exposed.** My round-1 review (notes/design_review_r1/statistics.md, committed in fe15541) published the seed-20261026, attempt-0 draw under the old 1+2+2 allocation. `design_numbers.py::draw` is the same algorithm: `default_rng(seed)`, then `rng.permutation` of S1, S2 and S3 in the same order with the same list order. The v1 draw only takes more elements from the same permutations, so the start of each permutation is already public:
  - S1 begins Pancreas, CNS/Brain;
  - S2 begins Soft Tissue, PNS, Ovary/Fallopian Tube;
  - S3 begins Eye, Prostate, Bladder/Urinary Tract, Cervix.
- **Consequence.** Under v1's 2+2+3 test allocation, the 7 test lineages are therefore Pancreas, CNS/Brain, Soft Tissue, PNS, Eye, Prostate and Bladder/Urinary Tract. Validation then contains Ovary/Fallopian Tube and Cervix.
  - **Attempt 0 almost certainly passes.** Test expression lines total 352 ≥ 250, and every lineage has ≥ 15 lines. Myeloid is not in validation, so (d) holds. (c) is very likely met.
  - I did not run the real seed; this follows from the published permutation prefixes.
- **The register is also wrong.** It says my round-1 statistics are "not line-level test information under the 7-lineage split". In fact my per-line positive counts and per-gene eligibility counts were computed on 5 of these 7 lineages.
- **Fix (before M1):**
  1. **Replace SPLIT_SEED with a value nobody can know before the freeze, via commit-reveal.**
     - **Preferred:** the user supplies a random integer *after* the M1 freeze commit, by message. It is recorded verbatim in the split README and the ledger.
     - **Fallback:** seed = int(SHA-256(M1 freeze commit ID)[:8], 16). The freeze commit is the one that contains `make_splits.py` and FROZEN.tsv. Log every commit made in the freeze window, so that grinding the hash would be visible.
     - Keep the `20261026 + attempt` redraw scheme relative to the new base.
  2. Amend §4.1 so it states that the old seed was exposed by the round-1 review and was retired, rather than "not evaluated".
  3. Correct the §11.5 register entry for the statistics reviewer. Attach a note to notes/design_review_r1/statistics.md, without editing its content, saying that its predicted tiers refer to a retired seed.
  4. Add a §12 check that `make_splits.py` refuses SPLIT_SEED = 20261026.
- **Effect on the design.** None beyond the seed. The simulated acceptance and size distributions in design_numbers.json remain valid.

No other new must-fix issues.

---

## 4. Should-fix (author's discretion; no re-review needed)

- **SF-a. C3 permutation threshold.** §6.3 uses the 95th percentile of the permutation null, which is one-sided α = 0.05. Every other confirmatory test is at 0.025. Use the **97.5th** percentile, with ≥ 2,000 permutations so the tail is stable.
- **SF-b. SESOI for C1b and C3.** §6.4 states SESOI only "in Δ̄_lin". Without common essentials, E[AP_random] is ≈ 0.018 (all-line round-1 aggregate), so an absolute 0.01 would mean a lift gain of about 0.5, which is effectively unattainable. Fix now, either:
  - **lift-scaled** (recommended): SESOI_cell = 0.087 × E[AP_random]_cell, i.e. 0.01 relative to the all-gene null, computed on fitting lines; or
  - C1b SESOI = 0.002 absolute.

  For C3, state that no SESOI applies, or set one in mΔAP_gene units (e.g. 0.01).
- **SF-c. Gatekeeping order.** In the single chain, a C3 failure blocks C4 and C5, which belong to a different pipeline, a different look and the PROMPT's "best predictor" goal. Context specificity is the claim most likely to fail. Pre-specify one of these:
  - (i) **C1 → C1b → C6 → C4 → C5 → C3** (C3 is still estimated and reported at T1, and confirmed only if reached); or
  - (ii) a graphical procedure (Bretz 2009). After C1b passes, split α into 0.0125 for the {C6 → C3} chain and 0.0125 for the {C4 → C5} chain, with full recycling to the other chain when one chain completes.

  Either is fine; leaving it as it is will likely forfeit confirmatory status for the W7 result.
- **SF-d. Failed-line rule details** (R2-2 above): symmetric for B7 and H1, a defined per-gene value, two sensitivities, and failure counts revealed at T1.
- **SF-e. Sign-flip statistic.** State that the C1 sign-flip permutes each line's ΔtAP sign within lineage and recomputes the *lineage-equal-weighted* Δ̄_lin, with the p-value computed as (1 + #{≥ observed}) / (1 + 10,000).

---

## 5. Items needing no change
- Lineage-equal-weighted primary with line-weighted secondary, the stratified bootstrap, and BCa triggered by |skew| > 1.
- CE_fit recomputed per fitting set, and the per-cell NA rule.
- The paired-SE one-SE rule, TG harmonisation, W7's runner-up rule, and the winner's-curse reporting.
- B2/B2m/B2c from the same label; B6 with medium features; MCC at a fixed cutoff as a secondary metric.
- Bootstrap calibration with lineage random effects (§12).
- Test-ledger mechanics and the rehearsal look.
