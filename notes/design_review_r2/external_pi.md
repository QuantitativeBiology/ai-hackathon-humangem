# Design review, round 2: external PI

Reviewer: external PI (systems biology / cancer metabolism), not an author of the design.
Reviewed: DESIGN.md v1 (commit d84a593), notes/design_review_r1/PI_decisions.md, and the change log and Appendix D.
For two of the new issues I also read the draw function in `scripts/design_numbers.py` (I did not run it), the round-1 statistics review §"Predicted split", and the time-limit comments in `data/Human-GEM-main/code/test/estimateEssentialGenes.py`.
I read no CRISPR label values and computed no draw of the real split seed.

## Verdict: revise, then freeze

All six of my round-1 must-fixes are adequately addressed (Part 1). v1 introduces **three new must-fixes** (Part 3):
- **N1:** the "unseen" test tier is in practice largely known, and label statistics on it sit in a committed review.
- **N2:** the PROMPT's "medium alone" number is contaminated by an expression-based gate.
- **N3:** the MILP failure rule handicaps the published-protocol baseline B7, which biases C6 in our favour.

All three are text-level fixes that need no new analyses. With them applied as specified, I do not need a third round; a PI check of the three edits suffices.

---

## Part 1. Round-1 must-fixes: verification

| Round-1 MF | Status in v1 | Verdict |
|---|---|---|
| MF1: generic + expression baseline | B8 (tie-break by higher expression) and B8′ (§8). "Adds beyond generic + expression requires P1 > B8" (§8). MCC at a fixed cutoff (§6.1). The B8 contrast is revealed at T1 (§10.2). tAP (§6.1) also removes the float-tie artefact behind my concern. | **Addressed.** |
| MF2: confirmatory structure, SESOI | Fixed sequence C1 → C1b → C3 → C4 → C5, one-sided α = 0.025; C2 is an estimate; SESOI 0.01 with three claim tiers; headroom closed; abstract sentence gives C1 + C1b + C6 (§6.3–6.4). | **Addressed.** See S1 (order) and S2 (SESOI scope). |
| MF3: transfer claim | Conditional wording; lineage t-interval (df = 6) plus exact sign-flip; ≥ 5 of 7 rule; "generalises" only if both hold; LOLO estimate labelled selection-contaminated (§6.4). The PI raised the test tier to 7 lineages against my advice, which I accept: it buys a usable lineage-level interval. | **Addressed.** |
| MF4: expressed pairs | Required secondary analysis (§5.4). C3 claim condition (§6.3). Part of the TG reports (W4). C3 is now paired against B3, which is stronger than I asked. | **Addressed.** |
| MF5: attribution on test | P1's full lever factorial on validation and test, with Shapley shares and "medium alone" on the generic parent, sealed until T2 (§10.1–10.2). | **Addressed in structure.** But the "medium alone" contrast is now contaminated by M1g (new must-fix N2). |
| MF6: Human1 full arm | P1(H1) on validation and test lines; within-parent contrasts at T1 (§3.1, §6.3, §13). | **Addressed.** |

My round-1 should-fix and optional dispositions (Appendix D) are acceptable. That includes the rejection of SF2: the library-class strata are adequate, and the leakage-surface argument is reasonable.

---

## Part 2. Answers to the round-2 open questions in my area

### Overall soundness
The design is now scientifically sound and proportionate, and it answers every PROMPT part:
- the baseline grid (W1);
- context models by tissue, background and their combination (W4);
- the method-vs-tissue test (W4);
- medium attribution (W2, §10);
- a priori tasks (W3);
- curation vs constraints (W5, §10);
- flux vs omics (W6);
- lineage hold-out (§4, §11);
- the auto-research loop with runtime (W7, §6.7);
- the SBML resource (W8).

Two things raise its value well beyond v0. B7 makes "better than the published per-line protocol" testable, and B8 makes "better than generic + expression" testable. Those, not "beats generic FBA", are the claims a top journal will judge.

The design is now very large (1,134 lines, ~60 frozen artefacts). The main residual risk is execution, not science. Holding the cut list (§13) and the milestone go/no-go rules is what matters.

### R2-1: should C6 (P1 vs B7) join the gatekeeping sequence?
**Yes, appended at the end: C1 → C1b → C3 → C4 → C5 → C6.**
- **Why add it at all.** C6 is the claim reviewers will care about most. Robinson 2020 and Gustafsson 2023 already showed per-line ftINIT models; "beats generic" alone is not novel. A confirmatory C6 lets the paper say "improves on the published protocol" with error control.
- **Why at the end.** Appending a hypothesis to a fixed sequence costs nothing for the earlier contrasts. Inserting C6 after C1b would let a C6 failure block C3–C5, and C3 is the PROMPT's context-specificity question. C6 can be confirmatory only if everything before it passes, so it may stay an estimate in practice. That is acceptable, and the abstract sentence already reports it.
- **Label it honestly.** If P1 contains F1 per line (which equals the B7 configuration), C6 measures "medium / tasks / curation levers on top of the published protocol". It should be labelled that way, not "a different extraction method beats ftINIT".

### R2-10: SESOI = 0.01 basis
**Accept 0.01 as an absolute value, fixed now, for C1 and C4.** A relative SESOI fixed from training at M2 would add a label-informed step for no gain.

Two corrections:
1. **Rationale (ii) conflates detectability with meaningfulness.** The MDE simulation shows that 0.01 is *detectable*. Whether it is *meaningful* rests on expert judgement. Say so.
   - Optional: REPORT.md also translates the observed Δ into an interpretable unit, such as the change in the mean number of true positives among each line's top-k genes, with k = that line's number of positives.
2. **Scope for C1b (should-fix S2).** Without common essentials, per-line prevalence is several-fold lower (round-1 statistics review), so an absolute 0.01 means a much larger relative effect. State that the SESOI tiers apply to C1 and C4 only. For C1b, either:
   - use "detectable" wording only; or
   - fix SESOI_C1b = 0.01 × [mean E[AP_random] without CE_fit ÷ mean E[AP_random] with all genes], computed on training lines at M2 before any model is scored.

### R2-3: compute of 2,940 CPU-h vs the ≈ 2,650 budget
- **The binding constraint is Gurobi wall-clock, not CPU-h.** About 1,525 MILP CPU-h on 14 sessions is ≈ 110 h, more than the whole ≈ 90 h budget. Cut MILP first.
- **B7/F1 per line:**
  - Run it on the **100-line panel, validation and test lines unconditionally**. That covers C6, H1.4 on the panel and validation, and B7 next to every validation candidate.
  - Run it on **the remaining training lines only if F1 is shortlisted after Stage A**. In that case the builds are needed for Stage B TG and are reused by the W8 resource (cache by config hash), so nothing is wasted.
  - This saves up to ~280 MILP CPU-h. With the H1 resource already cut-list item 2, the plan sits at ≈ 2,650 CPU-h and ≈ 90 h of MILP wall-clock.
- **H1.4 on training.** If F1 is not shortlisted, state that H1.4 is evaluated on panel lines only.
- The remaining blocks and the cut-list order are fine. Keep the 20% re-run reserve; it will be used.

### R2-11: do reviewers' all-line label views compromise the test tier?
**The aggregate exposures alone do not. One of them does, because of the seed (N1).** All-line prevalence distributions, the lineage ICC, the common-essential share of positives, per-library positive rates and ASNS-by-medium prevalence are:
- performance-blind;
- not gene × line information on a known test set;
- already published for ASNS (Lagziel).

Conditions:
- **No a priori lists from exposed agents.** Agents in the register must not author or approve a priori lists, candidate sets or decision thresholds. The attestation rule in §11.5 covers this; keep it.
- **SESOI simulation.** The SESOI simulation used all-line prevalence. That is power analysis, not tuning, and it is acceptable.

The real compromise is that the statistics reviewer **did evaluate the real split seed** under the v0 rules. In a committed file, they named the predicted test lineages (Pancreas, Soft Tissue, PNS, Prostate, Eye) and reported label statistics on them. Under v1's draw function, that information carries over to the new test tier. See N1.

### R2-2, R2-6, R2-7, R2-8 (briefly)
- **R2-2:** the no-information score for failed lines is right. It is conservative against the method that failed. Add a sensitivity row that excludes failed lines for all methods. N3 must hold, though, or the rule becomes a handicap for B7.
- **R2-6:** attribute the gate to context, not medium. This is fixed in N2.
- **R2-7:** "≥ 9 of 12 folds" is fine (P ≈ 0.07 under a fair coin). Stage A's panel TG has 10 lineages, so state its bar: **≥ 8 of 10 folds** (should-fix S3).
- **R2-8:** one direction suffices. B3 is already reported, and the reverse tie-break is the weaker baseline, so reporting it adds nothing.

---

## Part 3. New must-fixes introduced by v1

### N1. The test tier is not unseen: re-seed the split before freeze (§4.1, §7, §11.5)
- **Facts.**
  - `design_numbers.draw()` calls `rng = default_rng(seed)` and then `rng.permutation(strata[s])` for S1, S2 and S3 in that order. The test tier takes the first `ALLOC_TEST[s]` items of each permutation, and validation takes the next ones.
  - The strata are unchanged from v0 (same 23 lineages, same ranking).
  - With seed 20261026 + attempt 0, the permutations are therefore **identical** under v0 and v1. v1's test tier (2/2/3) is v0's test tier (1/2/2) **plus** perm_S1[1] and perm_S3[2], which were in v0's validation tier. It is a strict superset of v0's test tier.
  - The round-1 statistics review (committed in fe15541, read by the author and the PI) reimplemented §4.1 with the real seed. It named v0's test tier (Pancreas, Soft Tissue, PNS, Prostate, Eye) and validation tier (CNS/Brain, Ovary, Bladder, Cervix), and it reported label statistics on that test tier:
    - 11/192 lines with 0 positives under L-ge10 without CE;
    - per-gene eligibility counts;
    - "Pancreas and Bowel have higher prevalence and SD".
  - If attempt 0 passes v1's checks (simulated acceptance is 0.854), five of v1's seven test lineages are already named, and the other two are among four named v0-validation lineages.
  - §4.1 ("the real draw is unseen") and §11.5 ("None of these is line-level test information under the 7-lineage split, which has not been drawn") are therefore incorrect.
  - I deduced this from the code and did not run the draw. Nobody else should run it either.
- **Why it matters.** The data are aggregate, but the PI and the author would know the test tissues for the whole project. That undermines the locked-test claim and W7's "ideas must name a training lineage" guard. A top-journal reviewer, or our own red-team audit, would find it in git.
- **Change (label-blind, before freeze):**
  1. Retire 20261026. Set `SPLIT_SEED = int(SHA-256(frozen M1 commit of DESIGN.md + FROZEN.tsv)[:8], 16)`. This is deterministic and reproducible, and no one can know it before the freeze. Alternatively, the user supplies an integer at M1. Forbid computing the draw before the freeze commit, and forbid amending that commit after the split runs (the split script already refuses to overwrite).
  2. Correct §4.1 and §11.5. Add a register entry stating that the v0 seed was evaluated, which lineages were named, and which label aggregates were reported on them. Note that a named lineage may land in the new test tier by chance; that residual exposure is aggregate-only and is stated as a limitation.
  3. Record this as a pre-freeze change in §17, with the reason "exposure, label-blind".

### N2. "Medium alone" must not include the expression-based transporter gate (W2 ladder, §10.2, R2-6)
- **Problem.**
  - M1g gates Asp/Glu uptake on SLC1A1/2/3 expression, so it is a line-specific *expression* constraint.
  - M2(f) and M3 are built on M1g, so every serum and quantitative variant includes the gate.
  - The PROMPT asks for the error "attributable to medium misspecification alone". The T1 contrast "M0 → documented medium on the generic parent" will carry an expression effect if the documented medium is the shortlisted M2 or M3. The Shapley shares would then credit context information to the medium.
- **Change.**
  - Make the gate an orthogonal toggle: M1 → M2(f) → M3 without the gate, and M1g / M2g / M3g with it. This only adds generic sweeps per medium, at ≈ 0 extra CPU-h beyond the 2^3 gate states.
  - Define "medium alone" as M0 → the best **gate-free** documented variant.
  - In the factorial and Shapley attribution, count the gate in the context-constraint lever (or as its own player), never in medium.

### N3. The MILP failure rule handicaps B7 and biases C6 (§11.4, §8 B7, §12 Hart test)
- **Problem.**
  - §11.4 counts an MILP that stops at its time limit with gap > 1e-4 as a failed build, and the line then gets the no-information score.
  - The repository protocol that B7 is meant to reproduce exactly (`estimateEssentialGenes.py`: `MIP_GAP_ABS = 10`, `PROVE_ABS_GAP = 1`, `TIME_LIMIT = 900`) **accepts time-limited incumbents by design**. Its comments note that RPE1 step 1 hits the limit routinely, with the same model as at 1800 s, and that tie-break phases adopt unproven incumbents.
  - Applying our rule would turn normal B7 builds into "failed" random-level lines. That inflates C6, may push B7 over the 5% failing flag, and contradicts the one-shot Hart acceptance test, which assumes repository behaviour.
  - This is a control made weaker in our favour, which PROMPT forbids.
- **Change.**
  - B7, and F1 per line since it is the same configuration, use the **repository's acceptance semantics**: a time-limited incumbent is accepted, and its gap and status are recorded.
  - Only no incumbent, infeasibility, WT growth ≤ threshold or a solver error counts as a failure.
  - Report the share of time-limited builds per method.
  - Any stricter rule applied to our own MILP families must not be applied to B7 differently from F1.

---

## Should-fix (new or residual; author's discretion)

- **S1. C3's position blocks C4/C5.** With C1 → C1b → C3 → C4 → C5, failing C3 (P1 beating *expression* at per-gene ranking, a high bar) removes confirmatory status from C5 ("does the GEM add value over ML"). C5 is arguably the most important question for a top journal. Consider a Bretz-style graphical split after C1b: α/2 to C3, α/2 to C4 → C5, with α propagated on success, and C6 last. If the PI keeps the simple chain, state this consequence explicitly in §6.3.
- **S2. SESOI scope for C1b.** See R2-10.
- **S3. Fold-consistency bar for the 10-lineage Stage A panel.** See R2-7.
