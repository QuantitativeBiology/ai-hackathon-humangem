# Discrepancies, open questions and unsure points

For the user to review. **Main point:** none of the discrepancies changes a scientific result. Every headline number used in the explainer matches its result file, and the worked example reproduced the committed P1 scores exactly. The problems are stale status text, two numbering schemes for the same cells, an incomplete account of where wall time went, and wording that overstates what "documented medium" and "biomass artefact" mean. `scripts/report_build.py --check` at HEAD (101aed1) matches 565 of 570 marked numbers in REPORT.md; the 5 mismatches are listed in A3.

Ordered by how much each could mislead a reader.

## A. Documents vs code, results or logs

1. **"Biomass artefacts" overstates what the pool removal removes.** README.md ("What P1 is") lists BH4 and CoQ among "false essential calls created by the biomass reaction itself". The follow-up F2 classes the CoQ, BH4, heme, B6, B1 and B2 genes as *real* cofactor requirements and keeps them lethal (results/f2/README.md §4). Its decomposition shows that un-calling those 20 genes carries much of P1's gain (they explain 0.020 of the 0.023 gap between P1 and F2 on training; results/f2/decomposition.tsv). Shown: these genes are essential in 15% of test lines on average (../data/gene_groups.tsv of this explainer). The explainer therefore says the edit removes calls that are rarely essential, not that it removes only artefacts.
2. **Wall-time attribution omits about 13 h of stalls.** The M3–M4 and M5–M9 process reviews (notes/milestone_reviews/) found no record of rate-limit stalls and attributed the overrun to review and fix rounds; README.md, REPORT.md §8 and RESEARCH_TREE.md D9 repeat this. The main transcript shows usage-limit stalls on day 1 22:48 → ≈ 00:20, day 2 12:22 → 19:00 (weekly limit; the W3, W5 and W8 agents launched at 12:22 sat idle until 19:02), day 2 23:52 → 00:00 and day 3 03:54 → 05:00, plus a 4.4-h wait for a user answer (day 2 03:15 → 07:41).
3. **REPORT.md no longer matches its sources at HEAD.** `report_build.py --check` reports 5 mismatches, all caused by the post-hoc follow-ups adding look and manifest rows:

   | Marker | REPORT | Source at HEAD |
   |---|---|---|
   | Validation looks used | 6 | 7 |
   | Runs | 152 | 174 |
   | CPU-h | 201.2 | 210.4 |
   | Wall-h | 28.4 | 30.0 |
   | Elapsed h | 42.9 | 46.1 |

   Either pin the runtime sources to the "project complete" commit 8c80bdc or re-run `--update` and add a follow-up section (a PI decision).
4. **Stale status text after the follow-ups.**
   - RESEARCH_TREE.md's status summary says "F2 and F4 in progress"; its own Decisions section says all four follow-ups are done.
   - results/f2/README.md says "validation look prepared, not run"; the look ran as M10_F2 (look 7/40, commit b62bd80; results/m10/val/val_look_M10_F2_decision.json).
   - README.md says "non-author review in progress" (the review is done, RESEARCH_TREE.md item 10), cites "amendments A1–A14" (A15 exists) and "6 of 40" looks (7 used), and calls the artefact-only variant "exploratory, keeps 57% of the gain" (it is now validated at 60%: +0.031 vs +0.052).
   - REPORT.md §11 lists as owed the artefact-only variant's validation, the worst-cell traces and line-specific constraints; F2, F3 and F4 have now delivered them.
   - Tree items 4, 5, 6, 7, 9 and 11 remain marked "in progress".
   - Already fixed in commit c2034de: the M6 README title, the README's missing validation gene-prior value, the CLAUDE.md growth note (Human2 74.76 vs Human1 65.57) and the earlier status paragraphs.
5. **The factorial cells are numbered two ways.** RESEARCH_TREE.md D10/D11 and `CANDIDATE_LABEL` in scripts/m6_common.py use (1) medium only, (2) biomass only, (3) both. results/m6/train/summary.json (`tg.*.label`), REPORT.md §4.7 and the decision JSON's `rule` string use (1) biomass only, (2) medium only. Selection is unaffected (one candidate was within one SE). The explainer names cells by content.
6. **"Documented culture medium" overstates line specificity.** Of 1,361 screened lines, 485 (36%) receive classic RPMI + 10% FBS by default (medium sheet missing 210, unmappable 201, imputed 74), and 1,050 (77%) share one classic-RPMI model (results/m6/scores/line_assignment.tsv). With the shuffled (97% retained) and constant-medium (+0.0015) results, "realistic medium with serum, documented where available" is the accurate description.
7. **The README's human-written disclaimer names "metabolic synthetic lethal interaction" as a goal.** PROMPT.md, DESIGN.md and REPORT.md do not mention synthetic lethality or double knockouts, and no analysis addresses it. The explainer lists it under "not validated".
8. **The M5–M9 process review says it found no record of the wrapper exit-code event.** Commit 6c6cab9's message records it.
9. **The headline p value is written three ways**: "p < 0.001" (README.md), "p < 1e-4" (RESEARCH_TREE.md D15), 9.999e-05 = 1/10,001 in results/m6/test/test_look_T1_gates.json. The precise statement is "p ≤ 1e-4, the minimum attainable with 10,000 flips".
10. **Two CIs exist for the same contrast**: C1 [0.04822, 0.05125] (gates JSON) vs attribution (i) [0.04815, 0.05120] (results/m6/test/unsealed/T1/attribution.tsv), presumably different bootstrap draws (unverified). The explainer uses the gates JSON.
11. **The line-specificity CI excludes its own point estimate**: −0.01142 vs [−0.01830, −0.01149]. REPORT.md acknowledges it; the cause is not established; the verdict does not depend on it.
12. **The generic Human1 lethal count is 149 in RESEARCH_TREE.md 3.2 and 147 in the explainer.** The difference is the gene universe (2,887 model genes vs 2,866 with a CRISPR column), confirmed by the worked example's re-run.
13. **W3 runtime**: ≈ 24.5 CPU-h in results/w3/README.md (16.0 recorded + 8.5 killed) vs 20.9 recorded CPU-h in REPORT.md §8. Probably different scopes; not reconciled.
14. **Citation gaps in LITERATURE.md**: CORDA has no reference (the explainer cites Schultz & Qutub 2016 from memory, unverified); Machado et al. 2016 (GPRTransform) appears only in DESIGN.md A8; McSherry & Najork 2008 is flagged there as unverified against the publisher record.

## B. Open questions

1. **Open: should P1 be used with an expression tie-break?** The reference with an expression tie-break (0.418 on test) beats P1 (0.347). P1 with the same tie-break was never scored. It is the obvious practical recommendation and it is unvalidated.
2. **Open: why does removing real cofactor requirements help?** Genes that make heme, CoQ, BH4 and B-vitamin cofactors are lethal in the model but essential in only 15% of DepMap test lines on average. Possible reasons (cofactor carry-over or uptake from serum that the model omits, slow depletion within the screen's duration, or a too-rigid fixed-proportion biomass) are untested. The answer decides whether P1's step 2 is good biology or a good fit to the screen.
3. **Open: the biomass fix alone was never scored on test.** Only its validation effect (+0.032) and the descriptive difference P1 − medium-only on test exist.
4. **Open: lineage and genomic-background group models were never evaluated.** They were skipped after the context pre-check; F4 tested per-line genomic losses only. "Untested", not "refuted".
5. **Open: why are SPCS2 and SPCS3 un-called?** These signal-peptidase genes are essential in 93–94% of test lines (results/m6/test/unsealed/T1/per_gene.tsv) but leave the lethal set when the pool is removed or serum is added (../data/worked_example.json). Their link to the cofactor pool and carnitine in Human1 was not traced.
6. **Open: F3 traced only 20 CNS/brain cells, chosen by tie order.** Whether the worst errors in other lineages have the same causes is not shown.
7. **Open (process): was the user asked to accept the ≈ 2.3-day duration?** D9 projected 2.2–2.4 days against a 1–1.5-day target; the transcript's question prompts cover only F-4 and D7.
8. **Open (process): should stalls be recorded** in RESEARCH_TREE.md or the process reviews?
9. **Resolved: does P1 profit from its own mevalonate error?** No. The mevalonate core costs P1 −0.0014 on training; the super-additive gain comes from sterol genes rescued by serum cholesterol once the vitamin-D demand is gone (+0.0072 in both P1 and F2) and from an interaction residual (results/f2/decomposition.tsv).
10. **Resolved: do the follow-ups change the pipeline?** No. F1 and F4 fail; F2 is safer but costs 0.021 on validation; F3 is descriptive.

## C. Points where I am unsure

- **Subagent counts.** 56 before the explainer work (42 Opus, 14 Sonnet), from Agent tool calls in the main transcript, plus the four follow-up engineers. "About half were checkers" comes from task-description keywords, not from reading every brief.
- **Stall durations** are inferred from usage-limit messages and timestamp gaps. The weekly reset at "7pm Europe/Lisbon" equals 19:00 UTC after the 25 October DST change, consistent with the 19:01–19:02 resumes.
- **Gene-group essential fractions** (../data/gene_groups.tsv) use lethal sets from the classic-RPMI + serum model, which 1,050 of 1,361 lines share, and test-line fractions pooled over all 380 test lines whatever their medium. Lines on DMEM or IMDM have slightly different lethal sets.
- **The cut-off.** Lists use e ≥ 0.99 (ratio ≤ 0.01); the design's call rule is ratio < 0.01. No gene sits at exactly 0.01 in the four sweeps of the worked example; I did not check all 19 serum media.
- **"Subagents never committed"** rests on their READMEs and commit authorship, not a full audit.
- **The origin of the biomass idea** was checked: DESIGN v0 (b105f88) already had leave-one-out of the 13 biomass species as a diagnostic; the pool groups were added in the frozen v3 (61b0beb); the artefact rationale came from design review round 1 (notes/design_review_r1/algorithms_compute.md).
