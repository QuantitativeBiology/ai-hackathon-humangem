# M9 review of REPORT.md (non-author), 2026-10-28
Reviewer: non-author Claude Opus 5.5 subagent. Object: REPORT.md at 3c09d20 (working tree identical; README.md also unchanged since 3c09d20).
Probe: `head -c 50` on the raw CRISPR gene-effect file was BLOCKED by the label guard. I read no raw label file and nothing in data/locked/; only committed results/ tables, DESIGN, RESEARCH_TREE, LITERATURE, README and notes.

## Verdict: APPROVE WITH FIXES
I found no numerical error. The gatekeeping outcome (C1 and C1b confirmed, C3 failed, C4 and C5 estimates only) is reported correctly. Fixes F1–F8 change the claims and their wording and should be made before release. F9–F16 are completeness, consistency and clarity fixes.

## Checks run
- `.venv/bin/python scripts/report_build.py --check --require-filled`: "markers checked: 530; issues: 0; placeholders open: 0", exit 0.
- I spot-checked 22 numbers by opening their source files; all match. Abstract: val tAP 0.302/0.295 (results/w1/val_look_W1_M2_baseline.tsv); C1 +0.0497 [0.0482, 0.0513] p 9.999e-5; C1b +0.0374 [0.0351, 0.0397], SESOI_C1b 0.005723; t-interval [0.0474, 0.0521]; 7/7; exact p 0.0078 (test_look_T1_gates.json); C4 +0.602 [0.598, 0.605]; C5 +0.0138 [0.0116, 0.0167] (T2 gates); tAP 0.347/0.297 and headroom 8.3% [8.1, 8.6] (test_look_T1.tsv); 48 models, 2,114 lines, 753, 4,228, 93.3 MB (build_summary.json); split 1,361/394/380/170, 30 rare lines (split_meta.json, split_summary.tsv). Unsealed: B8+ 0.418, B2 0.895, H2−H1 −0.0084 [−0.0089, −0.0079] (w1_benchmark.tsv); Breast and Pleura rows (per_lineage.tsv); failure classes, e.g. complex 48,946 (48.8%), and upper bound 54.7% (failure_account.tsv); CN≥4 +0.0054 [−0.0129, 0.0235] (strata.tsv); (ii) +0.0127 (attribution.tsv); L-ge05/L-ge10/U∩ (sensitivities.tsv); the worst-cell Entrez IDs match the 20 symbols. Also: MCC/F1/precision/recall in §4.9; C3 −0.0114 and mΔAP −0.0002 vs 0.000117; W2-SHUF +0.0005 and 97%; W4 F0 −0.0080 [−0.0103, −0.0053] 0/10.
- Citations: all 13 numbered references match LITERATURE.md (authors, titles, DOIs). Support problems are in F7. Machado 2016 is declared as outside LITERATURE.md.

## Claim-language and gatekeeping fixes (must fix)
- **F1 · C3 overclaims an absence.** DESIGN §6.4 says C3 has no SESOI and takes only "detectable" wording. The report also contradicts itself: §4.10 says "25 above the null". Change:
  - Abstract l.19, "so it carries no line-specific information beyond expression" → "so P1 shows no detectable line-specific signal beyond expression".
  - §4.8 l.282, "P1 therefore carries **no line-specific information beyond expression**" → "P1 therefore shows **no detectable line-specific information beyond expression**".
  - §7 l.496: apply the same softening to "even that information is not line-specific".
- **F2 · C3 interval excludes its own estimate.** The estimate is −0.01142 and the CI is [−0.01830, −0.01149]; the expressed-pairs version does the same. §4.8 l.281 gives a "presumably" cause that has not been checked. Gene dropping is small: 0.69 genes on average out of 963. Duplicated lines under line resampling, which inflate B3's per-gene AP, are an equally likely cause. Replace the bullet with: "The percentile interval excludes the point estimate. The cause (gene dropping or line duplication in the bootstrap) is not established. Both C3 criteria fail independently (CI and permutation), so the verdict does not depend on it."
- **F3 · Flags do not travel with the confirmed claim.** W3-D1, W3-MEV, W2-SHUF and W2-MOPEN appear only in the abstract's Selection paragraph, §4.7 and §7. Add "(flags W3-D1, W3-MEV, W2-SHUF, W2-MOPEN)" in four places: the abstract C1 bullet (l.16), the §4.8 C1 and C1b rows in the status column, the Conclusion (l.29), and the §4.9 lead sentence.
- **F4 · W3-D1 understated.** Neither §4.5 l.226 nor §7 l.500 says that, under the control reading first used (symmetric), the W3 shortlist was empty (results/w3/README.md D-1). In other words, P1's main component qualified only under a reading fixed after training aggregates were seen. Add to §4.5 W3-D1: "Under the symmetric control used in the first scoring run, no W3 mode, this one included, qualified."
- **F5 · "Generalises" is justified by half its rule.** Abstract l.18 reads "Δ > 0 in 7/7 lineages, so P1 'generalises to new tissues'". Change it to "The lineage t-interval [0.0474, 0.0521] excludes 0 and Δ > 0 in 7/7 lineages, so P1 'generalises to new tissues', conditional on these 7 lineages." Also move "conditional on the 7 test lineages" into the C1 headline bullet (l.16), as §6.4(1) requires ("the headline says so").
- **F6 · C4/C5 read as tests.**
  - Abstract l.20, "P2 beats G\* by +0.602 … and the pre-registered B6 by +0.0138" → "P2 (non-GEM; it says nothing about GEMs, D12) exceeds G\* by an estimated +0.602 … and B6 by +0.0138".
  - §4.8 table: replace the C4/C5 p entries "< 0.001" with "— (not tested)". §6.3 says contrasts after a failure are point estimates and CIs only.
- **F7 · Citations.**
  - §3.6 l.114, "the failure mode reported by [7]": Opdam [7] attributes reaction-content variance to the expression threshold. The method-over-tissue clustering evidence is Richelle 2019a [8] (LITERATURE §3.2). Cite "[7, 8]" and add ref 8.
  - §4.4 l.213, "This direction agrees with the medium effects reported in DepMap [10]": [10] reports medium–dependency associations, with ASNS as the example. Hypoxanthine rescue of purine synthesis is Rossiter 2021 [11]. Rewrite as "Medium-dependent metabolic essentiality is reported in DepMap [10] and in physiological-medium screens, including purine synthesis rescued by hypoxanthine [11]", and add ref 11.
  - §4.1 l.170, "published 0.25–0.44 range (LITERATURE §2.1)": the §2.1 table gives 0.27–0.44, and 0.25–0.45 is LITERATURE §0 item 2. Make the range and its source agree.
  - Ref 92 is marked "verify against the Springer record" in LITERATURE.md; it needs that check before release.
  - CORDA (§3.6) has no reference.
- **F8 · The non-author recomputation claim is too broad.** §3.13 l.137 says "Every result relied on was recomputed by a non-author agent". notes/reviews/ contains no recomputation of the T1/T2 output numbers: M6_evaluate_test_review used synthetic labels, and the unseal SHA match compares the harness with itself. It also contains no review of W8. Either have an unexposed agent recompute C1, C1b, C3, C4 and C5 through labels.py, or write: "Training and validation results were recomputed by non-author agents. The test harness was reviewed on synthetic labels. The T1/T2 outputs and the W8 resource were not independently recomputed."

## PROMPT.md coverage (gaps)
- **F9 · No "follow-ups and remaining ideas" section,** although PROMPT.md requires one in REPORT.md. Add a §11 that collects the items now scattered through the report:
  - W7-03 artefact-only biomass (fixes W3-MEV);
  - genomic-background and lineage group models (Stage C);
  - W5 UniProt removal path and the pid30 curation;
  - the Human1 FROG report;
  - pFBA traces of the 20 worst cells;
  - CN, LOH and other line-specific inputs beyond medium;
  - W3 context claims checked on validation;
  - the "top-k" unit and the optimistic 16-lineage LOLO estimate promised in §6.4. Neither was produced, so list both as deviations.
- **F10 · Genomic backgrounds are missing.** PROMPT.md asks for models "of specific genomic backgrounds" (DESIGN §1, W4 Stage C: MTAP_del, ASS1, etc.). The report never says that these were not built. Add to §4.6: "A12 skipped Stage B/C, so no lineage, genomic-background or lineage × background models were built or evaluated. Genomic data entered only as the Stage A loss constraint (≤ 0.0002)."
- **F11 · Tasks "in which contexts".** results/w3/train/context_claims.tsv and the W3 README answer this: 1,413 of 1,416 q < 0.05 claims restate the global sign, and 3 H2 reversals are training-only. The report omits it; add one bullet to §4.5. Also state in §3.5 that the task set was fixed a priori (notes/apriori_lists_final.md, A2).
- **F12 · MCC/F1 and "overall" are incomplete.**
  - §4.9 gives F1 only for P57 G\*/P1. test_look_T1.tsv also has F1 for PGR (G\* 0.363 → P1 0.379), H2 and P1(H2); the papers' line-weighted means (P57 MCC 0.326 → 0.349); and accuracy and specificity. Add them as a row or a sentence.
  - unsealed/T1/binary_threshold_grid.tsv (the §6.8 L-ge06 grid) is never cited.
  - "Overall" AUPRC: add the pooled tAP column from w1_benchmark.tsv (G\* 0.295, P1 0.346, B2 0.890) to the §4.10 table.
- **F13 · SBML resource vs PROMPT.md.** The released models are P1 per medium signature, not context-specific models. Abstract l.27: "We release 48 SBML models covering 2,114 DepMap lines" → "We release 48 SBML models of P1, one per parent and medium signature (P1 has no line-specific parameters), mapped to 2,114 DepMap lines". §5 Deviations should state that no extraction-based context-specific models are released, and why.
- **Runtime alongside accuracy (minor).** DESIGN §6.7 wants runtime per result. §8 l.533 ("B2 costs seconds… hours") gives no numbers. Add CPU per line or in total for G\*, P1, B2, B6 and P2 next to their test tAP.

## Internal consistency (REPORT ↔ DESIGN ↔ README)
- **F14 · Amendment count.** DESIGN §17 has A1–A14 (A14: the train ∪ val refit clarification). REPORT l.11 and l.41 say "A1–A13". README's status table says A1–A14, but its pitfalls section says "A1–A10". Change REPORT to A1–A14 and README pitfalls to A1–A14.
- **F15 · README vs REPORT.**
  - README headline "(7 held-out lineages, 394 lines)": C1 used 380 expression lines. Write "380 lines with expression (394 in the tier)".
  - "about 96% of the medium gain" → 97% (REPORT, shuffled_control.tsv 0.01536/0.01588).
  - "MEMOTE run on 10 models; FROG run." → "FROG run on 1 Human2 model; Human1 not run".
  - "Compute budget: ≈1–1.5 days in total" reads as achieved. REPORT §8 has 41.9 h elapsed; say "target".
  - The README "generalises" sentence needs the F5 fix.
- **§7 l.495 uses the Human2 figure (42%) for a Human1 pipeline.** P1 and G\* are Human1 (45% on training; 45.3% on test, from 54.7%). Quote both.
- **§10 l.560, "The biomass artefact fix replicated from training to test".** B-loo alone was not scored on test. Mark the inference as descriptive: P1 − medium-alone ≈ 0.347 − 0.310.

## Clarity (F16)
- §4.4 l.213 and §4.4 W2-SHUF "[−0.0000, 0.0010]": write "[−0.00003, 0.0010]".
- §6 l.455: state the call rule and label behind "100,383 false-negative and 10,934 false-positive" calls.
- §4.10 per-lineage table: "(lines)" counts expression lines (CNS 126), while §2 has 133. Say so in the header.
- §4.8 C5 shrinkage (l.310): add that validation compared against W7's B6, and test against the B6 refit on train ∪ val (D16).
- Title and abstract Background: the headline pipeline is neither context-specific nor curated. Add to the abstract C1 bullet: "P1 is generic Human1 with a biomass edit and the documented medium; no extraction, no curation".
- §4.3: label the pid30 row as "sensitivity declared in code before scoring, not in DESIGN".
