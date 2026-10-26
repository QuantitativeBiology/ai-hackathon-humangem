# M1 process review (design frozen, 2026-10-26)

Fresh-process reviewer. Inputs: PROMPT.md (AI lessons, Research tree), LESSONS.md, RESEARCH_TREE.md, notes/label_exposure.md, git log. Restricted review files and raw data were not opened. Judgements come from the record, not from re-running anything.

## 1. What should we stop doing?

1. **Stop fixing leakage rules after the exposure has happened.** Evidence: the label guard hook arrived after the design author printed ~10 values (fe15541 -> 8fa72d5). The round-1 stats review published the split draw, which forced a re-seed. The rule "exposed agents must not author or approve lists" was then loosened by the PI (register note, 2026-10-26, F-4). The loosening is defensible, but it was made after the fact by the party it benefits. Do instead: write the exposure policy, with its disqualification threshold, before any agent touches labels. Have a non-PI agent or the user ratify any later change.
2. **Stop letting reviewers and builders write to shared, committed files without a content check.** Evidence: label copies (ge_U.parquet plus 11 derived files) and six label-reading scripts were left in scratch, and the PI swept them by hand. The seed slipped into a register entry. Aggregates were written into reviews and the register, so later "unexposed" agents became exposed (round-3 verifier, 073e12c, 936042d "redact aggregates"). Do instead: rely on the planned content-based scanner and a pre-commit grep for the seed and label patterns. Do not rely on agent discipline or human sweeps.
3. **Stop running a full-size review round when the previous one found little.** Evidence: round 2 had all 7 reviewers return "revise minimally, then freeze", 5 of 7 needed no round 3, and each took 1-7 min. Round 3 still added 17 fixes plus "v3 final check" fixes, so the check had value. The panel is cheap in time but costly in exposure and register bookkeeping. Do instead: for later milestones, scale the panel to the risk (2-3 reviewers, one of them unexposed), and run the design-numbers script before review rather than after.
4. **Stop hand-copying numbers and decisions between documents.** Evidence: v0 had wrong counts (375 vs 309). The register and RESEARCH_TREE carry overlapping, partly superseded text ("superseded in part" note; D3/D4 restate the seed history). RESEARCH_TREE.md is now ~10 KB with finding-level detail that belongs in notes/. Do instead: keep the tree to status, pointers and decisions, and move finding detail into notes/ with links. Generate quoted numbers from scripts.
5. **Stop putting tool/data facts in LESSONS.md.** Evidence: 3 of 13 current rules (medium, load-once, IsDefaultEntryForModel) are already in CLAUDE.md. They dilute the process rules.
6. **Stop starting work that depends on config not yet loaded.** Evidence: the hook was edited and then inactive until /hooks reload, and two probes got through. Treat "written" as not "enforced" until a live probe is blocked.

## 2. What worked well

1. **Independent, non-author recomputation.** Evidence: the baseline review recomputed all genes with its own code (max diff 2e-12) and found the float-tie AP instability, which led to ratio_r rounding and bitwise-identical reruns. Design review r1 found the design infeasible (non-growing media, ~12k CPU-h) before any compute was spent.
2. **Mechanical guards over promises.** Evidence: the PreToolUse label guard, hardened to fail closed, and the live probes that confirmed blocks after restart. The plan to add a label-copy scanner and tripwire follows the same pattern.
3. **Commit-reveal for the split seed, plus an unexposed verifier.** Evidence: the compromised seed was retired. The new seed was supplied by the user, kept out of files until freeze commit 61b0beb, and revealed after it (5a8c218). The round-3 unexposed Opus verifier ratified or amended the exposed author's protocol choices and re-approved the a priori lists.
4. **Recording decisions (D0-D4) with dates, and PI-decision files per review round.** Evidence: the 34 and 17 itemised PI decisions let v1 and v2 be built mechanically and let a reset context resume from the tree.
5. **Cheap label-free work before the design was frozen.** Evidence: the generic-model KO baseline and the literature survey ran in parallel with the design and found real model facts (UMPS/DHODH bypass, hypoxanthine/thymidine in Ham's).
6. **Batching blocked downloads into one user request** (notes/USER_DOWNLOADS.md) instead of retrying scripted fetches.

## 3. Pruning notes (LESSONS.md, 13 -> 11)

Dropped, now enforced or documented elsewhere:
- "Set a growth medium before FBA", "load SBML once / `with model:`" and "filter IsDefaultEntryForModel": all in CLAUDE.md, and the first is enforced in code by `set_hams_medium()`. These are tool facts, not process rules.
- "Round/clip KO ratios": enforced in the pipeline (ratio_r rounding) and by tie-aware expected AP in the metrics design. Kept only as part of the merged independent-recompute rule.
- "Hooks edited mid-session don't load": kept, but folded into the guard rule.
- "Reviewers must never run the real seed" and "Keep the secret out of every file": merged into one commit-reveal rule. Partly enforced after the reveal, since the seed is now public in D4. The reveal makes the grep part moot.
- "Reading a file that quotes label aggregates is an exposure" and the register-redaction practice: merged into one rule. Partly enforced once the scanner exists.
- The Model.csv swapped-columns note is a data fact. Move it to DATA.md or CLAUDE.md; it is dropped from LESSONS.

## 4. Proposed pruned LESSONS.md

```
# Lessons (process rules)

One line per rule, with evidence. Pruned at each milestone (M1: 13 -> 11; tool/data facts moved to CLAUDE.md/DATA.md, secret-seed rule retired after reveal).

- Have a non-author agent recompute key results with an independent code path before relying on them: it found the float-tie AUPRC instability (up to 0.018 AP between identical reruns) that the builder's 36 passing checks missed; round/clip ratios before ranking. (Baseline review, 2026-10-26)
- Run a specialised parallel review panel before building, scaled to risk: round 1 found the design infeasible (non-growing media, ~12k CPU-h) and two metric biases before any compute; round 2 found little, so size later panels down. (Design review r1/r2, 2026-10-26)
- Generate every number quoted in design/report docs from a script; hand-copied v0 counts were wrong (Nusinow∩CRISPR 375 vs 309). (Data review r1, 2026-10-26)
- Fix the leakage/exposure policy, including who is disqualified from authoring a priori lists, before any agent touches labels; do not loosen it afterwards without non-PI ratification. (Label guard and F-4, 2026-10-26)
- Enforce with mechanisms, not discipline: a hook edited mid-session is inactive until /hooks reload, so verify with a live probe that it blocks; it must fail closed. (Label guard, 2026-10-26)
- Keep label values and label aggregates out of committed process docs (registers, reviews, notes): reading them exposes later "unexposed" agents; scan content, not filenames. (Round-3 verifier, 2026-10-26)
- Reviewers get simulated seeds only and must not write label copies or label-derived scripts to shared scratch; the PI sweeps the scratchpad at each milestone until the scanner covers it. (Design review r2, AI review r2, 2026-10-26)
- Use commit-reveal for any secret that gates evaluation (split seed): keep it out of files, memory and prompts until the freeze commit, and grep before each commit. The PI's slip required amending a commit. (Seed handling, 2026-10-26)
- Batch blocked downloads (depmap.org, vendor pages, journal supplements block scripts) into one user request (notes/USER_DOWNLOADS.md) instead of retrying. (Data agent, 2026-10-26)
- Check column contents, not names: DepMap 26Q3 Model.csv OnboardedMedia/FormulationID are swapped. Record such quirks in DATA.md. (Data survey, 2026-10-26)
- Keep RESEARCH_TREE.md to status, decisions and pointers; put finding detail and review text in notes/ so the tree stays a lean resume point. (M1 process review, RESEARCH_TREE ~10 KB)
```
Count: 11 rules. (The swapped-columns rule is kept as a rule to check contents, with the fact itself moved to DATA.md.)
