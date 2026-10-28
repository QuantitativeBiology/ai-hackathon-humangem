# M5-M9 process review (2026-10-28, 00:46 to 08:24, 7.6 h; project about 2.3 days old)
Fresh-process reviewer. Read: CLAUDE.md, RESEARCH_TREE.md (D9-D16), LESSONS.md, the M3-M4 review, the M5/M6/M7/M9 reviews (headline and finding lines), the guard changelog, `git log` since 00:45 (84 commits). No data or label files opened, nothing re-run. Times come from commit stamps (+/- 30 min). I found no record of the wrapper exit-code event in the repo, so item (c)3 rests on the briefing.
## (a) Where the wall time went (7.6 h; the pre-registered looks took seconds each)
| Block | Wall | Note |
|---|---|---|
| M5 pre-check (A12 text, build, 3.2 CPU-h, recompute) 00:50-01:31 | 0.7 h | replaced Stage B/C (about 115 CPU-h budgeted) |
| M6 factorial: training, clean-tree rerun, spec review, A13, look 01:47-02:26 | 0.8 h | harness rule A13 plus F-A1..F-B2 fixes |
| W6 build, run (9.8 CPU-h), recompute 02:09-03:48 | 1.6 h, parallel | null result |
| evaluate_test.py: write, review, F1-F11 fixes, re-check 03:07-06:10 | about 3 h, the largest block | the only new harness; paid for itself (b) |
| W7 loop (2.4 h of the 4 h box, 8.9 CPU-h), review, 2 looks 02:40-05:31 | about 2.8 h, overlaps | ended with a non-GEM P2 |
| W8 resource: build, 40 then 48 models, MEMOTE, FROG, rebuild 05:03-08:24 | about 3 h, parallel | the longest tail; the final 2 h only polished (clean-commit rebuild, FROG) |
| Rehearsal, T1, T2, unseal 06:25-07:25 | 1 h | includes the aborted attempt |
| REPORT draft 05:59, fill 07:34, review, fixes 07:46 | about 2 h, parallel | 570 checked numbers |
Compute was about 25 CPU-h for the window; the time was review, fix rounds and waiting for the one-at-a-time ledger steps. Elapsed total is about 2.3 days against the M3-M4 forecast of 2.2-2.4, so that forecast was good. No rate-limit stall is recorded.
## (b) Checks that caught real problems
- **F-B1 (M6 spec review).** The D10 candidate order contradicted the DESIGN s9 simplest-first order. Fixed before any validation scoring. D11 shows it was decisive only in principle, since one candidate was within 1 SE, but it would have been an unregistered tie-break.
- **evaluate_test review F1/F11.** F1: a frozen tripwire test contradicted the frozen policy, and would have blocked the freeze. F11: CE_fit came out empty in every test-shaped dry run (synthetic validation labels at a flat rate), so C1b == C1 in dry runs and the C1b logic was never exercised. F2-F5: unseal bitwise re-derivation unpinned, refusals after the validation read, T3 not tied to T1's P1. All fixed before T1. This is the review to keep.
- **REPORT drafter: 39/40 gates.** The W8 commit said "all gates pass"; one Human2 model missed the random-KO gate by 8e-9. Caught while checking claims for the report, not by the build.
- **W7 step-7 runner-up look.** The PI noticed W7-05 beat G* but led B6 by only +0.0007 on validation (training +0.014). The second look (P2_select) applied the pre-registered B6 comparison and kept C5 honest (estimate +0.0138 at test).
- **Aborted rehearsal (PYTHONPATH leak).** The tripwire refused the worktree code at the ledger commit. No look was consumed, which is the mechanism working as designed.
- **M7 W7 review (F1/F4).** Flagged allow_dirty runs; forced the non-GEM P2 wording. **M9 review F1-F16:** no numeric error, but 8 claim fixes (C3 overclaim, flags missing from C1, "generalises" half-rule, recompute claim too broad).
- **Added cost, no result changed:** M5 recompute (kept, cheap), M7 W6 review (cheap), the 2-h REPORT polish and W8 clean rebuild.
## (c) PI process mistakes
1. Commit 6a71f7f said "all gates pass" for 39/40 models. It was fixed only in a later tree note and the 05:59 commit.
2. The M6 spec's `decision_rule` string kept the pre-correction D10 order after the order was fixed (D11 says so; REPORT.md l.580 repeats it). It was documented, not repaired, because re-pinning the frozen spec was costlier than the note.
3. A wrapper's exit code was misread (as briefed; I could not find it in the repo). Rule: read the run's own output file, never the wrapper's status.
4. PYTHONPATH exported into a ledger-committing process: aborted the first rehearsal and left a cleanup commit.
5. W7 ran with allow_dirty=True, against s11.4 and against the M3-M4 recommendation 6 (pinned clean worktree). The reviewer's bitwise reproduction rescued it.
6. Review findings came late again: evaluate_test was reviewed after being written, and F1 and F11 showed up at the end; a dry-run label-shape check would have found F11 on day one.
7. Process reviews recommended "no amendments for procedure"; A13 and A14 were still added. Both were mild.
Followed from the M3-M4 review: v1.3.0 cut; Stage B/C skipped; factorial cut to 3 cells; W7 box under 4 h; recompute-only for failed arms; W6/W8 in parallel with T1; LESSONS pruned to 14.
## (d) Proposed LESSONS.md (14 rule lines; exact final text)
```
- Have a non-author agent recompute key results with an independent code path, and dry-run every evaluate script label-free on real IDs before freezing it; before a frozen look, give tests stable names, commit every file provenance pins and run --preflight-only. Dry runs must use label shapes that exercise every branch (evaluate_test F11: empty CE_fit made C1b == C1). (Baseline, W1, harness and evaluate_test reviews)
- Decide review depth before the run: recompute only for arms that fail TG; full review for new harness code and TG passers; review once, then dry-run the fixes; run reviews in parallel with the next step. (Design r1/r2; M2 to M6 process reviews)
- Generate every number quoted in design/report docs and commit messages from a script; hand-copied counts were wrong (Nusinow 375 vs 309; "all gates pass" for 39/40). (Data review r1; W8 commit 6a71f7f)
- Fix the leakage/exposure policy before any agent touches labels; changes need user ratification (D5). Keep label values and aggregates out of committed process docs and scratch; use commit-reveal for secrets that gate evaluation. (F-4; M1 infra review)
- Enforce with mechanisms, not discipline; verify each with a live probe; fail closed. A hook edited mid-session is inactive until /hooks reload. (Label guard, 2026-10-26)
- Fix the control reading, the TG comparator, the eligible set and the candidate order (check it against DESIGN s9) before looking at results; commit scoring code before its first label-reading run. Fixing a spec string later leaves a stale text (D11). (M4 W3/W2 reviews; M6 F-B1)
- Every long job inherits the LP time limit (300 s + one retry), the 2x CPU kill and the pinned clean worktree; no allow_dirty. W3 lost 8.5 of 24.5 CPU-h to a GLPK stall; W7 ran dirty. (A10; M7 W7 review F1)
- Run looks and ledger commits from a clean main tree with no PYTHONPATH or other code-path overrides; read a run's own output, never a wrapper's exit status. (Rehearsal attempt 1)
- Record data quirks in DATA.md after checking column contents, not names (swapped Model.csv medium columns). (Data survey)
- Keep RESEARCH_TREE.md to status, decisions and pointers; amend DESIGN only when a reported number, label, split, endpoint, list or frozen SHA changes; procedure notes go to the tree or the guard changelog. (M1-M6 process reviews)
- Read results against the strongest label-free and label-prior baselines, not only the generic GEM (B2 0.90 vs G* 0.30; P1 closes 8% of the gap). Say so in the abstract. (W1; T1)
- Stub every PROMPT.md deliverable at the first milestone with a model; fill the report's claims table before any test look so wording is pre-declared. (M2; M9 review F1-F16)
- Biomass artefacts dominate generic false positives: look for structural false positives before adding data. Expression-pruning families all lost tAP. (W3 look 2; Stage A)
- Agents must never delete files they did not create. (M3)
```
## (e) Next time, differently
- Review the harness spec and dry-run label shapes (including a CE-like branch) before writing the harness, not after (F1/F11 cost about 3 h).
- Never let a commit message assert a gate result; generate it from the manifest.
- Start W8 the minute P1 freezes (it began 2.6 h later) and cap its polishing; the last 2 h bought only FROG and a clean-commit rebuild.
- Fix the spec text and the code in the same commit when a pre-registered order changes, so no stale string needs a note.
- Use one launcher script for looks that scrubs the environment (PYTHONPATH, worktree vars) and prints the run's own status.
- Pre-register the "challenger vs best baseline" step (W7 step 7, with the SE rule) as the default look plan, so it is one look, not two.
- Draft REPORT.md claims with placeholders earlier (it worked, 289 checked numbers pre-T1) and give the reviewer the claim-wording checklist (C3, flags, "generalises") up front.
- Keep forecasts: the M3-M4 schedule was accurate within 5%; do the same projection at M5 for any follow-up.
