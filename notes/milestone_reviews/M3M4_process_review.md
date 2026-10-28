# M3-M4 process review (2026-10-28, project ~1.4 days old)

Fresh-process reviewer. Read: CLAUDE.md, RESEARCH_TREE.md, LESSONS.md, the M2 review, DESIGN.md s13-14 and the A4-A11 rows of s17, the six M3/M4 reviews (word counts and verdict/blocking lines only, not the full text), the guard changelog and `git log` since 2026-10-27 09:00 (30 commits). No data or label files opened; nothing re-run. Wall times come from commit stamps and README runtimes, so treat them as +/- 1 h. Commits at 00:47 on 10-28 are batched and do not mark when work finished.

## (a) Where the wall time went (10-27 09:00 to 10-28 00:47, about 16 h)
| Block | Wall | Compute share |
|---|---|---|
| W4 code review (6.1k words), NO-GO, fixes B1-B3, A10, re-check (09:06-10:15, then to about 12:00 for fixes) | about 3 h | none |
| W2 library review (5.4k words), A9, re-check, freeze (to 12:21) | overlaps the above; about 2 h of agent time | 17 CPU-h of sweeps, minutes of wall |
| Stage A builds and sweeps (12:21-19:01) | about 6.2 h | 108.7 CPU-h; this is the only wall block that was real compute |
| Stage A scoring code, review (2.4k words), clean-tree rerun, W8 stub, W5 run and review (3k words) | about 2 h | under 1 CPU-h |
| W3 (19:06-23:31) | about 4.5 h | 24.5 CPU-h, of which 8.5 killed or stalled (H2 pseudo-run stalled on a GLPK LP for 5.7 CPU-h) |
| W3 review (2.2k words), clean-tree rerun, A11, harness extension + review + re-check (00:02-00:46), look 2 | about 2 h | the look took seconds |

Estimate: about 40% real compute, 40% review and fix rounds, 20% amendments, harness and bookkeeping. 21k words of review text were written for about 135 CPU-h of compute. No rate-limit stall is recorded in the repo, so I cannot attribute any wall time to one. The biggest avoidable waits were (1) the W3 review starting only after the run ended and (2) the full harness review for look 2. Stage A was GLPK/Gurobi-pool bound at 6.2 h against the 3.5 h projection (about 1.8x).

## (b) What caught real problems, and what did not
Caught real problems (keep):
- **W4 code review** (M3_W4_review.md B1-B3): unpooled Gurobi WLS sessions, builds that depend on worker history (F3 and the F0 GLPK repair), missing provenance. Any one of these would have invalidated or aborted the 6-h Stage A run.
- **W2 review** (M3_W2_review.md F1/B1-B3): the cysteine default contradicted the frozen text; the 2-ME bug (no cysteine uptake) was found in the A9 fix round; blind two-agent extraction agreed on 312/313 components.
- **M4 W3 review** (M4_W3_review.md): the control reading was fixed after the results were seen. It forced the W3-D1 and W3-MEV flags and withdrew the metabolite-pool candidate at the TG floor (winner's curse). This is what made the validation look honest.
- **Harness review** (M4_harness_shortlist_review.md F1): stale provenance pins would have blocked or tainted the look. It also confirmed the W1 toy look is byte-identical.
- **Non-author recomputation** in all of the above reproduced every headline number (bitwise for W3, 0 mismatches for W5). This is the core guarantee; it found no wrong headline in M3-M4.

Added cost, found nothing that changed a result:
- **Stage A scoring review** (2.4k words): numbers matched; findings were a provenance condition and a note on the H4.1 interval. The verdict was unchanged.
- **W5 edit review** (28 edits, 14 rejected): useful biology, but the arm failed TG regardless (-0.029 on 0/12 folds). Cheaper: a pre-declared negative result needs only the recompute, with the edit review run only if TG passes.
- **Harness re-check rounds and A11(g)**: the generic shortlist look needed config only, per D7, yet harness code changed, so a full review ran. The re-check also hit 3 tripwire tests that fail identically on HEAD (environment, not code).
- **Amendments A9, A10, A11** (about 3 of the 8 amendments in this window): A10 and A11(g) are procedure notes, not changes to a reportable number or frozen SHA. They break the M2 rule "amend only when a number, label path or SHA changes".
- **W5 as specified** was predictable: the POLR AND-to-OR merge was a known risk and the literal rule lost 0.026. A 1-CPU-h label-free check on lumped pseudo-reactions would have said so before the run.
- **Killed runs** (8.5 of 24.5 CPU-h in W3): the missing LP time limit was a known risk from the W4 review (B2); the W3 code did not inherit it.

## (c) Cuts for M5-M9 (leakage guarantees and non-author verification kept)
Apply now from the DESIGN s13 cut list (all are already effectively true):
1. **Item 1 (v1.3.0 sensitivity): cut.** Item 2 (M3 quantitative bounds): already cut. Item 3 (O-NGAM to one level): already moot, O-NGAM did not qualify. Item 4 (CORDA 15-line fallback): satisfied by the 10-line row (D8). Item 5 (W7 time box 8 h to 4 h): **apply**, with W7 scoped to GEM-feature plus gene-prior hybrids (B6 tuning and the 'artefact groups only' B-loo variant), since B2 (0.90) sets the bar.
2. **M5 Stage B/C is nearly empty.** No family passes TG and F2 R1 goes forward flagged "fails TG". Run Stage B only on F2 R1 and F4 on one LP configuration, and skip Stage C unless Stage B passes. This saves most of the budgeted 115 CPU-h and about 4 h of the 9 h M5 slot. Record the M5 verdict (H4.1: method > tissue) from Stage A directly.
3. **M6 factorial: reduce it.** Only W3 B-loo passes TG, and W2 medium moves about 1e-3 per line. The 24-combination panel factorial is mostly empty. Run: {G*, G*+B-loo} x {M0, M1/M2(f)} x {F2 R1, none}, which is 8 cells on validation, then freeze P1. Keep the pre-declared test contrasts (i)-(v) and T1. Drop the top-3-on-training step because there is nothing to rank.
4. **Review budget.** One non-author recompute per result the PI relies on, with a reviewer note under about 1.5k words. Run the review in parallel with the PI's next step, not in series, and pre-register stop/no-stop on its blocking list. Skip edit-level or code-level review of any arm that fails TG (W5, Stage A families); do the recompute only. Require a full code review only for new harness code (T1, T2, unseal), where one Opus review plus one dry-run re-check is the cap.
5. **W6 and W8 are the next risk, not W7.** Start the W8 SBML builds on H1+B-loo as soon as P1 is frozen. They need no labels, so run them in parallel with T1 and W6. W8 is budgeted at 150 CPU-h and 2 x 2,030 models; build the two P1 arms only and keep generic parents as demo exports unless time remains.
6. **Every long job gets an LP time limit, a CPU budget kill at 2x, and a pinned worktree** (s11.4). W3 skipped the worktree (A-8) and lost 8.5 CPU-h to a stall.
7. **Amendments only for reported numbers, labels, splits, endpoints, lists.** Use RESEARCH_TREE decision notes and the guard changelog otherwise. Do not run further hardening of guards (the M2 rule 6 still holds).

## (d) LESSONS.md (18 lines now; proposed 14)
- Merge lines 8, 9, 10, 11, 12 (exposure policy, commit-reveal, mechanisms, downloads, column contents) into two lines: exposure policy ratified before labels are touched, no aggregates in committed docs or scratch; secrets by commit-reveal, mechanisms verified by a live probe, fail closed. Move the blocked-downloads and column-contents lines out (USER_DOWNLOADS.md / DATA.md already hold them).
- Merge 2 and 13 (review once; dry-run fixes; keep the tree lean) and 3 into 7 ("every quoted number comes from a script").
- Drop 14 (pinning is now routine) and fold 16 into 17: the SBML stub exists (eb570af); keep only "read against B2/B8, not only G*".
- Merge 15 into line 1 (label-free dry run on real IDs plus the pre-freeze checklist).
- **Add three** (from this window): (i) every long job inherits the LP time limit and pinned-worktree rule from the previous review (W3 lost 8.5 CPU-h); (ii) fix the control reading and the TG comparator before looking at results, because the W3-D1 flag exists only because it was not; (iii) decide the review depth before the run: recompute only for arms that fail TG, full review for harness code and TG passers.
- Remove W5 AND-to-OR (line 18) from LESSONS and keep it in results/w5/README.md; it is a domain fact, not a process rule.

## (e) Realistic schedule from now (H = hours; wall includes review in parallel)
| Step | Wall |
|---|---|
| M5: Stage B on F2 R1/F4 with an LP time limit; verdict note; Stage C only if B passes | 3-4 h |
| M6: 8-cell validation factorial, P1 freeze, test contrasts, T1 (scalars), one harness review | 5-6 h |
| M7: W6 (40 CPU-h) and W8 P1 builds in parallel with T1; W7 4 h box; T2, unseal | 6-8 h |
| M8/M9: resource gates, manifest, REPORT.md, verify_report and citations, one non-author review | 4-5 h |
| Total | about 18-23 h, i.e. a finish at about 2.2-2.4 days elapsed |
- This is about 0.7 day over the user's 1-1.5 day target; the over-run comes from the M3-M4 review and fix rounds, not compute. The cuts above recover about 6-8 h against the DESIGN 36-42 h critical path (about 12 h of that is spent). Without them it ends at about 2.7 days.
- The first thing to settle with the user is whether 2.3 days is acceptable. If not, drop W7 entirely and report B2 as the strongest baseline with no hybrids; that ends at about 2 days.
- Risks: GLPK stalls (mitigated by the time limit), unprepared SBML export for F4 (needs its scaled bounds stored), and MEMOTE/FROG needing a separate venv (start it during T1).
