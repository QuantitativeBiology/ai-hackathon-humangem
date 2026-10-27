# Lessons (process rules)

One line per rule, with evidence. Pruned at each milestone (M1: 13 -> 11; M2: 21 -> 13, per notes/milestone_reviews/M2_process_review.md; tool facts moved to CLAUDE.md, data quirks live in DATA.md).

- Have a non-author agent recompute key results with an independent code path, and dry-run every evaluate script label-free on real IDs before freezing it. The baseline review found float-tie AUPRC instability (up to 0.018) that 36 passing checks missed. The M2 dry run found steps after the ledger 'started' row and frozen inputs that toy tests could not show. (Baseline review 2026-10-26; M2 W1 review 2026-10-27)
- Size review panels to risk; review once, then dry-run the fixes rather than writing a second review round. Round 1 found the design infeasible before any compute, and round 2 still caught the compromised seed. The second M2 round ended with three PI actions that the pre-flight already enforced. (Design r1/r2; M2 process review)
- Generate every number quoted in design/report docs from a script; hand-copied v0 counts were wrong (Nusinow∩CRISPR 375 vs 309). (Data review r1, 2026-10-26)
- Fix the leakage/exposure policy before any agent touches labels; changes need user ratification (D5). Keep label values and aggregates out of committed process docs and the shared scratch: reading them exposes later "unexposed" agents. Scan content, not filenames. (F-4; round-3 verifier; design review r2)
- Use commit-reveal for any secret that gates evaluation, and reveal it only when the reviewed code that uses it is committed. The PI's slip with seed 38 required amending a commit and left a window in which the tiers were computable. (Seed handling; M1 infra review)
- Enforce with mechanisms, not discipline, and verify each one with a live probe: a hook edited mid-session is inactive until /hooks reload, and it must fail closed. (Label guard, 2026-10-26)
- Batch blocked downloads into one user request (notes/USER_DOWNLOADS.md) instead of retrying. (Data agent, 2026-10-26)
- Check column contents, not names; record such quirks in DATA.md (e.g. swapped Model.csv medium columns). (Data survey, 2026-10-26)
- Keep RESEARCH_TREE.md to status, decisions and pointers; finding detail goes in notes/ or results/*/README.md. (M1, M2 process reviews)
- Pin documents by hash; never keep byte copies to preserve a pin. Example: results/w4/README_ratified_A4.md. (M2 process review)
- Pre-freeze checklist for a frozen look: give tests stable names before any full run, because the test-run guard remembers WIP names; commit or snapshot every workstream's code, because provenance refuses any untracked src/scripts/tests file; run --preflight-only. (A5; W1 look, 2026-10-27)
- Stub every PROMPT.md deliverable (SBML export, REPORT.md, auto-research loop) at the first milestone that has a model, so gaps show early. (M2 process review: none started after 11 h)
- Read results against the strongest label-free and label-prior baselines, not only the generic GEM: B2 (gene prior) is 0.90 tAP vs 0.30 for G*, and 42% of positive pairs are outside the addressable set. (W1, M2)
