# Lessons (process rules)

One line per rule, with evidence. Pruned at each milestone (M1: 13 -> 11; M2: 21 -> 13; M3-M4: 18 -> 14, per notes/milestone_reviews/M3M4_process_review.md; tool facts live in CLAUDE.md, data quirks in DATA.md, download requests in notes/USER_DOWNLOADS.md, domain findings in results/*/README.md).

- Have a non-author agent recompute key results with an independent code path, and dry-run every evaluate script label-free on real IDs before freezing it; before a frozen look, give tests stable names, commit every file provenance pins and run --preflight-only. The baseline review found float-tie AUPRC instability (up to 0.018) that 36 passing checks missed; the harness review found stale provenance pins. (Baseline review; M2 W1 review; M4 harness review)
- Decide review depth before the run: recompute only for arms that fail TG; full review for harness code and TG passers; review once, then dry-run the fixes; run reviews in parallel with the next step. The W5 edit review and the Stage A scoring review changed no result. (Design r1/r2; M2 and M3-M4 process reviews)
- Generate every number quoted in design/report docs from a script; hand-copied v0 counts were wrong (Nusinow∩CRISPR 375 vs 309). (Data review r1)
- Fix the leakage/exposure policy before any agent touches labels; changes need user ratification (D5). Keep label values and aggregates out of committed process docs and the shared scratch; scan content, not filenames. Use commit-reveal for secrets that gate evaluation (the seed-38 slip). (F-4; M1 infra review)
- Enforce with mechanisms, not discipline; verify each with a live probe; fail closed. A hook edited mid-session is inactive until /hooks reload. (Label guard, 2026-10-26)
- Fix the control reading and the TG comparator before looking at results: W3-D1 exists only because the matched-call-rate reading was settled after training aggregates were seen. (M4 W3 review)
- Every long job inherits the LP time limit (300 s + one retry), the 2× CPU kill and the pinned clean worktree; W3 skipped them and lost 8.5 of 24.5 CPU-h to a GLPK stall. (A10; M3-M4 process review)
- Check column contents, not names, and record such quirks in DATA.md (e.g. swapped Model.csv medium columns). (Data survey)
- Keep RESEARCH_TREE.md to status, decisions and pointers; amend DESIGN only when a reported number, label, split, endpoint, list or frozen SHA changes; procedure notes go to the tree or the guard changelog. (M1–M4 process reviews)
- Read results against the strongest label-free and label-prior baselines, not only the generic GEM: B2 (gene prior) is 0.90 tAP vs 0.30 for G*, and 42% of positive pairs are outside the addressable set. (W1, M2)
- Stub every PROMPT.md deliverable at the first milestone that has a model, so gaps show early. (M2 process review)
- Biomass artefacts dominate generic false positives: removing the cofactor pool gained +0.032 on validation, while every expression-pruning family lost tAP. Look for structural false positives before adding data. (W3 look 2; Stage A)
- Agents must never delete files they did not create; the W2 agent removed ~37 shared scratch files. (M3)
