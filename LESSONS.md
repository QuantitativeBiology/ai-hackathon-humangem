# Lessons (process rules)

One line per rule, with evidence. Pruned at each milestone (M1: 13 -> 11; M2: 21 -> 13; M3-M4: 18 -> 14; M5-M9: 17 -> 14, per notes/milestone_reviews/M5M9_process_review.md). Tool facts live in CLAUDE.md, data quirks in DATA.md, download requests in notes/USER_DOWNLOADS.md, domain findings in results/*/README.md.

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
