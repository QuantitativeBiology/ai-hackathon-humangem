# Lessons (process rules)

One line per rule, with evidence. Pruned at each milestone (M1: 13 -> 11; tool/data facts moved to CLAUDE.md/DATA.md, secret-seed rule retired after reveal).

- Have a non-author agent recompute key results with an independent code path before relying on them: it found the float-tie AUPRC instability (up to 0.018 AP between identical reruns) that the builder's 36 passing checks missed; round/clip ratios before ranking. (Baseline review, 2026-10-26)
- Run a specialised parallel review panel before building, scaled to risk: round 1 found the design infeasible (non-growing media, ~12k CPU-h) and two metric biases before any compute; round 2 was smaller and still caught the compromised split seed, so size later panels to risk rather than dropping them. (Design review r1/r2, 2026-10-26)
- Generate every number quoted in design/report docs from a script; hand-copied v0 counts were wrong (Nusinow∩CRISPR 375 vs 309). (Data review r1, 2026-10-26)
- Fix the leakage/exposure policy, including who is disqualified from authoring a priori lists, before any agent touches labels; do not loosen it afterwards without non-PI ratification. (Label guard and F-4, 2026-10-26)
- Enforce with mechanisms, not discipline: a hook edited mid-session is inactive until /hooks reload, so verify with a live probe that it blocks; it must fail closed. (Label guard, 2026-10-26)
- Keep label values and label aggregates out of committed process docs (registers, reviews, notes): reading them exposes later "unexposed" agents; scan content, not filenames. (Round-3 verifier, 2026-10-26)
- Reviewers get simulated seeds only and must not write label copies or label-derived scripts to shared scratch; the PI sweeps the scratchpad at each milestone until the scanner covers it. (Design review r2, AI review r2, 2026-10-26)
- Use commit-reveal for any secret that gates evaluation (split seed): keep it out of files, memory and prompts until the freeze commit, and grep before each commit. The PI's slip required amending a commit. (Seed handling, 2026-10-26)
- Batch blocked downloads (depmap.org, vendor pages, journal supplements block scripts) into one user request (notes/USER_DOWNLOADS.md) instead of retrying. (Data agent, 2026-10-26)
- Check column contents, not names: DepMap 26Q3 Model.csv OnboardedMedia/FormulationID are swapped. Record such quirks in DATA.md. (Data survey, 2026-10-26)
- Keep RESEARCH_TREE.md to status, decisions and pointers; put finding detail and review text in notes/ so the tree stays a lean resume point. (M1 process review, RESEARCH_TREE ~10 KB)
