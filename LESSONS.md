# Lessons (process rules)

One line per rule, with evidence. Pruned at each milestone.

- Set a growth medium before any FBA on Human-GEM: the raw model has all 1,654 exchanges open, so it grows without glucose or O2 and knockouts look non-essential. (Setup FBA test, 2026-10-26)
- Load Human-GEM SBML once (~15 s) and use `with model:` / `model.copy()` for perturbations. (Setup timing, 2026-10-26)
- DepMap expression has several rows per model; filter `IsDefaultEntryForModel == "Yes"` (2,314 → 2,130 rows, one per model). (Setup check, 2026-10-26)
- Round/clip FBA knockout ratios (6 dp) before ranking metrics: ~2,650 genes at ratio≈1 spread over ~60 floats changed AP by up to 0.018 between identical re-runs. (Baseline review, 2026-10-26)
- Have a non-author agent recompute results with an independent code path: it found the float-tie AUPRC instability that the builder's 36 passing checks missed. (Baseline review, 2026-10-26)
- DepMap 26Q3 `Model.csv` medium columns are swapped vs their names: `OnboardedMedia` = formulation ID (MF-xxx), `FormulationID` = text ("RPMI + 10% FBS"). Check column contents, not names. (Data survey, PI-verified 2026-10-26)
- depmap.org, Thermo/Sigma formulation pages and journal supplements block scripted downloads; batch these into one user request (notes/USER_DOWNLOADS.md) instead of retrying. (Data agent, 2026-10-26)
- Run a specialised parallel review panel before building: round 1 found the design infeasible (non-growing media/models, ~12k CPU-h matrix) and two metric biases (ties, common essentials) before any compute was spent. (Design review r1, 2026-10-26)
- Generate every number quoted in design/report docs from a script; the v0 design had several hand-copied counts wrong (e.g. Nusinow∩CRISPR 375 vs 309). (Data review r1, 2026-10-26)
- Reviewers must never run the real split seed or write label summaries into committed files: a round-1 review published the seed-20261026 draw, forcing a re-seed. Give reviewers simulated seeds only. (Design review r2, 2026-10-26)
- Hooks/settings edited mid-session don't load until /hooks or restart; verify with a live probe before relying on a guard (two probes got through). (Label guard, 2026-10-26)
- Agents leave large intermediates (incl. label copies) in the shared scratchpad; sweep it at each milestone. (AI review r2, 2026-10-26)
