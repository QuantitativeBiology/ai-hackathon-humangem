# Lessons (process rules)

One line per rule, with evidence. Pruned at each milestone.

- Set a growth medium before any FBA on Human-GEM: the raw model has all 1,654 exchanges open, so it grows without glucose or O2 and knockouts look non-essential. (Setup FBA test, 2026-10-26)
- Load Human-GEM SBML once (~15 s) and use `with model:` / `model.copy()` for perturbations. (Setup timing, 2026-10-26)
- DepMap expression has several rows per model; filter `IsDefaultEntryForModel == "Yes"` (2,314 → 2,130 rows, one per model). (Setup check, 2026-10-26)
