# Review round 1 — PI decisions (2026-10-26)

Reviews: `external_pi.md`, `statistics.md`, `ai_failure_modes.md`, `algorithms_compute.md`, `data.md`, `prior_work.md`, `domain_metabolism.md` (this folder). The PI read each reviewer's must-fix summary; the design author must read the full reviews. Default rule: every **must-fix** below is accepted unless listed under "Rejected / modified". **Should-fix and optional** items are at the author's discretion. Each one must be marked accepted / rejected (with a one-line reason) in the change log.

## Endpoint and statistics
1. **AP computation.** Use tie-aware exact expected AP (McSherry & Najork 2008). The null is E[AP_random] (= the tie-aware AP of a constant scorer), not π. sklearn `average_precision_score` is kept as a sensitivity row only. Add a toy verification with known answers. (stats MF1)
2. **Primary endpoint.** Keep per-line ΔAP (selected pipeline − best generic parent) on locked test lines. Summarise as the **lineage-equal-weighted** mean: the mean of each test lineage's mean line ΔAP.
   - Report the line-weighted mean as secondary.
   - Inference: line bootstrap stratified by lineage, giving a CI conditional on the test lineages.
   - Transfer claim, tiered as the external PI proposed: a lineage-level t-interval, plus a descriptive criterion that the gain holds in ≥ ⌈0.7·k⌉ of k test lineages. (stats MF5, ext-PI MF3)
3. **Test lineages: increase to 7** (validation 4, rest training). Keep the size-stratified, label-blind draw.
   - The "redraw until checks pass" rule may use only label-blind checks.
   - Record the number of redraws and every rejected draw.
4. **Common essentials.**
   - Add confirmatory **C1b**: the primary endpoint excluding common essentials, using a list derived from training ∪ validation lines only (never DepMap's all-line list, for evaluation).
   - Per-line AP requires ≥5 positives.
   - Cells that are not estimable are NA, not dropped silently. (stats MF3)
5. **Multiplicity.** Use a fixed-sequence gatekeeping order **C1 → C1b → C3 → C4 → C5**. C2 is an estimate only. Drop Holm. (stats MF4, ext-PI MF2)
6. **Minimum meaningful gain.** Fix it a priori, and report the fraction of the generic → gene-prior gap closed. The author proposes the value with a rationale that uses no labels.
7. **C3 (context specificity).** Compare against E[AP_random], add a within-lineage label-permutation control, and pair against B3 (expression-only). (stats MF2)
8. **Expressed-gene subset.** Context-specific gains must also be reported on expressed gene × line pairs; this is a required secondary analysis. (ext-PI MF4)

## Baselines
9. Add **B7: published Human-GEM protocol**.
   - Per-line ftINIT on RNA, at the exact configuration the repo uses (1+1 with its tie-breaking; pinned commit).
   - Scored by the growth ratio and by task essentiality.
   - Contrast C6 = P1 vs B7. (prior-work MF1, compute MF3)
10. Add **B8: generic + expression tie-break** (generic FBA score; ties broken by low expression). This guards against "gains" that are only tie-breaking. (ext-PI MF1)
11. **Gene prior.**
    - Add a **medium-aware gene prior** (prior within the same base medium, training lines only).
    - Give medium as a feature to the B6 ML baseline. (data MF6, prior-work MF3)
12. **Human1 is a full arm at T1.** P1 is reported on each parent against that parent's generic model. (ext-PI MF6)
13. Add the **MCC of binary calls at a pre-fixed cutoff** as a secondary metric (comparability with the literature).

## Leakage and agent guards
14. **Single label access point** `src/nextgem/labels.py`.
    - It writes `labels_{train,val,test}.parquet`; the test file is hash-recorded.
    - The API is fold-scoped: it returns only training-fold labels to fitting code, with nested CV for every LOLO fit and selection.
    - A Python audit-hook tripwire in all project scripts raises if any module other than `labels.py` opens the raw CRISPR files (`CRISPRGeneEffect`, `CRISPRGeneDependency`, `ScreenGeneEffect`, common-essential lists).
    - A leakage canary is required. (AI MF1, MF3)
15. **Rejected / modified (AI MF1):** encrypting raw label files with a user-held key is rejected, as too heavy for a multi-day agent run. A Claude Code PreToolUse hook that blocks direct reads of raw label paths is **proposed to the user** (it needs a settings change). The design must not depend on it.
16. **Denylist of label-derived sources** for fitting and evidence: De Kegel and Ito paralog SL tables (allowed only as post-hoc annotation), Project Score, any DepMap common-essential list, and DepMap co-essentiality outside training folds.
    - Drug screens are not features for the essentiality predictor.
    - State whether Human2 GPR edits used DepMap. If unknown, say so and treat it as a caveat. (AI MF2)
17. **Test ledger.**
    - Hash-chained, append-only, committed on every append.
    - The look cap is enforced in code; looks are idempotent.
    - A rehearsal look on validation precedes T1.
    - The split script refuses to overwrite an existing split. (AI MF4)
18. **Frozen manifest** of metric code, label files, split file and verification tests (hashes). CI-style check script. (AI MF5)
19. **Provenance.** Long jobs record the git commit and must run from a clean tree or a pinned worktree. Writes are atomic. (AI MF6)
20. **Solver and licence failures abort.** Never e=0, never a generic fallback, never a silent switch to GLPK; GLPK only where explicitly configured. A shared lock pool covers ≤16 Gurobi sessions. (AI MF7)
21. **T1 detail sealed.** At T1, only the pre-declared scalar endpoints are revealed. Per-lineage breakdowns and the failure account are computed and sealed (hash in ledger) until P2 is frozen. Validation looks during W2–W7 return scalars only. (stats MF6)
22. **Attribution.** The full factorial of P1's levers runs on validation (for decisions) and at T1 on test as sealed, pre-declared secondary estimates (Shapley shares; "medium alone" on the generic model). It is revealed with the rest at the end. (ext-PI MF5, reconciled with stats MF6)

## Data
23. **Medium tiers.** Screen medium (ScreenSequenceMap → ModelCondition) > onboarded medium (`Model.csv` text column; the columns are swapped) > imputed RPMI + 10% FBS (flagged). Report results by tier. Organoid and proprietary media form their own stratum. (data MF1)
24. **Screen library** (Avana / KY / Humagne-CD / Brunello / TKOv3) and organoid status go into the split file as covariates and reporting strata. (data MF2)
25. **Fallbacks.** Specify the evaluation set and fallbacks for lines lacking expression (exclude from context-specific arms; report n) or WGS CN (constraint omitted, flagged). (data MF3)
26. **Correct all numbers to the data reviewer's verified values.**
    - |U| = 2,866
    - H1∩H2 = 2,766
    - Nusinow ∩ CRISPR = 309
    - Metabolomics: 927 lines, 652 screened
    - The author must regenerate every count in DESIGN.md from a script (`scripts/design_numbers.py`), never by hand. (data MF5)

## Biology and algorithms
27. **Medium modelling.**
    - Add a fixed micronutrient/vitamin set required by biomass (linoleate, linolenate, tocopherols, retinoate, lipoate, B12, …), with a growth check per medium before any sweep. (compute MF1)
    - Bound uptakes relative to biomass requirements (log-fraction scheme), not 1000 or 10.
    - Remove nucleosides and hypoxanthine from the serum model.
    - Use cystine, not cysteine.
    - Gate Asp/Glu uptake on transporter expression.
    - Fix the sterol control.
    - Add a growth-rate calibration check. (domain MF1–3)
28. **Pruning repair rule.**
    - Add a fixed repair step so extracted models can grow (e.g. minimal re-addition to restore biomass and essential tasks, as in ftINIT), defined before any labels. (compute MF2)
    - Gene-loss constraints that kill a viable line's model: drop only that loss and log it. (domain MF4)
29. **Run matrix.**
    - Screen extraction configurations on a 100-line training panel first (label-blind gates plus training-fold labels).
    - Then run ≤4 configurations per line, with ≤1 ftINIT configuration.
    - Use measured runtimes (ftINIT 18–25 min/line). Budget ≈2,650 CPU-h, ≈90 h wall-clock. (compute MF3–4)
30. **MOMA:** build once with `add_moma` and use a unique-solution formulation, or drop it if too costly after a pilot. (compute MF5)
31. **W3 task modes.** Reformulate the ill-posed task-as-constraint and NGAM modes per the compute review. (compute MF6)
32. **Genomic backgrounds.** Use metabolism-linked groups with GEM mechanisms (e.g. ENO1/ENO2, ME2/ME3, MTAP, FH/SDH, IDH1/2), using ploidy-aware CN where available. Allow pooled within-lineage contrasts. Report n per group, with a minimum n defined a priori. (domain MF5, data MF4)

## Literature fixes (LITERATURE.md)
33. Correct Vieira 2022 (733 CCLE lines; parameters chosen on an MCF7 reference). Correct the claim about ftINIT at DepMap scale: Gustafsson 2023 ran 891 lines for essentiality, and the untested claim is specifically about method-vs-tissue clustering. Fix the reference numbering gap and the DeepMeta PMID status. Add the missed works the prior-work review lists.

## Process
34. **Prior art and reuse.** Add a section naming the Human-GEM `DepMapGeneEss`, `gradedEssentiality.py` and `raven_toolbox` (tasks, init) code, and stating what is new relative to them. (prior-work MF2)
