# Research tree

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-26, M1: design frozen)

DESIGN.md v3 is frozen after three review rounds (7 reviewers in rounds 1–2, an unexposed verifier in round 3) and user decision D2. The primary endpoint is the lineage-equal-weighted per-line ΔAP (tie-aware expected AP) of the selected pipeline versus the better generic parent (Human1 v1.19.0 / Human2 v2.1.1) on 7 locked test lineages, with confirmatory sequence C1→C1b→C3→C4→C5. MCC/F1/precision/recall/accuracy under the Human1/Human2 papers' own protocols are required secondary endpoints for comparability. ftINIT is excluded (user's own analysis); context-specific models use ≥3 LP families. Label leakage controls: a PreToolUse hook (active, verified), a labels.py-only access rule, a label-exposure register, sealed test outputs and a user-supplied split seed revealed after this freeze. The generic-model knockout baseline exists and was independently reviewed; nothing has yet been scored against labels.

Next (M1 implementation → M2): build and independently review the M1 infrastructure (labels.py with fold scoping, tripwire, label-copy scanner, ledgers, tie-aware AP and binary metrics with toy tests, make_splits.py), draw the split once with the user's seed, then benchmark the generic Human1/Human2 models (W1). Planned compute ≈1,020 CPU-h, all LP, ≈36–42 h wall-clock from the start of M2 once the code is verified; implementation time comes on top. Pending from the user: four medium sheets (MEM/EMEM, F12, McCoy's 5A, L-15); their lines are flagged as fallback until then.

## Decisions
- **D0 (2026-10-26)**: Human2 := Human-GEM v2.1.1 (local); Human1 := latest Human-GEM v1.x release (to fetch). Baseline infrastructure started before DESIGN.md is final because it uses only generic models and no DepMap labels for fitting. Git: commit locally; push only when the user asks. Time budget: multi-day effort (user, 2026-10-26); report runtime with every result.
- **D3 (2026-10-26, M1)**: DESIGN.md v3 frozen. Round 3 found the v3 draft not freeze-ready; all its text fixes were applied by the PI (verification.md "v3 final check"), and the unexposed verifier ratified/amended the exposed author's protocol choices. design_numbers.py made deterministic (tie order). Split seed: supplied by the user, revealed only after this freeze commit.
- **D5 (2026-10-26, user)**: F-4 — the PI's post-hoc relaxation of the list-authorship rule is not simply ratified; a blind re-derivation by a fully unexposed agent is reconciled with the approved lists. Changes to leakage/exposure rules need user ratification from now on.
- **D4 (2026-10-26, post-freeze reveal)**: split seed = **38**, supplied by the user. Timeline: user stated 38 before the freeze; the PI kept it out of all files, memory and agent prompts until freeze commit 61b0beb (one slip in an unpushed commit was amended away before any agent read it; LESSONS.md). The split is drawn once by make_splits.py after it is implemented and reviewed; the ledger records this entry.
- **D2 (2026-10-26, user)**: ftINIT dropped entirely (user's own analysis: ftINIT underperforms the generic model). Sole comparator = generic Human1/Human2 FBA (+ label-free baselines). Add MCC, F1, precision, recall, accuracy under the Human1/Human2 papers' protocols as required secondary endpoints for comparability; AUPRC stays primary. Method-vs-tissue test kept with ≥3 non-ftINIT families. Target compute ≈1–1.5 days (user: 5–6 days too much).
- **D1 (2026-10-26)**: results/ is committed to git (PROMPT.md requirement); files >50 MB go under results/**/tmp/ or *.large.* (ignored) and must be regenerable by a script.

## Tree
- [~] 0. Setup & record keeping (PROMPT.md, RESEARCH_TREE.md, LESSONS.md, DATA.md)
- [~] 1. Literature & datasets survey → LITERATURE.md
  - [x] 1.1 Literature survey (Opus, 16 min, ~8.6k words, 88 refs) → LITERATURE.md. PI spot-checks passed (v1.19.0 = last v1.x tag; 419/2,798 screened Human2 genes are DepMap common essentials).
    - [!] No leakage-free baseline exists: published GEM-vs-CRISPR is MCC ~0.25–0.44, best-of-settings, no lineage hold-out; pooled metrics inflated by gene priors (ROC 0.94 pooled vs 0.61 per gene) → per-gene AUPRC + gene-prior baselines must be primary.
    - [ ] Method-vs-tissue clustering never tested on Human1/2 + ftINIT at DepMap scale → variance partitioning + permuted-expression control.
    - [ ] Per-screen medium available locally (ScreenSequenceMap → ModelCondition.GrowthMedia); fetch CRISPRGeneDependency.csv.
    - [ ] Paywalled papers to request from user (list in LITERATURE.md, batched with data requests).
  - [x] 1.2 Dataset survey/download (Sonnet, 49 min): fetched Nusinow proteomics (375 lines), CCLE metabolomics (928 lines), PRISM 24Q2 + 20Q2 secondary, GDSC1/2 r8.5, CORUM 5.1 + Complex Portal, Ensembl 116 paralogs, Cell Model Passports list (~0.9 GB). Blocked items → notes/USER_DOWNLOADS.md (asked user).
    - [!] Medium annotation: 2,097/2,554 models (1,282/1,361 CRISPR lines); base media RPMI 907, DMEM 291, MEM/EMEM 157, DMEM/F12 104, IMDM 72, F12 56, McCoy's 31, L-15 23; ~300 organoid/proprietary (unmappable); 239 serum-free. Model.csv OnboardedMedia/FormulationID columns are swapped (PI-verified).
- [x] 2. DESIGN.md — v0 committed (b105f88, ~9.1k words, Opus author, 12 min)
    - v0 primary endpoint: paired Δ mean per-line AUPRC (selected pipeline − best generic model) on locked test lineages, line bootstrap stratified by lineage, success = 95% CI lower bound > 0; per-gene AUPRC (≥5 pos & ≥5 neg lines) confirms context specificity. Splits: scripts/make_splits.py, seed 20261026, 5 test / 4 validation lineages, rest training with LOLO-CV. Seeds [11,23,37,41,53,67,79,83,97,101].
  - [x] 2.1 Review round 1: 7 parallel reviewers (external PI, statistics, AI failure modes, algorithms/compute [Opus]; data, prior work, cancer-metabolism domain [Sonnet]; 3–33 min each) → notes/design_review_r1/. All verdicts: revise. PI decisions: notes/design_review_r1/PI_decisions.md (34 items).
    - [!] Common essentials = 86.5% of positives (L-ge05) → confirmatory C1b without them (list from train∪val only).
    - [!] sklearn step-wise AP rewards tie-breaking (+0.03–0.05/line under jitter) → tie-aware expected AP; B8 generic+expression tie-break baseline.
    - [!] Generic Ham's medium contains hypoxanthine + thymidine → all de novo purine genes, TYMS, DHFR predicted non-essential in H1/H2.
    - [!] Biomass needs micronutrients absent from RPMI/DMEM; naive expression pruning gave 12/12 non-growing models → growth checks + fixed repair rule before any sweep.
    - Compute after fixes ≈2,650 CPU-h, ≈90 h wall-clock (ftINIT-bound).
  - [x] 2.2 DESIGN v1 (design author, 19 min): all 34 PI decisions applied; ~13k words incl. appendices; numbers from scripts/design_numbers.py (label-free: IDs/headers only, PI-checked). Tie-aware expected AP; lineage-equal-weighted per-line ΔAP; fixed sequence C1→C1b→C3→C4→C5; SESOI 0.01; 7 test / 4 val / 12 train lineages (85% of simulated draws pass label-blind checks). Compute ≈2,940 CPU-h. 11 open questions for round 2.
    - Label exposure logged (DESIGN §11.5): author's `head` printed ~10 dependency values for one line (ACH-002730) before the label guard hook existed.
  - [x] 2.2b Label guard: PreToolUse hook (.claude/hooks/guard_labels.py) blocks agent reads of raw label files; user-approved 2026-10-26; active after /hooks reload.
  - [x] 2.3 Review round 2 (same 7 reviewers, 1–7 min each): all "revise minimally, then freeze"; 5/7 need no round 3. PI decisions: notes/design_review_r2/PI_decisions.md (17 items).
    - [!] Split seed 20261026 compromised: round-1 stats review published its draw → test tier known. Retired; new seed supplied by the user after the freeze commit (commit-reveal).
    - [!] Label hygiene: all-line label copy (ge_U.parquet + 11 derived files) left in scratch by a round-1 reviewer → deleted by PI; content-based scanner required. Guard hook hardened (fail-closed, labels_val, data/locked, recursive data searches) but inactive until user reloads /hooks.
    - Medium requirement redefined (non-limiting micronutrient floor; pFBA-based req for synthesisable serum lipids); transporter gate = context; time-limited MILP incumbents accepted in all arms (fair to B7); C6 last in gatekeeping.
    - Compute ≈2,950 CPU-h (≈1,550 MILP), wall-clock ≈125–150 h; 12 MILP + 16 LP workers (RAM-bound).
  - [x] 2.4 DESIGN v2 "candidate for freeze" (design author, 19 min): 17 decisions applied; medium-attribution set 835 lines (4 sheets), 9 two-media lineages; FH/SDH background dropped (n=2); Stage-A table ≈435 CPU-h; total ≈2,950–3,005 CPU-h. LOH via low-CN proxy (no LOH file).
  - [x] 2.5 Round 3 (unexposed Opus verifier, 12 + 7 min): v2 → 17 fixes; a priori lists re-approved from biology only (notes/apriori_lists_approved.md); v3 final check → text fixes (panel fold rule, workers 24+≤4 vs licence, Human1/Human2 protocol alignment: 21-step grid, NaN MCC, parent-model gene sets) applied by PI.
  - [x] 2.6 **M1: DESIGN.md v3 frozen** (61b0beb).
  - [x] 2.8 F-4 re-check (user decision D5): blind lists by a fully unexposed agent (57670ff, committed before seeing the approved lists), reconciled item by item → notes/apriori_lists_final.md (amendment A2). Corrections both ways (e.g. ATP5MC1/2/3 is an AND complex, not isozymes; ten wrong PMIDs fixed; FH/SDH dropped, n=3).
  - [x] 2.7 M1 implementation, independently reviewed (notes/reviews/M1_*_review.md): metrics (0c06695; 691 tAP + 17k binary cases, 0 mismatches), labels/splits/guards (code commit C = 639f3a9; reviewer GO after re-trying exploits). 156 tests pass. Amendments A1 (split procedure), A2 (final lists), A3 (scanner exempts pinned split file).
  - [x] 2.9 **Split drawn once** (S=38, 34575a9): test 7 lineages / 394 lines; validation 4 / 200; training 767. Labels built (c54e867; validation/test sealed). Residual guard bypasses documented as accepted risk (results/m1_infra/README.md).
- [~] 3. Baseline: generic Human1 & Human2 (pipeline agent, Opus; predictions only, no DepMap scoring until DESIGN.md fixes labels/splits)
  - [x] 3.1 Human1 = v1.19.0 (commit 022ed5c) fetched; Entrez mapping 2,874/2,887 (H1), 2,857/2,870 (H2) one-to-one; 2 unmapped (GK3P, SLC37A4 alt ENSG); join on Entrez (symbols renamed for 7–10 genes)
  - [x] 3.2 Single-gene KO, Ham's medium, Gurobi (16 WLS sessions OK): H1 149 / H2 154 lethal (ratio_r<0.01; none in 0.01–0.5); ~50 s per model on 15 workers. Verification 45/46 (only known H1 UMPS). Independent review (notes/review_baseline_ko.md): all genes recomputed with own code, max diff 2e-12 → accept with fixes; fixes applied (ratio_r rounding → bitwise-identical reruns; loud failure on non-optimal; per-release task files; SHA-pinned Human1 download).
    - [!] Human1 UMPS/DHODH non-lethal on Ham's: thymidylate synthase MAR04644 reversible in H1 (irreversible in H2) → UMP via thymidine bypass. Kept as released; a model difference to report.
    - Predicted ratios are near-binary → essentiality scores for AUPRC need a continuous score (e.g. ratio + tie-breaking / flux-based scores), to be fixed in DESIGN.md.
  - [ ] 3.3 Benchmark vs DepMap CRISPR (AUPRC overall / per gene / per lineage)
- [ ] 4. Medium (exchange bounds) — per-line culture media; error attributable to medium
- [ ] 5. Objectives / metabolic tasks (a-priori task set)
- [ ] 6. Context-specific models (omics integration; method-vs-tissue clustering test)
- [ ] 7. GPR curation (proteomics, paralogs); curation vs constraint contribution
- [ ] 8. Flux-state information content (subtypes, drug screens, metabolomics)
- [ ] 9. Resource: SBML models per DepMap 26Q3 line
- [ ] 10. REPORT.md
