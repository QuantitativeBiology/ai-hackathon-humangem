# Research tree

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-28, M4: first TG-passing lever)

M2 is complete. The generic models were scored on 686 training lines (12 lineages). H1 beats H2 on tie-aware AUPRC: tAP_line 0.300 vs 0.292. Their paper-protocol MCCs are similar: 0.329 vs 0.331 (P57 × L-ge06). An independent reviewer recomputed every headline number with its own code and matched it exactly. The single W1 validation look (1 of 40, scalars only, 198 lines in 4 lineages) set **G\* = Human1**: Δ̄_lin(H2−H1) = −0.0071 [−0.0079, −0.0063]. Label-only gene-prior baselines score far higher (B2 tAP 0.90). About 42% of positive gene × line pairs fall outside the H2 addressable set. These two facts set the realistic bar for context-specific GEMs.

Update (2026-10-28): Stage A context families (W4) and GPR curation (W5) fail the training gate. W3's cofactor-pool biomass leave-out is the first lever to pass it (+0.034 on H1, flags W3-D1/W3-MEV); it passed its validation look (+0.032 [0.031, 0.033], look 2/40). Next: M6 combination factorial, W6 flux information, W7 hybrids with gene priors.

## Decisions
- **D9 (2026-10-28, PI, after the M3–M4 process review)**: scope cuts recorded as A12: W7 is time-boxed at 4 h and limited to hybrids; v1.3.0 is dropped; Stage B covers only F2 R1/F4 on training; the factorial has ≤ 8 cells scored in one validation look. Projected finish ≈ 2.2–2.4 days elapsed, above the 1–1.5 day target; the overrun comes from review and fix rounds, not compute (≈ 160 of ≈ 1,020 CPU-h used). Review depth is now fixed before each run (LESSONS).
- **D8 (2026-10-27, user)**: CORDA is not dropped. It runs as a 10-line comparison (A7: one representative panel line per lineage, R1 and R3, real RNA, ≈ 33 CPU-h). It is descriptive only, not eligible for P1.
- **D7 (2026-10-27, user)**: two process simplifications from the M2 process review, recorded as amendment A6. (1) Guard-internal changes go to notes/guard_changelog.md while frozen label, metric, split and list SHAs are unchanged. (2) New validation looks reuse the reviewed harness: a non-author config review, a dry run and the pre-flight, with no full code review unless harness code changes.
- **D6 (2026-10-27, M2)**: G\* = Human1, fixed by the single pre-registered W1 validation look (09b69a2). The PI ratified the Stage-A specifications as amendment A4 (fixed label-free A5 grid, one solver per family, KO isozyme convention, A6 CORDA dropped by its frozen cost rule, Stage A on M0). Amendment A5 covers test renames and the staged check of git-ignored frozen inputs.
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
- [x] 3. Baseline: generic Human1 & Human2 (pipeline agent, Opus; predictions only, no DepMap scoring until DESIGN.md fixes labels/splits)
  - [x] 3.1 Human1 = v1.19.0 (commit 022ed5c) fetched; Entrez mapping 2,874/2,887 (H1), 2,857/2,870 (H2) one-to-one; 2 unmapped (GK3P, SLC37A4 alt ENSG); join on Entrez (symbols renamed for 7–10 genes)
  - [x] 3.2 Single-gene KO, Ham's medium, Gurobi (16 WLS sessions OK): H1 149 / H2 154 lethal (ratio_r<0.01; none in 0.01–0.5); ~50 s per model on 15 workers. Verification 45/46 (only known H1 UMPS). Independent review (notes/review_baseline_ko.md): all genes recomputed with own code, max diff 2e-12 → accept with fixes; fixes applied (ratio_r rounding → bitwise-identical reruns; loud failure on non-optimal; per-release task files; SHA-pinned Human1 download).
    - [!] Human1 UMPS/DHODH non-lethal on Ham's: thymidylate synthase MAR04644 reversible in H1 (irreversible in H2) → UMP via thymidine bypass. Kept as released; a model difference to report.
    - Predicted ratios are near-binary → essentiality scores for AUPRC need a continuous score (e.g. ratio + tie-breaking / flux-based scores), to be fixed in DESIGN.md.
  - [x] 3.3 **M2: W1 benchmark** (W1 agent, Opus, ≈41 min wall / 3.5 CPU-h; review notes/reviews/M2_W1_review.md, GO after re-check). Training: tAP_line H1 0.300 / H2 0.292 (lift 2.42 / 2.36), without CE_fit 0.145 / 0.142; MCC P57×L-ge06 0.329 / 0.331, PGR 0.334 / 0.331. Baselines: B2 0.900, B6 0.888, B8+ 0.42, B3 0.33, B0/B1 = null 0.124. SESOI_C1b = 0.00572 (frozen). Validation look (1/40, 191 s): H1 0.302, H2 0.295 → G\* = Human1.
    - [!] H1 > H2 in all 12 training lineages, and also on validation, so H1.2 (H2 improves on H1) is not supported. MCC is essentially tied.
    - [!] The gene-prior baseline (B2, 0.90) is about 3× the generic GEM. Any context-specific gain must be read against B2/B8, not only G\*.
    - [!] About 42% of positive pairs are outside the addressable set (GPR-mapped genes), which caps achievable recall for any GEM.
    - [-] The positive control B6 ≥ B2 fails with untuned defaults (0.888 vs 0.900). Tuning belongs to W7.
    - Task-augmented scores (P57) lower tAP on all genes but raise it without common essentials.
- [~] 4. Medium (exchange bounds) — per-line culture media; error attributable to medium
  - [x] 4.1 W2 pilots and media library (Opus; ≈ 17 CPU-h of sweeps). Blind two-agent extraction: 312/313 components agree and every difference is adjudicated. The library is frozen (63a51df; A9). 472 KO sweeps; M0 = baseline bitwise. Review notes/reviews/M3_W2_review.md (GO after re-check).
    - [!] Medium moves few genes (mean |Δratio_r| ≈ 1e-3 per line), so the shuffled-medium null will be tight.
    - [-] M3 calibration is infeasible as specified (growth ≈ 12.3 × the O2 bound; other uptakes stay at 1000), so M3 is reported as uncalibrated and excluded.
    - [!] Sterol branch: serum cholesterol rescues only DHCR7. Biomass needs vitamin D3 products reachable only via lathosterol → 7-dehydrocholesterol, a model-structure issue for W3/W5. CAD is rescued by orotate but not uridine (W5).
    - Bug found in the A9 fix round: 2-ME did not open cysteine uptake. It is fixed and changes 0 gene ratios.
- [~] 5. Objectives / metabolic tasks / biomass (W3)
  - [x] 5.1 W3 M3–M4 (efc4138, clean-tree rerun 603cb9a; review notes/reviews/M4_W3_review.md, numbers reproduced bitwise; ≈ 24.5 CPU-h + 2.1 gene-pFBA). P57 equals Human-GEM taskEssentialGenes (8/8).
    - [!] **B-loo:cofactor_pool_biomass passes TG**: H1 +0.0343 [0.0336, 0.0349], 12/12 folds (H2 +0.0369). It removes ~20 biomass-artefact calls (CYP26/24, carnitine, BH4, COQ) and ~15 real requirements (heme, PDXK, TPK1, RFK). Sole W3 candidate for the validation look W3_M4_shortlist, with flags W3-D1 (control reading fixed after seeing results; beats a post-hoc random-removal control by +0.062, loses to the B2-informed comparator) and W3-MEV (under serum it removes mevalonate-core essentiality) (A11).
    - [-] Metabolite-pool B-loo (+0.0056, CE-only) withdrawn (A11c). Task modes, O-MOMA, O-NGAM, O-cap, D-modes do not qualify; task modes help only without common essentials (as W1).
    - [!] With vitamin-D B-loo, serum cholesterol rescues post-squalene sterol genes dose-dependently; dihydroorotate B-loo lets uridine rescue CAD. ETF positive control fails in H1.
    - [-] A8 gene-pFBA tie-break: +0.0003 vs ratio_r (H1), breaks only 13–14 ties; far below B8+.
    - [x] **Validation look W3_M4_shortlist (look 2/40, 5537c83; harness 5e44fb6 reviewed + re-checked): gate PASSED.** Δ̄_lin vs G* = **+0.0322 [0.0308, 0.0334]** (198 lines, 4 lineages); tAP_lin 0.302 → 0.334; without CE_fit 0.149 → 0.168; expressed pairs 0.352 → 0.380. Training estimate +0.0343 replicated. H2 variant 0.329 (descriptive). Flags W3-D1, W3-MEV carried. This is the first lever to beat G* on held-out lineages; candidate component for P1 (M6).
    - [ ] W7 exploratory: 'artefact groups only' B-loo variant (reviewer suggestion, post hoc).
- [~] 6. Context-specific models (omics integration; method-vs-tissue clustering test)
  - [x] 6.1 W4 pilots (Opus): F0, F2 GIMME, F3 FASTCORE, F4 E-Flux, F5 LCP builder; repair rule restores 12/12 lines (0.62–0.88 of generic WT); loss rule order-independent; 1.7–3.3 CPU-min per build. CORDA 1.5–1.7 CPU-h per model → full A6 dropped (A4); replaced by a 10-line comparison row A6′ (A7, user D8).
  - [x] 6.1b GPRTransform (Machado 2016; user suggestion) assessed label-free (notes/methods/GPRTransform_assessment.md, ≈1.7 CPU-h). The gene-level KO is exactly equal to the Boolean KO on both parents. The paper's isozyme encoding inflates capacity (Human2 WT +25%). It has no regulatory network: "regulation" means only expression bounds on gene usage. The repo has no license, so the method was reimplemented.
    - [-] Gene-level E-Flux as family F6: no extra distinct KO values on the pilot → not adopted.
    - [x] Gene-pFBA KO-cost tie-break (A8): see 5.1; negligible.
  - [-] 6.4 **M5 panel pre-check (A12(b))**: context families on H1 + W3 B-loo, 30 panel lines, pinned 4729547, ≈ 3.2 CPU-h; recomputed independently (notes/reviews/M5_precheck_review.md). Δ̄ vs H1 + B-loo: F2 R1 +0.0031 [0.0017, 0.0045] 8/10; F4 +0.0044 [0.0035, 0.0052] 10/10 (fails viability gate). Both < +0.005, so the context lever is not selected: Stage B/C skipped; M6 factorial has no context axis. Note: on top of B-loo both families turn slightly positive (they were negative vs generic H1 in Stage A).
  - [x] 6.2 W4 code review (GO after B1–B3 fixes). Stage A run from a pinned worktree at c280413: 1,380 builds plus 600 re-sweeps, 0 failures, 108.7 CPU-h, 6.2 h wall.
  - [-] 6.3 **Stage A training scoring** (75c0805; review notes/reviews/M3_W4_stageA_review.md, independently recomputed). **No family passes the training gate vs G\* = H1.** Nested LOLO Δ̄: F0 −0.0080, F2 −0.0067, F3 −0.0102, F4 −0.0029 (gate-failed; GR scale artefact), F5 −0.0114, CORDA −0.0127 (10 lines, descriptive). Vs H2 nothing reaches +0.005 (best F4 +0.0038). F2 R1 is carried to M6 flagged "fails TG".
    - [!] Less pruning is always better (R1 > R2 > R3); the genomic constraint changes ≤ 0.0002.
    - [!] H4.1: method explains more variance than tissue (KO R² 0.199 vs 0.085). F0/F2/F3 carry tissue signal beyond the permuted control; F4 is method-dominated by construction.
    - Interpretation: expression-driven pruning removes the generic model's correct calls faster than it adds line-specific ones. The remaining levers are medium/tasks/GPR (W2/W3/W5) and W7 hybrids with gene priors.
- [~] 7. GPR curation (proteomics, paralogs); curation vs constraint contribution
  - [-] 7.1 **W5a/W5b as specified fail TG on both parents** (b280526; review notes/reviews/M3_W5_review.md, independently recomputed). Training LOLO Δ̄ vs parent: H1 −0.0286, H2 −0.0293 (0/12 folds). W5b = W5a (no co-essential pair reaches r ≥ 0.4). The ≥ 30% identity sensitivity gives −0.0021 / −0.0031. ≈ 0.7 CPU-h.
    - [!] Literal "Ensembl paralog" AND→OR merges Pol I/II/III subunits on the lumped MAR07161 and costs 0.026. All 12 POLR merges were rejected on biology in the edit review.
    - [!] The H1 TYMS irreversibility fix helps consistently (+0.0007, 9/12; +0.0021 without CE_fit) but stays below SESOI.
    - [!] CAD lethality under uridine is caused by (S)-dihydroorotate in the biomass pool MAR10064; removing it makes the CAD KO ratio 1. This goes to W3.
    - Open: the UniProt removal path (file requested from the user, USER_DOWNLOADS item 12).
- [ ] 8. Flux-state information content (subtypes, drug screens, metabolomics)
- [~] 9. Resource: SBML models per DepMap 26Q3 line
  - [x] 9.1 W8 stub (Sonnet, ≈ 17 min): src/nextgem/resource.py export and gates. The demo exports (2 generic parents and 2 F0 panel models) pass all 5 gates. REPORT.md skeleton.
  - PI decisions:
    - a KO gate floor of 1e-9 absolute next to 1e-6 relative (LP noise);
    - closed reactions kept with bounds (0,0);
    - `version` = the Human-GEM release.
  - Open: builders must emit the re-opened (repaired) reaction set; F4 needs its scaled bounds stored to be exportable; MEMOTE/FROG need a separate venv.
- [ ] 10. REPORT.md
