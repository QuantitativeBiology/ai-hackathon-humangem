# Design review round 3: final verification of DESIGN.md v2

Verifier: round-3 agent, not the design author, not in `notes/label_exposure.md`. Date: 2026-10-26.
Reviewed: DESIGN.md v2 (HEAD a8ba7b8, which is fbeb491 plus a RESEARCH_TREE line), PROMPT.md, CLAUDE.md, all seven `notes/design_review_r2/*.md`, `notes/design_review_r2/PI_decisions.md`, `notes/label_exposure.md`, `results/design/design_numbers.json`, `results/baseline/*_ko.parquet` (label-free KO ratios) and the Human-GEM 2.1.1 SBML (biomass structure).

**Label safety.**
- The hook probe `head -c 20 data/depmap/CRISPRGeneEffect.csv` was **blocked** ("Blocked by label guard"). A later `ls` of `data/locked` was also blocked.
- I did not open `notes/design_review_r1/statistics.md` or `notes/design_review_r1/data.md` (I read only the PI note in the first lines of the former, via `head -5 | grep`). I opened no label file and computed no split draw.
- **Disclosure (see V-17).** Two files I was told to read contain all-line *aggregate* label statistics, so I have now seen them:
  - `notes/label_exposure.md`, "Later entries": the common-essential share of positives and ASNS prevalence by medium;
  - `notes/design_review_r2/statistics.md` SF-b and R2-10: the null AP on all genes and without common essentials.

  This is the same aggregate exposure the PI declared for itself. I saw no gene × line values and no statistics for any line or lineage.

**Scope change applied (user decision, relayed by the coordinator; RESEARCH_TREE decision D2).**
- ftINIT is dropped. That removes B7, B7g, B7t, C6, per-line F1, the ftINIT rows of Stage A and Stage C, and the H1 ftINIT arm.
- Binary metrics (MCC, F1, precision, recall, accuracy) become required secondary endpoints, computed under the Human1 and Human2 papers' own protocols.
- Checks that concern only B7, C6 or ftINIT are marked **MOOT (v3)**. The sections that break when these are removed are listed in Part 4.

**Verdict: not freeze-ready as v2.**
- The v2 text has 10 failing items: V-1, V-2, V-4, V-5, V-6, V-7, V-8, V-9, V-13 and V-14.
- Four items are notes or housekeeping actions for the PI: V-10, V-15, V-16 and V-17.
- The v3 scope change also requires rewrites across 15 sections (Part 4).
- None of these is an architectural problem. Once the fixes and Part 4 are applied, a PI check against this file is enough; no fourth review round is needed.

---

## Part 1. The 17 PI round-2 decisions

| PD | Decision | Status | Where checked / notes |
|---|---|---|---|
| r2-1 | Retire seed; user-supplied seed after the freeze commit; `make_splits.py` refuses 20261026 and records the seed source, the freeze commit and every redraw; §4.1 corrected; register updated (4 entries) | **FAIL (V-1, V-2)**; the rest passes | §0, §4.1, §7, §11.5, §12, §15. All four register entries are present in §11.5 and in `notes/label_exposure.md`. The seed source goes to `SEED_SOURCE.md`, and the attempt log records every redraw. Two procedural defects remain (V-1, V-2). |
| r2-2 | C6 last: C1→C1b→C3→C4→C5→C6; B7g the only confirmatory comparator | PASS for v2; **MOOT (v3)** | §0, §6.3 (lines 26, 267, 269, 278), §16 R2-1 and Appendix D all agree. No stale v2 order remains (the v1 history row at line 1064 is historical). Under v3 the chain becomes C1→C1b→C3→C4→C5 (Part 4). |
| r2-3 | B7 renamed; validation + test; training lines deferred; gap stated | PASS (wording nit V-12); **MOOT (v3)** | §0, §2, §8, §15 |
| r2-4 | B8 both tie-break directions | PASS | §6.3, §8 (B8+/B8−/B8′; the claim requires all three), §10.2 |
| r2-5 | SESOI basis restated; C1b SESOI with its own basis | PASS, with the claim-language defect **V-4** | §6.4. The basis separates "convention" from "detectable", and SESOI_C1b is lift-scaled and frozen at M2 on training lines. But the claim tiers still hard-code 0.01 (V-4), and C5 has no SESOI statement (V-5). |
| r2-6 | Time-limited MILP incumbents accepted; failure = no incumbent / infeasible; symmetric no-information score; licence errors abort | PASS (note) | §11.4. DESIGN also counts "WT ≤ threshold after repair" and "persistent numerical LP status" as failures. These come from ext-PI N3 and v1, they apply symmetrically, and DESIGN reads PD6's self-contradictory "licence errors count as failed / still abort" correctly. Under v3 the MILP clause loses its object, but the failed-build rule stays for LP families. |
| r2-7 | Requirement: fixed non-limiting floor; serum req = pFBA ≥99% with c open; auxotrophs FVA-min; reference medium | PASS | §9 W2 "Biomass micronutrient floor", "Serum model", "Requirement req_c"; §12 W2 tests (non-binding assertion, req ≠ 0). |
| r2-8 | Gate a separate toggle; "medium alone" on gate-free variants; gate counts as context | **FAIL (V-6)**, otherwise PASS | Correct in the W2 ladder, W2 attribution (i) and §10.1 step 5. **§10.2 still says "medium alone (M0 → documented medium …)"** (V-6). |
| r2-9 | Attribution set = sheets on disk at M3; fallback flag; two-media lineages recomputed | PASS (note) | §9 W2, App. C: 835 / 557+180+60+38, 198 fallback, 9 lineages; all from dn keys. PD9 says "~870". That is the all-tier count (582+185+64+39); DESIGN restricts the set to the screen tier, as v1 already did (onboarded = sensitivity). This is consistent. Optionally add "(870 including the onboarded tier)". |
| r2-10 | F1 repair off-M0 | PASS; **MOOT (v3)** | W4 F1 row, §12 W4 |
| r2-11 | Loss rule; order-independent drop | **FAIL (V-7, minor)**, otherwise PASS | W4. "Low CN" is undefined and separate from the LOH proxy, and the step-4 order has no tie-break, so the §12 "shuffled input order gives the same result" test can fail on ties. |
| r2-12 | Frozen Stage-A table whose total matches §13 | PASS for v2 (V-9); **must be rebuilt (v3)** | Recomputed below (Part 3c). The arithmetic is exact, but the F3 row contradicts row A8 (V-9). |
| r2-13 | 12 MILP + 16 LP; final counts from ftINIT RSS at M3; ≥20 GB free; MILP lane from M2 | PASS; **MOOT/rewrite (v3)** | §11.4, §13, §14 |
| r2-14 | Content-based label-copy scanner, every milestone and pre-push | PASS | §11.1 (fingerprints and shape heuristic; scratch and `/tmp/claude-*`), §12 agent guards |
| r2-15 | `sitecustomize` tripwire (pandas, pyarrow, open); residual risk; hook fail-closed; liveness probe per session | PASS | §11.1, §12. The probe target differs between DESIGN (`data/locked/HOOK_PROBE`) and CLAUDE.md ("a label read is blocked"); both are acceptable. |
| r2-16 | Sealed T1 outputs in `data/locked/`; hashes only; pre-commit guard until T2 unseal | PASS | §10.2, §11.1 |
| r2-17 | openpyxl installed and pinned; guide maps low priority | PASS | `requirements.txt` has `openpyxl==3.1.5`; §3.2, §5.4 |

## Part 2. Round-2 reviewer must-fixes

| Reviewer item | Status | Note |
|---|---|---|
| ext-PI N1 / stats NM1 / AI N1: re-seed with a value no agent can know or grind | PASS in substance, but **FAIL V-1, V-2** | User integer after the freeze commit; a SHA is explicitly rejected. |
| ext-PI N2: "medium alone" without the gate | **FAIL V-6** (one sentence in §10.2) | |
| ext-PI N3 / compute NMF2: MILP rule vs B7 | PASS; **MOOT (v3)** | |
| AI N2: purge caches + content scanner + CLAUDE.md rule | PASS in DESIGN; **register inaccurate (V-10)** | `ge_U.parquet` and `splits.pkl` are gone. **`testtier.py` is still in the shared scratchpad**, together with `prev.py`, `stats1.py`, `apnull.py` and `splits_sim.py`. I did not open them; judging by their names, they are the round-1 label-reading scripts. The register and §11.5 say `testtier.py` was deleted. The CLAUDE.md rule is present. |
| AI N3: hook live and fail-closed; liveness probe | PASS | It was live for this session. |
| AI N4: sealed location + pre-commit window to T2 | PASS | |
| compute NMF1: floor not scaled by f; f only on synthesisable components | PASS | |
| compute NMF3: F1 repair | PASS; MOOT (v3) | |
| compute NMF4: Stage A enumerated | PASS (V-9 aside); rebuild for v3 | |
| compute NMF5: memory-based workers | PASS; rewrite for v3 (no MILP lane) | |
| data N1: media sheets actually on disk | PASS | |
| domain N1: req_c for substitutable components | PASS | |
| domain N2: heterozygous damaging alone ≠ loss; order-independent drop | PASS, with V-7 | |
| prior-work (no must-fix; B7 rename, C6 end) | PASS; MOOT (v3) | |

## Part 3. Internal consistency

### 3a. Gatekeeping order
- Every v2 statement reads C1 → C1b → C3 → C4 → C5 → C6: §0 line 26, §6.3 lines 267 and 269, §16, the §17 v2 table and Appendix D. **PASS.**
- §6.5 gives "B = 10,000 for C1–C6". This is consistent.
- **Inconsistency (moot under v3):** §6.3 calls C6 a "sealed scalar" at T1, but §10.2 lists C6 among the scalars *revealed* at T1. It disappears with C6.

### 3b. Seed procedure (§4.1, §7, §12, §14)
- §4.1 says the PI makes the freeze commit, the user then sends S, and `make_splits.py` runs once with refusals, the attempt log, the ledger and the reference-draw cross-check. §7 and §14 M1 agree.
- **V-1 (FAIL). The HEAD check cannot be satisfied.**
  - Step 1 records the freeze SHA in the test ledger, and `ledger.append()` commits on every append (§11.3). The PI may also commit `SEED_SOURCE.md`.
  - Either commit moves HEAD past the freeze commit, so step 3's refusal "HEAD is not the freeze commit" always triggers.
  - Fix in §4.1 step 3:
    - old: `     - HEAD is not the freeze commit, or `src/`, `scripts/` or `tests/` differ from it;`
    - new: `     - the freeze commit is not an ancestor of HEAD, or `git diff --quiet <freeze SHA> HEAD -- src scripts tests DESIGN.md FROZEN.tsv` fails, or those paths are dirty in the working tree (only ledger rows and `results/splits/SEED_SOURCE.md` may be committed between the freeze commit and the draw);`
  - Fix in the §12 Splits row:
    - old: `a non-freeze HEAD and dirty code`
    - new: `code, DESIGN.md or FROZEN.tsv differing from the freeze commit, and dirty code`
- **V-2 (FAIL). The retired range is checked on S, but the exposure is on S + attempt.**
  - A user value such as S = 20261000 reaches the retired draw `default_rng(20261026)` at attempt 26.
  - Fix in §4.1 step 3:
    - old: `     - S = 20261026, or S lies in 20261026…20261125 (the retired seed's attempt range);`
    - new: `     - any attempt seed S + attempt (attempt 0…99) falls in 20261026…20261125, i.e. S lies in 20260927…20261125; the check is repeated per attempt inside the draw loop;`
  - Fix in §12: change `refuses seed 20261026 (and its attempt range)` to `refuses any S whose attempt range overlaps the retired range 20261026…20261125`.
- **Note (pass).** Seeds 1–1000 were simulated, and seeds 1–20 are stored as reference draws. This does not expose the production draw, because S is unknown until after the freeze. Nothing needs to change.

### 3c. Compute: §13 vs the Stage-A table (v2)
- I recomputed every row from the stated unit costs (22 CPU-min per ftINIT, 1.5 CPU-min per LP):
  - A1 400 builds → 146.7 h; A2 60 → 22; A3 1,600 → 40; A4 400 → 10; A5 1,440 → 36;
  - A6: the real-expression builds of A1–A5 total 200 + 60 + 800 + 400 + 1,440 = 2,900 → 72.5 h;
  - A7 2,400 → 60; A8 ≤ 25; A9 60 → 22.
- Totals: 520 ftINIT, 8,840 LP, **435 CPU-h, 191 MILP**. These match the table and §13.
- §13 column sums are 1,985 CPU-h and 711 MILP. With +20% (400 / 140) they give **2,385 / 850**. The conditional block (+620) gives **3,005 / 1,470**. The MILP lane before M6 is 661 × 1.2 / 12 ≈ 66 h, and 3,005 − 2,650 = 355. All consistent. **PASS** for v2.
- **V-9 (FAIL).** The W4 F3 row says CORDA is included "only if the M3 pilot is ≤ 2 CPU-h per model", but row A8 and §13 cost it at ≤ 0.25 CPU-h per model. At 2 CPU-h × 100 lines, A8 would be 200 CPU-h, not 25.
  - old: `F3b = CORDA (`corda` 0.5.1), included only if the M3 pilot is ≤ 2 CPU-h per model`
  - new: `F3b = CORDA (`corda` 0.5.1), included only if the M3 pilot is ≤ 0.25 CPU-h per model (row A8)`
  - The v3 note on CORDA is in Part 4.
- **V-11 (minor).** The W4 "Stability" line adds an "expression-noise perturbation" that no Stage-A row costs. Under v3 it is replaced anyway (Part 4).

### 3d. `[dn:…]` keys
- All **81 references (77 distinct keys)** in DESIGN.md resolve in `design_numbers.json`, including the shorthand continuations (`.key`) and the wildcards (`backgrounds.*.…`, `baseline_ko.*.…`). Checked by a script on the JSON.
- I spot-checked every value quoted in the text: cohort, coverage, medium tiers, attribution set by base, fallback 198, unmappable 195 = 109 + 45 + 41, the 9/13 two-media lineages, backgrounds, feasibility shares, loss-rule counts and split simulation. All match.
- `_meta.label_values_read = false`. **PASS.**

## Part 4. Other failing items in v2 (with fixes)

- **V-4 (FAIL, PD5). Claim tiers hard-code 0.01, so C1b's SESOI is never used.**
  - §6.4 "Claim language", old:
    `  - "detectable improvement": C1 passes;`
    `  - "improves": C1 passes and the point estimate ≥ 0.01;`
    `  - "meaningfully improves": CI lower bound ≥ 0.01.`
  - New:
    `  - "detectable improvement": the contrast passes its test;`
    `  - "improves": it passes and its point estimate ≥ its SESOI (0.01; SESOI_C1b for C1b);`
    `  - "meaningfully improves": its CI lower bound ≥ that SESOI.`
- **V-5 (FAIL). C5 has no SESOI or wording rule.** It is in the confirmatory chain, but neither the SESOI list nor the claim-language list mentions it.
  - §6.4, old: `- **SESOI for C1, C4 and C6 = 0.01 absolute in Δ̄_lin (all genes).**`
  - New (v3, C6 gone): `- **SESOI for C1, C4 and C5 = 0.01 absolute in Δ̄_lin (all genes).**`
  - Also change `(C1, C1b, C4, C6, each with its own SESOI)` to `(C1, C1b, C4, C5, each with its own SESOI)`.
  - If the PI prefers, C5 can instead get "detectable" wording only, as C3 does. Either way, state it.
- **V-6 (FAIL, PD8).** §10.2 T1 sealed column:
  - old: `"medium alone" (M0 → documented medium on the generic parent)`
  - new: `"medium alone" (M0 → best gate-free documented variant on the generic parent; PD r2-8)`
- **V-7 (FAIL, minor, PD11).** In W4 "Loss rule":
  - old: `    - heterozygously damaged (code 1) **and** shows LOH, low CN or TPM < 1; or`
  - new: `    - heterozygously damaged (code 1) **and** shows LOH (proxied by relative CN < 0.5 until an LOH file exists) or TPM < 1; or`

  In "Order-independent drop" step 4:
  - old: `in descending order of the lost gene's expression, until growth is restored`
  - new: `in descending order of the lost gene's expression (ties broken by ascending Entrez ID), until growth is restored`
- **V-8 (FAIL, minor). Small-group flag misses NRF2_active (43 lines).** W4 notes:
  - old: `Groups below 50 lines (ME2, NAPRT, PHGDH) are flagged as small`
  - new: `Groups below 50 lines (NRF2, ME2, NAPRT, PHGDH) are flagged as small`
- **V-10 (PI action). The register does not match the disk.**
  - Delete `testtier.py`, `prev.py` and `stats1.py` (plus `apnull.py` and `splits_sim.py` if they read label files) from the shared scratchpad `/tmp/claude-10105/-data-tmpA-nextgem/c5c0f45b-d54c-4e57-bdef-3d0657510b12/scratchpad/`, and log the deletion.
  - Alternatively, change §11.5 and the register to "`ge_U.parquet` and `splits.pkl` deleted; label-reading scripts remain (code only)".
  - The remaining `*.parquet`/`*.pkl` files there are the baseline-review KO tables and a pFBA pickle, which are label-free by name and path. The scanner (PD14) should confirm this at M1.
- **V-12 (nit; moot in v3).** In the §8 B7 "Stated gap", "F1's training gate … evaluated on panel lines (and validation)" mixes up the training gate and validation.
- **V-13 (W3, model structure). The B-loo index is undefined.** "B-loo_j (j = 1…13)" does not match Human-GEM 2.1.1, whose MAR13082 has **9 direct reactants**: ATP, DNA, H2O, RNA, glycogen, and the cofactor, protein, lipid and metabolite pools.
  - old: `| B-loo_j (j = 1…13) | biomass | leave one biomass component out |`
  - new: `| B-loo_j | biomass | leave one component out; components = the direct reactants of MAR13082 except ATP/H2O, plus the cofactor-pool groups that are floor- or artefact-linked (retinoid and vitamin-D catabolites, tocopherol catabolites, carnitine, heme, ubiquinol, lipoate, cobamide); enumerated per parent by script at M3 and frozen |`
- **V-14 (W3, model structure). Part of the B-add list is already in the H2 biomass.** The H2 cofactor pool (MAR10065) already contains ubiquinol (CoQ), heme and lipoic acid, and the metabolite pool (MAR10064) contains CoA. Adding those demands changes only coefficients, which W3 itself rejects as uninformative.
  - Fix: "B-add components already present in the parent's biomass are skipped for that parent and logged". Revise the list as in `notes/apriori_lists_approved.md` §6.
- **V-15 (note, Appendix A vs biochemistry).** Several Appendix A expectations are wrong on biochemical grounds:
  - IMPDH2 and ADSL are listed as "relaxed by hypoxanthine";
  - TYMS and DHFR are listed as pan-medium essentials, but they are thymidine- or HT-rescued, and the Ham's baseline already shows ratio 1;
  - NAMPT is listed both as pan-medium and as medium-dependent.

  The corrected panels are in `notes/apriori_lists_approved.md`. The Appendix A addendum should cite that file.
- **V-16 (note).** §11.5 says the register is "to be created by the PI", and F-2 lists actions that are done (register, r1 note, CLAUDE.md rule). Update this to "created; see notes/label_exposure.md".
- **V-17 (PI action).** Two committed notes files contain label-derived aggregates:
  - `notes/label_exposure.md` "Later entries" contains label-derived numbers;
  - `notes/design_review_r2/statistics.md` (SF-b, R2-10) quotes all-line null-AP aggregates, but it has no warning header, unlike its round-1 counterpart.

  This makes every agent who reads the register exposed to the same degree as the PI, and that conflicts with CLAUDE.md ("never write label summaries into committed files"). Redact the values in the register (describe them, e.g. "the CE share of positives; ASNS prevalence by medium"), add the round-1-style PI note to `notes/design_review_r2/statistics.md`, and add a register line for this verifier (aggregate-only, after reading the assigned files).

## Part 5. Sections that break when ftINIT is removed (v3 checklist)

Each of these must be rewritten. I give the minimal change.

1. **§0 Summary.** Delete "and better than … B7"; the C6 bullet; and "C6" from the chain, which becomes **C1 → C1b → C3 → C4 → C5**. The abstract sentence gives **C1 and C1b**.
2. **§2 Prior art.** Keep the rows as "read, not run". Drop "B7 runs exactly this configuration" and the `raven_toolbox` `init.ftinit` use. Task checking via `raven_toolbox.tasks.check` and `taskEssentialGenes.py` can stay for W3. Remove "new" item 1's ftINIT clause.
3. **§3.1.** Remove "B7 (H2 only; §8)". In the Human1 "full arm" row, P1(H1) now uses LP families only.
4. **§6.3.**
   - Delete the C6 row.
   - Delete "P1 vs B7t" and the C6-labelling bullet.
   - The "Consequence of the simple chain" sentence becomes "a C3 failure removes confirmatory status from C4 and C5".
   - Add the binary endpoints (item 15).
5. **§6.4.** Remove C6 from the SESOI and claim lists (see V-5). "Abstract-level sentence must give C1, C1b and C6" becomes "C1 and C1b, plus C4/C5 once T2 has run".
6. **§7 Seeds.** Delete "ftINIT stability study with `resolve_ties` off: SEEDS[0:3]". Use SEEDS[0:3] for the replacement LP-family stability study (item 9).
7. **§8.** Delete the B7 row. The validation gate becomes "shown next to B2, B2m, B6 and B8±". The comparator set is now generic H1/H2 (G\*) + B0–B3, B6 and B8±/B8′, per the user.
8. **W1.** Delete H1.4. Add the binary-metric endpoints (item 15).
9. **W4.**
   - Delete the F1 family row, the "≤ 1 ftINIT per line" clause in Stage B and the F1 rows A1, A2 and A9.
   - Method-vs-tissue builds become **F0, F2, F3, F4 (+F3b)** × {R1, R3} × {real, permuted}. That is ≥ 3 established non-ftINIT families (GIMME, FASTCORE, E-Flux), plus F0 as the naive reference, which satisfies the user's "≥ 3" rule.
   - Recommendation: CORDA is now the only PROMPT-named established method left, so if its 0.25 CPU-h pilot fails, run it on a 30-line subset of the variance panel rather than dropping it. H4.1–H4.3 keep their wording.
   - The stability study becomes: GIMME (alternate optima) and FASTCORE with expression-noise perturbation, on 20 panel lines × SEEDS[0:3] (≈ 120 LP builds, ≈ 3 CPU-h).
   - Repair rule: delete the F1 sentence. GPR convention: delete "F1 uses RAVEN's own convention".
   - Positive controls stay (TYR/DCT, urea cycle).
10. **Stage-A table (frozen per PD12).** Rebuild it as LP-only:
    - A3 1,600 (40 h), A4 400 (10), A5 1,440 (36), A7 2,400 (60), A8 ≤ 100 (≤ 25) and stability 120 (3);
    - A6 re-sweeps over the real-expression LP builds: 800 + 400 + 1,440 = 2,640 → 66 h;
    - total ≈ **8,900 LP, ≈ 240 CPU-h, 0 MILP**.

    Recompute it with `design_numbers.py` or by hand and freeze.
11. **Stage C.** About 150 group models × ≤ 4 **LP** configurations: ≈ 15 CPU-h, 0 MILP.
12. **§10.1** step 6 (H1 pre-builds), **§10.2** (C6 at T1; B7 references), **§10.3** (no change), **§11.4**:
    - the "MILP time limits" paragraph is reduced to QP MOMA only;
    - the workers become LP workers plus ≤ 4 Gurobi sessions for QP MOMA, verification and interactive work;
    - the no-information rule stays for LP-family failures, and the "B7 and H1 arm" wording is removed.
13. **§12.** Delete the B7 Hart one-shot positive control. Replace it with a label-free W4 positive control (e.g. GIMME/FASTCORE toy networks with known answers) and the existing TYR/DCT check. `[dn:hart_reference]` then becomes unused; keep it only if the Hart 2015 lines are used as the Human1-paper binary-metric comparison (item 15).
14. **§13, §14, §15, §16, App. B.**
    - **§13:** rewrite the table. Remove the B7 block (200/200), the H1 arm (200/190, or ≈ 15 LP CPU-h for P1(H1) on LP families) and the conditional +620 and +545 blocks. Stage A drops to ≈ 240; MILP drops to ≈ 5–55 (QP MOMA plus the W7 cap). The new total is ≈ 1,900 CPU-h, well under budget. Recompute the wall clock: no MILP lane, so the LP lane becomes the critical path.
    - **§14:** delete "MILP lane starts" (M2), "B7 one-shot Hart test" and "ftINIT peak RSS → worker counts" (M3), and the H1 pre-builds (M5).
    - **§15:** delete the B7 limitation and add "ftINIT excluded by user decision D2 (it underperformed generic in the user's own analysis); the method-vs-tissue test covers LP families only".
    - **§16:** F-1 is resolved (no per-line F1).
    - **Appendix B:** delete it.
    - Record all of this in §17 as a pre-freeze v3 change with reason "user decision D2".
15. **New required secondary endpoints: binary metrics under each paper's protocol.** Specify these in §6.1/§6.3 before freeze. Open issues I flag:
    - **(a) Which call rule.** For each paper, state the KO call rule and the label it used, taken from the paper's methods by the prior-work agent and cited. Human1 is Robinson et al. 2020 (*Sci Signal*); Human2 is Luo et al. 2026. Candidate rules are a growth-ratio cutoff or task failure.
    - **(b) Protocol vs model.** Human1's essentiality protocol builds per-line tINIT/ftINIT models. With ftINIT dropped, "under the paper's protocol" can only mean its **call rule and metric definitions applied to our models**, not its extraction step. State this explicitly.
    - **(c) Labels.** If a paper's label definition differs from L-prob, it is an additional label. It is frozen now, reported only as a secondary analysis, and never used for selection.
    - **(d) Primary endpoint unchanged.** AUPRC (tAP) stays primary (PROMPT). The binary endpoints are estimates, outside the gatekeeping chain, reported for G\*, P1 and every baseline with bootstrap CIs (§6.5).
    - **(e) Accuracy.** It is dominated by negatives at this prevalence; report it, but MCC is the headline binary metric.
    - The existing fixed cutoff ratio_r < 0.01 stays as our own binary rule.

## Summary of fixes required before freeze
- V-1, V-2: seed procedure (§4.1 step 3, §12).
- V-4, V-5: claim language and C5 SESOI (§6.4).
- V-6: "medium alone" wording (§10.2).
- V-7: loss rule (W4).
- V-8: small-group flag (W4).
- V-9: CORDA threshold (W4 F3).
- V-13, V-14: B-loo/B-add definitions (W3).
- V-15: Appendix A expectations; use the approved lists.
- V-10, V-17: PI housekeeping (scratch scripts; redact the aggregates in the register and the r2 statistics review).
- Part 5: the v3 ftINIT removal and the binary endpoints.

---

# v3 final check (DESIGN.md v3, working tree, uncommitted; 2026-10-26)

**Inputs.** DESIGN.md v3 (1,439 lines), `notes/label_exposure.md` (updated), and `notes/design_review_r2/statistics.md` (PI note added). For §6.8: `papers/aaz1482_sm.pdf`, the Robinson main text (`papers/5191391c-….pdf`), `papers/pnas.2516511123.sapp.pdf`, the Luo main text, and the Human-GEM script `code/DepMapGeneEss/evaluateDepMapEssentialityPredictions.m`. Text was extracted with ghostscript; figure pages were rendered to PNG to read their axes.

**Label safety.** No label file was read; the listing of the locked directory was blocked.

**Verdict: not freeze-ready as written; freeze-ready once the edits below are applied.** All edits are text-level, and no further review round is needed. Three items block the freeze:
- **N-1:** the panel training gate cannot pass as written;
- **N-7:** §6.8 deviates from the documented paper and repository protocol (D1, D2, D6);
- **N-2:** worker counts versus the Gurobi slot pool (D5).

The rest are minor.

## A. V-items and Part-5 items in v3

| Item | Status | Note |
|---|---|---|
| V-1 | PASS | §4.1 step 3 and §12 |
| V-2 | PASS | §4.1 step 3 and §12 |
| V-4 | PASS | §6.4 |
| V-5 | PASS | §6.4 (C1, C4, C5; C5 is in the claim language) |
| V-6 | PASS | §10.1 (ii), §10.2 |
| V-7 | PASS | W4 loss rule and tie-break; §12 W4 |
| V-8 | PASS | W4 notes |
| V-9 | PASS | F3 row → A6 bands |
| V-10 | PASS | §11.5 and the register (deletion logged) |
| V-11 | PASS | Row A7 is costed; §7 seeds updated |
| V-12 | MOOT | |
| V-13 | PASS | W3 B-loo |
| V-14 | PASS | W3 B-add skip rule; the list is taken by reference |
| V-15 | PASS (one residual, N-5) | App. A points to the file; the v2 addendum is marked superseded |
| V-16 | PASS | §11.5 |
| V-17 | **Done in the tree but stale in DESIGN (N-8)** | The register is redacted and the PI note has been added to the r2 statistics review. DESIGN §16 F-2 and the §17 row still say "still open". |
| Part 5, item 1 (§0) | PASS | |
| Part 5, item 2 (§2) | PASS | |
| Part 5, item 3 (§3.1) | PASS | |
| Part 5, item 4 (§6.3 chain, C6 row, consequence) | PASS | |
| Part 5, item 5 (§6.4) | PASS | |
| Part 5, item 6 (§7) | PASS | |
| Part 5, item 7 (§8) | PASS | |
| Part 5, item 8 (W1) | PASS | Citation nit in N-6. |
| Part 5, item 9 (W4: F1 retired; F0, F2–F4 + CORDA; ≥ 3 non-ftINIT families; stability; repair; GPR convention) | PASS | |
| Part 5, item 10 (Stage A) | PASS (arithmetic) | Recomputed: A1 480/12, A2 120/3, A3 480/12, A4 840/21 (240 + 120 + 480 real builds), A5 720/18, A6 120/30 (band 2: 30 builds ≤ 30 h), A7 120/3. Total **2,880 builds / 99 ≈ 100 CPU-h**, matching §13. |
| Part 5, item 11 (Stage C) | PASS | |
| Part 5, item 12 (§10.1, §10.2, §11.4) | PASS, except N-2 (stale 16-slot pool vs workers) | |
| Part 5, item 13 (§12) | PASS | |
| Part 5, item 14 (§13–§16, App. B, §17 v3 row) | PASS, except N-3 (step-5 wall-clock) and N-8 | §13 column: 40+25+60+10+100+100+15+15+110+160+40+100+150 = 925, plus 95 = 1,020. Correct. |
| Part 5, item 15 (§6.8) | PASS in structure; **amend per N-7** | |
| `[dn:]` keys | PASS | 80/80 references resolve |
| Gatekeeping | PASS | C1→C1b→C3→C4→C5 everywhere in operative text (§0, §6.3, §6.5 "C1–C5", §6.4) |

## B. No ftINIT, B7 or C6 arm remains; lists are verbatim

- **B7, C6, ftINIT.** In operative text these appear only as exclusion statements: §0, §1, §2 (rows "read, not run" or "context only"), §3.2 (Hart inputs "not used"), W4 "(F1) retired", §11.2, §12 W4 ("replacing the removed ftINIT Hart test"), §15 and §16. Every other mention is in the v1/v2 revision history and the Appendix D dispositions, which are historical records. **PASS.**
- **Lists.** DESIGN no longer restates any list; it points to `notes/apriori_lists_approved.md` (§5.4, W2 controls, W3, W4 backgrounds, W5, W6, App. A). The file is unchanged since I wrote it. **PASS**, with two residual inline statements that contradict the approved file:
  - **N-4 (fail, minor).** W2 "Controls", mevalonate core:
    - old: `It is reported as a **diagnostic panel at every W2 step**, because serum cholesterol makes it non-essential in FBA until W3 B-add supplies the non-sterol demands.`
    - new: `It is reported as a **diagnostic panel at every W2 step**. Ubiquinol is already a biomass component in H2, so HMGCR, MVK, PMVK and MVD should stay essential when serum cholesterol is supplied; a loss of essentiality under serum flags a modelling error (HMGCS1, IDI1, FDPS and GGPS1 are non-essential already on M0, for isozyme/GPR reasons; approved lists §1.1, §5).`
    - §12 W3 positive control, old: `mevalonate core under B-add`; new: `no loss of mevalonate-core essentiality under serum + B-add`.
  - **N-5 (fail, minor).** §3.1 M0 bullet:
    - old: `As a result, de novo purine genes (PPAT, GART, PFAS, PAICS, ADSL, ATIC, IMPDH2) have ratio_r = 1 in both parents [dn:baseline_ko.human2.panel_ratio_r]. That is a systematic, medium-driven false-negative block, which W2 targets.`
    - new: `As a result, the hypoxanthine-rescued de novo purine genes (PPAT, GART, PFAS, PAICS, ATIC) have ratio_r = 1 in both parents [dn:baseline_ko.human2.panel_ratio_r], a systematic, medium-driven false-negative block that W2 targets. ADSL and IMPDH2 also have ratio_r = 1, but hypoxanthine cannot rescue them biochemically (IMPDH1 OR-logic; ADSL to be traced); these are model errors, not medium errors (approved lists §1.2).`

## C. Other defects found in v3

- **N-1 (FAIL; blocks the panel training gate).** §9 TG item 3 counts a fold with < 5 estimable lines as "not positive". The Stage-A panel has 3 lines per lineage fold, so every panel fold is "not positive" and no LP family or threshold can reach ≥ 8 of 10.
  - old: `For gates computed on the 10-lineage Stage-A panel, the bar is **≥ 8 of 10** folds (ext-PI r2 S3).`
  - new: `For gates computed on the 10-lineage Stage-A panel (3 lines per lineage), the bar is **≥ 8 of 10** folds, and a panel fold counts as "not positive" if fewer than 2 of its 3 lines are estimable; the < 5-line rule applies to the 12 training-LOLO folds only (ext-PI r2 S3; r3 N-1).`
- **N-2 (FAIL).** §11.4 keeps "Gurobi slots: 16 `flock` slots" but starts 28 LP workers plus ≤ 4 Gurobi sessions. Either the LP workers do not use Gurobi, which makes the unit costs (Gurobi-measured) invalid, or they exceed the slot pool. Fixed by the worker amendment in D5.
- **N-3 (FAIL, minor). §13 step 5 is mis-timed.** Step 5 is ≈ 290 CPU-h in 8 h, which needs ≈ 36 workers. At 24–28 workers it takes ≈ 10.5–12 h.
  - old: `| 5 | W8 ‖ W6 ‖ W7 (W7 time-boxed at **8 h wall**, ≤ 100 CPU-h); T2 | ≈ 290 | ≈ 8 h (H22 → H30) |`
  - new: `| 5 | W8 ‖ W6 ‖ W7 (W7 time-boxed at **8 h wall**, ≤ 100 CPU-h); T2 at the W7 stop, W8/W6 finish after | ≈ 290 | ≈ 11 h (H22 → H33) |`
  - Shift step 6 to H33 → H37 and the §14 M7–M9 hours by +3.
  - Expected-wall-clock bullet, old: `**Expected wall-clock: ≈ 30–36 h (≈ 1.3–1.5 days)**`; new: `**Expected wall-clock: ≈ 36–42 h (≈ 1.5–1.75 days) at 24–28 LP workers**; cut-list items 2 and 5 recover ≈ 4–6 h`.
  - At the upper end this exceeds the user's 1–1.5-day target. Flag it to the PI and user.
- **N-6 (citation nits).**
  - §6.8 "Model per line", Human1: change `main text p. 5` to `main text pp. 5–6`. The 621-line DepMap analysis is on p. 6.
  - W1 H1.2, old: `Luo 2026 claim this in the main text, Fig. 2D.`; new: `Luo 2026 claim this for DepMap 19Q3 in SI Fig. S2D (cited in the main text, p. 3); main-text Fig. 2D is the Hart-2015 version series.`
- **N-7 (FAIL).** §6.8 departs from the documented protocols in three places. See D1 (grid), D2 (undefined metrics) and D6 (gene set; Human2 threshold source).
- **N-8 (stale).** §16 F-2: delete the "V-17 … still open" bullet and add "V-17 done (register redacted; PI note on the r2 statistics review)". In the §17 V-17 row, change `Still a PI action` to `Done by the PI`. After this check, F-7 becomes "ratified/amended per verification.md v3 final check, D".
- **N-9 (PI action, register).** The opening paragraph of `notes/label_exposure.md` still says listed agents "must not author or approve a priori lists, candidate sets or decision thresholds". That contradicts the later PI decision (aggregate exposure does not disqualify). Add one line at the top: "Superseded in part by the PI decision of 2026-10-26 below (disqualifying exposure = gene × line values or post-freeze val/test summaries)".

## D. Ratification of the design author's v3 choices (F-4 / F-7)

I am unexposed by the PI's rule (aggregates only) and did not author DESIGN.md.

1. **Binary-metric robustness grid: AMEND.**
   - DESIGN chose −1.0 to 0.0 in steps of 0.2 as "our choice". The Human1 SM Fig. S10 (p. 28, axis read from the rendered page) shows **21 thresholds, −1.0 to 0.0 in steps of 0.05**. The Human-GEM evaluation script uses the same grid: `scan.threshold = linspace(-1,0,21)`, line 106.
   - §6.8, old: `Human1-paper robustness grid, **our choice of steps:** −1.0, −0.8, −0.6, −0.4, −0.2, 0.0.`
   - new: `Human1-paper robustness grid: 21 thresholds, −1.0 to 0.0 in steps of 0.05, as in SM Fig. S10 (p. 28) and `linspace(-1,0,21)` in `evaluateDepMapEssentialityPredictions.m` (line 106).`
2. **Undefined MCC = 0: AMEND.**
   - The repository script computes MCC = (TP·TN − FP·FN)/√(…) with no guard, so a zero denominator gives NaN. That is the papers' convention. For comparability, the paper-protocol summaries must use NaN-excluded means.
   - old: `A line with no predicted or no true essentials has undefined MCC. It is set to 0, counted and reported; F1 and precision are treated the same way.`
   - new: `A metric with a zero denominator is NaN, as in the repository script; NaN lines are excluded from the means and counted per method and protocol. A sensitivity row sets them to 0. F1 = 2TP/(2TP + FP + FN), as in the script.`
   - Note: generic P57 calls are never empty, so this matters only for P1 models.
3. **CORDA cost bands: RATIFY.** The bands are ≤ 0.25 CPU-h per model (full A6), (0.25, 1] (R1 × real/permuted on 15 lines) and > 1 (drop, reported).
   - The arithmetic matches A6 (≤ 30 CPU-h either way).
   - Each band is a cost rule, not a label-informed choice. The 15-line fallback keeps a PROMPT-named method in the variance partition, analysed within method (H4.2).
4. **30-line LP-only Stage-A panel (10 × 3): RATIFY, conditional on N-1.**
   - The power limitation is already stated (§15, W4 *Power*).
   - Without the N-1 fold rule the panel gate is unpassable, so N-1 is required.
   - Optional, not required: 10 × 5 = 50 lines costs ≈ +65 CPU-h (≈ +2.5 h wall) and doubles the within-lineage replication for H4.2.
5. **Worker counts: AMEND.**
   - 28 LP workers at 1.45 GB plus 4 Gurobi sessions is ≈ 47 GB. Keeping ≥ 20 GB free then needs ≥ 67 GB available on a shared host that had ≈ 69 GB available (compute r2).
   - The 16-slot Gurobi pool also cannot serve 32 sessions.
   - §11.4 Workers bullet, new text: `**Workers** (PI decision r2-13 as amended by D2 and r3): there is no MILP lane. Start at **24 LP workers + ≤ 4 Gurobi sessions** (≈ 41 GB). All 28 run on Gurobi only if an M3 check confirms that the WLS licence admits 28 concurrent sessions, in which case the `flock` pool is set to 28; otherwise LP workers beyond the licence pool run GLPK, only for families that pass the GLPK = Gurobi equality check and whose M3 GLPK per-line cost is ≤ 2× Gurobi. Raising to 28 LP workers at M3 requires measured RSS to keep ≥ 20 GB free. A job aborts below 10 GB.`
   - Replace the "16 `flock` slots" bullet with `Gurobi slots: a flock pool sized to the licence session limit verified at M3.`
   - §13 Workers line: `24–28 LP workers plus ≤ 4 Gurobi sessions (§11.4)`.
6. **The Human2-paper gaps.** All four gaps are now filled from documented sources rather than "our choice":
   - **−0.6: RATIFY, now with a documented source.** The Human2 paper and SI state no threshold (confirmed: SI pp. 2–3, lines 56–65). The Human-GEM evaluation script used for repository releases sets `threshold = -0.6` (line 117), and its header says it is copied from the Human1 paper's code. In §6.8, replace "Our choice: the same L-ge06, which inherits the Human-GEM lineage's documented choice" with `Not stated in the paper; the Human-GEM evaluation script uses gene effect < −0.6 (line 117), i.e. L-ge06.` In the table row "Label (DepMap)", Human2 column, add `(repository script: < −0.6)`.
   - **Gene set U: AMEND.** The script scores **the reference (parent) model's genes ∩ screened genes with a non-NaN effect**, and counts parent genes absent from a line model as non-essential (lines 121–146). U-wide scoring would penalise each parent for the other parent's genes and break comparability with the papers.
     - old: `**Gene set (not stated in either paper; our choice).** Genes in U with a non-NaN label in that line. A gene absent from the model is predicted non-essential, as in §3.3.`
     - new: `**Gene set (not stated in either paper; taken from the repository script, lines 121–146).** The parent model's genes (H1 for generic H1 and P1(H1); H2 for generic H2 and P1(H2)) that have a CRISPR column and a non-NaN label in that line. A parent gene absent from a line model is predicted non-essential. U-wide scoring is a sensitivity row.`
     - Table row "Gene set scored": `not stated in the papers; repository script: parent-model genes ∩ screened, non-NaN`.
   - **Per-line means: RATIFY.** Robinson reports per-line distributions plus the mean MCC ("2.8-fold mean increase", main p. 6). Luo shows per-line violins (SI Figs. S2D, S3A; the summary statistic is not stated). Keep the line-weighted mean as the paper-comparable summary and the lineage-equal mean as ours. Add the per-line median as a third descriptive summary, because Luo's violins show a central bar.
7. **Other F-7 items: RATIFY.**
   - The P57 fallback (> 0.2 CPU-h per line model → generic parents and P1(H2) only) is a cost rule.
   - The W7 8-hour time box is a stop rule, not label-tuned.
   - The non-GEM dichotomisation at the fitting-line positive rate is training-only.

## E. Check of the §6.8 extracted protocols against the cited pages

| Claim in §6.8 | Source checked | Result |
|---|---|---|
| H1: tINIT per line; 621 DepMap 19Q3 lines; 5 Hart lines | SM p. 14 ("Evaluation with 621 human cell lines", 19Q3); SM p. 13 (5 lines); main pp. 5–6 | Correct (page nit N-6) |
| H1: KO rule; reactions without a compensating isozyme set to 0; essential if ≥ 1 of 57 tasks impaired (data file S3) | SM p. 14 | Verbatim correct |
| H1: biomass-only variant | SM p. 14 "Gene essentiality using only biomass production"; Fig. S11 p. 29 | Correct |
| H1: gene effect < −0.6, from the dataset's publication; grid −1.0 to 0.0 | SM p. 14; Fig. S10 p. 28 | Correct. The grid is in 0.05 steps (D1). |
| H1: metrics are accuracy, sensitivity, specificity, MCC, hypergeometric + BH; MCC headline; accuracy uninformative | main p. 6; SM Fig. S9 p. 27 | Correct |
| H2: ftINIT "1+1"; 57 tasks; growth on Ham's | SI p. 2, lines 25–33 | Correct |
| H2: same KO call rule | SI pp. 2–3, lines 58–62 | Verbatim correct |
| H2: 625 DepMap 19Q3 lines; 908 Project Score lines | SI p. 3, lines 64–65; SI Fig. S3A p. 11 | Correct. The 625 is "RNA-Seq data from 625 cell lines" (SI line 64). |
| H2: threshold not stated | SI pp. 2–3; main Methods p. 9 | Correct for the paper. The repository script supplies it (D6). |
| H2: metrics MCC (Fig. 2D legend p. 3; SI S2D p. 10; S3A p. 11); "accuracy" (Methods p. 9) | as cited | Correct. Note that main Fig. 2D and the Methods "accuracy" refer to the Hart-line release checks, not DepMap (N-6). |
| Both: medium for task checks = each task's own inputs; GR = Ham's | SI line 30; Robinson data file S3 | Correct |
| Both: gene set "not stated" | papers | Correct for the papers, but the repository script states it (D6) |

## F. Freeze checklist
To freeze, apply the following, then re-run `scripts/design_numbers.py` (no data fact changes):
- N-1, N-2/D5, N-3 and N-7/D1, D2, D6 (blocking);
- N-4, N-5, N-6 and N-8 (minor);
- N-9 (register note, PI).
