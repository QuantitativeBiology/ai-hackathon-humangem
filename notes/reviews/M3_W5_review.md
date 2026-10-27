# M3 W5 (GPR curation) — independent non-author review

Reviewer: independent agent (Claude Opus), 2026-10-27. Not the W5 author. Scope: DESIGN.md v3 §9 W5 (+ §5–6, §10.3, §17 A1–A10), `results/w5/` and the uncommitted W5 code (`src/nextgem/w5.py`, `src/nextgem/w5_evidence.py`, `scripts/w5_*.py`, `tests/test_w5.py`). I did not edit DESIGN.md, FROZEN.tsv, the ledgers or the author's code, and I committed nothing.

**Hygiene.**
- Session-start probe `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **BLOCKED** by the label hook.
- Labels were read only through `nextgem.labels.train()`, in `results/w5/review/r2_train_score.py`. `labels.ids()` (ModelIDs only, no values) was also used, to define the screened lines for the expression check.
- I never called `validation()`/`test()`, never touched `data/locked/`, and never opened the De Kegel/Ito tables.
- GLPK only, one process (W3 was running).
- Outputs: `results/w5/review/` (148 KB) and this file.

Reviewer code (independent of `nextgem.w5`/`w5_evidence`):

| Script | What it does |
|---|---|
| `r1_labelfree.py` | Raw-file evidence check; own GPR edit application on cobra's GPR AST; KO sweeps |
| `r2_train_score.py` | Own tie-aware AP; LOLO Δ̄_lin |
| `r3_cad.py` | Label-free CAD check |

The edit-level table is in `r4_edit_review.tsv`.

## Verdict summary

| Item | Verdict |
|---|---|
| 1. Edit review (20% sample + ≥5-flip edits; 28 edits) | **Evidence mostly real, rule mostly followed, biology often wrong.** 6 accept / 14 reject / 8 questionable. All 12 POLR merges are rule-literal but biologically wrong (reject). CYP1A2 removal (both parents) fails its own evidence rule: the gene **is** detected in Nusinow. 8 removals of ER-membrane enzymes violate the "soluble" qualifier. |
| 2. Independent recompute | **Reproduced exactly.** The curated closure sets equal the author's (0 mismatches on 341/405 reactions). KO ratios are identical (max \|Δ\| = 0). Δ̄_lin matches to 5 decimals. TG fails for every arm. The POLR attribution is confirmed. |
| 3. The 13 deviations | All are defensible or immaterial to the conclusion. Two need README corrections (the soluble proxy; the detection test). The paralog reading is literal but exposes a DESIGN-rule flaw. |
| 4. CAD claim | **Confirmed** label-free in both parents. |
| 5. Leakage | **Clean.** W5b is fitted on `labels.fold(k)` inside `cv.lolo`. The scanner found 0 hits over `results/w5` (82 files, incl. this review). |
| **Overall** | **GO** for committing W5 as the M3/M4 record (negative result: shortlist empty), with the README corrections listed under "Required before commit". |

## 1. Edit review (DESIGN W5 "Review")

**Selection.**
- Random 20% per parent: `numpy.random.default_rng(0)`, one generator per parent, `choice(n, ceil(0.2 n), replace=False)`.
  - H2: 12 of 60 edits.
  - H1: 8 of 38 edits.
- **"Flips ≥ 5 genes".** I applied each edit alone to its parent (H1 on the TYMS base).
  - **No single edit flips more than 3 genes:** POLR merges flip 2–3, RRM 2, SEC11 2, and IDI2/LENEP/CYP7A1 removals 1 each (`r1_edit_flips.tsv`).
  - At the GPR-group level, the W5a POLR edits (group MAR07161/2) jointly flip 15 genes. I therefore added every POLR edit of both parents to the review set.
  - PAPOLB sits in the same group, so its "group flips" column shows 15, but it flips 0 alone.
- Reviewed set: **28 edits** (H2 16, H1 12).

**Evidence check** (`r1_evidence_check.tsv`, all 98 edits, raw files, own parsers):
- **Paralog pairs.** Every Ensembl 116 paralog pair cited exists. I give the type, min %id and LCA for each.
- **No shared complex.** No AND→OR pair shares a CORUM or Complex Portal complex. I matched by UniProt accession or symbol, with no organism filter, so this check is stricter than the author's.
- **Expression.** TPM < 1 in ≥ 95% of lines was recomputed from the raw DepMap 26Q3 CSV (default entries).
  - It holds for all 27 removal genes, both over the 1,292 screened lines and over all 2,130 lines.
  - Note: the author's validator only checks that the Entrez column exists, not the fraction. The fraction itself is correct.
- **Never detected in either proteomics set.** This leg fails for:
  - **CYP1A2.** Nusinow row `sp|P05177-2|CP1A2_HUMAN`: 9 of 378 channels quantified, up to 5 peptides. `EvidenceDB.protein()` and the validator match exact accessions only, so isoform rows (`P05177-2`) are missed. → 3 edits (W5a_H2_010, W5a_H2_019, W5a_H1_010) violate the rule.
  - **UGT2A1** (ambiguous; W5a_H2_029, W5a_H1_020). Nusinow row `Q9Y4X1-4`, 9 channels, 1 peptide. Q9Y4X1 is the pre-2020 accession that UniProt split into UGT2A1 (P0DTE4, the model's ID) and UGT2A2. The peptide may come from the shared UGT2A1/2A2 region. → Questionable.
- In total, 93 of 98 edits pass my evidence check; the 5 failures are the CYP1A2 and UGT2A1 removals. **None of these 5 edits changes any KO ratio**, so the scores are unaffected.

**Rule check: the "soluble enzyme" qualifier.** The proxy ("no 'Cell membrane' in genes.tsv, no transport reaction") accepts ER-membrane enzymes as soluble:
- CYP2A6/2A7/2A13/2C8/2C19/2F1/3A43/4Z1/26C1/1A2/7A1;
- UGT2B10/2B28/2A1/3A1;
- TECRL and ACER1.

Microsomal P450s have an N-terminal ER membrane anchor (Ruan et al. 2002, [doi:10.1042/BJ20021001](https://doi.org/10.1042/BJ20021001)). UGTs are ER-membrane enzymes (Meech et al. 2019, [doi:10.1152/physrev.00058.2017](https://doi.org/10.1152/physrev.00058.2017)).

These removals are biologically harmless, because TPM < 1 in ≥ 95% of lines is itself strong absence evidence. But they do not satisfy the DESIGN text, which requires "soluble" because proteomic non-detection is weak evidence for membrane proteins. Only one of them changes a KO ratio: **CYP7A1** (MAR01510, lathosterol oxidase `SC5D or CYP7A1`) makes SC5D lethal. That flip goes in the wrong direction for most lines; it is a sterol/medium artefact, see the W2 review.

**Edit-level table** (full notes and citations in `results/w5/review/r4_edit_review.tsv`):

| Edit(s) | Genes | Evidence real | Follows rule | Verdict | Reason (citations via PubMed) |
|---|---|---|---|---|---|
| H2_002–007, H1_002–007 (12) | POLR1/2/3 paralog sets | yes | yes, literally | **reject** | The subunits are paralogous across Pol I/II/III, which share a common ancestry (Werner & Grohmann 2011, [doi:10.1038/nrmicro2507](https://doi.org/10.1038/nrmicro2507)). But each polymerase is a distinct, non-redundant complex, and MAR07161 is a lumped pseudo-reaction that needs all three, so AND is correct. The merges make 13 U genes viable; 12 of them are near-universally essential. |
| H2_010, H1_010 | CYP1A2 | **no** (detected in Nusinow) | no | **reject** | Also ER-membrane; no KO effect |
| H2_015/017/027/032/042/048, H1_011/017 | CYP2A7, CYP3A43, UGT2B28, CYP26C1 (×3), CYP2C19, UGT2B10 | yes | no (soluble leg) | questionable | ER-membrane enzymes; biologically harmless; no KO effect |
| H2_039 | RRM2B / RRM2 | yes | yes | **accept** | RRM1 + (RRM2 or p53R2/RRM2B) is the established logic (Tanaka et al. 2000, [doi:10.1038/35003506](https://doi.org/10.1038/35003506)); a priori §1.3. It costs AP because RRM2 becomes viable. |
| H2_001, H1_001 | PAPOLB | yes | yes | accept | Testis-type cytoplasmic PAP (TPAP; Nakanishi et al. 2005, [doi:10.1016/j.ydbio.2005.10.017](https://doi.org/10.1016/j.ydbio.2005.10.017)); TPM < 1 in 100% of lines |
| H2_055 | FTMT | yes | yes | accept | Mitochondrial ferritin, soluble and tissue-restricted (Levi et al. 2001, [doi:10.1074/jbc.C100141200](https://doi.org/10.1074/jbc.C100141200)) |
| H1_021 | RBP3 | yes | yes (soluble, but a binding protein, not an enzyme) | accept | Retinoid-binding protein in a dehydrogenase OR; removal is justified on function too |
| H1_027 | LENEP | yes | yes | accept | MAR06508 is FMNAT (FLAD1). LENEP has no such activity. Removal is right; FLAD1 becomes lethal, which is a model consequence. |

Totals over the 28 reviewed edits: H2 3 accept / 7 reject / 6 questionable; H1 3 accept / 7 reject / 2 questionable.

**Score-changing edits outside the sample** (screened for completeness):

| Edit | Verdict | Reason |
|---|---|---|
| SEC11A/SEC11C → OR | accept | Two paralogous SPC complexes (Liaci et al. 2021, [doi:10.1016/j.molcel.2021.07.031](https://doi.org/10.1016/j.molcel.2021.07.031)) |
| COX7B/COX7B2 → OR | accept on biology | Tissue-specific COX isoforms (Sinkler et al. 2017, [doi:10.1155/2017/1534056](https://doi.org/10.1155/2017/1534056)). COX7B2 is unexpressed in cancer lines, so the edit costs AP. |
| HADHA/HADH → OR in `HSD17B10 and HADHA and HADH` | questionable | The direction is right: three 3-hydroxyacyl-CoA dehydrogenases are AND-ed, an isozyme error. The merge is partial and rests on a 13%-identity `other_paralog`. |
| HADHB/ACAA2 → OR (H1) | questionable | The correct logic is `(HADHA and HADHB) or ACAA2`; HADHA/HADHB form the membrane-bound trifunctional protein (Carpenter et al. 1992, [doi:10.1016/0006-291x(92)90501-b](https://doi.org/10.1016/0006-291x(92)90501-b)) |
| IDI2 removal | accept | Muscle-restricted (Clizbe et al. 2007, [doi:10.1074/jbc.M610922200](https://doi.org/10.1074/jbc.M610922200)); IDI1 becomes lethal |

## 2. Independent recompute

**KO ratios** (`r1_*`):
- I applied the W5a edit list myself to the released GPR strings (cobra GPR AST).
- The resulting single-gene closure sets equal the author's curated rules on all 341 (H2) and 405 (H1) reactions: **0 mismatches**.
- KO ratios (cobra semantics, GLPK, M0) for every gene of an edited GPR:
  - generic H2 equals `results/baseline` (Gurobi);
  - W5a equals the author's `ratio_W5a`, max |Δ| = 0, both parents.
- My full H1+TYMS sweep equals the author's `ratio_base` (max |Δ| = 0, 2,887 genes).
- Genes changed by W5a: 15 POLR, RRM2, RRM2B, SEC11A, SEC11C, COX7B, COX7B2, IDI1, SC5D, FLAD1, HADH, BAAT and AGXT (sub-threshold), plus HADHB and ACAA2 in H1.

**Training LOLO Δ̄_lin** (`r2_*`):
- Setup: L_prob, all U genes, 686 expression lines of the 12 training lineages, 95% CI from a lineage-stratified bootstrap (B = 2,000, SEEDS[0]).
- tAP is my own implementation. It was checked against brute-force enumeration of within-tie orderings (max |Δ| = 2e-16, 300 toys) and against `metrics.tap` (2e-16).

| Contrast | Reviewer Δ̄_lin [CI] | folds > 0 | Author | TG |
|---|---|---|---|---|
| H1_TYMS − H1 | +0.00069 [+0.00040, +0.00094] | 9/12 | +0.00069 [+0.00042, +0.00096], 9/12 | fail (< 0.005) |
| H1_W5a − H1 | −0.02855 [−0.02894, −0.02819] | 0/12 | −0.02855 | fail |
| H2_W5a − H2 | −0.02928 [−0.02950, −0.02906] | 0/12 | −0.02928 | fail |
| H1_W5a − H1_TYMS (GPR part) | −0.02924 | 0/12 | −0.0297 (GPR on released H1) | fail |
| *diagnostic:* H1_W5a without POLR merges − H1 | −0.00193 | 1/12 | (pid30 arm −0.0021) | — |
| *diagnostic:* H2_W5a without POLR merges − H2 | −0.00299 | 0/12 | (pid30 arm −0.0031) | — |

L_ge05 gives the same signs: W5a −0.0301 (H1) / −0.0304 (H2); TYMS +0.0004, 8/12.

The POLR merges account for −0.0266 of −0.0286 (H1) and −0.0263 of −0.0293 (H2). This confirms the author's attribution. The "without POLR" row is a reviewer diagnostic, not a candidate. W5b = W5a was confirmed from the code path and the empty `w5b_fold_edits.tsv`; I did not refit it.

## 3. The author's deviations (README §11; 13 items listed, the brief says 12)

| # | Deviation | Defensible reading? | Could it change "no arm passes TG"? |
|---|---|---|---|
| 1 | Edit unit = canonical GPR group | Yes. The same gene logic is edited in every reaction that carries it, and the per-reaction reading would score the 14-reaction RNR fix at 0 (I checked that MAR02352 alone is non-essential). | No. MAR07161 alone is essential (closing it gives ratio 0; MAR07162 gives 1), so the POLR edits would still be top targets per reaction. RNR would drop out, which removes both the RRM2 cost and the RRM2B gain. |
| 2 | "Top 200" = groups with any KO effect (109/92 groups = 928/1,064 reactions) | Yes. Below that, every group ties at 0, and only effect-bearing edits can move scores. | No. Groups without effect cannot change any score. |
| 3 | Ranking key (union of flips) | Yes; the DESIGN does not specify one | No: all effect groups are included anyway |
| 4 | AND→OR as pairwise merges applied as cliques | Yes. This is the minimal edit; the node-level alternative is more destructive. | No |
| 5 | Ensembl paralogs taken literally (incl. `other_paralog`, 13–26% id) | **A literal and defensible reading of the frozen text, but it exposes a DESIGN-rule flaw.** "Paralogs and not co-members of any complex" is not evidence of substitutability when the AND spans several distinct complexes in a lumped pseudo-reaction (§1.4d already flags MAR07161 as a pseudo-reaction). The ≥30%-identity variant is not the fix, and it was correctly kept as a non-eligible sensitivity: it also drops HADHA/HADH, and the POLR1B/2B/3B triplet is blocked only by the clique rule. | Not in the direction of passing: the pid30 / no-POLR arms are still −0.002 to −0.003. |
| 6 | "Soluble enzyme" proxy | **Partly defensible, but under-stated.** No UniProt topology file is on disk, so some proxy was needed. However, the proxy admits about 17 ER-membrane enzymes (CYPs, UGTs, TECRL, ACER1, CYP7A1). The README must say so. | No. Only CYP7A1 changes a score (SC5D → lethal). A strict reading would remove that one wrong-direction flip, which is worth far less than the +0.007–0.008 needed. |
| 7 | OR→AND requires all cross-branch pairs | Yes; this is the literal reading of "the branch genes" | No (0–1 eligible nodes) |
| 8 | W5b: co-essentiality replaces the co-abundance leg; "loose" variant as sensitivity | Yes. "Supports OR→AND" is ambiguous, and both readings were run. | No: 0 edits in all folds under both readings |
| 9 | Protected list verbatim; pan-cancer-absence removal not applied inside protected ORs (PRPS1L1) | Yes, the conservative choice. The §3 isoform note names specific genes, and PRPS1L1 is not among them. | No: PRPS1/PRPS2 still buffer, so no flip |
| 10 | H1 enumerated on the TYMS base | Yes (§3.1). The GPR-only row is reported separately. | No |
| 11 | `allow_dirty` provenance | Acceptable while commits are forbidden. I checked that every code SHA-256 in the final MANIFEST rows equals the current files (0 mismatches). The PI should commit exactly these files. | No |
| 12 | Negative-control pool; random AND→OR not paralog-restricted | Yes (matched on type and number) | No |
| 13 | Tests run with `NEXTGEM_NO_TEST_RECORD=1` | A process item; the PI's full run must record them | No |
| — | **Missing UniProt "unrelated function" removal path** | It is an admitted gap, not a reading. The DESIGN offers it as an alternative removal route. | It cannot rescue literal W5a (−0.029). Removing spurious OR branches only makes genes more essential (cf. LENEP→FLAD1, CYP7A1→SC5D). It could add true positives only where common-essential genes are buffered by mis-annotated partners. At ≈ 0.002 per gene (POLR: 0.026 over 13 genes), it would need ≳ 4 such genes just to bring a no-POLR arm to +0.005. Such an arm is not the pre-registered W5a. The record should say "W5a as specified (minus the UniProt path) fails"; it should not say "GPR curation cannot help". |

**Additional issue (not in the author's list).** Detection is tested by exact UniProt accession only (`EvidenceDB.protein`, `validate`). Isoform-suffixed rows (`P05177-2`) and superseded accessions (Q9Y4X1 for UGT2A1) are therefore missed. The validator also does not re-check the ≥ 95% TPM fraction. These change no score here, but they should be fixed or documented before any reuse of the code.

## 4. CAD claim (label-free, `r3_cad.json`)

Confirmed in both parents (GLPK, M0):
- MAR10064 ("Small metabolite pool for biomass reaction") consumes (S)-dihydroorotate MAM00180c at **0.0043**. Its product MAM10015c enters MAR13082 at 0.4835.
- Cytosolic dihydroorotate is made only by MAR04575 (DHODH) and MAR04608 (CAD), both reversible. H2 also has transport MAR20201 and exchange MAR20202, whose lower bound is 0 (closed).
- Ham's + uridine (H2 WT 92.13, H1 82.96):
  - CAD KO ratio = 0;
  - DHODH and UMPS KO ratios = 1;
  - **after deleting dihydroorotate from MAR10064, CAD KO = 1.0 in both parents.**
- Ham's + orotate: CAD KO = 1.0 (DHODH runs in reverse).
- Side finding: in released H1 on plain Ham's, CAD lethality is also due solely to the pool (pool removed → ratio 1), because DHODH/UMPS are non-essential there.

So the CAD lethality under uridine is a biomass-composition artefact, not a GPR or direction error, as the author states. The fix belongs to W3/biomass, not W5.

## 5. Leakage

- **W5b fitting.** `scripts/w5_score_train.py::fit_predict` receives `fit_long` from `cv.lolo`, which calls `labels.fold(k)`, and computes Chronos r only on it. Held-out labels reach only the harness `score` callback, which returns the predictions unchanged. `cv.lolo` asserts that held-out lines do not overlap the fitting lines.
- **Candidates are label-free.** `w5b_candidates.json` comes from `w5_curate.py`, which has no label import.
- **Dry run.** `w5_val_dryrun.py` uses `labels.train()` (the §5.3 fitting set for a validation look) and reports `validation_labels_read: false`.
- **What was written.** Only aggregate values: all-training pair r, and training positive fractions for about 100 changed genes. No gene × line values.
- **Scanner.** `scripts/scan_label_copies.py --roots results/w5`: **0 hits**, 82 files (W5 outputs and this review's files).

## Required before commit (non-blocking for the negative conclusion)

1. README §3/§11 must state that the soluble proxy admits ER-membrane enzymes (list them). The CYP1A2 removal (and UGT2A1) fails the "never detected" leg because of isoform-accession matching. No score is affected.
2. README "Why curation hurts" should add that the DESIGN AND→OR rule is unsafe for lumped pseudo-reactions. Recommend LESSONS.md or a §17 note before any W5 re-specification; any change to the rule is an amendment, not a W5 re-run.
3. README should say "W5a as specified, without the UniProt path". The PI commits the exact files whose SHA-256 values are in `MANIFEST.tsv`, and runs the W5 tests with recording on.

**GO**: commit W5 as the M3/M4 record with these corrections. The headline (no curated arm passes TG; shortlist empty; Δ_curation reported) is reproduced exactly and is robust to every deviation reviewed.

## Runtime and exposure

- **Runtime.** Reviewer wall ≈ 25 min. Compute: r1 121 s, r2 10 s, r3 48 s, scanner 1 s; ≈ 0.05 CPU-h. GLPK only, one process, no Gurobi session.
- **Label exposure.** Training-tier aggregates only:
  - LOLO Δ̄_lin, CIs and fold counts for the contrasts above (L_prob, L_ge05);
  - the author's training aggregates as written in `results/w5/README.md`.
- No gene × line label value was printed or written. No validation- or test-tier label or sealed output was read.
