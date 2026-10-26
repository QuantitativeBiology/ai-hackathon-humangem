# Design review, round 1: external PI

Reviewer: external PI (systems biology / cancer metabolism), not an author of the design.
Reviewed: DESIGN.md v0.1 (2026-10-26) against PROMPT.md. I also skimmed LITERATURE.md (§0, §2, §4.3–5, §8–9), DATA.md, RESEARCH_TREE.md, results/baseline/README.md and LESSONS.md.
Angle: whether the science is worth doing and sound; whether the design covers PROMPT.md; whether the endpoint and claims are proportionate; order of work; what would get the paper rejected.
I did not read or compute any CRISPR label values. I read only `ScreenSequenceMap.csv` metadata (library and screen-type counts) and the `data/depmap/` file listing.

## Overall verdict

This is a strong, unusually careful design. It has locked labels, a test ledger, label-free GPR targets, negative controls in every workstream and an honest gate fallback. It covers every part of PROMPT.md. My main concerns are not leakage. They are about **what the headline can and cannot show**:

1. The primary endpoint can be won by an artefact of score resolution. Generic FBA scores have 24–32 distinct values, and about 2,650 of ~2,850 genes are tied at e = 0. Any pipeline that breaks those ties with an expression-correlated quantity will gain mAP_line without any metabolic insight.
2. "CI lower bound > 0" on ~230 paired lines will declare trivially small gains a success.
3. "Transfer to unseen tissue" is claimed from 5 lineages.
4. The decomposition the PROMPT asks for (medium alone; curation vs constraints) is estimated on training or along an order-dependent ladder, not on test.

All four can be fixed with text changes and modest compute.

---

## Must-fix

### MF1. Add a trivial "generic + expression" baseline, so P1's win cannot be a tie-breaking artefact (§6.1, §7, §6.4, §9.3)
- **Problem.** results/baseline/README.md shows that generic FBA is near-binary: 154 lethal genes (H2), 2,625 genes at ratio 1, and only 32 distinct values.
  - mAP_line(G\*) is therefore essentially "precision of 154 calls, then a 2,600-gene tie".
  - Every context method that produces continuous scores, especially F4 EWB (bounds scaled by expression), will spread that tie by an expression-correlated ordering. Unexpressed genes also tend to have low dependency probability (see MF4).
  - C1 can then pass because "expression ranks the tie", not because the GEM got better. B3 (expression alone) does not rule this out, because it is not combined with the generic call.
- **Change.**
  - Add **B7 = lexicographic(e_G\*, B3)**: rank by the generic essentiality score, with ties broken by expression. It has no fitting.
  - Optionally add **B7′ = e_G\* with genes at TPM < 1 set to 0** (the "remove unexpressed genes, no network" reference; F0 is close but includes GPR closure).
  - Add B7 to the T1 set as a reported contrast (estimate + CI), next to C1 and C2.
  - Pre-declare the claim language: "the GEM optimisation adds information beyond simply combining the generic model with expression" requires P1 > B7. Otherwise the headline is "beats generic FBA, but no more than generic FBA + expression".
  - Also promote the fixed-operating-point binary contrast (MCC at ratio_r < 0.01, P1 vs G\*) from "secondary for FBA methods" to a reported T1 estimate. It is insensitive to tie resolution.

### MF2. Fix the confirmatory structure and add a smallest effect size of interest (§0, §6.3, §6.4 table, §9.1 gate)
- **Problems.**
  - Holm over C1–C5 mixes two looks taken days apart (T1, T2), so the T1 headline cannot be finalised until T2.
  - C2 is two-sided and *expected* to go against P1 (H1.3: the gene prior beats FBA), so it does not belong in a family of success claims.
  - "Lower bound > 0" with paired line resampling on ~230 lines will declare a Δ mAP_line of 0.003 a success. A reviewer will ask whether that is biologically meaningful.
- **Change.**
  - T1 family: fixed-sequence testing **C1 → C3**, each at α = 0.05 (C3 is tested only if C1 passes). This answers Q1 without splitting α.
  - C2 and B7 contrasts become estimates with CIs, not tests.
  - T2 family: fixed-sequence **C4 → C5**, as a separate family.
  - Pre-specify a SESOI before M2. It must be fixed a priori, so either set an absolute value (I suggest Δ mAP_line = 0.01) or define it as 10% of the training-tier gap between G\* and B2, measured at M2 on training labels only and frozen before any candidate is scored.
    - Report "**headroom closed**" = (P1 − G\*) / (B2 − G\*) with its CI.
    - Tiered language: "detectable improvement" if CI_low > 0; "meaningful improvement" if CI_low > SESOI.
  - The abstract-level sentence must also state the minus-common-essentials Δ (already computed in §6.5). That stops readers from attributing a gain on DepMap pan-essentials to context modelling.

### MF3. Make the transfer-to-unseen-tissue claim proportionate to 5 test lineages (§0, §4.1, §6.4, Q2)
- **Problem.** C1 resamples lines within lineage. That supports the claim "P1 improves on lines of these 5 held-out lineages", not "transfers to unseen tissues". The two-stage bootstrap over 5 lineages will be close to uninformative, and the design already says so.
- **Change. Pre-declare three claim tiers in §6.3:**
  1. "Improves on held-out lineages": C1, line bootstrap.
  2. "Consistent across held-out tissues": Δ_lineage > 0 in **≥ 4 of 5** test lineages, *and* in ≥ 75% of the 14 training LOLO folds. Flag the LOLO folds as selection-contaminated, but they are the only multi-tissue replication available.
  3. "Generalises to new tissues": only if the two-stage bootstrap CI_low > 0. State in advance that this is unlikely to be powered.
- Report a forest plot of per-lineage Δ for the 5 test lineages, with per-lineage CIs and prevalence.
- Keep 1 + 2 + 2 by stratum and constraint (d). Both are label-free and justified (RPMI/suspension biology). Do not enlarge the test tier at the cost of validation; with only 4 validation lineages, selection is already noisy.

### MF4. Context-specificity claims must hold on expressed gene × line pairs (§3 strata, §6.4 C3, W4 TG and decision, Q4)
- **Problem.** Chronos uses unexpressed genes as negative controls (LITERATURE §2.3 item 7), so unexpressed pairs carry low dependency partly *by construction of the label*.
  - Any extraction family that removes unexpressed genes collects these easy negatives.
  - mΔAP_gene can then be positive with no metabolic insight: the per-gene ranking across lines just tracks expression on/off.
  - Today "expressed pairs" is a reporting-only stratum.
- **Change.**
  - Keep U as the primary universe (answer to Q4).
  - Make the **expressed-pairs (TPM ≥ 1) estimate of C3 a required consistency condition**: the claim "context-specific models capture line-specific dependency" needs C3 > 0 overall *and* the same direction with CI_low > 0 on expressed pairs.
  - Add the expressed-pairs Δ to the W4 TG report, so shortlisting does not silently favour pure on/off pruning.
  - Apply the same condition to the W4 level comparison (per line vs lineage vs background).

### MF5. Measure the attribution the PROMPT asks for on test, order-independently (§9.3 T1 "ablation ladder", W5 ablation, W2 metrics)
- **Problem.** PROMPT asks for two quantities: "how much baseline error is attributable to medium misspecification alone" and "how much improvement comes from curation rather than from the constraints". The design answers them only partly.
  - W5's 2 × 2 factorial runs on *training*, where the candidates were selected, so the estimate is optimistic.
  - The T1 ladder (generic → +medium → +tasks → +curation → +extraction) is sequential. Its increments depend on the order, and the order is the one most favourable to the late levers.
- **Change.**
  - Include in T1, as part of the same pre-declared look (all cells are frozen, label-free builds), the **full factorial of P1's frozen levers**: ≤ 2^4 = 16 cells on test lines.
  - Report the main effects, the curation × constraints interaction, and **Shapley attributions** with paired-bootstrap CIs. Also report the sequential ladder, for readability.
  - The "medium alone" number must be the M0 → M_doc contrast on the **generic** parent on test (mappable-medium lines), with the W2 binary-flip attribution. This is cheap: KO runs per distinct medium, not per line.
  - Cost: ≤ 16 × ~230 lines of KO sweeps; most cells are generic (per medium), so this is well under 50 CPU-h.

### MF6. Human1 must be a full arm at T1, not a conditional replication (§2.1 "default parent", W8 "released if…", §9.3)
- **Problem.** PROMPT fixes both Human1 and Human2. The design optimises on H2 only, and the H1 replication of P1 is not in the T1 list.
  - If G\* turns out to be H1, C1 compares P1(H2) with generic H1. That confounds the parent model with the optimisation.
- **Change.** T1 must include **P1 applied to the H1 parent** (same frozen pipeline; W5 edits are already enumerated in H1). Report the within-parent contrasts P1(H1) − generic H1 and P1(H2) − generic H2, next to C1 (estimates).
  - Cost: about 230 test lines × ~6 CPU-min ≈ 25 CPU-h.
  - Releasing H1 SBML models is covered in SF8.

---

## Should-fix

### SF1. A per-gene baseline report that is not degenerate (W1, §6.2 AP_gene)
PROMPT asks for the baseline benchmark "per gene". For any generic model, AP_gene = π_g exactly, and the design notes this. Report it as required, but add a per-gene table that carries information: for each gene, the generic call (lethal / non-lethal), the fraction of lines in which it is essential (PPV_g for called genes, miss rate for uncalled genes), and the error class from §9.4. That table is what readers and curators will actually use.

### SF2. Library replication from data already on disk (§3, §9.3; LITERATURE §9(a) item 2 recommends Project Score and the design dropped it silently)
`ScreenSequenceMap.csv` lists Avana 2,280, KY 906 (Sanger Project Score library) and Humagne-CD 398 screens, and `ScreenGeneEffect.csv` holds per-screen effects.
- Add a pre-declared T1 sensitivity: C1 and C3 recomputed with labels from **KY-only** and **Avana-only** screens (fixed rule: screen gene effect < −0.5, mean over passing screens) for test lines that have them.
- This tests whether gains are specific to one library or assay, which is the main batch effect in this field (LITERATURE §2.2).
- It costs no downloads.

### SF3. Pre-registered known-answer cases for genomic-background models (W4 level 3, §9.4)
The background groups are well chosen, but nothing states what the models *should* predict. Before any scoring, record label-free predictions for literature cases:
- **expected to be captured by FBA:** ENO1 deletion → ENO2 (paralog collateral lethality); ME2 loss → ME3; ASS1 silenced → arginine-import dependency;
- **expected to be missed** (the mechanism is outside FBA): MTAP deletion → MAT2A/PRMT5, which goes through MTA inhibition of PRMT5; IDH1/2 hotspot (neomorphic 2-HG).

These cases then serve as interpretable case studies in the report and in the failure account. They also pre-empt the reviewer question "does your genomic constraint do anything biologically sensible?"

### SF4. W6 needs a "network mapping without simulation" comparator (W6 comparators)
Comparing X_KO/X_flux with X_RNA/X_prot confounds two things: (a) aggregating genes onto reactions through GPRs and the network, and (b) actually solving for fluxes.
- Add **X_rxnRNA** (and X_rxnProt): GPR-mapped reaction expression (AND = min, OR = sum, the same mapping as F4), with the same PCA-20 and the same learner.
- Claim "flux simulation carries information beyond its inputs" only if X_flux/X_KO beats X_rxnRNA, not just X_RNA. Otherwise a positive W6 result will be attributed to feature engineering.

### SF5. Define the medium attribution number explicitly, with a ceiling (W2 Metrics)
"Fraction of baseline error attributable to medium" needs a denominator. Add a **training-only medium oracle**: for each line, the best ladder variant chosen with its own labels. This bounds how much error *any* medium in the family could explain. Then report:
- (documented − M0) / (oracle − M0);
- the binary-flip fractions;
- ΔAP.

This gives the paper one interpretable number ("documented media recover X% of the medium-explicable error, which is Y% of total baseline error").

### SF6. Order of work and compute (§12, §13)
- Start the **W4 variance panel immediately after M1**, in the background on spare Gurobi sessions. It is label-free, it answers a PROMPT question independent of performance, and at ~8 h it is the longest single block before M5.
- Build W4 per-line models for **validation and test lines only for shortlisted configurations**. Builds are label-free, but building 4 configurations × 1,361 lines before TG wastes ~200 CPU-h.
- Build the **W8 resource right after T1**, in parallel with W6/W7, not after T2. The resource depends only on P1. Putting it last makes the deliverable hostage to a W7 overrun.
- **Cuts if time runs short**, in this order:
  1. Human1 v1.3.0 sensitivity (different gene namespace, low value; one sentence in Limitations suffices);
  2. M3 concentration-proportional bounds (no uptake-rate data; LITERATURE §4.3);
  3. O-NGAM to one level;
  4. per-line variance panel lines reduced from 100 to 60.

  Do not cut the lineage × background combinations; PROMPT names them.

### SF7. Report the "addressable" ceiling (§6.5, W1 failure tests)
About 1,750 H2 genes close no reaction when knocked out, so they are unpredictable by construction. Pre-declare a label-free "addressable set" (KO closes ≥ 1 reaction in the parent) and report:
- mAP_line on it;
- the fraction of positives outside it, i.e. the recall ceiling of any FBA method on that parent.

This explains where gains can come from (W5 curation is the lever that makes genes addressable), and it pre-empts "GEMs are just bad at this".

### SF8. Per-task breakdown at task level, and both parents in the resource (W3, W8)
- PROMPT asks for a "per-task" breakdown. The design gives 11 groups (T4_k) plus T5. Add an individual-task table for the 57 Essential tasks: genes whose KO fails the task (label-free), and the precision/AP of those calls on training and on test. Label it exploratory.
- W8: release H1-parent models as well (≈ 245 CPU-h, so schedule after the H2 release), or state in REPORT that they were not released and why. PROMPT says "use … Human1 and Human2" and asks for the generated models in SBML. A silent conditional is not enough.

### SF9. Housekeeping on labels (§2.2, §3, Q3)
- `data/depmap/CRISPRGeneDependency.csv` is now on disk (timestamp 2026-10-26 14:50), so the L-prob branch of the freeze rule applies. Remove "(pending)" in §2.2 and resolve Q3. Run the §3 verification (ID sets, SHA-256, training-only κ) at M1.
- Q3, common essentials: add a training-derived common-essential list (same DepMap rule, training lines only) as a sensitivity stratum. It is cheap and removes a reviewer objection.

---

## Optional

- **O1. Enzyme constraints.** Proteomics is used only as presence/absence (W4) and as GPR evidence (W5). A GEM-community reviewer will ask "why not GECKO/ecHuman?". Either add it as a logged W7 idea with a 10-line pilot under the existing cost cap, or state in §15/Limitations why it is out of scope (LITERATURE §3.1: no evidence that it improves binary essentiality). Note that F4 with the protein modality is a crude capacity-bound stand-in; say so.
- **O2. Validation looks.** 60 looks on 4 validation lineages (~170 lines) is generous. Consider 40 in total, with ≤ 20 for W7. The one-SE rule and shrinkage reporting partly mitigate this.
- **O3. Q7 (F5).** Keep F5 in W4, but tag it "label-calibrated" in every table and in P1's attribution. Exclude it from the method-vs-tissue panel; it already is. Moving it to W7 is acceptable but not necessary.
- **O4. Q6 (ftINIT).** Pooled AUPRC on 5 Hart lines is a weak acceptance test, because pooled metrics are dominated by gene identity. If it can be generated, add a reaction-content Jaccard against a MATLAB/RAVEN reference build for the same lines; otherwise keep the current test and flag it.
- **O5. Q8 (gate).** Keep the fallback ("gate not passed", still tested). A strict stop would make a negative result unreportable, which is worse science.
- **O6. Q5 (serum).** The a priori list is adequate. Keep the ×0.1/×10 sensitivity and do not add more serum variants.

---

## What would get the resulting paper rejected at a top journal (and which item addresses it)

1. "The gain over generic FBA is what you get from adding expression or gene identity to a near-binary score." → MF1, MF2 (headroom, minus-CE), MF4.
2. "Statistically significant, biologically trivial." → MF2 (SESOI, headroom closed).
3. "Transfer to unseen tissues is claimed from five lineages." → MF3.
4. "The decomposition (medium / tasks / curation / extraction) is order-dependent and estimated where it was selected." → MF5, SF5.
5. "Only one parent model was really optimised." → MF6, SF8.
6. "Results may be specific to the Avana library or the Chronos pipeline." → SF2.
7. "Flux features win because of GPR aggregation, not simulation." → SF4.
8. "No biologically interpretable successes or failures." → SF3, SF7, §9.4.

I have no changes to the leakage machinery (§4.2, §11), the seed plan (§5), the metric toy tests (§10), the W2 control set, or the W3 a priori task enumeration: none are needed from my angle.
