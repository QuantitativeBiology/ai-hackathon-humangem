# A priori lists: independent re-approval (open question F-4)

Approver: round-3 verifier agent. Date: 2026-10-26.

**Attestation.** Author not in notes/label_exposure.md; no label values viewed (no gene × line values, and no per-line, per-gene or per-lineage label statistics).

**Disclosure.** While reading files the task assigned, I saw two kinds of all-line *aggregate* label statistics, printed in `notes/label_exposure.md` and `notes/design_review_r2/statistics.md`:
- the common-essential share of positives and ASNS prevalence by medium;
- null-AP levels.

This is the same aggregate exposure the PI declared. **No decision below uses them.** The ASNS medium control rests on medium composition: DMEM lacks asparagine and RPMI contains it.

**Basis.** Every decision uses only:
- biochemistry, medium composition and published literature (cited by first author and year; DOIs are to be resolved by `verify_citations.py`);
- Human-GEM 2.1.1 / 1.19 structure;
- the label-free generic KO ratios on Ham's in `results/baseline/*_ko.parquet` (cited as "H1/H2 ratio").

**Scope (v3).** ftINIT is dropped (user decision D2). No list exists *only* for ftINIT. Two configurations go with it and are not lists: the B7t viability-task scoring and Appendix B. The Hart 2015 lines remain usable only as a binary-metric comparison if the Human1-paper protocol needs them.

Legend: **A** = approve, **M** = modify (reason given), **R** = remove.

---

## 1. Gene panels (`data/controls/gene_panels.tsv`)

### 1.1 Pan-medium essentials (expected essential in every standard medium)

| Gene(s) | Decision | Reason |
|---|---|---|
| CAD, DHODH, UMPS | **M** | Essential unless uridine or cytidine is supplied. HPLM contains uridine, so the expectation is "all media except uridine-containing ones (HPLM)". Thymidine does not rescue the pyrimidine ring. |
| CTPS1 | **M → buffering pair** | CTPS2 is a functional paralog. Move to §1.3 as CTPS1/CTPS2, an isozyme-logic test case. H1/H2 ratio 1. |
| TYMS | **M → medium-dependent** | Thymidine rescues it via TK1 salvage. Ham's/F12 and HPLM contain thymidine; H2 ratio 1 on Ham's. Expected essential in thymidine-free media (RPMI, DMEM). |
| DHFR | **M → medium-dependent** | The classic HT rescue of antifolates: hypoxanthine + thymidine (+ glycine) bypass DHFR. Expected essential in HT-free media. |
| RRM1, RRM2 | **A** | No dNDP salvage source in standard media. RRM2B (p53R2) is a weak paralog; the Note below covers the GPR. |
| GUK1 | **A** | Single GMP kinase. |
| GMPS | **A** (caveat) | Hypoxanthine does not bypass GMPS; only guanine or guanosine does. |
| PCYT1A | **A** (caveat) | Rate-limiting step of PC synthesis. PCYT1B is a low-expressed paralog (an expression-logic test); the PEMT route is hepatic. |
| CHKA | **M → buffering pair** | CHKB is a paralog with overlapping activity, and choline kinase loss is partially buffered. Move to §1.3 as CHKA/CHKB. |
| CS | **A** (caveat) | Single citrate synthase. Reductive carboxylation can supply cytosolic citrate, so the expectation is weaker in OXPHOS-low lines. Tag "moderate". |
| MDH2 | **A** (caveat) | Needed for TCA flux and aspartate synthesis. Pyruvate (an electron acceptor) partly relaxes it, so tag it "pyruvate-modulated" (Birsoy 2015; Sullivan 2015). |
| SDHA, SDHB, SDHC, SDHD | **M** | Move to a new class, **OXPHOS/TCA expected-essential, pyruvate- or uridine-modulated**. SDH-deficient tumours exist, so "pan-medium" overstates it. H1/H2 ratios 0.76/0.88 (partial). |
| NAMPT | **M → medium-dependent / pair** | It is listed in both pan-medium and medium-dependent (contradiction). NAMPT is bypassed when NAPRT is expressed **and** nicotinic acid is present; most media supply nicotinamide only. Keep it only under §1.2 and the NAPRT→NAMPT pair. |
| COASY | **A** | Single bifunctional enzyme for the last CoA steps. |
| Mevalonate core: HMGCS1, HMGCR, MVK, PMVK, MVD, IDI1, FDPS, GGPS1 | **A** | Non-sterol isoprenoids (CoQ, dolichol, prenyl groups) cannot be supplied by serum LDL (Goldstein & Brown 1990). |

Note: the generic models already disagree on several of these.
- HMGCS1, IDI1, FDPS and GGPS1 have H1/H2 ratio 1.
- HMGCS1 vs HMGCS2 is a **compartment GPR error**: HMGCS2 is mitochondrial and ketogenic. Flag it for W5.

New class, **OXPHOS/TCA expected-essential (pyruvate/uridine-modulated)**: SDHA–D, plus the complex I and III core subunits. Tag it "moderate". The expectation is essential in standard media. Uridine and pyruvate supplementation relaxes it (King & Attardi 1989).

### 1.2 Medium-dependent

| Gene(s) | Condition | Decision | Reason |
|---|---|---|---|
| ASNS | essential without Asn (DMEM), not with Asn (RPMI, F12, IMDM) | **A** | From medium composition (Lagziel 2020). |
| Sterol branch: FDFT1, SQLE, LSS, CYP51A1, MSMO1, NSDHL, HSD17B7, SC5D, EBP, DHCR7, DHCR24 | gain essentiality without serum lipoprotein | **A** (caveat on DHCR24) | Desmosterol can partly substitute for cholesterol in cultured cells, so DHCR24 is tagged "weaker". |
| De novo purine: PPAT, GART, PFAS, PAICS, ATIC | relaxed by hypoxanthine (F12-based media, HPLM) via HPRT1 | **A** | |
| ADSL | **M** | Not relaxed by hypoxanthine | ADSL also catalyses IMP→AMP (S-AMP → AMP). Only adenine (APRT) or adenosine (ADK) rescues it. The H1/H2 ratio of 1 on Ham's, which has no adenine, is a model artefact to trace. |
| IMPDH2 | **M** | Not relaxed by hypoxanthine | IMP→XMP is still required, and only guanine or guanosine rescues it. Buffered by IMPDH1 (pair in §1.3). Remove it from the hypoxanthine set in the Appendix A addendum. |
| TYMS, DHFR | thymidine (TYMS); hypoxanthine + thymidine (DHFR) | **A (moved here)** | See §1.1. |
| CAD, DHODH, UMPS | uridine (HPLM) | **A (moved here as an exception)** | See §1.1. |
| GLUL | essential only without glutamine | **A** | Most media carry Gln or GlutaMAX. |
| PYCR1, ALDH18A1 | essential without proline (DMEM) | **A** | PYCR2/PYCR3 paralogs are noted. ALDH18A1 also feeds ornithine. |
| SLC7A11 | cystine-only media vs 2-mercaptoethanol (cysteine available) | **A** | Add the transsulfuration alternative (CTH, CBS from methionine) as a note. |
| NAMPT / NAPRT | nicotinic acid + NAPRT expression | **A** | Hara 2007. |
| GOT1, GOT2 | relaxed by pyruvate; Asp uptake only with SLC1A3 (the gate) | **A** | Birsoy 2015; Sullivan 2015; Garcia-Bermudez 2018. |
| PHGDH, PSAT1, PSPH | essential without Ser/Gly | **A** | Standard media contain serine, so the expectation is non-essential except in PHGDH-amplified lines. |
| ASS1, ASL | **add** | essential only without arginine (none of the standard media) | Supports the ASS1_silenced group (Delage 2010). |

### 1.3 Collateral-lethality and buffering pairs

| Pair | Decision | Reason |
|---|---|---|
| ENO1 loss → ENO2 | **A** | Muller 2012 |
| ME2 loss → ME3 | **A** | Dey 2017 |
| NAPRT loss → NAMPT | **A** | Cerna 2012; Chowdhry 2019 |
| CTPS1/CTPS2, IMPDH1/IMPDH2, CHKA/CHKB, RRM2/RRM2B | **add** | Paralog-buffered isozymes. They test expression-aware OR logic. **RRM2B is lethal in both generic models (ratio 0)**, which signals a GPR AND-error, since RRM1·(RRM2 or RRM2B) is the known logic. Add it to the W5 positive-control candidates. |
| MTAP loss → MAT2A/PRMT5 | **A** (expected FBA failure) | Kryukov 2016; Mavrakis 2016; Marjon 2016 |

### 1.4 Known GEM false-positive classes

| Class | Decision | Reason |
|---|---|---|
| H2 β-oxidation (ACADL, HADHA, ECHS1) | **A** (caveat ECHS1) | Model fact: ratio 0 in H2 and 1 in H1. Fatty-acid oxidation is dispensable in glucose/glutamine culture. ECHS1 loss can be truly toxic through valine catabolism (methacrylyl-CoA), which FBA does not represent; tag it "uncertain". |
| Xenobiotic CYP/UGT/GST | **A** | No biomass role. |
| Urea cycle (CPS1, OTC, ARG1) | **A** | Inactive in most cancer lines. CPS1 is the exception in STK11-loss lines (Kim 2017); see §2. |
| Bile acid and steroid-hormone synthesis | **A** | Tissue-restricted. |
| Biomass-artefact lethals: CYP4F2, BBOX1, TMLHE, CYP26A1, CYP26B1, CYP24A1, CYP1A2 | **A** | Confirmed against the model. All are ratio 0 in both parents except CYP4F2 (H2 only). Mechanism: the H2 cofactor pool (MAR10065) demands retinoid and vitamin-D catabolites and L-carnitine. |

### 1.5 Known GEM false-negative classes

| Class | Decision | Reason |
|---|---|---|
| V-ATPase (ATP6V\*) | **A** | Generic ratio 1. The essentiality is non-stoichiometric (lysosomal pH, mTORC1). |
| PGD | **A** | Toxic accumulation of 6-phosphogluconate (Lin 2015); not stoichiometric. |
| GPX4 / ferroptosis | **A** | Yang 2014 |
| MTAP → MAT2A/PRMT5 | **A** | |
| Mitochondrial translation and assembly factors | **A** | Outside the GEM gene scope or GPR. |

### 1.6 CRISPR-artefact strata (§5.4)

**Multi-copy loci: A, with additions.** These are loci whose guides can hit near-identical copies (high-identity paralogs or pseudogenes).

| Locus | Decision | Note |
|---|---|---|
| UGT1A (shared exons 2–5) | A | |
| UGT2B | A | |
| GSTM1–5 | A | |
| GSTT1/2/2B | A | |
| AMY1A–C | A | |
| NAT1/2 | A | |
| GBA1/GBAP1 | A | GBA is now named GBA1. |
| AMY2A/2B | **add** | |
| SULT1A3/1A4 | **add** | |
| CYP2D6/CYP2D7 | **add** | |
| the CYP4F cluster (CYP4F2/3/11/12) | **add** | Relevant because CYP4F2 is an artefact lethal. |
| AKR1C1–3 | **add** | |
| ACSM2A/2B | **add** | |
| GSTA1/2 | **add** | |
| CES1/CES1P1 | **add** | |

**Rule (M).** To keep the list reproducible, complete it label-free at M1: add every gene in U with an Ensembl 116 paralog or pseudogene at ≥ 95% nucleotide identity. The seeds above are always included.

**Housekeeping-adjacent: A.** Define it by gene family so it is reproducible:
- aminoacyl-tRNA synthetases (cytosolic \*ARS and mitochondrial \*ARS2);
- V-ATPase (ATP6V\*);
- Fe-S biogenesis (NFS1, ISCU, FXN, HSCB, GLRX5, CIAO\*, NUBP1/2);
- GPI anchor (PIG\*, DPM1–3);
- N-glycosylation (ALG\*, DPAGT1, DDOST, RPN1/2, STT3A/B);
- OXPHOS assembly factors (NDUFAF\*, COA\*, SCO1/2, SURF1, BCS1L).

**Moonlighting list: M (it was unspecified, so it could not be approved as written).** Proposed:
- ACO1 (IRP1; Haile 1992);
- GAPDH (Sirover 2011);
- ENO1 (MBP-1; Feo 2000);
- PKM (nuclear PKM2; Yang 2012);
- ALDOA;
- HK2 (mitochondrial anti-apoptotic; Pastorino 2002);
- PGK1;
- ACLY and ACSS2 (nuclear acetyl-CoA; Wellen 2009; Li 2017);
- PDHA1, PDHB, DLAT (nuclear PDC; Sutendra 2014);
- FH (DNA repair; Jiang 2015).

## 2. Genomic-background groups (W4)

| Group | Decision | Reason / change |
|---|---|---|
| MTAP_del (relCN < 0.3) | **A** | Expected FBA failure for MAT2A/PRMT5. Also pre-register "increased reliance on de novo purine synthesis" as a weak GEM-addressable expectation (adenine salvage from MTA is lost). |
| ASS1_silenced (log2(TPM+1) < 1) | **A** | Expected: arginine import (SLC7A1) and pyrimidine synthesis via aspartate (CAD; Rabinovich 2015). Add CAD to the expected genes. |
| NRF2_active (KEAP1 loss or NFE2L2 hotspot) | **M** | Expected genes: add **GLS** (glutaminolysis; Romero 2017; Sayin 2017) to SLC7A11 and NADPH (G6PD, PGD, ME1, IDH1). CUL3 loss as a sensitivity definition only. |
| STK11_loss | **M** | Expected genes: add **DTYMK** (Liu 2013) to CPS1/pyrimidine (Kim 2017). DTYMK is generic-lethal (ratio 0), so it is uninformative for contrasts but kept for completeness. |
| NAPRT_silenced | **A** | NAMPT |
| FH_SDH_loss | **A (drop)** | 2 lines under the r2-11 rule. In-silico only. |
| PHGDH_amp (relCN ≥ 2) | **A** | PHGDH, PSAT1, PSPH (Possemato 2011; Locasale 2011) |
| ME2_loss | **A** | ME3 (Dey 2017) |
| ENO1_del | **A (in-silico)** | ENO2 (Muller 2012) |
| IDH_mut | **A (drop)** | Too few lines; the 2-HG reaction would have to be added. |
| KRAS_hotspot (negative control) | **M** | KRAS is not metabolically neutral: KRAS-driven PDAC uses GOT1-dependent glutamine metabolism (Son 2013) and macropinocytosis. Keep it as a secondary "non-metabolic-driver" control. Add the **primary** negative control: size-matched random within-lineage member sets, SEEDS[0:5], with the same contrast. |

## 3. Isozyme exceptions (W5 "never edited automatically")

All of LDHA/B, PFKL/M/P, ENO1/2/3, ALDOA/B/C, CKB/CKM, COX4I1/2, COX6A1/2, COX7A1/2, ATP5MC1/2/3, ATP6V1B1/2 and ATP6V0A1–4 are **A**. Each family forms mixed or alternative oligomers, or alternative subunit isoforms of one complex, so OR is biochemically correct. The CORUM co-membership and co-abundance rule could otherwise turn them into AND.
- ATP5MC1/2/3 encode an identical mature c-subunit.
- The ATP6V0A isoforms are compartment-targeted; that is a note only.

**M, add** the following, where OR logic is established and the co-abundance rule could misfire:
- PGAM1/PGAM2;
- PRPS1/PRPS2;
- CTPS1/CTPS2;
- IMPDH1/IMPDH2;
- RRM2/RRM2B: protected only against OR→AND. The current RRM2B lethality is an AND-error, and a manual fix is allowed.

Tissue-restricted isoforms (PKLR, PDHA2, MAT1A, ENO3, ALDOB, CKM) remain eligible for the "pan-cancer absent" branch removal. That edit is legitimate curation, not an isozyme-logic error.

## 4. Metabolic-drug classes and compounds (W6)

**Rules (A):**
- the primary mechanism targets a Human-GEM enzyme or transporter;
- kinase inhibitors are excluded;
- each compound carries a prodrug flag and its activating gene(s);
- a compound needs ≥ 100 lines with features in PRISM, GDSC2 or CTRP (checked by script).

Compounds not found in any screen are dropped by that script, not by hand.

| Class | Compounds (target; prodrug → activating genes) | Decision |
|---|---|---|
| Antifolates | methotrexate (DHFR; uptake SLC19A1, polyglutamation FPGS); pemetrexed (TYMS, DHFR, GART; SLC19A1, FPGS); pralatrexate (DHFR; SLC19A1, FPGS) | **A** (compounds specified) |
| TYMS inhibitors | 5-FU (TYMS; prodrug → UMPS, or UPP1/UCK2, TYMP/TK1; catabolism DPYD); floxuridine (TYMS; prodrug → TK1); raltitrexed (TYMS; SLC19A1, FPGS); capecitabine (TYMS; prodrug → CES1/2, CDA, TYMP) | **A** (capecitabine added) |
| RNR inhibitors | hydroxyurea (RRM2); gemcitabine (RRM1; prodrug → DCK); triapine (RRM2) | **A** |
| DHODH inhibitors | brequinar; teriflunomide; leflunomide (prodrug, non-enzymatic conversion: no GEM activator); BAY-2402234 | **A** (BAY-2402234 added) |
| IMPDH inhibitors | mycophenolic acid / mycophenolate mofetil (IMPDH1/2; mofetil → CES1); ribavirin (IMPDH; prodrug → ADK); mizoribine (IMPDH; prodrug → ADK) | **A** (mizoribine added) |
| Purine analogues | **M (specified)**: 6-mercaptopurine and 6-thioguanine (prodrug → HPRT1; targets PPAT / de novo purine); fludarabine, cladribine, clofarabine (prodrug → DCK; target RRM1 and DNA synthesis) | M |
| Statins and FDPS inhibitors | simvastatin and lovastatin (HMGCR; lactone prodrugs, esterase); atorvastatin, fluvastatin, pitavastatin, rosuvastatin (HMGCR); zoledronate, alendronate, ibandronate, risedronate (FDPS) | **A** |
| FASN, ACLY, SCD inhibitors | TVB-3166 and TVB-2640/denifanstat, GSK2194069, orlistat (FASN thioesterase; note its lipase off-targets); BMS-303141, SB-204990, NDI-091143 (ACLY); bempedoic acid (ACLY; prodrug → SLC27A2/ACSVL1, liver-restricted: flagged); A939572, CAY10566, MF-438 (SCD) | **A** |
| GLS inhibitors | CB-839/telaglenastat, BPTES (GLS) | **A** |
| Mutant-IDH inhibitors | ivosidenib, AGI-5198 (IDH1); enasidenib (IDH2); vorasidenib (IDH1/2) | **A (caveat)**: the 2-HG reaction is absent from the GEM, so the response has no GEM handle except IDH expression. Report it as a negative-expectation class. |
| ODC1 inhibitor | eflornithine/DFMO | **A** |
| NAMPT inhibitors | FK866, GNE-617, GMX1778/CHS-828 | **A** (GMX1778 added; KPT-9274 excluded: dual PAK4) |
| OXPHOS inhibitors | **M (broadened from "complex I")**: IACS-010759, phenformin, metformin (weak, mechanism debated: flagged), rotenone (complex I); antimycin A, atovaquone (complex III, MT-CYB); oligomycin (ATP synthase) | M |
| **Add:** MCT1 | AZD3965 (SLC16A1) | **add** (transporter in Human-GEM; Polański 2014) |
| **Add:** GLUT1 | BAY-876 (SLC2A1) | **add** |
| **Add:** LDHA | GSK2837808A, GNE-140 (LDHA/LDHB) | **add** |
| **Add:** CPT1 | etomoxir (CPT1A; off-target at high dose: flagged) | **add** |
| **Add:** SHMT | SHIN1/SHIN2 (SHMT1/2) | **add** (Ducker 2017) |
| **Add:** ACC | ND-646, firsocostat (ACACA/ACACB) | **add** |
| **Add:** system xc⁻ | erastin, sulfasalazine (SLC7A11; erastin also hits VDAC: flagged) | **add** |

Positive-control pairs (W6): DHODH–brequinar, DHFR–methotrexate, HMGCR–statins, TYMS–5-FU and IMPDH–mycophenolate are **A**.

## 5. Candidate task sets (W3)

| Set / mode | Decision | Reason |
|---|---|---|
| Essential task list: 57 tasks in 11 groups (`metabolicTasks_Essential.txt`) | **A** | Fixed external list. It is also the source of the GR medium. |
| CellFie consensus list: 195 tasks | **A** | |
| T0–T3, T4_k (k = 1…11), T5 (CellFie-inspired, 25th/75th percentiles) | **A** | |
| D-BS, D-ER (1% of max) | **A** | |
| O-cap, O-NGAM k ∈ {0, 0.1, 0.3}, O-MOMA | **A** | |
| B-loo_j (j = 1…13) | **M** | "13" does not match H2: MAR13082 has 9 direct reactants (ATP, DNA, H2O, RNA, glycogen, and the cofactor, protein, lipid and metabolite pools). New definition: the direct reactants except ATP and H2O, plus the cofactor-pool groups linked to the floor or to artefacts. Those groups are retinoid and vitamin-D catabolites, tocopherol catabolites, carnitine, heme, ubiquinol, lipoate and cobamide. The set is enumerated per parent by script at M3 and frozen. |
| B-add list: CoQ10, heme and heme A, dolichol-P, GPI-anchor precursor, lipoate, CoA, Fe-S cluster, SAM sink, spermidine | **M** | In H2, **ubiquinol (CoQ), heme and lipoic acid are already in the cofactor pool, and CoA is in the metabolite pool**. Adding them changes only coefficients, which W3 rejects. New rule: skip and log any component already in the parent's biomass. Keep heme A, dolichol-P, GPI-anchor precursor, Fe-S cluster, SAM sink and spermidine. **Add an NAD⁺ turnover sink** (consumption by PARPs and sirtuins releases nicotinamide). NAD⁺ is in the pool only statically, which makes NAMPT non-essential in FBA (ratio 1) against the well-established NAD turnover biology. **Add a glutathione turnover sink** (GSH export and conjugation; GCLC/GSS). |
| Negative controls: random pseudo-tasks; matched call rate | **A** | |
| Positive controls: ETF capability FPs under T3; mevalonate core under B-add | **A** (mevalonate) | HMGCR, MVK, PMVK and MVD are already lethal on M0 because ubiquinol is in biomass. The B-add control is informative only for genes whose M0 ratio is 1 (HMGCS1, IDI1, FDPS, GGPS1), where the cause is isozyme or GPR, not demand. Restate it as "**no loss of mevalonate-core essentiality** when serum cholesterol is added under B-add". |

## 6. Summary of changes
- **Moved:**
  - TYMS and DHFR → medium-dependent;
  - CTPS1 and CHKA → buffering pairs;
  - NAMPT out of pan-medium;
  - SDHx → a new OXPHOS/TCA "moderate" class.
- **Corrected:** ADSL and IMPDH2 are not hypoxanthine-rescued.
- **Added:**
  - ASS1/ASL (Arg);
  - buffering pairs, including the RRM2B GPR error;
  - multi-copy loci plus an identity rule;
  - housekeeping families defined;
  - the moonlighting list specified;
  - GLS for NRF2, DTYMK for STK11, CAD for ASS1;
  - a random-group negative control for backgrounds;
  - five isozyme exceptions;
  - drugs specified per class, plus seven new classes;
  - NAD⁺ and GSH turnover sinks.
- **Redefined:** B-loo (the set) and B-add (skip components already in biomass). The mevalonate B-add control is restated.
- **Removed:** none. ftINIT-only items (B7t scoring, Appendix B) go with D2.
