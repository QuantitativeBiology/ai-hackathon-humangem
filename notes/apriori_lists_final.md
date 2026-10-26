# A priori lists: final reconciled version

Reconciles `notes/apriori_lists_blind.md` (written blind, commit 57670ff) with `notes/apriori_lists_approved.md` (written by an agent with declared *aggregate-only* label exposure). Every item carried over from the approved file was re-judged on biochemistry, medium composition, published literature and Human-GEM structure; nothing was accepted because it was already approved.

**What I read of the approved file.** Its own attestation and disclosure (lines 5-11) *name* two kinds of all-line aggregate statistics its author saw in files I have not opened. The values are not in that file and I did not look for them. No decision below uses them; the ASNS medium item rests on the vendor formulation sheets (DMEM high-glucose contains no asparagine; RPMI contains 0.38 mM).

**Evidence tags.** `[MS]` Human-GEM structure verified by me (reaction ID / GPR / stoichiometry) · `[PRED]` generic-model KO ratio from `results/baseline` (predictions, no labels) · `[OMICS]` computed by me from non-label DepMap omics · `[LIT]` literature. All PMIDs/DOIs below were resolved against PubMed in this session (see the citation ledger, §B of this file).

According to PubMed, the citations used here resolve as listed in **§C. Citation ledger**, with DOI links given there.

---

# 1. Final gene panels

## 1.1 Pan-medium essentials — expected essential in every standard medium

Tagged `strong` unless marked. 60 genes + 20 ARSs.

**A. Aminoacyl-tRNA charging (20, `strong`, non-metabolic mechanism)** — AARS1, CARS1, DARS1, EPRS1, FARSA, FARSB, GARS1, HARS1, IARS1, KARS1, LARS1, MARS1, NARS1, QARS1, RARS1, SARS1, TARS1, VARS1, WARS1, YARS1. Substrates of `MAR10062` → `protein_pool_biomass` [MS]; lethal in both models [PRED]. Scored separately (§1.6g).

**B. Pyrimidine core (5, `strong`)** — CAD, DTYMK, RRM1, RRM2, GUK1. *Exception:* CAD, DHODH and UMPS move to §1.2 (uridine/HPLM). GUK1 is the only GMP kinase in the model (`MAR04020`, `MAR04135`, sole gene) [MS], ratio 0.999 [PRED].

**C. Purine core not rescued by hypoxanthine (4, `strong`)** — ADSL, ADSS2, IMPDH2, GMPS. Hypoxanthine+HPRT1 yields IMP only; AMP still needs ADSS→ADSL and GMP still needs IMPDH→GMPS. PPAT, GART, PFAS, PAICS and ATIC move to §1.2.

**D. Folate / 1C (3, `strong`)** — SHMT2, MTHFD2, plus MTHFD1 (`moderate`). DHFR moves to §1.2 (HT rescue). [LIT Ducker & Rabinowitz 2017, PMID 27641100]

**E. Cofactor activation and transport (9, `strong`)** — COASY, PPCS, NADK, TPK1, RFK, FLAD1, PDXK, SLC25A26, SLC5A6. NAMPT and PNPO move to §1.2.

**F. Fe-S assembly and lipoylation (6, `strong`; Human2 only)** — NFS1, ISCU, LYRM4, NFU1, GLRX5, HSPA9. Absent from Human1 [MS]; ratio 1 in Human2 [PRED]. (FXN is present in both at ratio 1; CIAO1, CIAO2A, NUBP1, NUBP2 are in neither model [MS].)

**G. Lipid / membrane core (8, `strong`)** — CDS2, CDIPT, PCYT1A, SPTLC1, SPTLC2, KDSR, TECR, HSD17B12.

**H. Mevalonate trunk / isoprenoid branch (8, `strong`)** — HMGCS1, HMGCR, MVK, PMVK, MVD, IDI1, FDPS, GGPS1. Serum LDL rescues the sterol branch but not FPP/GGPP for prenylation and dolichol [LIT Goldstein & Brown 1990, PMID 1967820]. **Curation flag:** HMGCS1/HMGCS2 are a compartment error — HMGCS2 is mitochondrial and ketogenic, so any `OR` between them is wrong; HMGCS1, IDI1, FDPS, GGPS1 are all ratio 1 [PRED], i.e. currently false negatives.

**I. Hexosamine / UDP-sugar core (5, `strong`)** — GFPT1, GNPNAT1, PGM3, UAP1, DPM1.

**J. Glycolysis core (7, `strong` biologically, `model-masked`)** — GAPDH, PGK1, TPI1, ALDOA, ENO1, PKM, GPI. All ratio 1 [PRED] because of isozyme `OR` logic (§1.5c); kept in the panel and tagged so the miss is attributable.

**K. OXPHOS / TCA, pyruvate- or uridine-modulated (new class, `moderate`)** — SDHA, SDHB, SDHC, SDHD; CS; MDH2; the core complex I and complex III subunits; ATP5F1A/B/C and the ATP-synthase subunits. Respiration's essential proliferative role is aspartate synthesis, and pyruvate or uridine supplementation relaxes it [LIT Birsoy 2015, PMID 26232224; Sullivan 2015, PMID 26232225; King & Attardi 1989, PMID 2814477]. SDH-deficient tumours exist, so "pan-medium" would overstate it. Model ratios 0.66-0.88 [PRED].

**Deliberately excluded:** the NDUF* subunits as a *strong* class (glucose-grown cells tolerate complex I loss), GLS, PHGDH, FASN, ACACA, SQLE, GCLC, SLC7A11, GPX4, ASNS, ASS1, NAMPT, TYMS, DHFR.

## 1.2 Medium-dependent essentials (rescuing component)

Medium distribution [OMICS, `Model.csv`, noting the DATA.md column swap]: RPMI 931, DMEM 298, EMEM/MEM 163, DMEM:F12 104, IMDM 73, F12 57, αMEM 31, McCoy's 31, L-15 23, Waymouth 5, organoid/stem/proprietary 347, other 34, not recorded 457; 1,730 formulations name serum, 824 do not (331 of those are the serum-free organoid/stem class).

| Gene(s) | Rescued by | Where | Evidence |
|---|---|---|---|
| CAD, DHODH, UMPS | uridine (and uric acid inhibits UMPS) | none of the vendor sheets; serum; **HPLM and Plasmax contain uridine** (both in `data/media/`) | [LIT Cantor 2017, PMID 28388410] |
| TYMS, TK1 | thymidine | F12 and DMEM/F12 only (0.365 mg/L) | exchange `MAR09423` [MS]; Human1 `MAR04644` reversibility lets thymidine also rescue UMPS [MS] |
| DHFR | hypoxanthine **+** thymidine (HT) | F12, DMEM/F12 | [LIT Urlaub & Chasin 1980, PMID 6933469] |
| PPAT, GART, PFAS, PAICS, ATIC | hypoxanthine → HPRT1 | F12, DMEM/F12; serum | **ADSL, IMPDH2 and GMPS are not in this set** (see §1.1C) |
| ASNS | asparagine | RPMI 0.38 mM, IMDM 0.17, DMEM/F12 0.05; **DMEM high-glucose none** | `MAR09062` [MS]; [LIT Krall 2016, PMID 27126896; Lagziel 2020, PMID 33046912] |
| GLUL | glutamine | 2-4 mM in most media; "no-glutamine"/GlutaMAX formulations differ | glutamine exchange `MAR09063` [MS]; ratio 1 [PRED] |
| ALDH18A1 (P5CS), PYCR1/2/3 | proline | RPMI 0.17 mM, IMDM 0.35, DMEM/F12 0.15; **DMEM high-glucose none** | all ratio 1 [PRED] |
| GOT1, GOT2 (+ SLC1A3 as the aspartate-uptake gate) | pyruvate; exogenous aspartate only where SLC1A3 is expressed | pyruvate in RPMI/IMDM/DMEM-F12, absent from DMEM w/o pyruvate | [LIT Birsoy 2015, PMID 26232224; Sullivan 2015, PMID 26232225; Garcia-Bermudez 2018, PMID 29941933] |
| PC, PCK2, ETC subunits | pyruvate (and uridine) | as above | `MAR09133` [MS]; [LIT PMID 2814477] |
| sterol branch FDFT1, SQLE, LSS, CYP51A1, TM7SF2, MSMO1, NSDHL, HSD17B7, SC5D, EBP, DHCR7, **DHCR24 (`weaker`: desmosterol partly substitutes for cholesterol)** | serum LDL cholesterol | serum only → the 331 serum-free lines are not rescued | `MAR09285` closed under Ham's; 11 lethal in Human2 [MS/PRED]; [LIT Garcia-Bermudez 2019, PMID 30760928] |
| FASN, ACACA, SCD | serum palmitate/oleate | serum; DMEM/F12 adds 0.042 mg/L linoleate | `MAR00611`, `MAR00650` [MS] |
| GCLC, GCLM, GSS | reduced glutathione; cystine | **RPMI uniquely contains GSH (1 mg/L)** | `MAR09351` [MS]; lethal in both models [PRED] |
| SLC7A11 | cystine level; 2-mercaptoethanol or free cysteine; transsulfuration (CBS, CTH) from methionine | cystine in all sheets; free cysteine in DMEM/F12 | `MAR09363`, `MAR09065` [MS] |
| NAMPT (modifier NAPRT) | nicotinic acid **plus** NAPRT expression | nicotinamide in all sheets; **nicotinic acid in none** | `MAR09378` vs `MAR09142` [MS]; [LIT Hara 2007, PMID 17604275] |
| PNPO | pyridoxal instead of pyridoxine | **IMDM supplies pyridoxal**; RPMI/DMEM/F12 supply pyridoxine | `MAR09400` vs `MAR09144` [MS] — sharpest single medium contrast; PDXK stays essential either way |
| BBOX1, TMLHE, ALDH9A1, SHMT1(H2)/SHMT2(H1) | serum carnitine | serum only | carnitine in `MAR10065`; `MAR09292` closed [MS] |
| MMAA, MMAB, MMUT | vitamin B12 | RPMI, DMEM/F12, IMDM; **DMEM high-glucose none** | `MAR08614-16` [MS] |
| ACACA, PC, PCCA/PCCB, SLC5A6 | biotin | RPMI 0.2, IMDM 0.013, DMEM/F12 0.0035 mg/L; **DMEM high-glucose none** | [MS] |
| LIAS, LIPT1, LIPT2 | lipoic acid | **DMEM/F12 only** | `MAR09167` [MS] |
| ODC1, SRM, SMS, AMD1 | putrescine | **DMEM/F12 only** | `MAR11426` [MS]; all ratio 1 [PRED] |
| ALAD, HMBS, UROS, UROD, CPOX | serum hemin | serum | `MAR09107` closed; all lethal [PRED] |
| ISYNA1 | myo-inositol | all sheets | `MAR09361` [MS] |
| SEPHS2, selenoproteins (GPX4) | selenium | **IMDM only**; otherwise serum | no selenite exchange in either model — representability gap [MS] |
| PHGDH, PSAT1, PSPH | serine + glycine | present in every sheet | expectation: non-essential except in PHGDH-amplified lines; **PSAT1 is nonetheless model-lethal for a pool artefact** (§1.4b) |
| ASS1, ASL | arginine (citrulline bypass) | arginine in all sheets | [LIT Delage 2010, PMID 20104527; Rabinovich 2015, PMID 26560030] |

## 1.3 Collateral lethality and paralog buffering

**Collateral lethality.** (1) MTAP loss → PRMT5 / MAT2A / RIOK1 [LIT Mavrakis 2016, PMID 26912361; Kryukov 2016, PMID 26912360; Marjon 2016, PMID 27068473]; PRMT5 is in **neither** model, MAT2A is in both → only the SAM-flux half is representable; additional weak GEM-addressable expectation: loss of adenine salvage from MTA raises reliance on de novo purine synthesis. (2) ENO1 deletion → ENO2 [LIT Muller 2012, PMID 22895339]. (3) ME2 deletion (SMAD4 locus) → ME3 [LIT Dey 2017, PMID 28099419]. (4) FH loss → HMOX1 [LIT Frezza 2011, PMID 21849978]. (5) ASS1 silencing → arginine auxotrophy and aspartate diversion into pyrimidines via CAD [LIT PMID 26560030]. (6) NAPRT loss or amplification → NAMPT dependence or resistance [LIT Shames 2013, PMID 24097869; Chowdhry 2019, PMID 31019297].

**Paralog-buffering pairs** (the model's `OR` is right here; expression-based pruning would break it) [LIT De Kegel & Ryan 2019, PMID 31652272; Thompson 2021, PMID 33637726]: ENO1/ENO2/ENO3, ALDOA/ALDOC/ALDOB, LDHA/LDHB, HK1/HK2/HK3/HKDC1, PFKL/PFKM/PFKP, PKM/PKLR, PRPS1/PRPS2, IMPDH1/IMPDH2, **CTPS1/CTPS2**, **CHKA/CHKB**, **RRM2/RRM2B**, NME1/NME2, ADSS1/ADSS2, AK1-AK4, GLS/GLS2, GLUD1/GLUD2, ME1/ME2/ME3, MTHFD1/MTHFD1L/MTHFD2/MTHFD2L, SHMT1/SHMT2, ACSS1/ACSS2/ACSS3, PTDSS1/PTDSS2, CDS1/CDS2, GYG1/GYG2, FADS1/FADS2, UGT1A*, SULT1A1/SULT1A2, AKR1C1/AKR1C2, GOT1/GOT2 (compartment-separated, partially redundant).

**Model GPR error to fix (verified here):** every RNR reaction (`MAR02352`, `MAR02357`, `MAR04611`, `MAR04612`, `MAR04614`, 14 reactions) carries `(RRM2B and TXN2 and RRM1 and RRM2) or (RRM2B and TXN and RRM1 and RRM2)` [MS], which makes RRM2B lethal in both models (ratio 0) [PRED]. The established logic is `RRM1 and (RRM2 or RRM2B)`. This is both a false positive and a W5 positive-control candidate.

## 1.4 GEM false-positive classes

(a) **Cofactor/vitamin pool dead-ends** (`MAR10065`, coefficient 0.0012, 36 forced species): CYP26A1, CYP26B1, CYP24A1, CYP1A2, CYP4F2 (H2), KL, TMLHE, BBOX1, ALDH9A1, SHMT1(H2)/SHMT2(H1), MMAA, GCH1, SPR, PTS(H1), SLC5A6, COQ2/3/6/7. The tocopherol chain `MAR06433`/`MAR06435`/`MAR06436` (AND-chains of β-oxidation enzymes) is why ACADL, ACADSB, ECHS1, HADHA, HADHB, ACAA2, AMACR are lethal in Human2 and not Human1 [MS/PRED]. **Caveat (`uncertain`) on ECHS1:** ECHS1 loss can be genuinely toxic through valine catabolism (methacrylyl-CoA), a toxicity FBA cannot represent, so a "false positive" call for ECHS1 may be right for the wrong reason.
(b) **Small-metabolite pool dead-ends** (`MAR10064`): PSAT1 is the sole producer of 3-phosphoserine (`MAR03841`), demanded at 0.0025, so PSAT1 is lethal whatever the medium serine; also UGDH, PGM3, GNPNAT1, CAD [MS].
(c) **Glycogen in biomass** (0.4062): UGP2, GYG1/GYG2, GYS1 [MS].
(d) **Pseudo-reaction GPRs** (true essentials, non-metabolic mechanism; 91 of Human2's 154 lethal genes): `MAR13086` 31 nucleoporins; `MAR07161` all POLR1/2/3*, POLRMT, MTPAP; `MAR06978` SPCS1-3 + SEC11A/C; `MAR10062` the 20 ARSs [MS]. Report AUPRC with and without.
(e) **Medium misspecification**: sterol branch, heme pathway, FASN/ACACA, GCLC/GCLM/GSS, carnitine synthesis, COQ* [MS].
(f) **Single-gene transport annotation**: Human1 NPC1L1 (`MAR00018`), SCARB1 (`MAR00019`), ABCB6; both models SLC25A26, CDIPT, SLC5A6 [MS/PRED].
(g) **Over-specified AND complexes**: `SEC11A and SEC11C`; `MAR06435`; the 31-gene nucleoporin AND; the RNR AND above [MS].
(h) **Unbounded-uptake artefacts**: TYMP, SUCLG1, DERA, PGM2, SLC25A3, SLC16A1 and the OXPHOS partials, where thymidine acts as a bulk carbon source [PRED].
(i) **Xenobiotic metabolism** (CYP, UGT, GST, SULT, NAT families): no biomass role; model lethality is always an artefact [MS].
(j) **Urea cycle** (CPS1, OTC, ARG1, ASS1, ASL in the hepatic sense): inactive in most cancer lines — with the documented exception of CPS1 in STK11/LKB1-loss lung lines [LIT Kim 2017, PMID 28538732]. All three are ratio 1 [PRED], so this class is currently a *latent* false-positive class that medium or task changes could switch on.
(k) **Bile-acid and steroid-hormone synthesis**: tissue-restricted (AMACR already appears via the bile-acid subsystem) [MS].

## 1.5 GEM false-negative classes

(a) **Structural non-predictions**: 1,801/2,887 (H1) and 1,750/2,870 (H2) genes close no reaction → ratio 1 by construction [PRED].
(b) **Genes missing from the model**: PRMT5, OGDHL, CIAO1, CIAO2A, NUBP1, NUBP2, GBA1, CYP2D7 (both models); the Fe-S module only in Human2; 59 H1-only and 42 H2-only genes [MS].
(c) **Isozyme `OR` masking the only expressed paralog**: GAPDH, PGK1, TPI1, ALDOA, ENO1, PKM, GPI, LDHA, HK1/HK2, PFKP — all ratio 1 [PRED]. Largest fixable class; §3 bounds the fix.
(d) **Compartment-blind redundancy**: GOT1/GOT2, MDH1/MDH2, IDH1/IDH2, SHMT1/SHMT2, ACSS1/ACSS2 [MS].
(e) **Demands absent from biomass**: polyamines (ODC1, SRM, SMS, AMD1 — all ratio 1), lipid-peroxide defence (GPX4, SLC7A11, ACSL4) [LIT Yang 2014, PMID 24439385], selenoprotein synthesis (SEPHS2), nucleotide-pool turnover, mitochondrial translation and OXPHOS assembly factors (outside GEM scope), NAD+ and GSH turnover.
(f) **Non-stoichiometric toxicity**: PGD loss is deleterious through 6-phosphogluconate accumulation and LKB1-AMPK signalling, not through flux [LIT Lin 2015, PMID 26479318]; ratio 1 [PRED].
(g) **Non-stoichiometric, compartment-physiological essentiality**: the V-ATPase (ATP6V1A, ATP6V1B2, ATP6V0A1-4 …), essential through lysosomal pH and mTORC1, all ratio 1 [PRED].
(h) **Alternative de novo routes culture cells do not use**: NAMPT bypassed by tryptophan→quinolinate→NaMN; PPAT/GART/ATIC bypassed by hypoxanthine (a Ham's input) [MS].
(i) **No non-growth ATP maintenance**: `MAR03964` has bounds [0, 1000] [MS], so OXPHOS and glycolysis genes come out at 0.57-1.00 instead of lethal. Fixed by the NGAM variant in §5.

## 1.6 CRISPR-artefact strata

(a) **mtDNA-encoded (13 per model)**: MT-ND1-6, MT-ND4L, MT-CYB, MT-CO1-3, MT-ATP6, MT-ATP8. Not editable by Cas9 [LIT Gammage 2018, PMID 29179920]; 6 of them are model-partial-lethal [PRED]. Exclude from scoring, report separately.
(b) **High-identity paralogs and pseudogenes.** Rule (reproducible, label-free): include every model gene with an Ensembl 116 paralog **or pseudogene** at ≥95% nucleotide identity; my protein-coding-only pass gives 67 model genes at ≥95% and 111 at ≥90% [OMICS]. Seed loci always included: UGT1A cluster (shared exons 2-5), UGT2B, UGT2B17, GSTM1-GSTM5, GSTT1/GSTT2/GSTT2B, GSTA1/GSTA2, AMY1A-C, AMY2A/AMY2B, NAT1/NAT2, SULT1A1/1A2, SULT1A3/1A4, CYP2D6/CYP2D7 (CYP2D7 not in the models), the CYP4F cluster (CYP4F2/3/11/12 — relevant because CYP4F2 is an artefact lethal), AKR1C1-AKR1C3, ACSM2A/ACSM2B, CES1/CES1P1, GBA1/GBAP1 (GBA1 not in the models), CALM1/2/3, UBE2D2/UBE2D3, PPIAL4*, USP17L*, NUDT4/NUDT4B, NUDT10/NUDT11, ALG10/ALG10B, CKMT1A/CKMT1B, UQCRH/UQCRHL, DDX19A/DDX19B, PGAM1/PGAM4, GK/GK3/GK3P, RPE/RPEL1, PRPS1/PRPS2, GLUD1/GLUD2.
(c) **Amplicon / copy-number effect** [LIT Meyers 2017, PMID 29083409; Dempster 2021, PMID 34930405]. Model genes at relative CN ≥ 2 in ≥8% of 1,763 WGS lines [OMICS]: PDPR 33.4, CLPS 22.9, SQLE 11.3, NDUFB9 10.8, RNF139 10.3, SLC39A4 9.6, GPT 9.5, NTAQ1 9.5, ADCY8 9.3, ATP5F1E 9.1, SLC6A19 9.0, DGAT1 9.0, PCK1 8.8, PYCR3 8.7, GPAA1 8.7, EXT1 8.6, GFUS 8.6, HAS2 8.4, RAE1 8.4, CYP24A1 8.3, SDHA 8.3 (%). Three of these (SQLE, RAE1, CYP24A1) are model-lethal — a correlated artefact in the model's own positives.
(d) **Recurrently deleted / unscreenable**: USP9Y 72.4%, ACOT1 40.2%, UGT2B17 32.6%, MTAP 19.1%, HACD4 4.7%, STS 3.4%, GSTM1 1.1% of lines at CN < 0.3 [OMICS].
(e) **Gene-ID hazards** [PRED]: KLHL23, USP17L2 (no CRISPR column); GK3P and SLC37A4 as `ENSG00000281500` (no Entrez); 11 multi-Entrez genes; 7-10 renamed symbols. Match on Entrez only; exclude these.
(f) **Moonlighting / structural roles** (flux-independent essentiality): ACO1/IRP1 [LIT Haile 1992, PMID 1502165]; GAPDH [LIT Sirover 2011, PMID 21640161]; ENO1/MBP-1 [LIT Feo 2000, PMID 10802057]; PKM (nuclear PKM2) [LIT Yang 2011, PMID 22056988]; HK2 (mitochondrial, anti-apoptotic) [LIT Pastorino 2001, PMID 11751859]; ALDOA; PGK1; ACLY [LIT Wellen 2009, PMID 19461003] and ACSS2 [LIT Li 2017, PMID 28552616]; PDHA1/PDHB/DLAT (nuclear PDC) [LIT Sutendra 2014, PMID 24995980]; FH (DNA repair) [LIT Jiang 2015, PMID 26237645]; GPI (neuroleukin).
(g) **Housekeeping-adjacent, defined by family so it is reproducible**: cytosolic `*ARS` and mitochondrial `*ARS2`; V-ATPase `ATP6V*`; Fe-S biogenesis (NFS1, ISCU, FXN, HSCB, GLRX5, and CIAO*/NUBP1/2 where present — these four are absent from both models); GPI anchor (`PIG*`, DPM1-3); N-glycosylation (`ALG*`, DPAGT1, DDOST, RPN1/2, STT3A/B); OXPHOS assembly factors (`NDUFAF*`, `COA*`, SCO1/2, SURF1, BCS1L); plus the §1.4d pseudo-reaction genes (nucleoporins, RNA polymerases, signal peptidase, ARSs).

---

# 2. Final genomic-background groups

Denominators [OMICS, default entry per model]: CN 1,763 lines; expression 2,130; mutation matrix 2,395; all three 1,724. Thresholds fixed a priori: deep deletion CN < 0.3; amplification CN ≥ 2.0; expression-null log2(TPM+1) < 1.0.

| Group | Definition | n | Expected genes (GEM-representable) |
|---|---|---|---|
| MTAP_loss | CN < 0.3 **∪** TPM-null | **400** (336 CN, 341 expr) | MAT2A, de novo purine (PPAT, GART, PAICS, ATIC); PRMT5 not in model → expected FBA failure |
| **NEG-1** CDKN2A_del_MTAP_intact | CDKN2A < 0.3, MTAP CN ≥ 0.5, MTAP TPM ≥ 2 | **157** | specificity control for MTAP_loss |
| ASS1_silenced | ASS1 TPM-null | **254** | SLC7A1 (arginine import), **CAD** (aspartate → pyrimidines) |
| ASNS_low | ASNS TPM < 4 | **81** | ASNS; interacts with medium asparagine |
| NRF2_active | KEAP1 damaging ∪ KEAP1 CN<0.3 ∪ NFE2L2 hotspot | **68** | **GLS**, SLC7A11, and NADPH supply G6PD, PGD, ME1, IDH1 [LIT Romero 2017, PMID 28967920; Sayin 2017, PMID 28967864] |
| STK11_loss | damaging ∪ CN < 0.3 | **82** (Lung+KRAS 21; CPS1-expressed 35) | CPS1 [PMID 28538732]; **DTYMK** [LIT Liu 2013, PMID 23715154] — DTYMK is generic-lethal (ratio 0) so it is uninformative for contrasts, kept for completeness |
| NAPRT_silenced | TPM-null | **59** | NAMPT |
| PHGDH_amp | CN ≥ 2 (or CN ≥ 2 ∪ TPM ≥ 9: 115) | **42** | PHGDH, PSAT1, PSPH [LIT Possemato 2011, PMID 21760589; Locasale 2011, PMID 21804546] |
| ME2_loss | CN < 0.3 (27) ∪ TPM-null (26) | **~35** | ME3 [PMID 28099419] |
| SQLE_null | TPM-null | **13** | sterol branch, LDLR-dependent cholesterol uptake [PMID 30760928] |
| HK1low_HK2high | HK1 TPM<2 ∧ HK2 TPM≥4 | **38** | HK2 (isozyme switch at `MAR04319`) |
| LDHB_null | TPM-null | **58** | LDHA (`MAR04281` `LDHB or LDHA`) |
| SLC7A11_null | TPM-null | **90** | cystine import, transsulfuration |
| CDS1_null | TPM-null | **808** | CDS2 |
| VHL_loss | damaging ∪ CN<0.3 (38); Kidney-only 30 | **38 / 30** | glycolysis, reductive carboxylation [LIT Metallo 2011, PMID 22101433] |
| ENO1_del | CN < 0.3 | **0** (hemizygous-low proxy: 39) | ENO2 — **in-silico only**, the 26Q3 panel has no homozygous ENO1 deletion |
| IDH_mut | IDH1 R132 ∪ IDH2 R140/R172 (long format) | **18** (hotspot matrix: 15) | **secondary / in-silico only**: 2-HG is not in the GEM, so the mechanism is not representable |
| **DROPPED** FH_SDH_loss | biallelic-like SDHA-D or FH | **3** | untestable; the permissive definition (67) is almost all monoallelic/passenger |
| **NEG-2** random within-lineage sets | size-matched random member sets per lineage, same contrast, fixed seed list | by design | **primary** negative control for every background contrast |
| **NEG-3** TTN_damaging_nonhyper | TTN damaging, burden ≤ 90th pct (≤ 88) | **162** | no metabolic mechanism |
| **NEG-4** SMARCA4_loss | damaging ∪ CN<0.3 | **97** | frequent driver, no representable enzyme change |
| **NEG-5 (secondary)** KRAS_hotspot | hotspot matrix | — | **not** metabolically neutral: KRAS-driven PDAC uses GOT1-dependent glutamine metabolism and macropinocytosis [LIT Son 2013, PMID 23535601]; keep only as a "non-metabolic-driver" secondary control |
| Zero/near-zero, recorded as not attempted | deep deletion of PGD, UMPS, PTDSS2 (0 each); HPRT1-null (1); FH deep loss (1) | — | — |

---

# 3. Final isozyme exceptions (`OR` relations never edited automatically)

Protected against per-line expression- or proteomics-driven pruning and against `OR → AND` conversion by a co-abundance/complex-membership rule:

- Glycolysis and its branches: `MAR04363` ENO1/ENO3/ENO2; `MAR04375`, `MAR04356`, `MAR08762` ALDOC/ALDOB/ALDOA; `MAR04281`, `MAR04280`, `MAR01453` LDHB/LDHA, LDHB/LDHAL6B; `MAR04319` HKDC1/HK1/HK2/HK3; `MAR04379`, `MAR04301`, `MAR04567` PFKP/PFKL/PFKM; `MAR04171`, `MAR04193`, `MAR04210` PKM/PKLR; `MAR04371` PGAM2/PGAM1/BPGM/PGAM4 (**PGAM1/PGAM2 added**).
- Transaminases and 1C: `MAR03829` GOT1/GOT1L1; `MAR03912`, `MAR06725` GOT1/GOT2; `MAR04340` MTHFD2/MTHFD1L/MTHFD2L.
- Nucleotide terminal steps: IMPDH1/IMPDH2, CTPS1/CTPS2, NME1/NME2, ADSS1/ADSS2, PRPS1/PRPS2.
- Energy-transfer and complex-subunit isoforms (**added**): CKB/CKM (and CKMT1A/CKMT1B); inside the complex IV GPR, `(COX4I1 or COX4I2)`, `(COX6A1 or COX6A2)`, `(COX6B1 or COX6B2)`, `(COX7A1 or COX7A2)` [MS, verified in the complex IV AND]; V-ATPase ATP6V1B1/ATP6V1B2 and ATP6V0A1-A4.
- Others: GLUD1/GLUD2, ME1/ME2/ME3 (keep as separate reactions, never merge), ACSS3/ACSS2/ACSS1, SULT2B1/SULT1A3 (`MAR01919`).

**Modified from the approved list — ATP5MC1/2/3.** These three do encode an identical mature c-subunit, but in the released models they are **AND-ed** inside the single ATP-synthase GPR (`ATP5F1D and … and ATP5MC1 and ATP5MC2 and ATP5MC3 and …`) [MS]. So there is no `OR` to protect; the needed curation is `AND → OR` for this triplet, and the same applies to the RRM2B case in §1.3.

**Modified — tissue-restricted isoforms (PKLR, PDHA2, MAT1A, ENO3, ALDOB, CKM).** Removing such a branch is legitimate curation **only** as a pan-cancer-absence operation: the gene must be shown unexpressed across the whole panel, the test logged, and the edit applied identically to every model. Per-line pruning of these branches remains forbidden, because that is exactly what converts §1.3 buffering pairs into false positives.

**Curation targets (not exceptions):** `TECR or TECRL` (`MAR02194`); `PTGR1 or HSD17B12 or TECRL` (`MAR01074`); the 27-reaction `ALDH3A2 or ALDH2 or ALDH1B1 or ALDH9A1 or ALDH7A1` block; the large CYP `OR` blocks in retinol/arachidonate metabolism; `SEC11A and SEC11C`; `MAR06435`/`MAR06436`; the RNR AND; the ATP5MC triplet; any HMGCS1/HMGCS2 relation.

---

# 4. Final metabolic-drug list

**Inclusion rules** (adopted): the primary mechanism must target a Human-GEM enzyme or transporter; kinase inhibitors are excluded (this drops dichloroacetate/PDK1, OSU-03012 and BX-912); every compound carries a prodrug flag with its activating gene(s); a compound needs ≥100 lines with features in PRISM, GDSC2 or CTRP, checked by script — compounds absent from the local screens are dropped by that script, not by hand.

| Class | Compounds (target; prodrug → activating genes) |
|---|---|
| Antifolates | methotrexate (DHFR; uptake SLC19A1, FPGS), pemetrexed (TYMS/DHFR/GART), pralatrexate, raltitrexed, pyrimethamine, trimethoprim |
| TYMS | 5-FU (prodrug → UMPS or UPP1/UCK2, TYMP/TK1; catabolism DPYD), floxuridine (→ TK1), capecitabine (→ CES1/2, CDA, TYMP), trifluridine+tipiracil (TYMP) |
| RNR | hydroxyurea (RRM2), gemcitabine (RRM1; → DCK), triapine (RRM2) |
| DHODH | brequinar, teriflunomide, leflunomide (non-enzymatic conversion), BAY-2402234 |
| IMPDH | mycophenolic acid / mofetil (→ CES1), ribavirin (→ ADK), mizoribine (→ ADK) |
| Purine analogues | 6-mercaptopurine, 6-thioguanine (→ HPRT1; de novo purine), fludarabine, cladribine, clofarabine (→ DCK), azaserine, allopurinol (XDH) |
| Pyrimidine analogues | cytarabine (→ DCK), decitabine, azacitidine |
| Glutaminase / amidotransferase | CB-839/telaglenastat, BPTES (GLS), acivicin (CTPS1) |
| Statins and FDPS | simvastatin, lovastatin (lactone prodrugs), atorvastatin, fluvastatin, pitavastatin, rosuvastatin (HMGCR); zoledronate, alendronate, ibandronate, risedronate (FDPS) |
| SQLE | terbinafine |
| FASN / ACLY / SCD / ACC | TVB-3166, TVB-2640/denifanstat, GSK2194069, orlistat, cerulenin (FASN); BMS-303141, SB-204990, NDI-091143 (ACLY); bempedoic acid (ACLY; → SLC27A2, liver-restricted, flagged); A939572, CAY10566, MF-438 (SCD); ND-646, firsocostat, BAY ACCi (ACACA/ACACB) |
| CPT1 / FAO | etomoxir (CPT1A; off-target at high dose, flagged), perhexiline |
| Glycolysis and transporters | 2-deoxyglucose (HK1/HK2), lonidamine, 3-bromopyruvate, TT-232/TLN-232, shikonin (PKM2); AZD3965 (SLC16A1) [LIT Polański 2014, PMID 24277449]; BAY-876 (SLC2A1); GSK2837808A, GNE-140 (LDHA/LDHB) |
| OXPHOS | IACS-010759, rotenone, phenformin, metformin (weak, mechanism debated, flagged) (complex I); antimycin A, atovaquone (complex III, MT-CYB); oligomycin-A (ATP synthase); CPI-613/devimistat (PDHA1, DLD, DLST) |
| Redox / ferroptosis | buthionine sulfoximine (GCLC/GCLM), sulfasalazine and erastin (SLC7A11; erastin also hits VDAC, flagged), disulfiram (ALDH2), elesclomol, N-acetyl cysteine, α-lipoic acid |
| NAD | FK866/daporinad, GNE-617, GMX1778/CHS-828, STF-118804 (NAMPT; modifier NAPRT). KPT-9274 excluded (dual PAK4) |
| 1C / SHMT | SHIN1, SHIN2 (SHMT1/2) [LIT Ducker 2017, PMID 29073064] |
| Polyamines | eflornithine/DFMO (ODC1) |
| Mutant IDH | ivosidenib, AGI-5198 (IDH1), enasidenib, AGI-6780 (IDH2), vorasidenib. **Negative-expectation class**: 2-HG is absent from the GEM, so the only handle is IDH expression |
| Off-model controls | tigecycline, doxycycline, chloramphenicol, linezolid (mitoribosome); bortezomib, belinostat, alisertib, mitoxantrone (no model target) |

**Positive-control pairs**: DHODH–brequinar; DHFR–methotrexate; TYMS–5-FU/floxuridine/pemetrexed; IMPDH2–mycophenolate; RRM1/RRM2–hydroxyurea/gemcitabine; HMGCR–statins; SQLE–terbinafine; FASN–GSK2194069/TVB-2640; ACACA–ND-646; GLS–CB-839/BPTES; SLC16A1–AZD3965; SLC2A1–BAY-876; CPT1A–etomoxir; NAMPT–FK866 (with NAPRT as modifier); GCLC/GCLM–BSO; SLC7A11–sulfasalazine/erastin; ATP5F1A–oligomycin; PDHA1/DLD/DLST–CPI-613; HK1/HK2–2-DG; ODC1–DFMO; SHMT1/2–SHIN1/2; IDH1–ivosidenib/AGI-5198 (restricted to IDH_mut).

**Availability note (verified in the local files)** [OMICS]: present in PRISM and/or GDSC — methotrexate, pemetrexed, pralatrexate, raltitrexed, 5-FU, floxuridine, capecitabine, trifluridine, tipiracil, gemcitabine, cytarabine, hydroxyurea, brequinar, teriflunomide, leflunomide, mycophenolic acid, mercaptopurine, fludarabine, cladribine, clofarabine, azaserine, allopurinol, acivicin, CB-839, BPTES (GDSC 1425), AZD3965, UK-5099, 2-deoxyglucose, lonidamine, 3-bromopyruvate, TT-232 (GDSC 436), shikonin (GDSC 170), oligomycin-A, atovaquone, phenformin, metformin, CPI-613 (GDSC 415), simvastatin, lovastatin, atorvastatin, terbinafine, orlistat, cerulenin, GSK2194069, BAY ACCi (GDSC 1362), etomoxir, perhexiline, BSO, sulfasalazine, erastin, disulfiram, elesclomol, daporinad/FK866 (GDSC 1248), STF-118804, ivosidenib, enasidenib, vorasidenib, AGI-5198 (GDSC 1913), AGI-6780 (GDSC 562/1634), EPZ015666, GSK3326595. **Not found locally** (the ≥100-line script will drop them): IACS-010759, antimycin A, rotenone (only dihydrorotenone, GDSC 1827), ML210, RSL3, FIN56, TVB-2640, TVB-3166, ND-646, firsocostat, BAY-876, GSK2837808A, GNE-140, SHIN1/2, DFMO, triapine, BAY-2402234, mizoribine, ribavirin, bisphosphonates, ACLY and SCD compounds, AG-270, olutasidenib, ADI-PEG20, asparaginase, syrosingopine, KPT-9274.

---

# 5. Final candidate task sets and objective/constraint variants

Fixed a priori; any choice among them must be made on held-out lineages, never against the benchmark. Task-file inventory [MS]: Essential 57 tasks (ER 9, BS 18, SU 17, IC 12, GR 1); Full 256; CellFie consensus 195; VerifyModel 21.

| Set | Content | Size | Rationale |
|---|---|---|---|
| T0 | `GR` growth on Ham's only | 1 | baseline reference |
| T1 | all 57 Essential tasks as constraints, `GR` objective | 57 | viability, not only growth [LIT Robinson 2020, PMID 32209698] |
| T2 | the 9 Essential ER tasks | 9 | fixes the no-maintenance artefact (§1.5i) |
| T3 | 8 Essential BS nucleotide tasks | 8 | dominant class of true proliferative dependencies |
| T4_k | the 12 Essential IC precursor tasks, enforced singly or in groups (k = 1…12) | 12 | central-carbon genes become load-bearing without editing biomass |
| T5 | Essential BS PC/PE/PS/PI + cholesterol + SU choline/inositol (+ Full 104-111) | 7 (+8) | membranes are the largest non-protein biomass sink; serum lipids are the main medium variable |
| T6 | Essential BS TPP/CoA/FAD/heme + Full 112-118 | 4 (+7) | replaces the dead-end cofactor pool (§1.4a) with defensible demands |
| T7 | CellFie consensus, **scored** as a task-activity vector, not enforced | 195 | context-specific model building and flux features [LIT Richelle 2021, PMID 34761247] |
| T8 | Full 42-50, 52-68, 77-99 (amino-acid synthesis, uptake, degradation) | ~49 | amino-acid composition is where the media differ most |
| T9 | every `SHOULD FAIL = TRUE` task (Full 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 72, 74, 158, 170, 257) + the 21 VerifyModel leak tests | 37 | model-quality gate and negative control; never an objective |
| D-BS, D-ER | demand variants: force BS or ER task flux at 1% of its maximum | 2 | graded, not binary, task enforcement |
| O-cap, O-NGAM, O-MOMA | objective variants: capped uptake; non-growth ATP maintenance on `MAR03964` at k ∈ {0, 0.1, 0.3} of WT ATP turnover; MOMA instead of FBA for the KO | 5 | `MAR03964` currently has bounds [0, 1000] [MS], which is the direct cause of §1.5i |
| B-loo_j | leave-one-out on the biomass reaction. **Definition corrected:** `MAR13082` has **9** direct reactants (ATP, DNA, H2O, RNA, glycogen, cofactor pool, protein pool, lipid pool, metabolite pool), not 13 [MS]. The set is the direct reactants except ATP and H2O, plus the artefact-linked cofactor-pool groups (retinoid and vitamin-D catabolites, tocopherol catabolites, carnitine, heme, ubiquinol, lipoate, cobamide); enumerated per parent by script and frozen | ~14 | isolates each §1.4a-c false-positive driver |
| B-add | biomass-addition variants: heme A, dolichol-P, GPI-anchor precursor, Fe-S cluster, SAM sink, spermidine, **NAD+ turnover sink**, **GSH turnover sink**. **Rule:** skip and log any component already in the parent's biomass — ubiquinol, heme and lipoate are already in `MAR10065` and CoA is in `MAR10064` [MS], so adding them would only change coefficients | 8 | the NAD+ sink is the principled fix for the NAMPT false negative (§1.5h), the GSH sink for GCLC/GSS |
| T10 | de-artefacted biomass: `MAR13082` with the cofactor-pool and glycogen terms zeroed, plus T3 as constraints | 1 + 8 | one pre-registered variant that removes the two documented FP drivers; declared here precisely because it edits biomass and must not be tuned |
| Controls | negative: random pseudo-tasks, and a call-rate-matched comparator. Positive: **no loss of mevalonate-core essentiality** when serum cholesterol is supplied under B-add — informative only for the genes whose M0 ratio is 1 (HMGCS1, IDI1, FDPS, GGPS1), where the cause is GPR/isozyme rather than demand; HMGCR, MVK, PMVK, MVD are already lethal on M0 because ubiquinol is in biomass [MS/PRED] | — | guards against task-set overfitting |

---

# A. Difference table

Resolution key: **K** keep (in both, or blind-only and retained) · **A** add (approved-only, independently judged justified) · **D** drop · **M** modify.

| # | Item | Blind | Approved | Res. | Biological reason / citation |
|---|---|---|---|---|---|
| **1.1 Pan-medium essentials** | | | | | |
| 1 | CAD, DHODH, UMPS | 1.1 | 1.1→1.2 | **M** | Uridine rescues the whole ring; HPLM and Plasmax contain uridine, and uric acid inhibits UMPS [PMID 28388410]. Blind text already said "serum uridine"; moved for consistency |
| 2 | TYMS | 1.1 **and** 1.2 | →1.2 | **M** | Thymidine + TK1 fully rescues; the blind file listed it twice — a real inconsistency, resolved to 1.2 |
| 3 | DHFR | 1.1 and 1.2 | →1.2 | **M** | HT rescue [PMID 6933469]; same duplication resolved |
| 4 | NAMPT | 1.1 and 1.2 | →1.2 only | **M** | Bypassed only with NAPRT **and** nicotinic acid, which no medium supplies [PMID 17604275]; duplication resolved |
| 5 | PNPO | 1.1 and 1.2 | — | **M** | IMDM supplies pyridoxal, bypassing PNPO but not PDXK [MS]; moved to 1.2 |
| 6 | Purine block PPAT/GART/PFAS/PAICS/ATIC | 1.1 | 1.2 | **M** | Hypoxanthine+HPRT1 rescue; retained in 1.2 |
| 7 | ADSL, IMPDH2, GMPS | 1.1 (ADSL also in 1.2 purine block) | "not hypoxanthine-rescued" | **M** | Correct and a genuine error in the blind list: hypoxanthine gives IMP only; AMP needs ADSS2→ADSL, GMP needs IMPDH→GMPS. Kept in 1.1 |
| 8 | 20 ARSs | 1.1A | — | **K** | Charged tRNAs feed `MAR10062` [MS]; scored separately as non-metabolic |
| 9 | RRM1, RRM2 | 1.1 | A | **K** | No dNDP salvage in any standard medium |
| 10 | GUK1 | — | A | **A** | Sole GMP kinase; single-gene GPR on `MAR04020`/`MAR04135` [MS], ratio 0.999 [PRED] — justified independently |
| 11 | CS | — | A (moderate) | **A** | Single citrate synthase, no isozyme; weaker in OXPHOS-low lines because reductive carboxylation can supply citrate. Tagged `moderate` |
| 12 | MDH2 | — | A (moderate) | **A** | Needed for TCA flux and aspartate synthesis; pyruvate relaxes it [PMID 26232224, 26232225] |
| 13 | SDHA-D as a new OXPHOS/TCA `moderate` class | partial predictions only | M | **A** | SDH-deficient tumours exist; pyruvate/uridine relax the requirement [PMID 2814477]. Model ratios 0.76/0.88 [PRED] |
| 14 | IDI1, GGPS1 | — | A | **A** | IDI1 is the only non-peroxisomal IPP isomerase of consequence; GGPS1 is the only GGPP source for Rab/Rho prenylation [PMID 1967820] |
| 15 | HMGCS1 vs HMGCS2 compartment GPR flag | — | A | **A** | HMGCS2 is mitochondrial/ketogenic; any OR with HMGCS1 is a compartment error [MS] |
| 16 | CTPS1 | 1.3 pair | 1.1→pair | **K** | Agreement: CTPS2 buffers; belongs in 1.3 |
| 17 | CHKA | — | →pair | **A** | CHKB has overlapping choline-kinase activity; both ratio 1 [PRED]. Added to 1.3, not 1.1 |
| 18 | PCYT1A, COASY | 1.1 | A | **K** | Rate-limiting PC step; last CoA steps. PCYT1B noted as a low-expressed paralog |
| 19 | Fe-S module (H2 only); lipid core; hexosamine core | 1.1F/G/I | — | **K** | Blind-only, retained; CIAO1/CIAO2A/NUBP1/NUBP2 are absent from both models [MS] |
| 20 | Glycolysis core | 1.1H | — | **M** | Retained but tagged `model-masked`: all ratio 1 because of isozyme OR [PRED] |
| **1.2 Medium-dependent** | | | | | |
| 21 | ASNS | 1.2 | A | **K** | DMEM high-glucose has no asparagine, RPMI has 0.38 mM [vendor sheets; PMID 27126896, 33046912] |
| 22 | Sterol branch | 1.2 | A + DHCR24 caveat | **M** | Desmosterol partly substitutes for cholesterol, so DHCR24 is tagged `weaker` |
| 23 | GLUL | — | A | **A** | Glutamine synthetase is a dependency only in Gln-free/low-Gln medium; ratio 1 [PRED] |
| 24 | ALDH18A1, PYCR1/2/3 | — | A | **A** | DMEM high-glucose contains no proline [vendor sheet]; all ratio 1 [PRED] |
| 25 | GOT1, GOT2 + SLC1A3 gate | 1.5d only (as FN) | A | **A** | Aspartate is the limiting output of respiration; pyruvate relaxes it and uptake needs SLC1A3 [PMID 26232224, 26232225, 29941933] |
| 26 | SLC7A11 | 1.2 | A + transsulfuration note | **M** | CBS/CTH from methionine, and 2-ME-reduced cysteine, are the alternatives |
| 27 | HPLM / Plasmax as medium variants | — | implied by uridine | **A** | Both are in `data/media/`; they change the pyrimidine and purine expectations |
| 28 | PHGDH/PSAT1/PSPH; ASS1/ASL; NAMPT/NAPRT | 1.2 | A | **K** | Agreement [PMID 20104527, 26560030, 17604275] |
| 29 | 15 further medium rows (thymidine, HT, pyridoxal, lipoate, putrescine, biotin, B12, carnitine, heme, selenium, serum FA, pyruvate, RPMI GSH, inositol, cystine) | 1.2 | — | **K** | Blind-only, each tied to a vendor formulation sheet and an exchange ID [MS] |
| **1.3 Collateral lethality / buffering** | | | | | |
| 30 | ENO1→ENO2; ME2→ME3; NAPRT→NAMPT; MTAP→MAT2A/PRMT5 | 1.3 | A | **K** | [PMID 22895339, 28099419, 24097869/31019297, 26912361/26912360/27068473] |
| 31 | "Cerna 2012" as the NAPRT citation | — | A | **D** | Could not be resolved; replaced by Shames 2013 (PMID 24097869) and Chowdhry 2019 (PMID 31019297) |
| 32 | CHKA/CHKB, CTPS1/CTPS2, IMPDH1/IMPDH2 | 2 of 3 | add | **A** | Paralog-buffered; they test expression-aware OR logic [PMID 31652272, 33637726] |
| 33 | RRM2/RRM2B and the RNR GPR error | — | add | **A** | Verified independently: all 14 RNR reactions use `(RRM2B and TXN2 and RRM1 and RRM2) or (…TXN…)`, so RRM2B is lethal in both models (ratio 0) [MS/PRED]. Correct logic `RRM1 and (RRM2 or RRM2B)` |
| 34 | MTA → de novo purine reliance in MTAP-null | — | A | **A** | MTAP salvages adenine from MTA; its loss raises de novo purine demand |
| 35 | 24-pair buffering list; FH→HMOX1; ASS1 aspartate diversion | 1.3 | — | **K** | Blind-only, retained [PMID 21849978, 26560030] |
| **1.4 False positives** | | | | | |
| 36 | H2 β-oxidation lethals (ACADL, HADHA, HADHB, ECHS1, ACAA2, AMACR, ACADSB) | 1.4a | A + ECHS1 caveat | **M** | Mechanism verified (`MAR06433/06435/06436` AND-chains feeding `MAR10065`) [MS]; ECHS1 tagged `uncertain` because methacrylyl-CoA toxicity in valine catabolism is real but non-stoichiometric |
| 37 | Xenobiotic CYP/UGT/GST/SULT/NAT | 1.4a (CYP only) | A | **M** | Broadened to the whole xenobiotic class; no biomass role |
| 38 | Urea cycle (CPS1, OTC, ARG1) | — | A | **A** | Inactive in most lines; all ratio 1 [PRED], so it is a *latent* FP class that medium/task changes could switch on. CPS1 exception in STK11-loss [PMID 28538732] |
| 39 | Bile-acid and steroid-hormone synthesis | implicit (AMACR) | A | **A** | Tissue-restricted pathways with no biomass role |
| 40 | Biomass-artefact lethals (CYP26A1/B1, CYP24A1, CYP1A2, CYP4F2, BBOX1, TMLHE, KL …) | 1.4a | A | **K** | Agreement; blind version adds the mechanism and 10 further genes |
| 41 | Pool dead-ends (PSAT1, UGDH, PGM3, GNPNAT1), glycogen/UGP2, pseudo-reaction GPRs, single-gene transport, unbounded-uptake partials | 1.4b-h | — | **K** | Blind-only, each verified in the model [MS] |
| **1.5 False negatives** | | | | | |
| 42 | V-ATPase (ATP6V*) | — | A | **A** | Essentiality is non-stoichiometric (lysosomal pH, mTORC1); all ratio 1 [PRED] |
| 43 | PGD | — | A | **A** | 6-phosphogluconate accumulation and LKB1-AMPK inhibition, not flux [PMID 26479318]; ratio 1 [PRED] |
| 44 | GPX4/ferroptosis; mitochondrial translation and assembly | 1.5e | A | **K** | [PMID 24439385] |
| 45 | Structural non-predictions (61% of genes), missing genes, isozyme masking, compartment blindness, missing biomass demands, alternative routes, no NGAM | 1.5a-g | — | **K** | Blind-only, all verified [MS/PRED] |
| **1.6 CRISPR artefacts** | | | | | |
| 46 | Multi-copy seed loci (UGT1A, UGT2B, GSTM1-5, GSTT1/2/2B, AMY1A-C, NAT1/2, GBA1/GBAP1, AMY2A/2B, SULT1A3/1A4, CYP2D6/2D7, CYP4F cluster, AKR1C1-3, ACSM2A/2B, GSTA1/2, CES1/CES1P1) | partial (AMY, SULT, GSTT2, AKR1C1/2, ACSM2, GSTA) | A + adds | **A** | Justified: these are shared-exon or pseudogene-adjacent loci where guides cross-target. GBA1 and CYP2D7 are not in the models [MS] |
| 47 | ≥95% identity rule | protein-coding paralogs only (67 genes) | ≥95% paralog **or pseudogene** | **M** | Extended to pseudogenes — that is what catches UGT1A, GBAP1, CES1P1, CYP2D7 |
| 48 | Housekeeping-adjacent families | implicit in 1.4d | A | **A** | Family-level definition is reproducible; noted that CIAO*/NUBP1/2 are absent from the models |
| 49 | Moonlighting list | 11 genes, no citations | 11 genes with citations | **M** | Merged; every citation resolved [PMID 1502165, 21640161, 10802057, 22056988, 11751859, 19461003, 28552616, 24995980, 26237645] |
| 50 | mtDNA stratum; amplicon frequencies; recurrently deleted genes; gene-ID hazards | 1.6a, c, d, e | — | **K** | Blind-only, computed from non-label omics [OMICS] |
| **2. Genomic groups** | | | | | |
| 51 | MTAP_del | 336/341/400 | A | **M** | Added the de novo purine expectation |
| 52 | ASS1_silenced (n=254) | A | A + CAD | **M** | Aspartate diversion to pyrimidines [PMID 26560030] |
| 53 | NRF2_active (n=68) | A | + GLS, G6PD, PGD, ME1, IDH1 | **M** | Glutaminolysis dependence and NADPH imbalance [PMID 28967920, 28967864] |
| 54 | STK11_loss (n=82) | A | + DTYMK | **M** | [PMID 23715154]; DTYMK is generic-lethal so it cannot discriminate — kept with that caveat |
| 55 | FH_SDH_loss | permissive 67, biallelic-like **3** | drop (2 lines) | **D** | Independent agreement on dropping: the permissive set is almost entirely monoallelic/passenger |
| 56 | IDH_mut | 15 (matrix) / **18** (long format) | drop | **M** | Demoted rather than dropped: 18 lines exist, but 2-HG is not in the GEM, so the mechanism is unrepresentable — secondary/in-silico only |
| 57 | ENO1_del | **0** deep deletions (39 hemizygous-low) | A (in-silico) | **M** | Agreement, with my count as the reason: in-silico only |
| 58 | KRAS_hotspot as negative control | not used as a control | M | **M** | Accepted: KRAS is not metabolically neutral (GOT1-dependent glutamine metabolism, macropinocytosis) [PMID 23535601] → secondary control only |
| 59 | Random within-lineage size-matched sets | — | add | **A** | The only control that holds lineage composition fixed; becomes the primary negative control |
| 60 | TTN_damaging_nonhyper (162), SMARCA4_loss (97) | NEG-3, NEG-4 | — | **K** | Blind-only, retained as event-matched negative controls |
| 61 | CDKN2A_del_MTAP_intact (157); ASNS_low; SQLE_null; HK1/HK2; LDHB_null; SLC7A11_null; CDS1_null; VHL; PHGDH amp-or-high; STK11 subgroups | blind | — | **K** | Blind-only, each with a counted definition [OMICS] |
| **3. Isozyme exceptions** | | | | | |
| 62 | LDHA/B, PFKL/M/P, ENO1/2/3, ALDOA/B/C | 3 | A | **K** | Mixed or alternative oligomers; OR is biochemically correct |
| 63 | CKB/CKM, COX4I1/2, COX6A1/2, COX7A1/2, ATP6V1B1/2, ATP6V0A1-4, PGAM1/2, PRPS1/2, CTPS1/2, IMPDH1/2 | 4 of 10 | A/add | **A** | Verified that the COX pairs are already ORs inside the complex IV AND [MS]; subunit-isoform ORs must survive a co-abundance rule |
| 64 | ATP5MC1/2/3 | — | A (protect OR) | **M** | Independently checked: the released models **AND** them inside the ATP-synthase GPR [MS]. There is no OR to protect; the fix is AND→OR |
| 65 | RRM2/RRM2B | — | protect vs OR→AND | **A** | Accepted, modified: the current state is already an AND error, so a manual AND→OR fix is required, not merely protection |
| 66 | Tissue-restricted isoforms (PKLR, PDHA2, MAT1A, ENO3, ALDOB, CKM) | protected (PKLR, ENO3, ALDOB) | eligible for removal | **M** | Compromise, stated as a rule: removal allowed only as a logged pan-cancer-absence edit applied to every model; per-line expression pruning stays forbidden |
| 67 | Curation targets (TECR/TECRL, ALDH OR block, CYP ORs, SEC11A/C, MAR06435) | 3 | — | **K** | Blind-only, retained |
| **4. Drugs** | | | | | |
| 68 | Inclusion rules (GEM target, no kinase inhibitors, prodrug flag, ≥100 lines) | — | A | **A** | Makes the list reproducible and drops look-alike targets |
| 69 | Dichloroacetate (PDK1) | in list | excluded by rule | **D** | PDK1 is a kinase, not a GEM reaction |
| 70 | Prodrug activators for 5-FU, capecitabine, gemcitabine, mycophenolate, thiopurines, statin lactones, bempedoic acid | — | A | **A** | The activating gene is often the real dependency (DCK, UMPS, UPP1/UCK2, HPRT1, CES1) |
| 71 | Added compounds: triapine, BAY-2402234, mizoribine, ribavirin, capecitabine, bisphosphonates, TVB-3166/2640, ACLY and SCD compounds, ND-646/firsocostat, GNE-617, GMX1778, IACS-010759, rotenone, antimycin A, BAY-876, GSK2837808A/GNE-140, SHIN1/2, DFMO | partial | A/add | **A** | All hit a Human-GEM enzyme or transporter; availability is decided by the ≥100-line script [PMID 24277449, 29073064] |
| 72 | IDH inhibitors | flagged (no 2-HG) | A (negative-expectation) | **K** | Agreement |
| 73 | Verified PRISM/GDSC names and IDs; "absent from local screens" list; off-model negative-control drugs | blind | — | **K** | Blind-only; it is what makes the ≥100-line rule executable now [OMICS] |
| **5. Task sets** | | | | | |
| 74 | Essential-57 and CellFie-195 | T1, T7 | A | **K** | Fixed external lists [PMID 32209698, 34761247] |
| 75 | B-loo set size | 9 direct reactants [MS] | 13 → corrected to 9 | **A** | Independent agreement; the blind file lists the same 9 |
| 76 | B-add "skip if already in biomass" | — | M | **A** | Verified: ubiquinol, heme, lipoate are in `MAR10065`; CoA is in `MAR10064` [MS] |
| 77 | NAD+ and GSH turnover sinks | — | add | **A** | The NAD+ sink is the mechanistic fix for the NAMPT false negative (§1.5h); the GSH sink for GCLC/GSS |
| 78 | O-NGAM k ∈ {0, 0.1, 0.3}; O-cap; O-MOMA; D-BS; D-ER | NGAM proposed in prose | A | **A** | `MAR03964` has bounds [0, 1000] [MS], the direct cause of the OXPHOS partials |
| 79 | Mevalonate positive control restated | — | A | **A** | HMGCR/MVK/PMVK/MVD are already lethal on M0 via ubiquinol; only HMGCS1, IDI1, FDPS, GGPS1 are informative [PRED] |
| 80 | Random pseudo-tasks + call-rate-matched comparator | — | A | **A** | Guards against task-set overfitting |
| 81 | T0, T2-T6, T8-T10 | blind | — | **K** | Blind-only, retained; T10 is the one pre-registered biomass edit |

---

# B. Citation ledger (all resolved against PubMed in this session)

Every PMID below was checked by retrieving the record and matching first author, journal, year, volume, pages and title. According to PubMed:

**Corrected (the blind file's PMID was wrong):** Cantor 2017 Cell 169:258 → **28388410** ([DOI](https://doi.org/10.1016/j.cell.2017.03.023)); Garcia-Bermudez 2019 Nature 567:118 → **30760928** ([DOI](https://doi.org/10.1038/s41586-019-0945-5)); Krall 2016 Nat Commun 7:11457 → **27126896** ([DOI](https://doi.org/10.1038/ncomms11457)); Dey 2017 Nature 542:119 → **28099419** ([DOI](https://doi.org/10.1038/nature21052)); Romero 2017 Nat Med 23:1362 → **28967920** ([DOI](https://doi.org/10.1038/nm.4407)); Chowdhry 2019 Nature 569:570 → **31019297** ([DOI](https://doi.org/10.1038/s41586-019-1150-2)); Richelle 2021 Cell Rep Methods 1:100040 → **34761247** ([DOI](https://doi.org/10.1016/j.crmeth.2021.100040)); Muller 2012 Nature 488:337 → **22895339**; Mavrakis 2016 Science 351:1208 → **26912361** and Kryukov 2016 Science 351:1214 → **26912360** (the blind file had these two swapped).

**Confirmed as written:** Ducker & Rabinowitz 2017 Cell Metab 25:27 — 27641100 ([DOI](https://doi.org/10.1016/j.cmet.2016.08.009)); Urlaub & Chasin 1980 PNAS 77:4216 — 6933469 ([DOI](https://doi.org/10.1073/pnas.77.7.4216)); Rabinovich 2015 Nature 527:379 — 26560030 ([DOI](https://doi.org/10.1038/nature15529)); King & Attardi 1989 Science 246:500 — 2814477 ([DOI](https://doi.org/10.1126/science.2814477)); Shames 2013 Clin Cancer Res — 24097869 ([DOI](https://doi.org/10.1158/1078-0432.CCR-13-1186)); De Kegel & Ryan 2019 PLoS Genet — 31652272 ([DOI](https://doi.org/10.1371/journal.pgen.1008466)); Thompson 2021 Nat Commun 12:1302 — 33637726 ([DOI](https://doi.org/10.1038/s41467-021-21478-9)); Gammage 2018 Trends Genet 34:101 — 29179920 ([DOI](https://doi.org/10.1016/j.tig.2017.11.001)); Kim 2017 Nature 546:168 — 28538732 ([DOI](https://doi.org/10.1038/nature22359)); Gross 2014 Mol Cancer Ther 13:890 — 24523301 ([DOI](https://doi.org/10.1158/1535-7163.MCT-13-0870)); Dixon 2012 Cell 149:1060 — 22632970 ([DOI](https://doi.org/10.1016/j.cell.2012.03.042)); Viswanathan 2017 Nature 547:453 — 28678785 ([DOI](https://doi.org/10.1038/nature23007)); Frezza 2011 Nature 477:225 — 21849978; Dang 2009 Nature 462:739 — 19935646; Metallo 2011 Nature 481:380 — 22101433; Meyers 2017 Nat Genet 49:1779 — 29083409; Robinson 2020 Sci Signal 13:eaaz1482 — 32209698; Possemato 2011 Nature 476:346 — 21760589; Locasale 2011 Nat Genet 43:869 — 21804546.

**Added for the reconciled items:** Birsoy 2015 Cell 162:540 — 26232224; Sullivan 2015 Cell 162:552 — 26232225; Garcia-Bermudez 2018 Nat Cell Biol 20:775 — 29941933; Lagziel 2020 Nat Metab 2:1369 ("Mind your media") — 33046912; Hara 2007 J Biol Chem 282:24574 — 17604275; Delage 2010 Int J Cancer 126:2762 — 20104527; Liu 2013 Cancer Discov 3:870 — 23715154; Sayin 2017 eLife 6 — 28967864; Goldstein & Brown 1990 Nature 343:425 — 1967820; Lin 2015 Nat Cell Biol 17:1484 — 26479318; Yang 2014 Cell 156:317 — 24439385; Ducker 2017 PNAS 114:11404 (SHIN inhibitors) — 29073064; Son 2013 Nature 496:101 — 23535601; Dempster 2021 Genome Biol 22:343 (Chronos) — 34930405; Marjon 2016 Cell Rep 15:574 — 27068473; Polański 2014 Clin Cancer Res 20:926 — 24277449; Wellen 2009 Science 324:1076 — 19461003; Li 2017 Mol Cell 66:684 — 28552616; Sutendra 2014 Cell 158:84 — 24995980; Jiang 2015 Nat Cell Biol 17:1158 — 26237645; Haile 1992 PNAS 89:7536 — 1502165; Sirover 2011 Biochim Biophys Acta 1810:741 — 21640161; Feo 2000 FEBS Lett 473:47 — 10802057; Yang 2011 Nature 480:118 — 22056988; Pastorino 2001 J Biol Chem 277:7610 — 11751859.

**Dropped as unresolvable:** "Cerna 2012" (NAPRT), replaced by PMIDs 24097869 and 31019297. No citation in this file carries a "(verify)" marker.

---

# C. Attestations

**"Author has never viewed DepMap label values or label aggregates."** — True. The probe `head -c 20 data/depmap/CRISPRGeneEffect.csv` was blocked by the label guard before any work began and was never retried by any route. I read no file under `data/labels` or `data/locked`, and no `CRISPRGeneEffect.csv`, `CRISPRGeneDependency.csv`, `ScreenGeneEffect.csv`, `CRISPRInferredCommonEssentials.csv`, `AchillesCommonEssentialControls.csv` or `AchillesNonessentialControls.csv`. The only DepMap files I opened are `Model.csv`, `OmicsCNGeneWGS.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, `OmicsSomaticMutations.csv`, `OmicsSomaticMutationsMatrixDamaging.csv` and `OmicsSomaticMutationsMatrixHotspot.csv`, plus the PRISM, GDSC and Ensembl files. I viewed no label value and no label statistic of any kind, including aggregates; where the approved file discloses that its author saw two categories of all-line aggregates, it does not state their values and I did not seek them.

**"Author has never computed, run or inspected the project's train/validation/test split draw for any real seed."** — True. I never executed, imported, read or inspected `scripts/make_splits.py` or `src/nextgem/splits.py`, and I ran nothing that calls them. The only project code I used was `src/nextgem/model.py` (read, to confirm the medium/objective definitions) and `results/baseline/*` (read). Everything else I ran was my own throwaway scripts in the session scratchpad: an SBML/GPR dump via cobra, pandas reads of the non-label omics files, and an Ensembl-paralog identity count. The "random within-lineage size-matched sets" control in §2 and the seeded random pseudo-tasks in §5 are *specified* here and were not drawn.
