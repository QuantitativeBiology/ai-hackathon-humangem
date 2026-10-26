# A priori lists for the essentiality benchmark (blind)

Author has not read DESIGN.md, LITERATURE.md, RESEARCH_TREE.md, notes/ or any label data; lists derived blind.

Probe `head -c 20 data/depmap/CRISPRGeneEffect.csv` was **BLOCKED** by the label guard before any work began.

## Sources used
- `PROMPT.md`, `CLAUDE.md`, `DATA.md`, `data/README.md`.
- Models: `data/Human1/model/` (Human1 = v1.19.0) and `data/Human-GEM-main/model/` (Human2 = v2.1.1): SBML, `genes.tsv`, `reactions.tsv`, `Human-GEM.txt`; task files in `data/Human-GEM-main/data/metabolicTasks/`.
- `results/baseline/` (generic-model KO predictions; **predictions only, no labels**).
- Non-label DepMap omics/annotation: `Model.csv`, `OmicsCNGeneWGS.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, `OmicsSomaticMutations*.csv`.
- `data/media/` vendor formulation sheets (RPMI-1640 ATCC, DMEM high-glucose, DMEM/F12, IMDM); `data/prism/`, `data/gdsc/` compound metadata; `data/paralogs/` (Ensembl 116 Compara).

## Evidence tags
`[MS]` model structure (reaction ID / GPR / stoichiometry, verified in this session) · `[PRED]` generic-model KO prediction from `results/baseline` · `[OMICS]` computed from non-label DepMap omics in this session · `[LIT]` literature. Identifiers marked *(verify)* are cited from memory and must be checked before any public write-up; the biology they support does not depend on the identifier.

Gene IDs: Ensembl IDs are identical in Human1 and Human2 for every gene below unless stated; Human-GEM gene IDs are Ensembl IDs (`model/genes.tsv`).

---

# 1. Gene panels

## 1.1 Pan-medium essential metabolic genes (expected essential in any standard medium)

Criterion used: the product/function is not supplied at usable concentration by RPMI, DMEM, DMEM/F12, IMDM, MEM/EMEM, F12, McCoy's 5A or L-15 (± 5-20% FBS), **and** no isozyme or salvage bypass exists in human cells. 58 genes in 9 classes.

**A. Aminoacyl-tRNA charging (20)** — `AARS1, CARS1, DARS1, EPRS1, FARSA, FARSB, GARS1, HARS1, IARS1, KARS1, LARS1, MARS1, NARS1, QARS1, RARS1, SARS1, TARS1, VARS1, WARS1, YARS1` (e.g. AARS1 ENSG00000090861, KARS1 ENSG00000065427). Rationale: charged-tRNA species are the substrates of the protein pseudo-reaction `MAR10062` → `protein_pool_biomass` → biomass `MAR13082` [MS]; all are lethal in both models [PRED]. Biologically pan-essential (translation). **Caveat: mechanism is non-metabolic — see 1.6g; these are "easy positives" that inflate AUPRC.**

**B. De novo pyrimidines (7)** — `CAD` (ENSG00000084774), `DHODH` (ENSG00000102967), `UMPS` (ENSG00000114491), `DTYMK` (ENSG00000168393), `TYMS` (ENSG00000176890), `RRM1` (ENSG00000167325), `RRM2` (ENSG00000171848). No medium above contains uridine, cytidine or orotate; serum uridine is low single-digit µM [LIT Cantor et al. 2017 Cell 169:258, PMID 28388411]. CAD/DTYMK/RRM1/RRM2 lethal in both models; UMPS/DHODH lethal only in Human2 — in Human1 `MAR04644` (TYMS) is reversible so medium thymidine routes dTMP→dUMP→uracil→UMP via `UPRT` (`MAR04343`) [MS, baseline README].

**C. De novo purines (9)** — `PPAT` (ENSG00000128059), `GART` (ENSG00000159131), `PFAS`, `PAICS` (ENSG00000128050), `ADSL` (ENSG00000239900), `ATIC` (ENSG00000138363), `IMPDH2` (ENSG00000178035), `GMPS` (ENSG00000163655), `ADSS2`. RPMI/DMEM/IMDM/MEM/McCoy/L-15 contain no purine source; only F12 and DMEM/F12 contain hypoxanthine (2.39 mg/L) [media sheets] → see 1.2. GMPS lethal in both models; the rest are `ratio_r = 1` [PRED] — a false-negative class (1.5a).

**D. Folate / one-carbon (3)** — `DHFR` (ENSG00000228716), `SHMT2` (ENSG00000182199), `MTHFD2` (ENSG00000065911). Folate (but not reduced folate) is in all media; thymidine+hypoxanthine rescue of DHFR is available only in F12/DMEM-F12 [media sheets]. [LIT Ducker & Rabinowitz, Cell Metab 2017;25:27, PMID 27641100 *(verify)*].

**E. Cofactor activation and transport (11)** — `COASY` (ENSG00000068120), `PPCS` (ENSG00000127125), `NADK` (ENSG00000008130), `NAMPT` (ENSG00000105835), `TPK1` (ENSG00000196511), `RFK` (ENSG00000135002), `FLAD1` (ENSG00000160688), `PDXK` (ENSG00000160209), `PNPO` (ENSG00000108439), `SLC25A26` (ENSG00000144741, SAM carrier), `SLC5A6` (ENSG00000138074, pantothenate/biotin/lipoate transporter). Media supply the vitamins, not the active cofactors, so the activating kinases/synthetases cannot be bypassed. All except NAMPT and FLAD1 are lethal in both models [PRED]; NAMPT is a model false negative via the tryptophan→quinolinate→NaMN route (1.5f).

**F. Fe-S cluster assembly and lipoylation (6, Human2 only)** — `NFS1` (ENSG00000244005), `ISCU` (ENSG00000136003), `LYRM4` (ENSG00000214113), `NFU1` (ENSG00000169599), `GLRX5` (ENSG00000182512), `HSPA9` (ENSG00000113013). **Absent from Human1** (new in Human2) [MS]; `ratio_r = 1` in Human2 [PRED] → false negatives (1.5b). Pan-essential because every aconitase/respiratory/radical-SAM enzyme needs clusters.

**G. Lipid/membrane core (8)** — `CDS2` (ENSG00000101290), `CDIPT` (ENSG00000103502), `PCYT1A` (ENSG00000161217), `SPTLC1` (ENSG00000090054), `SPTLC2` (ENSG00000100596), `KDSR` (ENSG00000119537), `TECR` (ENSG00000099797), `HSD17B12` (ENSG00000149084). Serum supplies fatty acids and cholesterol but not PI, PC head-group activation or sphingoid bases. CDIPT/TECR/HSD17B12 lethal in both models; CDS2/PCYT1A/SPTLC1/2/KDSR are ~1 [PRED] (1.5c/1.5a).

**H. Mevalonate trunk, isoprenoid branch (6)** — `HMGCS1` (ENSG00000112972), `HMGCR` (ENSG00000113161), `MVK` (ENSG00000110921), `PMVK` (ENSG00000163344), `MVD` (ENSG00000167508), `FDPS`. Rationale and the key distinction from 1.2: serum LDL rescues the **cholesterol** branch but not FPP/GGPP for prenylation and dolichol, so statin-level blockade is rescued by mevalonate, not by LDL [LIT Garcia-Bermudez et al. Nature 2019;567:118, PMID 31527799 (cholesterol auxotrophy via the post-squalene branch)].

**I. Hexosamine / UDP-sugar core (5)** — `GFPT1` (ENSG00000198380), `GNPNAT1` (ENSG00000100522), `PGM3` (ENSG00000013375), `UAP1` (ENSG00000117143), `DPM1` (ENSG00000000419). N-glycosylation and GPI anchors require UDP-GlcNAc and dolichol-P-mannose; no medium supplies them. GNPNAT1/PGM3 lethal in both models [PRED].

**Deliberately excluded from 1.1** (common misclassification): all NDUF* complex I subunits (glucose-grown cells tolerate complex I loss; complex IV/V are the dependency), GLS, PHGDH, FASN, ACACA, SQLE, GCLC, SLC7A11, GPX4, ASNS, ASS1 — all medium- or context-dependent (1.2, 2).

## 1.2 Medium-dependent essentials (rescuing component named)

Medium distribution in `Model.csv` (note the DATA.md quirk: `FormulationID` holds the free text) [OMICS]: RPMI 931, DMEM 298, EMEM/MEM 163, DMEM:F12 104, IMDM 73, F12 57, αMEM 31, McCoy's 31, L-15 23, Waymouth 5, organoid/stem/proprietary 347, other 34, **not recorded 457**; 1,730 formulations mention serum, 824 do not (457 of those are blank, 331 are the serum-free organoid/stem class).

| Gene(s) (model ID) | Rescuing medium component | Which media contain it | Model/structure or literature evidence |
|---|---|---|---|
| `TYMS` (ENSG00000176890), `TK1` | thymidine | F12 and DMEM/F12 only (0.365 mg/L); absent from RPMI, DMEM, IMDM, MEM | exchange `MAR09423`; Human1 reversibility of `MAR04644` makes thymidine rescue UMPS too [MS] |
| `DHFR` (ENSG00000228716) | hypoxanthine **+** thymidine (HAT-type rescue) | F12, DMEM/F12 (hypoxanthine 2.39 mg/L) | `MAR09358`, `MAR09423` [MS]; classic DHFR-null CHO genetics [LIT Urlaub & Chasin PNAS 1980, PMID 6933469 *(verify)*] |
| de novo purine block `PPAT, GART, PFAS, PAICS, ADSL, ATIC` | hypoxanthine → `HPRT1` salvage | F12, DMEM/F12; serum hypoxanthine elsewhere | hypoxanthine is one of the 44 Ham's inputs of task `GR` [MS] |
| `CAD, DHODH, UMPS` | uridine (and uric acid inhibits UMPS) | none of the sheets; serum only | [LIT Cantor 2017, PMID 28388411] |
| `ASNS` (ENSG00000070669) | asparagine | RPMI 0.38 mM, IMDM 0.17 mM, DMEM/F12 0.05 mM; **DMEM high-glucose has none** | exchange `MAR09062` [MS]; asparagine as exchange factor [LIT Krall et al. Nat Commun 2016;7:11457, PMID 27161942] |
| `ASS1, ASL` | arginine (and citrulline bypass) | arginine in all media | `MAR09066`, `MAR09201` [MS]; ASS1-null arginine auxotrophy [LIT Rabinovich et al. Nature 2015;527:379, PMID 26560030 *(verify)*] |
| `GLS` (ENSG00000115419) | glutamine level, plus glutamate/α-KG | IMDM 4 mM, RPMI 2 mM, "no-glutamine" DMEM/F12 formulations 0 mM | `MAR09063`, `MAR09071` [MS]; [LIT Gross et al. Mol Cancer Ther 2014;13:890, PMID 24523301] |
| `PHGDH, PSAT1, PSPH` | serine + glycine | serine 0.25-0.4 mM and glycine in every sheet | `MAR09069`, `MAR09067`; **but** `PSAT1` is lethal in both models for a pool artefact, not for serine (1.4b) [MS/PRED] |
| `GCLC, GCLM, GSS` | reduced glutathione and cystine | **RPMI uniquely contains GSH (1 mg/L)**; cystine in all | `MAR09351`, `MAR09363` [MS]; all three lethal in both models [PRED] |
| `SLC7A11, GPX4` | cystine, plus serum lipids/α-tocopherol (ferroptosis buffering) | cystine in all; tocopherol only from serum | [LIT Dixon et al. Cell 2012;149:1060, PMID 22632970; Viswanathan et al. Nature 2017;547:453, PMID 28678785] |
| `PC, PCK2`, OXPHOS subunits | pyruvate (and uridine) | RPMI 1 mM, IMDM 1 mM, DMEM/F12 0.5 mM; **DMEM high-glucose w/o pyruvate has none** | `MAR09133` [MS]; pyruvate/uridine rescue of respiration-deficient cells (ρ0 genetics) [LIT King & Attardi Science 1989;246:500, PMID 2814477 *(verify)*] |
| cholesterol branch `FDFT1, SQLE, LSS, CYP51A1, TM7SF2, MSMO1, NSDHL, HSD17B7, EBP, SC5D, DHCR7, DHCR24` | LDL/serum cholesterol | serum only → the 331 serum-free lines cannot be rescued | exchange `MAR09285` is closed under Ham's; 11 of these are lethal in Human2 [MS/PRED]; [LIT PMID 31527799] |
| `FASN` (ENSG00000169710), `ACACA` (ENSG00000278540), `SCD` | serum palmitate / oleate | serum; linoleate+linolenate are Ham's inputs; DMEM/F12 adds 0.042 mg/L linoleic acid | `MAR00611`, `MAR00650`, `MAR09035/36` [MS]; FASN/ACACA lethal in both models [PRED] |
| `BBOX1, TMLHE, ALDH9A1, SHMT1`(H2)/`SHMT2`(H1) | serum carnitine | serum only (no sheet lists carnitine) | carnitine is a `cofactor_pool_biomass` component (`MAR10065`, coeff 0.0345) and `MAR06979/06982/06980` are the only route; exchange `MAR09292` closed [MS]; all lethal [PRED] |
| `MMAA` (+`MMAB`, `MMUT`) | vitamin B12 | RPMI 0.005, DMEM/F12 0.68, IMDM 0.013 mg/L; **DMEM high-glucose has no B12** | `MAR08614-08616`; cobamide-coenzyme in `MAR10065` [MS] |
| `ACACA, PC, PCCA/PCCB` (biotin-dependent carboxylases), `SLC5A6` | biotin | RPMI 0.2, IMDM 0.013, DMEM/F12 0.0035 mg/L; **DMEM high-glucose has no biotin** | biotin in `MAR10065`; transport `MAR09609` [MS] |
| `PNPO` (ENSG00000108439) | pyridoxal (vs pyridoxine) | **IMDM supplies pyridoxal**; RPMI/DMEM/F12 supply pyridoxine | `MAR09400` (pyridoxal) vs `MAR09144` (pyridoxine) [MS] → PNPO should be dispensable in IMDM-grown lines but `PDXK` essential in both cases. Sharpest single medium-contrast prediction in this list. |
| `LIAS, LIPT1, LIPT2` | lipoic acid | **DMEM/F12 only (0.105 mg/L)**; lipoate is a Ham's input | `MAR09167`; lipoate in `MAR10065` [MS] |
| `ODC1, SRM, SMS, AMD1` | putrescine | **DMEM/F12 only (0.081 mg/L)** | exchange `MAR11426` [MS] |
| `ISYNA1` | myo-inositol | all sheets (RPMI 35 mg/L) | `MAR09361`; inositol is a Ham's input and a `CDIPT` substrate [MS] |
| `CBS, CTH` (transsulfuration) | cystine/cysteine | cystine in all sheets; DMEM/F12 also free cysteine | `MAR09363`, `MAR09065` [MS] |
| `SEPHS2`, selenoproteins (`GPX4`) | selenium | **IMDM only (selenite 0.0173 mg/L)**; otherwise serum | no selenite exchange in the models — a representability gap [MS] |
| heme pathway `ALAD, HMBS, UROS, UROD, CPOX` | serum hemin / heme | serum | exchange `MAR09107` closed under Ham's; all five lethal in both models [PRED] |
| `NAMPT` vs `NAPRT` | nicotinamide (all media) vs **nicotinic acid (none of the sheets)** | — | `MAR09378` vs `MAR09142` [MS]; NAPRT rescue of NAMPT inhibition needs exogenous NA [LIT Chowdhry et al. Nature 2019;569:570, PMID 31086336 *(verify)*] |

## 1.3 Collateral lethality and paralog buffering

**Collateral lethality (passenger deletion → paralog/pathway dependency)**
1. `MTAP` deletion (9p21, with `CDKN2A`) → `PRMT5` / `MAT2A` / `RIOK1` dependency via MTA accumulation [LIT Mavrakis et al. Science 2016;351:1208, PMID 26912360; Kryukov et al. Science 2016;351:1214, PMID 26912361]. **PRMT5 is absent from both models**; `MAT2A` (ENSG00000168906) is present → only the SAM-flux half is representable.
2. `ENO1` deletion (1p36, glioblastoma) → `ENO2` dependency [LIT Muller et al. Nature 2012;488:337, PMID 22895344]. Model-representable: `MAR04363` GPR `ENO1 or ENO3 or ENO2` [MS].
3. `ME2` deletion (18q, co-deleted with `SMAD4`, pancreas) → `ME3` dependency [LIT Dey et al. Nature 2017;542:119, PMID 28178238]. Representable: ME1/ME2/ME3 appear in separate compartment reactions [MS].
4. `FH` loss → `HMOX1` dependency (haem oxygenase) [LIT Frezza et al. Nature 2011;477:225, PMID 21849978]. `HMOX1` is in both models.
5. `ASS1` silencing → arginine auxotrophy and rerouting of aspartate into pyrimidines [LIT PMID 26560030 *(verify)*].
6. `NAPRT` amplification/`NAPRT` methylation → NAD-route switch and NAMPT-inhibitor resistance/sensitivity [LIT PMID 31086336 *(verify)*; Shames et al. Clin Cancer Res 2013;19:6912, PMID 24097869 *(verify)*].

**Paralog-buffering pairs** (knock-out of one is buffered; a GEM with an `OR` GPR reproduces this automatically, so these are the pairs where the model is *right* and expression-based pruning would break it) [LIT De Kegel & Ryan, PLoS Genet 2019;15:e1008466, PMID 31652272 *(verify)*; Thompson et al. Nat Commun 2021;12:1302, PMID 33637726 *(verify)*]:
`ENO1/ENO2(/ENO3)`, `ALDOA/ALDOC(/ALDOB)`, `LDHA/LDHB`, `HK1/HK2(/HK3/HKDC1)`, `PFKL/PFKM/PFKP`, `PKM/PKLR`, `PRPS1/PRPS2`, `IMPDH1/IMPDH2`, `CTPS1/CTPS2`, `NME1/NME2`, `ADSS1/ADSS2`, `AK1/AK2/AK3/AK4`, `GLS/GLS2`, `GLUD1/GLUD2`, `ME1/ME2/ME3`, `MTHFD1/MTHFD1L/MTHFD2/MTHFD2L`, `SHMT1/SHMT2`, `ACSS1/ACSS2/ACSS3`, `PTDSS1/PTDSS2`, `CDS1/CDS2`, `GYG1/GYG2`, `FADS1/FADS2`, `UGT1A*` cluster, `SULT1A1/SULT1A2`, `AKR1C1/AKR1C2`, `GOT1/GOT2` (compartment-separated, partially redundant through the malate-aspartate shuttle).

## 1.4 Known GEM false-positive classes (model calls essential, cells do not need it)

**(a) Dead-end species in the cofactor/vitamin biomass pool `MAR10065`** (coefficient 0.0012 in `MAR13082`, yet it forces 36 species including retinoid, vitamin D and tocopherol catabolites) [MS]. Forced genes, all lethal [PRED]: `CYP26A1` (ENSG00000095596), `CYP26B1` (ENSG00000003137), `CYP24A1` (ENSG00000019186), `CYP1A2` (ENSG00000140505), `CYP4F2` (H2 only), `KL`, `TMLHE`, `BBOX1`, `ALDH9A1`, `SHMT1` (H2) / `SHMT2` (H1), `MMAA`, `GCH1`, `SPR`, `PTS` (H1), `SLC5A6`, `COQ2/3/6/7`. Example chain: `MAR10065` needs `3-carboxy-alpha-chromanol`, produced only by `MAR06436`/`MAR06435` whose GPRs are AND-chains of β-oxidation enzymes — this is why `ACADL, ACADSB, ECHS1, HADHA, HADHB, ACAA2, AMACR` are lethal in Human2 and not in Human1 (Human1 lacks `MAR06433/06435`) [MS/PRED]. These enzymes are not required for tocopherol catabolism in a cultured cancer cell.

**(b) Dead-end species in the small-metabolite pool `MAR10064`** [MS]: `PSAT1` (ENSG00000135069) is the **only** producer of 3-phosphoserine (`MAR03841`), which `MAR10064` demands at 0.0025 — so PSAT1 is lethal irrespective of medium serine; `UGDH` (sole UDP-glucuronate producer `MAR04122`); `PGM3`/`GNPNAT1` (UDP-GlcNAc, coeff 0.052); `CAD` via dihydroorotate.

**(c) Glycogen in biomass (coefficient 0.4062, the largest non-protein term)** [MS] forces `UGP2` (lethal in Human2) and the glycogenin route `GYG1/GYG2`, `GYS1`. Cultured cancer cells do not accumulate glycogen at that rate.

**(d) Pseudo-reaction GPRs (true CRISPR essentials, but for non-metabolic reasons → "easy positives")** [MS]: `MAR13086` (RNA nuclear transport) is a 31-gene AND of nucleoporins → `NUP35/37/42/43/50/54/58/62/85/88/93/98/107/133/153/155/160/188/205/210/214, TPR, NDC1, SEH1L, RAE1, GLE1, SEC13, RANBP2, AHCTF1, POM121, DDX19B` (+`AAAS` in Human1); `MAR07161` (RNA polymerisation) is a 30-gene AND → all `POLR1*/POLR2*/POLR3*`, `POLRMT`, `MTPAP`; `MAR06978` → `SPCS1 and SPCS2 and SPCS3 and SEC11A and SEC11C`; `MAR10062` → the 20 ARSs. 91 of Human2's 154 lethal genes are in these four GPRs. **Recommendation: report AUPRC with and without this stratum.**

**(e) Medium misspecification** (uptake closed although culture supplies it): cholesterol branch, heme pathway, `FASN`/`ACACA`, `GCLC/GCLM/GSS` (RPMI GSH), carnitine synthesis, `COQ*` — all with open exchanges in the released model that Ham's closes [MS].

**(f) Single-gene annotation of a transport step** (one wrong annotation = one lethal gene): Human1 `NPC1L1` (ENSG00000015520, `MAR00018` sole cholesterol transporter), `SCARB1` (ENSG00000073060, `MAR00019`), `ABCB6`, and in both `SLC25A26`, `CDIPT`, `SLC5A6` [MS/PRED].

**(g) Over-specified AND complexes**: `SEC11A and SEC11C` (paralogous catalytic subunits, not both obligate); `MAR06435`'s `AMACR and ACADSB and ECHS1 and (HADH or HSD17B10) and ACAA2`; the 31-gene nucleoporin AND [MS].

**(h) Unbounded-uptake artefacts** (uptake bound 1000): the 0.50-0.88 partial-growth genes `TYMP, SUCLG1, DERA, PGM2, SLC25A3, SLC16A1` and all OXPHOS/ATP-synthase subunits, where thymidine is consumed as a bulk carbon source [PRED, baseline README]. Treat as medium-dependent, not as predictions.

## 1.5 Known GEM false-negative classes

**(a) Structural non-predictions**: 1,801/2,887 (Human1) and 1,750/2,870 (Human2) genes close **no** reaction when knocked out, so `ratio_r = 1` by construction — ~61% of the gene space carries no information [PRED]. Any gene whose only reactions are blocked, or that sits in an `OR` with an unexpressed paralog, is a structural false negative.

**(b) Genes missing from the model**: `PRMT5`, `OGDHL` (absent from both); the Fe-S/lipoylation module (`NFS1, ISCU, FDX1/2, NFU1, GLRX5, LYRM4, HSPA9, IBA57, ISCA1/2, HSCB`) exists **only in Human2**, 59 genes are Human1-only and 42 Human2-only [MS]. No mechanism → guaranteed miss.

**(c) Isozyme `OR` masking the only expressed paralog**: `ENO1, ALDOA, GAPDH, PGK1, TPI1, PKM, GPI, LDHA, HK1/HK2, PFKP` are all `ratio_r = 1` [PRED] although glycolysis is a real dependency, because each sits in an `OR` with a tissue-restricted paralog (`ENO3` muscle, `ALDOB` liver, `PKLR` liver/erythroid, `PGAM4` retrogene) [MS]. This is the single largest fixable false-negative class — and the reason section 3 exists.

**(d) Compartment-blind redundancy**: `GOT1/GOT2`, `MDH1/MDH2`, `IDH1/IDH2`, `SHMT1/SHMT2`, `ACSS1/ACSS2` share reactions or mirrored reactions across compartments [MS], so the malate-aspartate shuttle and cytosolic NADPH supply are never essential.

**(e) Demands absent from biomass**: polyamines (`ODC1, SRM, AMD1`), lipid-peroxide defence (`GPX4, SLC7A11, ACSL4`), selenoprotein synthesis (`SEPHS2`), nucleotide-pool turnover, dolichol/GPI load beyond `DPM1`, mitochondrial translation, redox maintenance under oxidative stress. None of these can be essential in FBA because no flux must carry them.

**(f) Alternative de novo routes that culture cells do not use**: `NAMPT` is non-lethal because NAD can be made from tryptophan via quinolinate/`QPRT` [MS]; `PPAT/GART/ATIC` are non-lethal because the purine pool is also reachable from hypoxanthine (a Ham's input) [MS]; cells in these media rely on salvage, so the model's "either route" logic erases the dependency.

**(g) No non-growth ATP maintenance**: `MAR03964` (ATP + H2O → ADP + Pi) has bounds [0, 1000], i.e. no enforced maintenance flux [MS]. With ATPM = 0 the model can satisfy biomass with minimal ATP turnover, which is why OXPHOS and glycolysis genes come out at 0.57-1.00 rather than lethal. A pre-registered ATPM variant is the obvious fix.

## 1.6 CRISPR-artefact strata (essentiality unreliable for non-metabolic reasons)

**(a) Mitochondrially encoded genes (13 in each model)** — `MT-ND1 (ENSG00000198888), MT-ND2, MT-ND3, MT-ND4, MT-ND4L, MT-ND5, MT-ND6, MT-CYB (ENSG00000198727), MT-CO1 (ENSG00000198804), MT-CO2, MT-CO3, MT-ATP6 (ENSG00000198899), MT-ATP8`. Cas9/sgRNA complexes are not imported into mitochondria, so mtDNA is not edited and the measured effect is ~0 regardless of biology [LIT Gammage et al. Trends Genet 2018;34:101, PMID 29179920 *(verify)*]. The models nevertheless predict partial lethality (0.57-0.88) for 6 of them [PRED] → guaranteed false positives. **Exclude from scoring; report separately.**

**(b) High-identity paralogs (sgRNA cross-targeting + functional buffering)** — 67 model genes have a ≥95%-identity within-species paralog, 111 at ≥90% (Ensembl 116) [OMICS]. Families to stratify: `CALM1/CALM2/CALM3`, `UBE2D2/UBE2D3`, `AMY1A/AMY1B/AMY1C/AMY2A/AMY2B`, `USP17L*` (20+ members), `PPIAL4*`, `SULT1A1/SULT1A2`, `SULT1A3/SULT1A4`, `GSTT2/GSTT2B`, `NUDT4/NUDT4B`, `NUDT10/NUDT11`, `AKR1C1/AKR1C2`, `ALG10/ALG10B`, `ACSM2A/ACSM2B`, `CKMT1A/CKMT1B`, `UQCRH/UQCRHL`, `DDX19A/DDX19B`, `PGAM1/PGAM4`, `GK/GK3/GK3P`, `RPE/RPEL1`, `PRPS1/PRPS2`, `GLUD1/GLUD2`, `GSTA1/GSTA2`, `ALPG/ALPP`, `UBE2A/UBE2B`, `FKBP1A/FKBP1C`, `TPSAB1/TPSB2`, `CDY1/CDY1B/CDY2A/CDY2B`.

**(c) Amplicon / copy-number effect** [LIT Meyers et al. Nat Genet 2017;49:1779, PMID 29083409]. Model genes most often at relative CN ≥ 2 (fraction of 1,763 WGS default-entry lines) [OMICS]: `PDPR` 33.4%, `CLPS` 22.9%, `SQLE` 11.3%, `NDUFB9` 10.8%, `RNF139` 10.3%, `SLC39A4` 9.6%, `GPT` 9.5%, `NTAQ1` 9.5%, `ADCY8` 9.3%, `ATP5F1E` 9.1%, `SLC6A19` 9.0%, `DGAT1` 9.0%, `PCK1` 8.8%, `PYCR3` 8.7%, `GPAA1` 8.7%, `EXT1` 8.6%, `GFUS` 8.6%, `HAS2` 8.4%, `RAE1` 8.4%, `CYP24A1` 8.3%, `SDHA` 8.3% (chiefly 8q22-8q24 and 17q25 amplicons). Chronos corrects CN, but residual bias remains, and three of these (`SQLE`, `RAE1`, `CYP24A1`) are model-lethal [PRED] — i.e. a correlated artefact in exactly the genes the model calls positive.

**(d) Recurrently deleted / unscreenable genes** (effect undefined or all-zero where the gene is absent): `USP9Y` 72.4% of lines at CN < 0.3 (Y loss / female lines), `ACOT1` 40.2%, `UGT2B17` 32.6% (common germline CNV), `MTAP` 19.1%, `HACD4` 4.7%, `STS` 3.4%, `GSTM1` 1.1% [OMICS].

**(e) Gene-ID hazards from the baseline mapping** [PRED]: `KLHL23` and `USP17L2` have no CRISPR column; `GK3P` (pseudogene) and `SLC37A4` as alt-locus `ENSG00000281500` have no Entrez and no column; 11 model genes carry several Entrez IDs; 7-10 symbols differ between the model and DepMap at equal Entrez. Match on Entrez only, and exclude these from scoring.

**(f) Moonlighting / structural roles** (CRISPR effect not attributable to flux, so even a perfect GEM must miss or mis-attribute it): `ENO1` (plasminogen receptor, nuclear MBP-1), `GAPDH` (mRNA binding, apoptosis), `PKM` (nuclear PKM2 kinase), `ALDOA` (actin binding, V-ATPase), `FH` (nuclear DNA-repair signalling), `ACO1/IREB2` (iron-responsive element binding), `GPI` (neuroleukin/cytokine), `HK2` (VDAC/apoptosis), `IDH1`, `ACLY` (histone acetylation), `SHMT2` (BRISC deubiquitinase complex).

**(g) Non-metabolic machinery inside metabolic pseudo-reactions**: the 91 genes of 1.4d (nucleoporins, RNA polymerases, signal peptidase, ARSs). Essentiality is real but its driver is transcription/translation, so including them in the metric measures pathway annotation, not metabolic modelling.

---

# 2. Genomic-background groups with a GEM-representable metabolic mechanism

Denominators (default entry per model) [OMICS]: WGS copy number **1,763** lines; expression **2,130**; damaging-mutation matrix **2,395**; all three available **1,724**; at least one **2,476**. Thresholds fixed a priori: deep deletion = relative CN < 0.3; amplification = relative CN ≥ 2.0; expression-null = log2(TPM+1) < 1.0; "low" and "high" thresholds are stated per group. Counts are over the lines where the needed data type exists.

| # | Group | Definition (data source) | n lines | Top lineages | Representable mechanism |
|---|---|---|---|---|---|
| 1 | MTAP loss (CN) | `OmicsCNGeneWGS` MTAP < 0.3 | **336** | CNS 72, Lung 32, Lymphoid 31, Pancreas 26 | MTA/SAM flux; `MTAP` reactions removed. PRMT5 not in model → partial |
| 2 | MTAP loss (expression) | TPM log2 < 1 | **341** | CNS 69, Lymphoid 41, Lung 32 | as above |
| 3 | MTAP loss (union) | 1 ∪ 2 | **400** | CNS 79, Lymphoid 53, Lung 39 | primary MTAP group |
| 4 | **NEG-1** CDKN2A-del, MTAP intact | CDKN2A < 0.3 **and** MTAP CN ≥ 0.5 **and** MTAP TPM ≥ 2 | **157** | CNS 19, Skin 17, Pancreas 14 | specificity control for group 3 (same 9p21 event, no MTAP mechanism) |
| 5 | ENO1 hemizygous-low | ENO1 CN < 0.6 **and** TPM < 9 | **39** | PNS 6, CNS 6, Breast 5 | ENO1→ENO2 collateral lethality. **Deep ENO1 deletion: 0 lines** — the canonical test is unavailable in 26Q3 |
| 6 | ME2 deep deletion | ME2 CN < 0.3 | **27** | Pancreas 8, Eso/Stomach 5 | ME2→ME3 collateral lethality |
| 7 | ME2 expression-null | ME2 TPM log2 < 1 | **26** | Pancreas 10, Eso/Stomach 5 | as above |
| 8 | **NEG-2** SMAD4-del, ME2 intact | SMAD4 < 0.3, ME2 CN ≥ 0.5, ME2 TPM ≥ 2 | **14** | Pancreas 7 | specificity control for 6/7; **underpowered, report as descriptive only** |
| 9 | PHGDH amplification | PHGDH CN ≥ 2 | **42** | Breast 8, Bone 6, Eso/Stomach 6 | serine-synthesis flux; PHGDH/PSAT1/PSPH |
| 10 | PHGDH amp or high | CN ≥ 2 **or** TPM ≥ 9 | **115** | Breast 13, Bone 13, Skin 11 | as above, better powered |
| 11 | ASS1 expression-null | ASS1 TPM log2 < 1 | **254** | Lymphoid 66, Skin 55, Myeloid 25 | urea-cycle break → arginine auxotrophy, aspartate rerouting |
| 12 | ASNS low | ASNS TPM log2 < 4 | **81** | Lymphoid 15, Skin 11, Myeloid 10 | asparagine auxotrophy; interacts with medium asparagine (1.2) |
| 13 | KEAP1/NFE2L2 axis | KEAP1 damaging **or** KEAP1 CN < 0.3 **or** NFE2L2 hotspot | **68** | Lung 24, Eso/Stomach 9, H&N 7 | NRF2-driven glutamine/glutathione flux [LIT Romero et al. Nat Med 2017;23:1362, PMID 29035354 *(verify)*] |
| 14 | IDH1/IDH2 hotspot (matrix) | `OmicsSomaticMutationsMatrixHotspot` > 0 | **15** | CNS 5, Biliary 3, Bone 2 | 2-HG production; `IDH1`/`IDH2` reaction reversal [LIT Dang et al. Nature 2009;462:739, PMID 19935646] |
| 15 | IDH1 R132 / IDH2 R140,R172 (long format) | `OmicsSomaticMutations` ProteinChange | **18** | CNS 5, Bone 5, Biliary 3 | preferred definition (3 lines missed by the hotspot matrix) |
| 16 | STK11 loss | damaging **or** CN < 0.3 | **82** | Lung 41, Cervix 9 | CPS1-dependent pyrimidine supply [LIT Kim et al. Nature 2017;546:168, PMID 28538732 *(verify)*] |
| 17 | STK11 loss + KRAS hotspot + Lung | intersection | **21** | Lung 21 | the published context for 16 |
| 18 | STK11 loss + CPS1 expressed | STK11 loss and CPS1 TPM ≥ 3 | **35** | Lung 19 | mechanism-matched subgroup |
| 19 | VHL loss | damaging **or** CN < 0.3 | **38** | Kidney 30 | HIF-driven glycolysis / reductive carboxylation [LIT Metallo et al. Nature 2011;481:380, PMID 22101433] |
| 20 | VHL loss, Kidney | intersection | **30** | Kidney 30 | lineage-clean version of 19 |
| 21 | SDHx/FH loss (permissive) | any of SDHA-D, FH damaging or CN < 0.3 | **67** | Lymphoid 11, Lung 6 | TCA break, HMOX1 synthetic lethality |
| 22 | SDHx/FH **biallelic-like** | deep deletion, or damaging + CN < 0.7, or ≥2 damaging alleles, in non-hypermutant lines | **3** | — | **Not testable.** Group 21 is almost entirely monoallelic/passenger: hard judgement call recorded here |
| 23 | NAPRT expression-null | NAPRT TPM log2 < 1 | **59** | Bone 12, PNS 10, CNS 10 | NAD route switch; NAMPT dependency |
| 24 | SQLE expression-null | SQLE TPM log2 < 1 | **13** | Lymphoid 7, Uterus 3 | cholesterol auxotrophy [LIT PMID 31527799]; small but mechanistically crisp |
| 25 | HK1-low / HK2-high | HK1 TPM < 2 and HK2 TPM ≥ 4 | **38** | Lymphoid 17, Liver 6 | isozyme switch at the glycolysis entry `MAR04319` |
| 26 | LDHB expression-null | LDHB TPM log2 < 1 | **58** | Eso/Stomach 14, Pancreas 11 | `MAR04281` GPR `LDHB or LDHA` → LDHA dependency |
| 27 | SLC7A11 expression-null | TPM log2 < 1 | **90** | Lymphoid 13, Soft tissue 11 | cystine import / GSH supply |
| 28 | CDS1 expression-null | TPM log2 < 1 | **808** | Lymphoid 155, CNS 113, Skin 77 | `CDS1/CDS2` paralog pair → CDS2 dependency; largest well-powered group |
| 29 | **NEG-3** TTN damaging, non-hypermutant | TTN damaging and burden ≤ 90th percentile (≤ 88 damaging genes) | **162** | Lung 32, Skin 27 | **negative control: no metabolic mechanism at all** |
| 30 | **NEG-4** SMARCA4 loss | damaging or CN < 0.3 | **97** | Lung 38, Ovary 9 | second negative control: strong, frequent driver event with no GEM-representable enzyme change |

Groups with **0-1 lines** in 26Q3, recorded so they are not attempted: deep deletion of `ENO1`, `PGD`, `UMPS`, `PTDSS2` (0 each); `HPRT1` expression-null (1); `FH` deep loss (1).

---

# 3. Isozyme exceptions: `OR` relations that must never be edited automatically

These GPRs encode genuine functional redundancy in cancer cells. Expression- or proteomics-based pruning (dropping the "unexpressed" isozyme and so collapsing `OR` to a single gene, or converting `OR` to `AND`) would manufacture false positives in exactly the genes that CRISPR shows are buffered. Any curation pipeline must carry an explicit allow-list.

| Reaction(s) | GPR as released | Why it must stay an `OR` |
|---|---|---|
| `MAR04363` (2-PG ⇄ PEP) | `ENO1 or ENO3 or ENO2` | ENO2 rescues ENO1 loss; this is the textbook collateral-lethality pair [LIT PMID 22895344] |
| `MAR04375`, `MAR04356`, `MAR08762` | `ALDOC or ALDOB or ALDOA` | ALDOA/ALDOC co-expression buffers either single loss |
| `MAR04281`, `MAR04280`, `MAR01453` | `LDHB or LDHA`, `LDHB or LDHAL6B`, `LDHA or LDHB` | LDHA/LDHB double-knock-out is needed to block lactate flux |
| `MAR04319` and the glucose hexokinase reactions | `HKDC1 or HK1 or HK2 or HK3` | HK1/HK2 buffer each other in most lineages (see group 25) |
| `MAR04379`, `MAR04301`, `MAR04567` | `PFKP or PFKL or PFKM` | PFK isoforms hetero-oligomerise; single losses are tolerated |
| `MAR04171`, `MAR04193`, `MAR04210` | `PKM or PKLR` | PKM2 loss is buffered by PKM1 splicing (one gene) and PKLR where expressed; pruning makes PKM falsely essential in every line |
| `MAR03829`, `MAR03912`, `MAR06725` | `GOT1 or GOT1L1`, `GOT1 or GOT2`, `GOT1 or GOT2 or TAT` | malate-aspartate shuttle redundancy; GOT1 dependency is conditional, not constitutive |
| `MAR04340` | `MTHFD2 or MTHFD1L or MTHFD2L` | mitochondrial 1C redundancy |
| purine/pyrimidine terminal steps | `IMPDH1 or IMPDH2`, `CTPS1 or CTPS2`, `NME1 or NME2`, `ADSS1 or ADSS2` | paralog buffering, confirmed by the need for combinatorial screens [LIT PMID 33637726 *(verify)*] |
| `MAR00153`, `MAR03796`, `MAR03797` | `ACSS3 or ACSS2 or ACSS1` | acetate activation is redundant across compartments |
| PRPP synthesis | `PRPS1 or PRPS2` | 97% identical paralogs; also an sgRNA cross-targeting pair (1.6b) |
| glutamate dehydrogenase | `GLUD1 or GLUD2` | ≥95% identity, both expressed in many lines |
| malic enzymes | `ME1`/`ME2`/`ME3` separate reactions | do **not** merge into one AND; ME2→ME3 buffering is the published mechanism [LIT PMID 28178238] |
| `MAR01919` | `SULT2B1 or SULT1A3` | SULT1A3/1A4 are a duplicated locus |
| `MAR04371` | `PGAM2 or PGAM1 or BPGM or PGAM4` | PGAM4 is a retrogene, but PGAM1/BPGM redundancy is real |

**Not exceptions — these are the curation targets** (the `OR` is annotation error, and they are already the models' top false positives/negatives): `TECR or TECRL` (`MAR02194`; TECRL is heart/muscle-restricted), `PTGR1 or HSD17B12 or TECRL` (`MAR01074`), the 27-reaction `ALDH3A2 or ALDH2 or ALDH1B1 or ALDH9A1 or ALDH7A1` group, the large CYP `OR` blocks in retinol/arachidonate metabolism, `SEC11A and SEC11C` (paralogs wrongly ANDed), and the AND-chains `MAR06435`/`MAR06436`.

---

# 4. Metabolic drug classes, compounds and positive-control pairs

Names verified against `data/prism/primary-screen-replicate-collapsed-treatment-info.csv`, `data/prism/repurposing_24Q2/..._Extended_Primary_Compound_List.csv` and `data/gdsc/screened_compounds_rel_8.5.csv` [OMICS]. "P" = present in PRISM, "G(id)" = GDSC drug id.

| Class | Compounds (screen) | Target gene(s) in both models |
|---|---|---|
| Antifolates | METHOTREXATE (P, G1008), PEMETREXED (P, G428), PRALATREXATE (P), RALTITREXED (P), PYRIMETHAMINE (P), TRIMETHOPRIM (P) | `DHFR` ENSG00000228716, `TYMS` ENSG00000176890, `GART`, `ATIC` |
| Thymidylate / pyrimidine antimetabolites | 5-FLUOROURACIL (P, G179/G1073), FLOXURIDINE (P), CAPECITABINE (P), TRIFLURIDINE + TIPIRACIL (P), GEMCITABINE (P, G135/1190/1393), CYTARABINE (P, G1006), DECITABINE, AZACITIDINE (P) | `TYMS`, `TYMP` ENSG00000025708, `RRM1`, `RRM2`, `DCK`, `UCK2`, `UPP1` |
| Purine antimetabolites | MERCAPTOPURINE (P), CLADRIBINE / CLOFARABINE / FLUDARABINE (P, G1813), MYCOPHENOLIC-ACID (P, G2055), AZASERINE (P), ALLOPURINOL (P) | `HPRT1`, `DCK`, `IMPDH1/IMPDH2`, `PPAT`, `XDH` |
| DHODH | BREQUINAR (P), TERIFLUNOMIDE / LEFLUNOMIDE (P), ATOVAQUONE (P; annotated DHODH but acts on complex III) | `DHODH` ENSG00000102967 (and `MT-CYB` for atovaquone) |
| Ribonucleotide reductase | HYDROXYUREA (P) | `RRM1`, `RRM2`, `RRM2B` |
| Glutaminase / amidotransferase | CB-839 / telaglenastat (P), BPTES (G1425), ACIVICIN (P) | `GLS` ENSG00000115419, `CTPS1` |
| Lactate / pyruvate transport | AZD3965 (P), UK-5099 (P; PRISM annotates SLC16A1, true target MPC1/MPC2) | `SLC16A1` ENSG00000155380, `MPC1/MPC2` |
| Glycolysis | 2-DEOXYGLUCOSE (P), LONIDAMINE (P), 3-BROMOPYRUVATE (P), TT-232/TLN-232 (G436), SHIKONIN (G170) | `HK1/HK2`, `SLC2A1`, `PKM`, `GAPDH` |
| OXPHOS / TCA | OLIGOMYCIN-A (P), ATOVAQUONE (P), PHENFORMIN (P, G196), METFORMIN (P), CPI-613/devimistat (P, G415), DICHLOROACETATE (P) | `ATP5F1A` and subunits, `MT-CYB`, complex I `NDUF*`, `PDHA1`, `DLD`, `DLST` |
| Fatty-acid synthesis / oxidation | GSK2194069 (P), ORLISTAT (P), CERULENIN (P), BAY ACCi (G1362), ETOMOXIR (P), PERHEXILINE (P) | `FASN` ENSG00000169710, `ACACA` ENSG00000278540, `CPT1A` |
| Cholesterol / isoprenoid | SIMVASTATIN, LOVASTATIN, ATORVASTATIN (P), TERBINAFINE (P) | `HMGCR` ENSG00000113161, `SQLE` ENSG00000104549 |
| Redox / ferroptosis | BUTHIONINE-SULFOXIMINE (P), SULFASALAZINE (P), ERASTIN (P; PRISM annotates VDAC2), DISULFIRAM (P), ELESCLOMOL (P), N-acetyl cysteine (G2499), alpha-lipoic acid (G2498) | `GCLC`, `GCLM`, `GSS`, `SLC7A11`, `ALDH2`, `GPX4` |
| NAD | DAPORINAD/FK866 (G1248), STF-118804 (P) | `NAMPT` (modifier: `NAPRT`, group 23) |
| SAM / methylation | EPZ015666 (P), GSK3326595 (P) | **`PRMT5` is in neither model** — representability gap; use `MAT2A`, `MTAP` |
| Mutant IDH | IVOSIDENIB, ENASIDENIB, VORASIDENIB (P), AGI-5198 (G1913), AGI-6780 (G562/G1634) | `IDH1`, `IDH2` (groups 14/15) |
| Mitochondrial translation (off-model) | TIGECYCLINE, DOXYCYCLINE, CHLORAMPHENICOL, LINEZOLID (P) | no model gene — use as "metabolic-looking but unrepresentable" controls |

**Positive-control drug–gene pairs** (sensitivity should track the model's predicted essentiality of the named gene; fix these before looking at any correlation): methotrexate↔`DHFR`; pemetrexed & 5-FU & floxuridine↔`TYMS`; brequinar & teriflunomide↔`DHODH`; mycophenolic acid↔`IMPDH2`; hydroxyurea & gemcitabine↔`RRM1`/`RRM2`; CB-839 & BPTES↔`GLS`; AZD3965↔`SLC16A1`; statins↔`HMGCR`; terbinafine↔`SQLE`; GSK2194069 & orlistat & cerulenin↔`FASN`; BAY ACCi↔`ACACA`; etomoxir↔`CPT1A`; FK866/daporinad↔`NAMPT` (with `NAPRT` expression as the modifier); BSO↔`GCLC`/`GCLM`; sulfasalazine & erastin↔`SLC7A11`; oligomycin-A↔`ATP5F1A`; CPI-613↔`DLD`/`DLST`/`PDHA1`; 2-DG↔`HK1`/`HK2`; AGI-5198 & ivosidenib↔`IDH1` (restricted to groups 14/15).

**Negative-control drug pairs** (no model target; expect no association with any predicted dependency): BORTEZOMIB, BELINOSTAT, ALISERTIB, MITOXANTRONE, SGC-CBP30, GSK2606414 (G1618).

**Compounds searched for and absent from the local screens** (do not plan analyses on them): IACS-010759, ML210, RSL3, FIN56, TVB-2640, ND-646, GSK2194069 in GDSC, AG-270, olutasidenib, ADI-PEG20, L-asparaginase/pegaspargase, BAY-876, syrosingopine, antimycin A, KPT-9274.

---

# 5. Candidate metabolic task sets for objectives and constraints

Chosen from the released task files **by category only** (`metabolicTasks_Essential.txt`: 57 tasks = ER 9, BS 18, SU 17, IC 12, GR 1; `metabolicTasks_Full.txt`: 256 tasks; `metabolicTasks_CellfieConsensus.txt`: 195; `metabolicTasks_VerifyModel.txt`: 21 leak tests). Fixed a priori; any selection among sets must be made on held-out lineages, never against the benchmark.

| Set | Content | Size | Rationale |
|---|---|---|---|
| **T0** baseline | `GR` "Growth on Ham's media (biomass production)" as the sole objective | 1 | current baseline; reference point for every comparison |
| **T1** essential-all | all 57 Essential tasks as feasibility constraints, `GR` as objective | 57 | tests whether enforcing viability (not just growth) changes KO calls; the task set the Human-GEM authors curated for exactly this purpose [LIT Robinson et al. Sci Signal 2020;13:eaaz1482, PMID 32209698] |
| **T2** energy/redox | Essential ER tasks (aerobic rephosphorylation of ATP/GTP/CTP/UTP, oxidative phosphorylation, oxidative decarboxylation, Krebs-cycle NADH, ubiquinol-to-proton, ubiquinol-to-ATP) | 9 | the models' OXPHOS genes sit at 0.57-0.88 because there is no ATP maintenance (1.5g); forcing ER flux is the principled fix |
| **T3** nucleotide demand | Essential BS de novo ATP, CTP, GTP, UTP, dATP, dCTP, dGTP, dTTP | 8 | nucleotide synthesis is the dominant class of true proliferative metabolic dependencies (1.1B/1.1C) |
| **T4** central precursors | the 12 Essential IC tasks (G3P, mito acetyl-CoA, mito AKG, E4P, F6P, GAP, G6P, mito OAA, PEP, pyruvate, R5P, mito succinyl-CoA) | 12 | makes central-carbon genes load-bearing without touching biomass coefficients |
| **T5** membrane/lipid | Essential BS PC/PE/PS/PI + cholesterol de novo + SU choline/inositol uptake, optionally Full 104-111 (PC, PE, PS, PI, CL, SM, ceramide, LacCer at physiological substrates) | 7 (+8) | membranes are the largest non-protein biomass sink (`MAR10063`), and serum lipid supply is the main medium variable (1.2) |
| **T6** cofactors | Essential BS thiamin→TPP, CoA from pantothenate, FAD from riboflavin, heme + Full 112-118 (CoA, NAD, NADP, FAD, thioredoxin, THF, pyridoxal-P) | 4 (+7) | replaces the dead-end cofactor pool `MAR10065` (1.4a) with explicit, defensible cofactor demands |
| **T7** CellFie consensus | all 195 CellFie tasks, scored as a task-activity vector (not hard constraints) | 195 | context-specific model building and flux-state features; the standard scored task set [LIT Richelle et al. Cell Rep Methods 2021;1:100040, PMID 35475159 *(verify)*] |
| **T8** amino-acid handling | Full 42-50 (NEAA de novo), 52-68 (EAA uptake), 77-99 (degradation) | ~49 | amino-acid composition is where RPMI, DMEM, F12 and IMDM differ most (1.2) |
| **T9** should-fail gate | every task with `SHOULD FAIL = TRUE` (Full 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 72, 74, 158, 170, 257) **+** all 21 VerifyModel leak tests | 37 | negative control and model-quality gate: a context-specific model that passes these is broken. Never an objective |
| **T10** de-artefacted biomass (pre-registered variant) | `MAR13082` with the `cofactor_pool_biomass` term and the glycogen term set to 0, plus T3 as constraints | 1 objective + 8 | removes the two documented false-positive drivers (1.4a, 1.4c) while keeping proliferative demand. **Declared here as a variant precisely because it edits biomass; it must not be tuned after seeing labels** |

Hard choices recorded for the task sets: (i) T10 edits the biomass reaction, which is a bigger intervention than adding constraints — it is listed as one fixed variant, not a family; (ii) T7 is scored, not enforced, because enforcing 195 tasks in a context-specific model is likely infeasible; (iii) the Essential `GR` task is identical in Human1 and Human2 and all 44 inputs map to the same exchange IDs, so T0-T6 are directly comparable across the two releases [MS, baseline README].
