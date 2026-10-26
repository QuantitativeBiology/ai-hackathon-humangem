# Design review r1: cancer metabolism / CRISPR-screen domain (biological realism)

Reviewer angle: medium and serum modelling, control genes, biomass and tasks, CRISPR artefacts, genomic backgrounds, "metabolic drugs".
Read: PROMPT.md, DESIGN.md v0.1, LITERATURE §4–7, results/baseline/README.md, the GR task input list in `data/Human-GEM-main/data/metabolicTasks/metabolicTasks_Essential.txt`, and the baseline `ratio_r` of about 50 named genes in `results/baseline/human{1,2}_ko.parquet`. No labels were opened and no performance was computed. The gene lists below are my expectations from the DepMap/medium literature and are to be written into the design as pre-registered, not tuned.

## Verdict
The statistical architecture is sound. The biology of W2 (medium/serum) is the weakest part. As written, the ladder will probably erase the best-predicted pathways, nucleotide and mevalonate, rather than repair medium error. Three things drive this: the units of the bounds, the nucleoside/hypoxanthine content of the serum model, and the lack of a transporter gate. There are also fixable gaps in the controls, the biomass additions, the background definitions and the drug list.

## What the baseline already shows (label-free, from the parquet files)

M0 = the 44 GR-task inputs of Human-GEM, each with uptake 1000. Notable inputs: **hypoxanthine, thymidine, aspartate, glutamate, asparagine, proline, alanine, glycine, cysteine (not cystine), linoleate, linolenate, alpha/gamma-tocopherol, retinoate**. Absent: cholesterol, palmitate, oleate, uridine, pyruvate, lactate, NH4, cystine.

Predicted `ratio_r` for named genes (H1 / H2):
- **De novo purine genes** PPAT, GART, PFAS, PAICS, ADSL, ATIC and IMPDH2 all give ratio 1 in both models. Hypoxanthine at 1000 bypasses the whole pathway. GMPS is lethal. In DepMap RPMI/DMEM + 10% FBS these genes are essential in most lines, so M0 produces a systematic false-negative block. This is not mentioned in the design or README.
- **Thymidine at 1000 and TYMS.** TYMS = 1 in both models.
  - DHFR = 1 in both models (cause not diagnosed; likely a folate-route or GPR artefact, worth a pFBA trace).
  - DHODH and UMPS = 1 in H1 and 0 in H2. This is the known TYMS reversibility bypass.
- **Sterol branch.** SQLE, LSS and DHCR7 are lethal at M0, because there is no cholesterol. This is the serum-free state.
- **Upstream mevalonate.** HMGCR and MVK are lethal. HMGCS1 and FDPS are 1, probably isozyme or GPR artefacts.
- **Other lipid genes.** FASN and ACACA are lethal, because palmitate/oleate are absent. SCD = 1 (SCD5 OR).
- **Aspartate/asparagine/serine.** GOT1, GOT2 and ASNS = 1, with Asp and Asn at 1000. PHGDH = 1, PSAT1 = 0 (an inconsistent GPR or serine-route artefact; check).
- **Other genes.** NAMPT = 1 (nicotinamide at 1000 plus NAPRT/NMNAT redundancy; check). SLC7A11 = 1 (cysteine, not cystine, is in the medium). NDUFA1 and COX4I1 = 1. UQCRC1 is about 0.77. H2 predicts the β-oxidation genes ACADL, HADHA and ECHS1 lethal; H1 does not.

So the generic baseline already carries a large medium-driven error. The attribution in W2 is meaningful, but it is dominated by nucleotide/lipid availability.

## MUST-FIX

**M1. Serum and supplement bounds are in the wrong units and the serum model will break the best-known true essentials (§W2 M2a/M2b, S2, "serum uptake bound = 10").**
- Reason 1: bound 10 is not 1% of anything meaningful.
  - Biomass is 1 g/gDW per unit flux, and Human-GEM demands for individual nucleotides, cofactors and lipids are 1e-4 to 1e-1 mmol/gDW.
  - A bound of 10 mmol/gDW/h for hypoxanthine, uridine, cholesterol or oleate is therefore effectively unlimited. It is about 100–1000 times what is needed.
  - At the same time, a glucose bound of 1000 is about 100 times a real cancer-cell glucose uptake (about 5–20 mmol/gDW/h from NCI-60 CORE / Jain 2012).
  - WT growth is 65 a.u. (§2.1, README), against a real value of about 0.02–0.05 /h. The system has no quantitative anchoring, so "serum bound ×0.1 / ×10" reaches nowhere near a regime where supply could be limiting.
- Reason 2: nucleosides and hypoxanthine do not belong in S2 at any bound that rescues.
  - FBS contributes only low-µM hypoxanthine, uridine and thymidine at 10%.
  - DepMap shows de novo purine and pyrimidine genes as strongly essential in RPMI/DMEM + 10% FBS. Salvage does not buffer them there.
  - Nucleoside rescue (uridine for DHODH/brequinar, hypoxanthine for purine inhibitors) needs 10–100 µM of added nucleoside.
  - Including nucleosides in S2 therefore erases correct calls: PPAT/GART/PFAS/PAICS/ADSL/ATIC, CAD/DHODH/UMPS, CTPS1, TYMS/DHFR.
  - It also removes NDUF/UQCR/COX essentiality. Under ETC loss, DHODH loses its CoQ re-oxidation and pyrimidine synthesis stops unless a nucleoside is supplied. This is the main reason ETC genes are essential in standard media.
- Required change:
  - (a) Express every non-core uptake bound as a **fraction of the biomass-stoichiometric requirement** of the product it can replace, at WT growth. Sweep that fraction on a log grid {0, 0.001, 0.01, 0.1, 1, ∞}, instead of using a single absolute 10 with ×0.1/×10.
  - (b) **Remove nucleosides and hypoxanthine from S2.** Make a separate, labelled supplement variant ("S-nuc") used only for media that actually contain them (Ham's F-12 and DMEM/F-12 contain hypoxanthine and thymidine; HPLM contains hypoxanthine; check the user-supplied compositions for each).
  - (c) For M3 ("bound_i = 1000·c_i/c_glucose"), replace the single linear rule.
    - Uptake saturates at transporter Km (about 0.05–1 mM for amino acids) and is not proportional to c.
    - It is also demand-limited for trace nutrients.
    - Use bound_i = min(Vmax_class, k·c_i) with a few measured anchors, or abandon M3 and report only the availability ladder plus the log-fraction sweep.
  - (d) Add a **growth-calibration check** (label-free, M3 pilot): with glucose about 10 mmol/gDW/h, glutamine about 2, O2 about 10–15 and lactate/NH4 secretion free, the WT growth rate should land in 0.01–0.06 /h. Without a calibrated regime, the OXPHOS genes (which show the partial ratios in the README) and the energy-related essentiality have no physiological meaning.

**M2. Add a transporter/permeability gate to "availability" (§W2 M1).**
- Some nutrients present in the medium are not taken up by most cells. Open uptake of 1000 via any transporter is wrong for them.
- Aspartate (RPMI 0.15 mM) is not usefully imported unless SLC1A3 (or SLC1A1/2) is expressed. This is why GOT2 and aspartate synthesis remain essential in RPMI.
- Glutamate is mostly exchanged against cystine through xCT and gives little net uptake.
- Asparagine is imported by SLC1A5/SLC38A-type transporters and is limited by low uptake. ASNS dependency tracks DMEM (no Asn) and ASNS/Asn availability.
- Required change:
  - Define M1g as M1 plus a **pre-registered** list of uptake gates, fixed before labels: Asp and Glu = 0 unless SLC1A3/SLC1A2/SLC1A1 TPM ≥ threshold; cystine versus cysteine handled explicitly (RPMI/DMEM contain **cystine only**; M0 has cysteine only).
  - Add cystine and make cysteine uptake unavailable by default, so SLC7A11 is meaningful.
  - Test M1 versus M1g as separate ladder steps.
- Positive controls: GOT2, GOT1, ASNS (in DMEM), PYCR/ALDH18A1 (DMEM lacks Pro), SLC7A11.

**M3. Positive control for the sterol branch is mis-specified (§W2 controls; also LITERATURE §4.4).**
- "Cholesterol-synthesis genes gain essentiality under serum-free media" is true only for the **post-FPP sterol branch**: FDFT1, SQLE, LSS, CYP51A1, MSMO1, NSDHL, HSD17B7, SC5D, EBP, DHCR7, DHCR24.
- The upstream mevalonate genes HMGCS1, HMGCR, MVK, PMVK, MVD, IDI1, FDPS and GGPS1 are essential in most DepMap lines **even with serum cholesterol**. They supply non-sterol isoprenoids: CoQ10 side-chain, heme A, dolichol, protein prenylation.
- Human-GEM biomass (13 species) has none of these as demands. Adding serum cholesterol (S1) will therefore make HMGCR/MVK/FDPS non-essential in FBA, which creates new false negatives that were not present at M0.
- Required change (§W2, §W3, §9.4):
  - Split the control list into sterol-branch (expect gain of essentiality when serum-free) and mevalonate-core (expect essential in all media).
  - Add the cofactor/prosthetic demands as a pre-registered biomass addition (see S1) and report the mevalonate-core genes as a check on M2a/M2b.

**M4. Genomic constraint "gene-loss KO" can make a viable line non-viable (§W4 genomic constraint, §6.1 failed context models).**
- A line that has lost a gene (homozygous deletion) is alive. If the forced KO makes the line's model infeasible or zero-growth, the model, not the line, is wrong (bypass missing, or GPR error).
- "WT growth ≤ 1e-6 gives the generic parent's scores, flagged" discards all other line-specific information, and >5% such failures marks the method as failing.
- Required change: pre-register a **viability-preserving rule**. Apply gene-loss KOs one at a time, and drop any individual loss that makes the WT model non-growing (log it as a "model-contradicted loss", with the gene list). Do not fall back for the whole line. Report the list of contradicted losses; it is itself an error-analysis result.

**M5. Backgrounds (§W4 levels) are mostly infeasible under whole-lineage hold-out, and several have no GEM mechanism.**
- The group rule is ≥10 training lines and ≥5 test lines. BRAF V600, STK11, KEAP1 and IDH1/2 are lineage-restricted (skin/thyroid/colorectal; lung; glioma/myeloid/soft tissue). With five test lineages, most groups will have none or fewer than 5 test lines. Expect only KRAS hotspot and possibly MTAP deletion to qualify.
- KRAS, BRAF, PIK3CA are signalling drivers. A stoichiometric GEM has no handle on them except through expression. A "KRAS GEM" is a null by construction, not a test.
- IDH1/2 hotspot is neomorphic (2-HG production from αKG). It is not represented in Human-GEM, so a background model cannot encode it unless a 2-HG-producing reaction is added.
- FH/SDH-loss lines are too few for any group.
- Required change:
  - Define **metabolism-linked** backgrounds, each with a mechanism, an expected dependency and a count (computed label-free from mutation/CN/expression):
    - enzyme loss with a paralog or bypass: ENO1 deletion (1p36, expect ENO2), ME2 loss (18q with SMAD4, expect ME3), NAPRT silenced (expect NAMPT), ASS1 silenced (arginine route), PHGDH/PSAT1 high or amplified, SLC7A11/NRF2-active (KEAP1, NFE2L2), CPS1-high (STK11 loss), UCK/DCK/other salvage-pathway losses;
    - MTAP deletion kept as a deliberate failure case (see below);
    - signalling drivers kept only as an optional negative-control group.
  - For group evaluation, allow a **within-lineage mutant versus WT contrast pooled across lineages** (lineage as stratum), not a requirement that the group is split by whole lineages. Report the number of qualifying groups before any scoring.
  - If IDH1/2 is kept, add the neomorphic 2-HG reaction as a pre-registered edit.

## SHOULD-FIX

**S1. Biomass additions (§W3 B-loo_j only removes components).**
- DepMap-essential metabolic genes frequently serve demands that are not in MAR13082: CoQ10, heme A and heme, dolichol-P (N-glycosylation), GPI anchors, lipoate, CoA, FeS cluster, SAM (methylation sink), tRNA modification, polyamines.
- Add a pre-registered **B-add** family, one component at a time at a fixed small coefficient, plus the full set.
- Check each gene group's lethality against the co-factor and prosthetic demands (COQ2–9, HMBS/UROD/FECH, DPAGT1/ALG*, PGAP*, LIAS/LIPT1, DHFR/MTHFD1, MAT2A).
- Also add cell-type composition variants (lipid pool ±50%, DNA/RNA scaled with doubling time, glycogen on/off). These are cheap, and LITERATURE §5 says composition, not coefficients, drives essentiality.

**S2. Serum-free and defined supplements (§W2 parser).**
- "OPAC", "Tumor Stem Media" and other serum-free formulations contain B27/N2/ITS, which are not "no serum". B27 supplies linoleic and linolenic acid, vitamin E, carnitine, ethanolamine, selenium, transferrin and albumin-bound lipids.
- The supplement table must cover: B27, N2, ITS, NEAA (adds Ala, Asn, Asp, Gly, Glu, Pro, Ser, which directly affects ASNS, PYCR, GOT), sodium pyruvate, GlutaMAX (dipeptide), 2-mercaptoethanol (bypasses xCT/cystine), HEPES, hydrocortisone, matrigel/BME.
- DMEM with/without pyruvate and with high/low glucose is ambiguous ("DMEM + 10% FBS"). Run both and flag the ambiguity.
- Pyruvate rescues ETC/aspartate-dependent essentials.

**S3. HPLM/Plasmax variants (§W2).**
- Audit exchange-ID coverage in Human1/Human2 for each component: uric acid, 3-hydroxybutyrate, citrulline, creatinine, taurine, carnitine, hypoxanthine and others. Report unmapped components.
- Do not use these for selection (already stated). Keep the expected-direction checks, but list also the expected direction for pyruvate-containing media: GOT/MDH and OXPHOS essentiality is relaxed by pyruvate.

**S4. Negative-control and positive-control gene panels should be fixed in advance as files (§W2, §10).**
Create `data/controls/gene_panels.tsv` with columns gene, class, medium-condition expectation, evidence reference, and use it in W1, W2, W3, W4 and W5. Suggested panels, to be justified by citation before the freeze:
- **Pan-medium essentials** (expected essential in FBA and DepMap): CAD, DHODH, UMPS, CTPS1, TYMS, DHFR, RRM1/2, GUK1, PPAT, GART, PFAS, PAICS, ADSL, ATIC, IMPDH2, GMPS, PCYT1A, CHKA, CS, MDH2, SDH subunits, HMGCR/MVK/PMVK/MVD/FDPS (mevalonate core), NAMPT, COASY.
- **Medium-dependent** (expected to move with medium):
  - ASNS (DMEM versus RPMI/Asn);
  - sterol branch (serum-free);
  - purine and pyrimidine de novo (Ham's F-12 / DMEM/F-12 / HPLM versus RPMI/DMEM);
  - GLUL (Gln-free);
  - PYCR1/ALDH18A1 (DMEM, Pro-free);
  - SLC7A11 (cystine versus 2-ME);
  - NAMPT/NAPRT (nicotinic acid);
  - GOT1/GOT2 (pyruvate, Asp transporter);
  - PHGDH/PSAT1/PSPH (Ser/Gly).
- **Collateral-lethality pairs** (positive controls for the genomic constraint and for GPR OR logic): ENO1 loss → ENO2, ME2 loss → ME3.
- **Known GEM false-positive classes** (expected not essential in DepMap): β-oxidation (ACADL, HADHA/B, ECHS1 in H2), xenobiotic CYP/UGT/GST, urea cycle (OTC, ARG1), bile acid and steroid hormone genes, drug-metabolising enzymes, tissue-specific synthases for unexpressed pathways.
- **Known GEM false-negative classes** (expected essential, no model mechanism): V-ATPase subunits (ATP6V*), aminoacyl-tRNA synthetases if included, proteasome-adjacent, PGD (6-phosphogluconate toxicity, not NADPH), MTAP-loss → MAT2A/PRMT5 (MTA inhibition of PRMT5 is not stoichiometric), GPX4/ferroptosis genes.

**S5. CRISPR-side artefacts (§3 labels, §9.4 failure account).** Add these label-side strata, computed from annotation files before scoring:
- **Multi-target guides / multi-copy loci:** UGT1A/UGT2B, GSTM/GSTT, AMY1A/B/C, NAT1/2, GBA/GBAP1, SLC25A5/6-like duplicated loci, HIST clusters. Flag genes whose guides have more than one perfect genomic match (from the Avana/KY guide map) and report with and without them.
- **Amplicon effects:** flag genes within a recurrent local amplicon (relative CN ≥ 4 or a segment-level focal amplification). Chronos corrects CN only partly, so amplified loci appear falsely essential. §W1 failure test 3 stratifies by CN; make "CN ≥ 4" and "multi-target" explicit strata at the label level, not only in the failure account.
- **Non-metabolic reason for essentiality:** pre-register a "housekeeping-adjacent" stratum of metabolic-annotated genes with core cellular roles: aminoacyl-tRNA synthetases and tRNA-modifying enzymes, V-ATPase, proteasome-linked, Fe-S biogenesis, GPI/N-glycosylation, nuclear-pore-linked, complex-assembly factors. Moonlighting enzymes (GAPDH, PKM, ENO1, PGK1, TPI1, ACLY) go in a separate small list. Report AP with and without them.
- **Mitochondrial genes:** the 13 MT- genes are absent from DepMap CRISPR (already handled in §2.3). Note that nuclear-encoded mitochondrial translation and assembly factors are not in the GEM, which is a source of OXPHOS false negatives.
- **Doubling time:** growth rate shifts effect-size distributions per line (slower lines have weaker Chronos effects). If DoublingTime is available in Model.csv or Passports, add it as a covariate in the per-line prevalence report and in the shuffled-medium control.

**S6. GPR curation rules (§W5 evidence table) misfire for isozyme heteromers and compartment-split paralogs.**
- "OR → AND if all branch genes are co-members of a CORUM/Complex Portal complex" would wrongly AND true isozymes that co-assemble: LDHA/LDHB, PFKL/M/P, ENO1/2/3, ALDOA/B/C, CKB/CKM, plus tissue-specific subunit isoforms (COX4I1/I2, COX6A1/2, COX7A1/2, ATP5MC1/2/3, ATP6V1B1/2, ATP6V0A1–4).
- Required change: OR → AND only if the branch genes are **not Ensembl paralogs of each other** (distinct subunit families). AND → OR only if paralogs. Add a pre-registered exception list for heteromeric isozymes.
- "Paralog-buffered" (§9.4) should require that the OR partner acts in the **same compartment and same reaction** (cytosolic GOT1 cannot buffer mitochondrial GOT2; MTHFD1 vs MTHFD2/1L; IDH1 vs IDH2; ME1/2/3). Use UniProt or HPA localisation.
- "Never detected in Nusinow" has a detection bias against membrane transporters and low-abundance proteins. Restrict the protein-absent evidence to soluble enzymes, and use RNA only for transporters.

**S7. "Metabolic drugs" (§W6.2) needs a curated definition, and the assay medium differs from the CRISPR medium.**
- "Annotated target is in U" is too broad (kinase inhibitors hitting off-target metabolic kinases, mislabelled repurposing annotations).
- Pre-register a curated list: compounds whose **primary mechanism** is inhibition of, or a substrate for, a Human-GEM metabolic enzyme or transporter.
  - Expected included classes: antifolates (methotrexate, pemetrexed, pralatrexate, trimetrexate: DHFR/TYMS), TYMS (5-FU, floxuridine, raltitrexed), ribonucleotide reductase (hydroxyurea, gemcitabine, triapine), DHODH (brequinar, leflunomide/teriflunomide), IMPDH (mycophenolate, ribavirin), purine analogues (6-MP, 6-TG, fludarabine, cladribine, clofarabine), statins and FDPS inhibitors (HMGCR, FDPS), FASN/ACLY/SCD inhibitors, GLS (CB-839), IDH-mutant inhibitors (ivosidenib, enasidenib; mutant-specific, with caveat), ODC1 (eflornithine), NAMPT (FK866, GNE-617), complex I (IACS-010759, phenformin).
  - Excluded: kinase inhibitors, cytotoxic non-metabolic agents.
- Each compound gets a **prodrug/activation flag** (5-FU needs UPRT/TK, gemcitabine needs DCK). The GEM target-KO score ignores activation, so predictions for activated drugs need the activating gene as a second feature.
- PRISM pools are screened in **PRISM's own medium** (RPMI-based, uniform across lines), and GDSC lines use their own medium. W6 must set the medium for each drug dataset to the assay medium. The W2 documented-medium is for CRISPR only. Record this in §W6 and §2.2.
- Drug responses at one dose (PRISM) track proliferation rate; add doubling time or baseline growth as a covariate, as well as lineage.
- Expected positive controls: DHODH–brequinar, DHFR–methotrexate, HMGCR–statins, TYMS–5-FU, IMPDH–mycophenolate. Report uridine/hypoxanthine-dependent shifts as a medium check.

**S8. H1 TYMS reversibility (§2.1).**
Keeping Human1 "as released" is fine for baseline G\*. But any **optimised** H1 pipeline must apply the irreversible TYMS (MAR04644 lb = 0) as a fixed, label-free biology correction (humans do not run TYMS backwards at relevant rates, and UPRT activity on uracil is minimal). Without it, the H1 replication of P1 (§W8) will be penalised on pyrimidine genes for a reason unrelated to the optimisation. Pre-register this before M1 as a "parent patch", not as a W5 edit.

## OPTIONAL

**O1. Oracle ceiling for "error attributable to medium".**
- Besides the documented-medium ladder, report an upper bound: for each false call at M0, whether any single nutrient add/remove (from the medium component list) flips it to correct, against a null of random single-nutrient flips at the same rate.
- Expect the attributable share to be small (tens of genes: nucleotide, lipid, Asn/Asp, cystine) compared with the FAO/tRNA/V-ATPase/non-metabolic error blocks. Say so in the pre-registration, so a small share is not read as a failure.

**O2. MTAP-deletion as a labelled failure case.** MAT2A, PRMT5 and RIOK1 dependency in MTAP-deleted lines acts through MTA accumulation inhibiting PRMT5, which FBA cannot represent. Also, the 9p21 deletion removes CDKN2A. Include MTAP (and PRMT5/MAT2A) in §9.4 as an expected unfixable mechanism.

**O3. Tumour-physiology medium (TIFM).** Optionally add one tumour-interstitial-fluid-like variant (Sullivan 2019) next to HPLM/Plasmax, as in the prompt's "worse for tumour physiology".

**O4. Hypoxia as a GEM variant.** Optionally add a low-O2 uptake bound variant (a few mmol/gDW/h). Fixed uptake bounds expected to affect only OXPHOS and heme/desaturase genes; label-free check.

## Expected failure modes (for the pre-registered failure account)
1. W2 will improve DMEM/RPMI-specific genes (ASNS, GOT, PYCR, sterol branch) but break nucleotide and mevalonate-core genes unless M1 above is applied. A "documented medium" gain at the first attempt would be suspicious.
2. HPLM/Plasmax will lose essentiality of de novo purine genes. This is correct biology, but the DepMap labels were measured in RPMI/DMEM, so apparent agreement loss is not a model error.
3. The shuffled-medium control will have a non-trivial gain because the media are not independent of lineage (suspension versus adherent, DMEM-adherent lines). Within-lineage permutation is right. Also permute within lineage × serum-free status.
4. Context-specific families will gain mostly on unexpressed-gene easy negatives (stratum in §3 handles this). Genomic constraint adds a few positive controls (ENO1, ME2) but will mostly add ratio-1 genes.
5. GPR curation will mostly shift isozyme/complex cases (OXPHOS subunit isoforms, ATP6V, COX4I2); see S6 for false curation.
6. Per-gene AP will be dominated by tissue-restricted pathways (heme, melanin, urea cycle, steroid) that merely track lineage, not cell-line metabolism.
7. The performance ceiling for any GEM is bounded by non-metabolic essentiality in U. State a rough upper bound (fraction of DepMap-essential U genes with a metabolic mechanism) in the report.
