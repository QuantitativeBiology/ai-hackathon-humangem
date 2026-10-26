# LITERATURE.md — critical survey for the nextgem design

Literature scientist report, 2026-10-26. Scope: the eight topics set by the PI for DESIGN.md. Sources were PubMed/PMC full texts, bioRxiv, the Human-GEM GitHub repository and release notes, and our local copies of Human-GEM v2.1.1 and DepMap 26Q3. Every citation was checked against PubMed or the publisher record. Numbered references are at the end. Papers read only as abstracts are marked *(abstract only)*. Inaccessible papers that matter are listed in §10.

**Evidence tags.** **[R]** = replicated by ≥2 independent groups or datasets; **[S]** = a single study (or one group's closely related studies); **[P]** = preprint or not peer-reviewed; **[L]** = verified by us in local files or the repository.

---

## 0. Key takeaways (one screen)

1. **"Human1" needs an exact definition.** Published Human1 is v1.3.0: 13,417 reactions, 10,138 metabolites, 3,625 genes [1]. The **last v1.x release is v1.19.0** (2024-06-24): 12,971 reactions, 8,455 metabolites, 2,887 genes [L, GitHub]. The counts in the Human2 paper (2,848 genes, 12,931 reactions, 8,461 metabolites) [2] are **exactly v2.0.0**. Our local v2.1.1 has 12,748 reactions, 8,374 metabolites and 2,870 genes, with 130 GPRs changed since v2.0.1 [L]. All three are "Human2" by name, so we should benchmark the exact versions.
2. **The reported GEM essentiality baselines are low and hard to compare.** Most report MCC ≈ 0.25–0.45 on binarised labels, after selection among many parameter sets [1, 5, 7, 17]. The only threshold-free numbers we found come from the Human-GEM repository (v2.1.1, ftINIT, Hart 2015 lines): **AUPRC 0.25–0.33 against base rates of 0.14–0.18, pooled AUPRC 0.289 against a base rate of 0.158** [L]. No published study reports per-gene or per-lineage AUPRC on DepMap at full scale with lineage hold-out. That gap is ours to fill.
3. **Pooled metrics are inflated by gene identity.** One graph-learning predictor reports ROC 0.94 on random gene–line pairs but a mean per-gene AUROC of only 0.61 [6, P]. About 15% of Human2 metabolic genes screened in DepMap 26Q3 are DepMap "inferred common essentials" (419/2,798) [L]. Simple gene-prior baselines must be in the benchmark.
4. **"Context-specific models cluster by method, not tissue"** is well supported for reaction content [7, 8, 9]. But almost all of that evidence comes from one lab lineage (the Lewis lab), with Recon1/Recon2.2/iHsa and NCI-60 or ≤4 cell lines. ftINIT itself *has* been run at DepMap scale: Gustafsson et al. built ftINIT/tINIT models for 891 DepMap lines and scored essentiality by MCC [20]. What has **not** been tested is the method-vs-tissue question: variance partitioning of model content and predictions by extraction method vs lineage, on Human1/Human2, at DepMap scale. Testing that is a contribution. *(Corrected 2026-10-26, design review r1, PI decision 33.)*
5. **Medium is a large, replicated confounder for metabolic essentiality** [R: 10, 11, 12]. In DepMap, culture medium associates with CRISPR dependency for 26% of variable metabolic genes versus 3% of other genes [10]. DepMap 26Q3 has per-condition medium annotation linked to each screen (`ScreenSequenceMap.ModelConditionID → ModelCondition.GrowthMedia`) [L]. Serum is the main unmodelled component.
6. **GPR errors are a primary error source.** This is shown for E. coli with AUPRC [13], for Human1 versions [1] and for Human2 (wrong "OR" isozymes created false bypasses) [2]. Paralog buffering explains part of the variable essentiality [R: 14, 15].
7. **Hybrid ML/FBA beats plain FBA in microbes and CHO cells** [16], but no peer-reviewed human-cancer study has shown this with a leakage-free split. The one cancer GEM+ML study (MCC 0.27 → 0.33) used random cross-validation in 50 breast lines [17].
8. **Flux states vs omics:** we found no rigorous head-to-head showing that GEM flux profiles separate subtypes or predict drug response better than the expression or proteomics they were built from. Proteomics predicts drug response about as well as transcriptomics [18]. This remains an open, testable question with a high bar.

---

## 1. Human1 and Human2

### 1.1 Human1 (Robinson et al. 2020) [1]
- Human1 merged HMR2, iHsa and Recon3D into one consensus model. The curation removed 8,185 duplicate reactions and 3,215 duplicate metabolites and re-balanced 3,226 reactions. The published model had 13,417 reactions, 10,138 metabolites (4,164 unique) and 3,625 genes, with 100% stoichiometric consistency, 99.4% mass balance and 98.2% charge balance (MEMOTE). The paper's own version was **v1.3.0** [1]. GPRs were built by combining HMR2, Recon3D and iHsa rules with complex information from CORUM [1].
- **Essentiality benchmark.**
  - Data: tINIT context-specific models from RNA-seq for the 5 Hart 2015 TKO lines (GBM, RPE1, HCT116, DLD1, HeLa) and for **621 DepMap lines**.
  - Essentiality definition: a gene was essential if its knockout made any of **57 metabolic tasks** fail, biomass production included.
  - Results: accuracy was 86–88% for all reference models, "the same … as predicting all genes non-essential". MCC was low overall, but Human1-derived models had **>2.5× the MCC** of HMR2/Recon3D-derived models on Hart lines and **2.8× the MCC** of Recon3D on DepMap.
  - Robustness: results held across a range of CRISPR thresholds (Fig. S10). Task-based essentiality increased sensitivity over biomass-only essentiality (Fig. S11) [1].
  - Exact MCC values are only in figures and the SI, which we could not access (§10) [S].
- Across the 27 pre-release versions, gene and GPR updates gave the largest gains in essentiality prediction. Mass balancing and ETC reversibility fixes also helped [1, S].
- ecGEMs (GECKO) were built for 11 NCI-60 lines. Only 11 were used because exchange data for the other lines were "deemed unreliable due to early depletion of one or more nutrients" [1].

### 1.2 Human2 (Luo et al., PNAS 2026) [2]
- **LLM-assisted curation.**
  - GPT-4 prompts, using UniProt and HPA function text, assessed all **26,246 gene–reaction pairs** in Human1. 7,774 were judged consistent, 2,195 inconsistent and 16,277 inconclusive (excluded).
  - Manual review of the 2,195 found 210 false flags, leaving **1,985 confirmed wrong pairs**. Fixing them updated 1,135 GPRs and removed 203 non-metabolic genes [2].
- **Manual curation** revised 775 reactions and refined 717 lipid GPRs. It added 210 reactions, 207 metabolites and 55 genes, and removed 76 duplicate reactions and 50 duplicate metabolites. Localisation came from UniProt, HPA and DeepLoc2. Overall, from v1 to v2 there were 1,864 GPR and 775 reaction changes [2].
- **Final counts reported:** 2,848 genes, 12,931 reactions and 8,461 metabolites. These are exactly the v2.0.0 README numbers [L].
- **Reported benchmarks** [2, S]:
  - MCC of gene essentiality across Human-GEM versions against the Hart 2015 CRISPR data (Fig. 2D, values not in the text).
  - ftINIT cell-line models from Hart 2015 and DepMap RNA-seq: "Human2 outperforming Human1 in MCC".
  - "Consistent performance" on Project Score (Sanger) as an independent CRISPR benchmark.
  - 112 inborn errors of metabolism: ec-Human2 accuracy 65.4%.
  - Flux consistency on NCI-60: ec-Human2 ~81% vs ec-Human1 79%. MEMOTE score 81%.
- **Mechanism claimed:** correcting inaccurate "OR" relationships removed false alternative pathways. Those pathways had made genes look non-essential and made metabolic tasks pass wrongly [2].
- **Critical notes.**
  - Neither paper reports AUPRC or per-lineage performance. Both use binarised labels and MCC on context-specific (tINIT/ftINIT) models, not the generic model. So no published baseline exists for "generic Human1 vs generic Human2 FBA vs DepMap" [S].
  - The LLM step excluded 62% of pairs as inconclusive. GPR errors therefore remain in the parts of the model that are poorly annotated [2].
  - One 2026 study chose Recon 2M.2 over Human1 for essentiality. It claims the "inflated GPR annotations" of newer models create bypasses [17, S, secondary claim]. That claim was not independently tested there.

### 1.3 Repository facts relevant to our pipeline [L]
- **Releases.** v1.19.0 (2024-06-24) is the last v1.x release. v2.0.0 followed on 2026-03-30, then v2.0.1, v2.1.0 and v2.1.1 (2026-10-07).
- **Gene counts over releases.** 3,625 (v1.3–v1.10) → 3,067 (v1.11, 2022-02) → 2,920 (v1.14) → 2,887 (v1.19) → 2,848 (v2.0.0/2.0.1) → 2,870 (v2.1.x). Much of the "Human1→Human2" gene cleanup therefore already happened inside v1.x.
- **File defects.** The v2.0.1 and v2.1.0 SBML/MAT files are incomplete; the errata say to use v2.1.1 or the YAML. From v2.0.1, SBML IDs carry `M_/R_/G_` prefixes and `fbc:label` holds Ensembl IDs. Joins to DepMap must go through `genes.tsv` (Entrez).
- **Spontaneous flag.** v2.1.0 marks diffusion and non-enzymatic reactions as `spontaneous`, so extraction keeps them without a gene. This matters for any extraction method we reimplement.
- **Repository essentiality test (v2.1.1).**
  - Script: `code/test/gradedEssentiality.py`. Labels are Hart 2015 Bayes factors at 5% FDR per line. Prediction is ftINIT plus the single-deletion growth ratio on Ham's medium.
  - Results: viability-task MCC 0.30–0.37; growth-ratio AUROC ≈ 0.66; **AUPRC 0.25–0.33 against base rates of 0.136–0.179** (`data/testResults/gene-essential_summary.md`).
  - The script splits tasks into *viability* (GR, ER) and *capability* (SU, BS, IC). It counts capability-only essentiality as a false positive. Example: ETFA/ETFB/ETFDH are essential only for β-oxidation and phospholipid tasks, while their growth ratio is 1.0 [L].

---

## 2. Gene-essentiality prediction with GEMs vs CRISPR screens

### 2.1 Reported baselines (all on binarised labels; none lineage-held-out)

| Study | Model / extraction | Data | Metric | Result | Tag |
|---|---|---|---|---|---|
| Robinson 2020 [1] | Human1, tINIT, 57-task essentiality | 621 DepMap lines; 5 Hart lines | MCC | Human1 2.8× Recon3D; values in SI | S |
| Opdam 2017 [7] | Recon1/2.2; 6 MEMs × 4 thresholds × 3 constraints | 4 lines (A375, HL60, K562, KBM7) | 1-tailed Wilcoxon p (CRISPR LFC of predicted essentials) | extracted > generic; MEM choice explains most variance in accuracy | S |
| Richelle 2019a [8] | Recon2.2/iHsa, 6 MEMs ± task protection | 20 NCI-60 lines | % predicted essentials with negative score | task protection +≤5% | S |
| Tobalina 2016 [19] | fast context-specific reconstruction | 174 CCLE lines vs Achilles | agreement | conflicting; many essentials also appear with *random* expression; biomass strongly influences results | S |
| Gustafsson 2023 [20] | Human1, ftINIT vs tINIT | 891 DepMap lines | MCC | ftINIT ≈ tINIT, >10× faster | S |
| Vieira 2022 (Troppo) [5] | Human-GEM template, >6,000 models, many threshold × MEM combos | **733 CCLE lines**; parameter combinations chosen using **MCF7 reference models** (fluxomics), then applied to all lines; CERES thresholds −1.5…−0.5 | MCC | best combination 0.44 vs <0.4 for Robinson; medium-constrained uptake improved prediction | S (configuration selected on one reference line, not a lineage hold-out; exact grid/selection details to re-check in full text) |
| Kim 2026 [17] | Recon 2M.2, tINIT, MOMA ± RF on KO-flux features | 50 breast lines | MCC, sens. | MOMA MCC 0.27 (sens 0.37) → RF MCC 0.33 (sens 0.55, precision 0.31) | S, random 5-fold CV |
| Human-GEM repo v2.1.1 | Human2, ftINIT, growth ratio | 5 Hart lines | AUPRC / base rate | 0.25–0.33 / 0.14–0.18 | L |

- **Comparability caveat.** Label definitions differ: BAGEL Bayes factor, CERES/Chronos cut-offs, negative score. Metrics differ too: MCC, p-value, % agreement. The model scope also differs: only genes in the model, and sometimes only genes whose knockout changes anything. None of these numbers can serve directly as our baseline. We must recompute generic Human1 and Human2 on one label set.
- **Older evidence** [S]: FBA predicted essential genes beyond chance in ccRCC but not in prostate adenocarcinoma, compared with in-vitro screens [21]. Performance is context-dependent even within one pipeline.

### 2.2 Label definitions (DepMap)
- **Chronos gene effect** [22]. It models read counts with a negative binomial and per-sgRNA and per-line knockout efficacy, and includes a copy-number correction. The authors recommend a *global* shift and scale so that the median of nonessential genes is 0 and the median of essential genes is −1 [22]. Thresholds like −0.5 or −1 are therefore relative to control sets, not absolute.
- **Probability of dependency (`CRISPRGeneDependency`).** DepMap staff recommend binarising with P(dependency) > 0.5. They describe it as the probability that the score comes from the essential rather than the non-essential distribution, "not intended as a value used to determine significance". This source is the DepMap forum, not peer-reviewed, so the DepMap release notes must be cited in DESIGN [P].
  - Our local 26Q3 folder has `CRISPRGeneEffect.csv` and `ScreenGeneEffect.csv`, but **no `CRISPRGeneDependency.csv`** [L]. The data agent should fetch it.
- **Common essentials.** DepMap ships `CRISPRInferredCommonEssentials` (2,263 genes; **419 overlap Human2 genes screened in 26Q3**, ≈15% of 2,798) [L]. It also ships `AchillesCommonEssentialControls` (1,247 genes; 209 in Human2), which is the positive-control set from Hart 2015 and Blomen 2015, not derived from the gene-effect matrix (DepMap forum) [P].
- **Independent core-fitness definitions:** CEG/NEG [23, 24] and CoRe [25]. The Broad (Achilles) and Sanger (Project Score) screens are highly concordant. Batch effects come mainly from the library and the assay length [26, R]. Integrated Broad+Sanger datasets exist [27]. Project Score [28] is the natural external test set.

### 2.3 Pitfalls
1. **Class imbalance.** Accuracy is meaningless here; Robinson et al. say so explicitly [1]. AUPRC depends on prevalence, so it must always be reported with the base rate (§8).
2. **Common-essential inflation.**
   - Many true positives are pan-essential genes (nucleotide synthesis, CoA, heme, OXPHOS). A model that predicts "essential everywhere" for these scores well in pooled metrics while explaining nothing context-specific.
   - DeepMeta shows the effect: ROC 0.94 on a random 10% of gene–line pairs, but a 10-fold mean per-gene AUROC of 0.61 (0.49 with scrambled labels) [6, P].
   - Lagziel et al. found that metabolic genes' cell-line-specific dependencies are *harder* to predict from omics than other genes'. Predictive models existed for 32% of variable metabolic genes (CRISPR) [10, S].
3. **Copy-number and cutting-toxicity artefacts.**
   - CRISPR cutting in amplified regions causes gene-independent depletion [29, 30, R]. Structural rearrangements cause this too [31]. CRISPRcleanR [32] and Chronos [22] correct for it.
   - Chronos notes that for essential genes the CN–effect relation is **non-monotonic** (more copies, less complete knockout) [22]. CN-driven residuals should be checked in our false positives and false negatives (e.g. stratify errors by CN).
4. **Genes absent from screens.**
   - Chronos does not accept sgRNAs that target multiple genes [22]. Highly similar paralogs or duplicated metabolic genes may be missing or weakly measured. 70 Human2 genes with Entrez IDs are absent from 26Q3 CRISPR (2,868 − 2,798) [L].
   - Missing genes must be excluded from evaluation, not counted as non-essential.
5. **Paralog masking.** About half of constitutively expressed genes are never essential in any DepMap screen, and these are enriched for paralogs [15, R with 14]. GEM "OR" rules can make a gene look non-essential when the isozyme is absent or non-functional in that line (§6).
6. **Medium confounding** (§4) and **growth-format/serum confounding** (organoids vs 2D) [12].
7. **Unexpressed genes as negatives.** Chronos uses unexpressed genes as negative controls and counts them as false positives if they fall in the bottom 15% [22]. A GEM that removes unexpressed genes gets an easy "non-essential" call on them. Report performance with and without trivially-unexpressed genes.

### 2.4 Non-FBA and hybrid ML/FBA predictors
- **FlowGAT** (graph attention network on FBA wild-type flux graphs) [33]: in E. coli its predictions are "close to those of FBA", not better [S].
- **Flux Cone Learning** [16, S]: Monte Carlo sampling of each deletion's flux cone feeds a random forest. It beats FBA for essentiality in E. coli, *S. cerevisiae* and CHO cells, and argues that "optimality assumptions are not required". Its split is random by gene with 80/20 hold-out. It has **not been tested on human cancer lines**. Sampling is expensive (the yeast dataset was 4.4 GB).
- **FluxGAT** (arXiv, NeurIPS 2024 workshop) reports about 2× FBA sensitivity [P, not peer-reviewed; not read in full].
- **MOMA + random forest, breast cancer** [17, S]: MCC 0.27 → 0.33 on 50 lines. Classifiers were trained on knockout-flux vectors with SMOTE/undersampling and 5-fold CV. Genes recur across lines, so random CV likely leaks gene identity. **No lineage hold-out.**
- **DeepMeta** (Cell Rep 2025; preprint read) [6, P/S]:
  - Inputs: a GNN on INIT tissue-model enzyme graphs, sample expression relative to GTEx, and MSigDB features.
  - Labels: DepMap P > 0.5.
  - Results: per-gene AUROC 0.61. External validation on Sanger CRISPR, DEMETER2, a pediatric CRISPR set and PRISM.
  - **No FBA baseline.**
- **Network-informed omics regression** [10, S]: random forests on neighbour-enzyme expression and copy number, oncogene mutations, medium and lineage. Neighbour-enzyme expression was the top feature for 34% of CRISPR models. **Medium was the top feature for 11%**; the gene's own expression was top for only 8%.
- **How much hybrids beat FBA in human cells:** published gains are modest (ΔMCC ≈ 0.06 [17]). They come from designs that likely overstate generalisation. No study shows hybrid > FBA under lineage hold-out [inference from the above].

---

## 3. Context-specific model extraction

### 3.1 Methods (families as in [8])
- **GIMME-like** (minimise flux through low-expression reactions while keeping an objective): GIMME [34].
- **iMAT-like** (MILP trade-off between keeping high-expression reactions and removing low ones): iMAT [35], INIT/tINIT [36, 37], ftINIT [20]. ftINIT has two steps and includes most reactions without GPRs in mode "1+0". Its default is 1 TPM, and it is >10× faster than tINIT, which makes bootstraps and ensembles feasible [20].
- **MBA-like** (protect a core set, prune the rest): MBA [38], mCADRE [39], FASTCORE [40], FASTCORMICS/rFASTCORMICS [41, 42].
- **Flux-weighting / enzyme constraints:** E-Flux (expression as bounds) [43]. GECKO/ecModels [44] and GECKO 3.0 (deep-learning k_cat) [45] were used for ec-Human1/2 [1, 2]. Enzyme constraints mainly bound fluxes and growth. We found no study showing they improve *binary* knockout essentiality in human cells [no evidence found].
- **Phenotype-informed:** PRIME (growth-rate-informed models, >280 lines; predicted drug targets, MLYCD validated) [46, S].
- **Newer (2020–2026):**
  - Threshold and alternate-optima methods: StanDep, with cluster-specific thresholds that keep housekeeping reactions [47]. DEXOM enumerates the alternative optimal networks [48].
  - Pipelines and benchmarks: Troppo (Python pipeline; global+local thresholds best; medium constraints help) [5]. Gopalakrishnan 2023 (alternate optima; mCADRE most reproducible in mammalian models; tasks must be "explicitly and quantitatively protected") [9]. The 2021 GEMbench benchmark (4 MEMs on Recon3D, cell lines and patients) [49].
  - Task inference: CellFie (task activity inferred from omics, 195-task consensus) [50].
  - Flux-inference tools, not knockout predictors: Compass [51], scFEA [52], METAFlux [53].

### 3.2 Evidence that models cluster by method rather than tissue
- **Richelle 2019a** [8]: 44 NCI-60 lines × 6 MEMs × 2 reference models. PC1 explained 38% of reaction-content variance, and **>60% of PC1 was explained by the extraction method**. Inter-line similarity within a method varied widely, while the input transcriptomes were highly correlated. Protecting data-inferred tasks reduced the method effect, made line identity explain more variance, and slightly improved essentiality [S].
- **Opdam 2017** [7]: the *expression threshold* drives most reaction-content variance; the *MEM* drives most variance in essentiality accuracy; uptake constraints matter least for essentiality but affect functions and growth. Stringent thresholds (top 10%, mean) gave better essentiality but fewer functions [S].
- **Gopalakrishnan 2023** [9]: the spread of alternate optima depends strongly on the algorithm and on parent-model topology (fatty-acid metabolism, transport) [S].
- **Counter-evidence (single method).** tINIT models of 33 tumours cluster with their tissue of origin [1]. Human2 organ models from GTEx and from independent GEO data have similarity coefficients > 0.9 [2]. ftINIT models group by tissue *and* by technology (bulk vs single-cell) [20]. So within one method, biology is visible; the problem is *between* methods and thresholds.
- **Assessment.** The finding is consistent across these papers [R within the Lewis-lab lineage, plus [49]]. But it is (i) on Recon1/2.2/iHsa, (ii) on ≤44 lines, (iii) about reaction content more than predictions, and (iv) not tested as a method-vs-tissue variance partition on Human1/Human2 at DepMap scale (ftINIT models of 891 DepMap lines exist [20], but were scored only for essentiality MCC).

### 3.3 How to test it (design-ready)
- **Full-factorial build.** Lines × lineages × methods (≥3 families: ftINIT, iMAT-like, MBA-like such as FASTCORE/mCADRE, GIMME) × threshold rules (global percentile, ftINIT 1 TPM, StanDep/local [47]) × with and without task protection.
- **Variance partitioning** (PERMANOVA on Jaccard distances, or PC-wise R² as in [7, 8]) at three levels: reaction content, task capability, and knockout-essentiality profiles. Compare with the same partition on the input expression.
- **Negative control:** permute expression across lines within each gene. If models from permuted data separate lines as well as real models, the separation is method noise. Tobalina et al. found many predicted essentials recur with random expression [19].
- **Positive controls:** known lineage-specific metabolic features. Examples: melanocyte pigment (TYR/DCT) present in melanoma models; urea cycle and bile acid synthesis absent from most lines; ASS1-low lines needing arginine import.
- **Report:** the fraction of variance explained by method, threshold and lineage, with bootstrap CIs over lines.

### 3.4 Thresholding
- Thresholding is the main determinant of which reactions are "active" [54, S; consistent with 7]. Global thresholds remove low-expressed housekeeping enzymes; gene-specific thresholds help [47].
- In Troppo, *global* thresholds affected essentiality more than local ones, and the best combinations used global + local [5].
- **Rule for us:** fix the threshold family a priori, or tune it on training lineages only. Any configuration chosen with knowledge of the evaluation labels is optimistic. Troppo chose its parameter combinations using MCF7 reference models and then applied them to 733 lines [5]; its MCC 0.44 is therefore the score of a selected configuration, not a lineage-held-out estimate.

---

## 4. Medium and exchange bounds

### 4.1 Evidence that medium changes essentiality
- **DepMap analysis** [10, S but large]:
  - Medium associated with CRISPR dependency for 26% of variable metabolic genes (15/58) versus 3% of other genes (19/716). For RNAi the figures were 37% vs 18%.
  - ASNS was the strongest example: more essential in DMEM, which lacks asparagine.
  - The associations largely survived controlling for adherent vs suspension culture (~90%) and for lineage (64%).
- **Physiological-medium CRISPR screens** [11, S]:
  - K562 genome-wide screens in HPLM vs RPMI found 525 HPLM-essential and 427 RPMI-essential hits (|z| > 2), not enriched for metabolism libraries. Medium-specific hits did not correlate with expression changes.
  - In a 4-line follow-up, GPT2, MPC and GLS phenotypes were traced to HPLM-specific metabolites. Purine-synthesis genes were RPMI-essential, likely because HPLM supplies hypoxanthine. ACLY was RPMI-essential, consistent with HPLM acetate.
  - Of >800 DepMap screens, **76% used RPMI or DMEM and 82% used 10% FBS** [11].
- **DepMap NextGen (2026)** [12, S]:
  - 3D organoids in serum-free OPAC vs 2D RPMI–FBS. Medium effects concentrated in lipid-metabolism dependencies (cholesterol synthesis genes).
  - SCD dependency was stronger in NextGen models, probably masked by lipid-rich serum. Growth format mainly affected integrin and actin genes.
  - Organoid and spheroid screens are in DepMap; our 26Q3 `ModelCondition` lists OPAC conditions [L]. These models need their own medium definition or must be stratified.
- **Mechanistic medium effects:** HPLM rewires metabolism. Uric acid inhibits UMPS and reduces 5-FU sensitivity [55]. Plasmax: supraphysiological pyruvate in commercial media stabilises HIF1α; arginine reverses ASL; selenium prevents ferroptosis [56]. Review [57]. Together these are [R] for "medium changes metabolic phenotype".
- **In silico:** constraining uptake to medium-matched metabolites improved Troppo's flux and essentiality predictions [5, S]. In Opdam's grid, uptake constraints mattered least for essentiality but mattered for functions and growth [7, S]. In E. coli, many false essentiality calls came from vitamins and cofactors "available to mutants despite being absent from the medium" [13]. The human analogue is serum-derived metabolites (lipids, cholesterol, carnitine, nucleosides) [inference]. Cheng et al. 2023 found that cholesterol-synthesis genes are predicted essential in colorectal-cancer GEMs only when no cholesterol uptake reaction is allowed, and compared predictions with DepMap colorectal lines (abstract read) [4, S].

### 4.2 Sources of per-line medium annotation
- **DepMap 26Q3** [L]:
  - `Model.csv` has `OnboardedMedia`/`FormulationID`.
  - `ModelCondition.csv` has `GrowthMedia` (195 formulation IDs; blank for 1,214 of 3,385 conditions), `SerumFreeMedia` and `PlateCoating`.
  - `ScreenSequenceMap.csv` links each CRISPR screen to its `ModelConditionID`.
  - Common formulations: "RPMI + 10% FBS" (585), "DMEM + 10% FBS" (215), "OPAC:Ad+++" (119), "RPMI + 20% FBS" (96), "Tumor Stem Media" (84), "MEM + 10% FBS" (60).
  - **Use the screen's condition**, not the onboarding medium.
- **Formulation sources:** the composition behind each formulation ID is not in these files. It must be mapped to vendor formulations (RPMI-1640, DMEM, MEM, F-12, DMEM/F12, McCoy's 5A), with explicit HEPES, glutamine and pyruvate additions.
- **Fallbacks:** ATCC product sheets give the recommended base medium and serum per line (vendor documentation, not peer-reviewed). Cellosaurus (RRIDs, identity) [58] is useful for ID mapping. We could not confirm that Cellosaurus systematically annotates culture medium [unverified].
- **Serum:** FBS is undefined. HPLM studies used dialysed FBS to remove undefined polar metabolites [11]. Any "medium" model must make serum assumptions explicit (e.g. open uptake of lipids, cholesterol and trace metabolites at low bounds) and test their sensitivity.

### 4.3 Uptake-rate data
- **NCI-60 CORE** exometabolomics: 219 metabolites; glycine consumption correlates with proliferation [59].
- **Dolfi 2013:** fluxes scale with cell volume and protein synthesis, so per-cell rates need size normalisation [60].
- **Zielinski 2017:** NCI-60 flux states; glucose uptake correlates with chemoresistance [61].
- **Caveat:** Robinson et al. kept only 11 NCI-60 lines with reliable exchange data [1]. Quantitative uptake bounds exist for very few DepMap lines. Medium composition gives *qualitative* availability, not rates. FBA essentiality under uptake bounds of 1000 is mostly an availability question [inference]. The most recent Human-GEM precedent for quantitative exchange bounds is Huijer et al. 2026 [3]: regression-based exchange fluxes during exponential growth (REGP) constrained Human-GEM cell-line models much better than the common consumption-release (CORE) method; it does not address essentiality (abstract read) [S].

### 4.4 Implication for "error attributable to medium"
The design should run a nested comparison on identical genes and lines:
1. generic Ham's medium (the GR-task 44 inputs) [L];
2. a "documented medium" model (basal formulation + glutamine/pyruvate + an explicit serum set);
3. an "HPLM/plasma-like" variant;
4. an "all exchanges open" null.

Positive controls are genes with literature-documented medium dependence: ASNS (DMEM) [10]; GPT2, MPC1/2 and GLS (HPLM) [11]; purine synthesis (hypoxanthine) [11]; cholesterol synthesis (serum-free media) [12].

---

## 5. Biomass and metabolic tasks

- **Biomass.**
  - Human1 has one generic human biomass reaction built from tissue and cell composition data [1].
  - Seven human biomass formulations exist. Metabolite *composition* changes essentiality predictions, while *coefficients* mainly change growth-rate predictions [62, S]. Consistent with this, randomising all coefficients by ±20% did not change essentiality accuracy [7, S].
  - Cell-size and protein-content variation across NCI-60 [60] means quantitative growth predictions need line-specific scaling.
  - Human2's WBM work used organ-specific biomass variants: non-regenerating organs drop nuclear dNTPs, and RBCs drop amino acids [2]. These are examples of *a-priori* biology-driven variants.
- **Task lists available locally** [L]:
  - `metabolicTasks_Essential.txt`: 57 tasks in five categories (ER, IC, SU, BS, GR), including "Growth on Ham's media".
  - `metabolicTasks_Full.txt`: 257 tasks.
  - `metabolicTasks_CellfieConsensus.txt`: 195 tasks, from the Richelle/Lewis curation [8, 50]; the original curation had 210 tasks covering 7 systems [8].
- **Task-based essentiality.**
  - Robinson et al. found that 57-task essentiality increased sensitivity over biomass-only essentiality, at some cost in specificity [1].
  - The Human-GEM repository now separates viability from capability tasks, because capability-only essentiality produces false positives against proliferation screens (e.g. ETF complex) [L].
  - Richelle 2019a protected data-inferred tasks during extraction; essentiality was still scored on biomass [8].
- **Tasks beyond biomass for target finding.** Guil et al. 2024 added the essential metabolic tasks to genetic-minimal-cut-set analysis of Human GEMs; generic toxicities detected rose from 106 to 281, and new targets appeared in two myeloma lines (abstract read) [89, S]. Barrena et al. 2023 integrated Human1 with regulatory networks for gMCS-based essentiality and compared with large-scale silencing data [90, S].
- **Other objectives.**
  - pFBA [63]; MOMA [64] (used for knockouts in [17]).
  - Multi-objective / Pareto optimisation of biomass vs other demands in NCI-60 [65, S].
  - Objective-free approaches (flux sampling + ML) argue that optimality is unnecessary [16, S].
  - ATP maintenance: we found no human-cancer study that isolates the effect of the ATP maintenance (NGAM) setting on essentiality [no evidence found]; it should be tested as a sensitivity parameter.
- **Design implication.** Fix the candidate task set a priori from these lists before looking at labels, e.g. Human-GEM Essential tasks by category plus CellFie consensus tasks grouped by system. Evaluate tasks as (a) extra knockout criteria (viability-only vs +capability) and (b) protected functions during extraction [8].

---

## 6. GPR quality

### 6.1 Known GPR error classes (with evidence)
1. **Spurious "OR" isozymes.** These create false bypasses and false non-essential calls. This was the main mechanism behind Human2's gains [2]. Isoenzyme GPR mapping was also the key error source for E. coli iML1515 [13] [R].
2. **AND/OR confusion for complexes.** Subunits are coded as isozymes or vice versa. Human1 used CORUM to fix complexes [1]. Repository releases keep correcting isoenzymes (e.g. HADHA vs ECHS1 assignments, v1.19.0) [L].
3. **Non-metabolic or wrong genes attached to reactions.** Human2 removed 203 genes [2]. v2.0.0 removed non-transporter genes from dozens of transport reactions [L].
4. **Missing genes and gap-filled gene-less reactions.** Gene-less reactions cannot be knocked out and can carry bypass flux. v2.1.0 adds a `spontaneous` flag [L].
5. **Localisation errors** create or remove compartment-specific routes [2].
6. **Reference-model differences** in GPRs drive differences between models extracted from different references [8].

### 6.2 Paralog buffering
- **De Kegel & Ryan 2019** [14]: across 558 lines, paralogs are less often essential than singletons. Buffering (synthetic lethality) explains **13–17%** of variable paralog essentiality, more for whole-genome-duplication pairs and pairs in complexes.
- **Dede 2020** [15]: half of constitutively expressed genes are never essential, enriched for paralogs. enCas12a dual knockout found 24 masked synthetic-lethal pairs.
- **Combinatorial screens:**
  - Thompson 2021: 1,191 pairs → 105 combinations with a fitness effect [66].
  - Parrish 2021: 12% of >2,000 paralogs are synthetic lethal in at least one context [67].
  - Ito 2021: 5,065 paralog pairs [68, abstract only].
  - Flister 2025: synthetic lethality is infrequent and varies with background [69, abstract only].
  - Predictive classifier: De Kegel 2021 [70, abstract only].
- **Lagziel 2019:** isozyme and neighbour expression correlates with metabolic gene dependency [10].
- **Implication.** For a GEM "A OR B" rule, gene A should be predicted essential in a line where B is absent or non-functional *in that line*. That requires line-specific evidence of B (protein, not just RNA). Paralog pairs measured in dual-knockout screens give direct tests of OR rules.

### 6.3 Data and tools for curation
- **Protein complexes:** CORUM 4.0 (5,204 complexes, 5,299 genes) [71] and CORUM 3.0 [72]; Complex Portal 2022 (>1,100 human complexes, CC0) [73].
- **Proteomics:**
  - CCLE proteomics: 375 lines; within-pathway correlations absent from RNA [74].
  - ProCan-DepMapSanger: 949 lines, 8,498 proteins. Proteome ≈ transcriptome for drug-response prediction; protein biomarkers of CRISPR dependencies not seen at RNA level [18].
  - Complex subunits are co-regulated post-transcriptionally, and some subunits are rate-limiting [75]. Subunit co-abundance is evidence for AND rules; absence of a paralog protein is evidence to prune an OR branch.
- **Co-essentiality modules** in DepMap assign function to genes (e.g. TMEM189 as plasmanylethanolamine desaturase) [76, S]. They are a source of GPR hypotheses, but **circular if evaluated on the same DepMap labels** (§8).
- **Tools:**
  - GPRuler (automatic GPR reconstruction from 9 databases; often more accurate than the original models after review) [77].
  - GeTPRA transcript-level rules (Recon 2M.2) [78].
  - Human2's LLM/GPT checker [2].
  - RAVEN/Human-GEM `code/GPRs` [L].
- **Curation vs constraints.** No paper we found separates the essentiality gain due to GPR curation from the gain due to new constraints. Robinson attributes the largest version-to-version gains to GPR updates [1], but medium and extraction were held fixed only implicitly. Our ablation design is new.

---

## 7. Beyond essentiality: flux states, drugs, metabolomics

- **Flux vs omics for subtypes.**
  - Flux patterns differ across cancers and can act as a "fingerprint" (3,740 models, Recon3D, 4 MEMs) [49, S]. That study did **not** compare against clustering the input omics.
  - Knockout-flux UMAP + k-means on proteome-built models identified sensitisation targets (GOT1, GPI, SLC1A5) in drug-resistant MCF7, with experimental validation [79, S, abstract only].
  - **We found no rigorous test that flux states separate subtypes better than their own inputs.** Flux states are a deterministic function of expression plus model plus medium. Any gain must exceed what equally flexible models on raw expression or proteomics achieve, under the same CV [inference].
- **Drug response.**
  - PRISM (repurposing library, 578 lines in DeepMeta's use) [80]; GDSC [81]; CTRP [82, 83].
  - PRIME models predicted responses to metabolic-enzyme inhibition [46].
  - Glucose-uptake-dominated flux states correlate with chemoresistance in NCI-60 [61].
  - DeepMeta reported that compounds targeting predicted dependencies had lower PRISM log-fold changes [6, P].
  - Expression outperforms genomics for predicting vulnerabilities [84, P]. Proteome ≈ transcriptome [18]. These are the omics baselines a flux-based predictor must beat.
  - Medium also changes drug response (uric acid → 5-FU) [55], so drug-screen media should be recorded.
- **Metabolomics.** CCLE: 225 metabolites × 928 lines [85]; asparaginase sensitivity linked to ASNS hypermethylation. We found no peer-reviewed study validating GEM-predicted flux or metabolite states against CCLE metabolomics at scale [no evidence found]. GEMCAT is a bioRxiv method that validates predicted metabolite changes [P, not read].
- **Statistical trap.** Lineage drives expression, drug response and metabolite levels. Correlations between flux and phenotype must be computed **within lineage** or with lineage as a covariate, or they will reproduce the lineage signal [inference; cf. 10].

---

## 8. Benchmarking practice

- **AUPRC under imbalance.**
  - Precision–recall is more informative than ROC when positives are rare [86]. For GEM essentiality specifically, AUPRC was recommended because it "demonstrated utility … relative to alternative accuracy metrics" [13].
  - Counterpoint: AUPRC is not universally superior and favours subgroups with higher prevalence [87, NeurIPS 2024, abstract only]. Per-lineage AUPRC therefore has to be reported as AUPRC/base rate (or ΔAUPRC vs a prevalence-matched null) when lineages differ in prevalence.
  - Report AUROC as a secondary metric.
- **Units of evaluation:**
  - **per gene** (across lines): tests context-specificity and is immune to gene-prior inflation;
  - **per line / per lineage** (across genes);
  - **pooled** (gene × line).
  - DeepMeta's 0.94 vs 0.61 shows why the per-gene number must be primary for "context-specific" claims [6].
  - Report separately on (i) all genes, (ii) genes excluding DepMap common essentials, and (iii) "selective" genes (variable across lines).
- **Baselines** (minimum set, non-mechanistic):
  - prevalence (random);
  - gene prior (fraction of *training* lines where the gene is essential);
  - gene expression (low expression → non-essential);
  - generic-model FBA.
  - Lagziel-style network-neighbour expression models [10] are the strongest published non-FBA baseline for metabolic genes.
- **Leakage and splits.**
  - Leakage is pervasive in ML-based science; Kapoor & Narayanan give eight types [88].
  - For us, the main channels are: the same gene in train and test via random gene–line splits [17]; thresholds, tasks or extraction parameters tuned on test labels [5]; co-essentiality/GPR curation derived from the evaluation screens [76]; and common-essential lists derived from the full DepMap matrix.
  - Use **leave-one-lineage-out** (or grouped K-fold by lineage), with all fitting inside training folds and a fixed seed list.
  - External test: Project Score (Sanger) lines and screens [28], noting the library and assay-length batch effects [26].
- **Uncertainty.**
  - Bootstrap over *cell lines* (and over genes for per-line metrics). Gene–line pairs are not independent.
  - Paired bootstrap or permutation for model-vs-model differences on identical gene × line sets.
  - Report CIs for ΔAUPRC.
- **Multiple testing.** Use Benjamini–Hochberg across genes, lineages and tasks, as in [10]. Pre-register the comparisons in DESIGN.md to limit the garden of forking paths.
- **Runtime.** ftINIT is >10× faster than tINIT [20]. Flux-cone sampling is heavy (GB-scale per organism) [16]. Report runtime per model and per line alongside accuracy.

---

### 8.1 Searched, not found (2026-10-26)
Absence claims in this document are absence of evidence from our searches, not proof of absence. The searches were made by the literature agent (this document) and by the prior-work reviewer (notes/design_review_r1/prior_work.md) on 2026-10-26, in PubMed, bioRxiv and general web search. **The exact query strings were not logged**; the literature agent must add them here before M1. Neither search found any of:
- per-gene or per-lineage AUPRC of GEM essentiality predictions on DepMap at full scale with lineage hold-out;
- an ablation that separates the gain from GPR curation from the gain from new constraints;
- a head-to-head of flux-state features against their own input omics, with matched model capacity;
- a released, DepMap-keyed resource of per-line Human1/Human2 SBML models other than the Human-GEM repository's own `code/DepMapGeneEss` pipeline (which produces models but does not distribute a versioned per-line resource). Metabolic Atlas (metabolicatlas.org) was **not** checked for hosted cell-line models; to do before REPORT.md.

Flagged by the prior-work reviewer but **not verified** (not cited): bioRxiv 10.1101/593277, "On the inconsistent treatment of gene-protein-reaction rules in context-specific metabolic models" (authors and status unknown).

## 9. Implications for our design (actionable)

### (a) Models and labels
1. **Pin three model versions and report each:** Human1 = **v1.19.0** (last v1.x; decision D0), with **v1.3.0** (the published Human1) as a sensitivity check; Human2 = **v2.1.1** (local), noting that the Human2 paper describes v2.0.0. Version differences inside "Human1" are large (3,625 → 2,887 genes).
2. **Primary label:** Chronos gene effect from `CRISPRGeneEffect` (continuous, for threshold-free AUPRC over ranks of predicted growth ratio), plus a binary label from `CRISPRGeneDependency` P > 0.5 (fetch it; missing locally). Sensitivity labels: gene effect < −0.5 and < −1.
   - Exclude genes not screened.
   - Report results with and without the 419 common essentials.
   - External validation on Project Score.

### (b) Metrics and baselines
3. **Primary metric: per-gene AUPRC** (context-specificity), plus per-lineage AUPRC/base rate and pooled AUPRC, all with line-bootstrap CIs.
4. **Mandatory baselines:** random/prevalence, gene prior (from training lineages only), expression-only, generic Human1 and Human2 FBA. The repository's ftINIT pooled AUPRC of 0.29 against a base rate of 0.16 is a sanity anchor, not a target.

### (c) Medium, extraction, tasks, GPRs
5. **Medium:** derive each screen's medium from `ScreenSequenceMap → ModelCondition.GrowthMedia`. Map formulations to compositions with an explicit, versioned serum model. Stratify or exclude serum-free/organoid conditions. Quantify medium attribution by nested variants on identical genes and lines, using medium-sensitive genes as positive controls.
6. **Extraction:** run ≥3 MEM families × ≥2 threshold rules × ± task protection, and test "method > lineage" by variance partitioning with a permuted-expression negative control. Use ftINIT bootstraps for stability. Tune thresholds only inside training lineages.
7. **Tasks:** fix the candidate set a priori (Human-GEM Essential by category, CellFie consensus by system) before any label is seen. Score viability tasks and capability tasks separately, since capability-only calls are known false positives.
8. **GPRs:** target curation at reactions whose OR/AND logic flips knockout predictions. Use evidence independent of the evaluation labels: CORUM/Complex Portal, ProCan/CCLE protein presence and subunit co-abundance, and paralog dual-knockout screens.
   - Treat DepMap co-essentiality as training-fold-only evidence.
   - Report gains from curation and gains from constraints as separate ablations.

### (d) Hybrids, flux-information tests, reporting
9. **Hybrid ML/FBA** (e.g. knockout-flux or flux-cone features) is promising but unproven in human cells. Evaluate it only with lineage hold-out and gene-grouped CV, against the same features without the GEM.
10. **Flux-information tests:** compare flux-derived features against their own input omics with matched model capacity, within-lineage, using PRISM/GDSC/CTRP restricted to metabolic-target drugs and CCLE metabolomics. A null result is a legitimate outcome.
11. **Report runtime** per model build and per knockout sweep. Fix seeds for any MILP or sampling step, since alternate optima are a known source of variance [9, 48].

---

## 10. Papers we could NOT access in full that matter (for the PI to request)

| Title | DOI | Why it matters |
|---|---|---|
| Luo et al., *Reconstruction of human metabolic models with LLMs* — **SI Appendix and Datasets S01–S09; numeric values of Fig. 2D/E** | 10.1073/pnas.2516511123 | Exact MCC per Human-GEM version, the curated GPR change lists, DepMap/Project Score benchmark details. Main text was read via PMC; figure values and SI were not available. |
| Robinson et al., *An atlas of human metabolism* — **Supplementary Figs S8–S11, Table S3** | 10.1126/scisignal.aaz1482 | MCC values for 621 DepMap lines, threshold sensitivity, biomass-only vs 57 tasks, per-version effects. Main text read via PMC. |
| Wu et al. (DeepMeta), *Cell Reports* 2025 final version (PDF now in `papers/`, not yet read) | 10.1016/j.celrep.2025.115945 | Only the bioRxiv preprint was read (10.1101/2025.02.26.640327); the final evaluation and AUPRC numbers may differ. |
| Behan et al., *Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens*, Nature 2019 | 10.1038/s41586-019-1103-9 | Project Score methods and labels (BAGEL, CRISPRcleanR) for our external test set; not in PMC. |
| Hart et al., *High-resolution CRISPR screens reveal fitness genes…*, Cell 2015 | 10.1016/j.cell.2015.11.015 | Label definitions behind the Human-GEM repository benchmark (BF thresholds); not in PMC. |
| Ito et al., *Paralog knockout profiling…*, Nat Genet 2021 | 10.1038/s41588-021-00967-z | 5,065 paralog pairs in a dual-knockout screen: direct test data for OR rules. Abstract only. |
| De Kegel et al., *Comprehensive prediction of robust synthetic lethality between paralog pairs*, Cell Syst 2021 | 10.1016/j.cels.2021.08.006 | Paralog-SL predictor and features usable for GPR prioritisation. Abstract only. |
| Chen et al., *GECKO Toolbox 3.0*, Nat Protoc 2024 | 10.1038/s41596-023-00931-7 | Needed if we build ec-models from ProCan proteomics. Abstract only. |
| Lim et al., *Genome-scale knockout simulation and clustering … drug sensitization targets*, PNAS 2025 | 10.1073/pnas.2425384122 | Closest precedent for "knockout flux profiles" analysis. Abstract only (PMC12232641 exists; not read). |
| Agren et al., *tINIT*, Mol Syst Biol 2014; Shlomi et al., *iMAT*, Nat Biotechnol 2008 | 10.1002/msb.145122; 10.1038/nbt.1487 | Original algorithm details for faithful reimplementation. Not read in full. |
| McDermott et al., *A closer look at AUROC and AUPRC under class imbalance*, NeurIPS 2024 (arXiv 2401.06091) | — (arXiv:2401.06091) | Formal caveats on per-subgroup AUPRC; abstract only. |

Also unread but available on PMC if needed: Bernstein 2023 [13], Moscardó García 2021 [62], Dai 2019 [65], Gonçalves 2022 [18] full texts. We only used their abstracts.

---

## References
(PMID where available; all DOIs checked against PubMed or the publisher record.)

1. Robinson JL et al. An atlas of human metabolism. *Sci Signal* 2020;13:eaaz1482. doi:10.1126/scisignal.aaz1482. PMID 32209698.
2. Luo J et al. Reconstruction of human metabolic models with large language models. *PNAS* 2026;123(15):e2516511123. doi:10.1073/pnas.2516511123. PMID 41950094.
3. Huijer CAM, Jiao X, Chen Y, Yu R. Improved flux profiling in genome-scale modeling of human cell metabolism. *Cell Rep Methods* 2026;6(1):101275. doi:10.1016/j.crmeth.2025.101275. PMID 41529686. (Added 2026-10-26, design review r1; abstract read.)
4. Cheng CT, Lai JM, Chang PM, Hong YR, Huang CF, Wang FS. Identifying essential genes in genome-scale metabolic models of consensus molecular subtypes of colorectal cancer. *PLoS One* 2023;18:e0286032. doi:10.1371/journal.pone.0286032. PMID 37205704. (Added 2026-10-26; abstract read.)
5. Vieira V, Ferreira J, Rocha M. A pipeline for the reconstruction and evaluation of context-specific human metabolic models at a large-scale (Troppo). *PLoS Comput Biol* 2022;18:e1009294. doi:10.1371/journal.pcbi.1009294. PMID 35749559. (733 CCLE lines; MCF7 reference used for parameter choice — per PubMed abstract, checked 2026-10-26.)
6. Wu T, Zhao X, …, Liu XS. Precise metabolic dependencies of cancer through deep learning and validations (DeepMeta). *Cell Rep* 2025;44(7):115945. doi:10.1016/j.celrep.2025.115945. PMID 40616842 (confirmed in PubMed 2026-10-26; first author corrected from "Tao W"). Preprint: bioRxiv 2025, doi:10.1101/2025.02.26.640327 (the version read).
7. Opdam S et al. A systematic evaluation of methods for tailoring genome-scale metabolic models. *Cell Syst* 2017;4:318–329. doi:10.1016/j.cels.2017.01.010. PMID 28215528.
8. Richelle A et al. Increasing consensus of context-specific metabolic models by integrating data-inferred cell functions. *PLoS Comput Biol* 2019;15:e1006867. doi:10.1371/journal.pcbi.1006867. PMID 30986217.
9. Gopalakrishnan S et al. Guidelines for extracting biologically relevant context-specific metabolic models using gene expression data. *Metab Eng* 2023;75:181. doi:10.1016/j.ymben.2022.12.003. PMID 36566974.
10. Lagziel S, Lee WD, Shlomi T. Inferring cancer dependencies on metabolic genes from large-scale genetic screens. *BMC Biol* 2019;17:37. doi:10.1186/s12915-019-0654-4. PMID 31039782.
11. Rossiter NJ et al. CRISPR screens in physiologic medium reveal conditionally essential genes in human cells. *Cell Metab* 2021;33:1248–1263. doi:10.1016/j.cmet.2021.02.005. PMID 33651980.
12. Neiswender JV et al. A dependency map enhanced with next-generation 3D cancer models. *Nature* 2026;657:789–799. doi:10.1038/s41586-026-10843-7. PMID 42557315.
13. Bernstein DB, Akkas B, Price MN, Arkin AP. Evaluating E. coli genome-scale metabolic model accuracy with high-throughput mutant fitness data. *Mol Syst Biol* 2023. doi:10.15252/msb.202311566. PMID 37888487.
14. De Kegel B, Ryan CJ. Paralog buffering contributes to the variable essentiality of genes in cancer cell lines. *PLoS Genet* 2019;15:e1008466. doi:10.1371/journal.pgen.1008466. PMID 31652272.
15. Dede M et al. Multiplex enCas12a screens detect functional buffering among paralogs otherwise masked in monogenic Cas9 knockout screens. *Genome Biol* 2020;21:262. doi:10.1186/s13059-020-02173-2. PMID 33059726.
16. Merzbacher C, Mac Aodha O, Oyarzún DA. Accurate prediction of gene deletion phenotypes with Flux Cone Learning. *Nat Commun* 2025;16:8492. doi:10.1038/s41467-025-63436-9. PMID 41006282.
17. Kim et al. Integrating genome-scale metabolic modeling with machine learning improves gene essentiality prediction in triple-negative breast cancer. *Int J Mol Sci* 2026;27:5059. doi:10.3390/ijms27115059. PMID 42278581.
18. Gonçalves E et al. Pan-cancer proteomic map of 949 human cell lines. *Cancer Cell* 2022. doi:10.1016/j.ccell.2022.06.010. PMID 35839778.
19. Tobalina L et al. Assessment of FBA based gene essentiality analysis in cancer with a fast context-specific network reconstruction method. *PLoS One* 2016;11:e0154583. doi:10.1371/journal.pone.0154583. PMID 27145226.
20. Gustafsson J et al. Generation and analysis of context-specific genome-scale metabolic models derived from single-cell RNA-Seq data (ftINIT). *PNAS* 2023;120:e2217868120. doi:10.1073/pnas.2217868120. PMID 36719923.
21. Gatto F, Miess H, Schulze A et al. Flux balance analysis predicts essential genes in clear cell renal cell carcinoma metabolism. *Sci Rep* 2015;5:10738. doi:10.1038/srep10738. PMID 26040780 (abstract-level use only).
22. Dempster JM et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biol* 2021. doi:10.1186/s13059-021-02540-7. PMID 34930405.
23. Hart T et al. Measuring error rates in genomic perturbation screens: gold standards for human functional genomics. *Mol Syst Biol* 2014;10:733. doi:10.15252/msb.20145216. PMID 24987113.
24. Hart T et al. Evaluation and design of genome-wide CRISPR/SpCas9 knockout screens. *G3* 2017;7:2719. doi:10.1534/g3.117.041277. PMID 28655737. (Hart 2015 *Cell*, the TKO screen source, is now ref 91.)
25. Vinceti A et al. CoRe: a robustly benchmarked R package for identifying core-fitness genes in genome-wide pooled CRISPR-Cas9 screens. *BMC Genomics* 2021. doi:10.1186/s12864-021-08129-5. PMID 34789150.
26. Dempster JM et al. Agreement between two large pan-cancer CRISPR-Cas9 gene dependency data sets. *Nat Commun* 2019. doi:10.1038/s41467-019-13805-y. PMID 31862961.
27. Pacini C et al. Integrated cross-study datasets of genetic dependencies in cancer. *Nat Commun* 2021. doi:10.1038/s41467-021-21898-7. PMID 33712601.
28. Behan FM et al. Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens. *Nature* 2019;568:511. doi:10.1038/s41586-019-1103-9. PMID 30971826.
29. Aguirre AJ et al. Genomic copy number dictates a gene-independent cell response to CRISPR/Cas9 targeting. *Cancer Discov* 2016;6:914. doi:10.1158/2159-8290.CD-16-0154. PMID 27260156.
30. Munoz DM et al. CRISPR screens provide a comprehensive assessment of cancer vulnerabilities but generate false-positive hits for highly amplified genomic regions. *Cancer Discov* 2016;6:900. doi:10.1158/2159-8290.CD-16-0178. PMID 27260157.
31. Gonçalves E et al. Structural rearrangements generate cell-specific, gene-independent CRISPR-Cas9 loss of fitness effects. *Genome Biol* 2019;20:27. doi:10.1186/s13059-019-1637-z. PMID 30722791.
32. Iorio F et al. Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting (CRISPRcleanR). *BMC Genomics* 2018;19:604. doi:10.1186/s12864-018-4989-y. PMID 30103702.
33. Hasibi R, Michoel T, Oyarzún DA. Integration of graph neural networks and genome-scale metabolic models for predicting gene essentiality (FlowGAT). *NPJ Syst Biol Appl* 2024. doi:10.1038/s41540-024-00348-2. PMID 38448436.
34. Becker SA, Palsson BO. Context-specific metabolic networks are consistent with experiments (GIMME). *PLoS Comput Biol* 2008;4:e1000082. doi:10.1371/journal.pcbi.1000082. PMID 18483554.
35. Shlomi T et al. Network-based prediction of human tissue-specific metabolism (iMAT). *Nat Biotechnol* 2008;26:1003. doi:10.1038/nbt.1487. PMID 18711341.
36. Agren R et al. Reconstruction of genome-scale active metabolic networks for 69 human cell types and 16 cancer types using INIT. *PLoS Comput Biol* 2012;8:e1002518. doi:10.1371/journal.pcbi.1002518. PMID 22615553.
37. Agren R et al. Identification of anticancer drugs for hepatocellular carcinoma through personalized genome-scale metabolic modeling (tINIT). *Mol Syst Biol* 2014;10:721. doi:10.1002/msb.145122. PMID 24646661.
38. Jerby L, Shlomi T, Ruppin E. Computational reconstruction of tissue-specific metabolic models: application to human liver metabolism (MBA). *Mol Syst Biol* 2010;6:401. doi:10.1038/msb.2010.56. PMID 20823844.
39. Wang Y, Eddy JA, Price ND. Reconstruction of genome-scale metabolic models for 126 human tissues using mCADRE. *BMC Syst Biol* 2012;6:153. doi:10.1186/1752-0509-6-153. PMID 23234303.
40. Vlassis N, Pacheco MP, Sauter T. Fast reconstruction of compact context-specific metabolic network models (FASTCORE). *PLoS Comput Biol* 2014;10:e1003424. doi:10.1371/journal.pcbi.1003424. PMID 24453953.
41. Pacheco MP et al. Integrated metabolic modelling reveals cell-type specific epigenetic control points of the macrophage metabolic network (FASTCORMICS). *BMC Genomics* 2015;16:809. doi:10.1186/s12864-015-1984-4. PMID 26480823.
42. Pacheco MP et al. Identifying and targeting cancer-specific metabolism with network-based drug target prediction (rFASTCORMICS). *EBioMedicine* 2019;43:98. doi:10.1016/j.ebiom.2019.04.046. PMID 31126892.
43. Colijn C et al. Interpreting expression data with metabolic flux models: predicting Mycobacterium tuberculosis mycolic acid production (E-Flux). *PLoS Comput Biol* 2009;5:e1000489. doi:10.1371/journal.pcbi.1000489. PMID 19714220.
44. Sánchez BJ et al. Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic constraints (GECKO). *Mol Syst Biol* 2017;13:935. doi:10.15252/msb.20167411. PMID 28779005. GECKO 2.0: Domenzain I et al. *Nat Commun* 2022;13:3766. doi:10.1038/s41467-022-31421-1. PMID 35773252.
45. Chen Y et al. Reconstruction, simulation and analysis of enzyme-constrained metabolic models using GECKO Toolbox 3.0. *Nat Protoc* 2024;19:629. doi:10.1038/s41596-023-00931-7. PMID 38238583.
46. Yizhak K et al. Phenotype-based cell-specific metabolic modeling reveals metabolic liabilities of cancer (PRIME). *eLife* 2014;3:e03641. doi:10.7554/eLife.03641. PMID 25415239.
47. Joshi CJ et al. StanDep: capturing transcriptomic variability improves context-specific metabolic models. *PLoS Comput Biol* 2020;16:e1007764. doi:10.1371/journal.pcbi.1007764. PMID 32396573.
48. Rodríguez-Mier P et al. DEXOM: diversity-based enumeration of optimal context-specific metabolic networks. *PLoS Comput Biol* 2021;17:e1008730. doi:10.1371/journal.pcbi.1008730. PMID 33571201.
49. Jalili et al. Exploring the metabolic heterogeneity of cancers: a benchmark study of context-specific models (GEMbench). *J Pers Med* 2021;11:496. doi:10.3390/jpm11060496. PMID 34205912.
50. Richelle A et al. Model-based assessment of mammalian cell metabolic functionalities using omics data (CellFie). *Cell Rep Methods* 2021;1:100040. doi:10.1016/j.crmeth.2021.100040. PMID 34761247.
51. Wagner A et al. Metabolic modeling of single Th17 cells reveals regulators of autoimmunity (Compass). *Cell* 2021. doi:10.1016/j.cell.2021.05.045. PMID 34216539.
52. Alghamdi N et al. A graph neural network model to estimate cell-wise metabolic flux using single-cell RNA-seq data (scFEA). *Genome Res* 2021. doi:10.1101/gr.271205.120. PMID 34301623.
53. Huang Y et al. Characterizing cancer metabolism from bulk and single-cell RNA-seq data using METAFlux. *Nat Commun* 2023. doi:10.1038/s41467-023-40457-w. PMID 37573313.
54. Richelle A et al. Assessing key decisions for transcriptomic data integration in biochemical networks. *PLoS Comput Biol* 2019;15:e1007185. doi:10.1371/journal.pcbi.1007185. PMID 31323017.
55. Cantor JR et al. Physiologic medium rewires cellular metabolism and reveals uric acid as an endogenous inhibitor of UMP synthase (HPLM). *Cell* 2017. doi:10.1016/j.cell.2017.03.023. PMID 28388410.
56. Vande Voorde J et al. Improving the metabolic fidelity of cancer models with a physiological cell culture medium (Plasmax). *Sci Adv* 2019;5:eaau7314. doi:10.1126/sciadv.aau7314. PMID 30613774.
57. Ackermann T, Tardito S. Cell culture medium formulation and its implications in cancer metabolism. *Trends Cancer* 2019;5:329. doi:10.1016/j.trecan.2019.05.004. PMID 31208694.
58. Bairoch A. The Cellosaurus, a cell-line knowledge resource. *J Biomol Tech* 2018;29:25. doi:10.7171/jbt.18-2902-002. PMID 29805321.
59. Jain M et al. Metabolite profiling identifies a key role for glycine in rapid cancer cell proliferation. *Science* 2012. doi:10.1126/science.1218595. PMID 22628656.
60. Dolfi SC et al. The metabolic demands of cancer cells are coupled to their size and protein synthesis rates. *Cancer Metab* 2013;1:20. doi:10.1186/2049-3002-1-20. PMID 24279929.
61. Zielinski DC et al. Systems biology analysis of drivers underlying hallmarks of cancer cell metabolism. *Sci Rep* 2017;7:41241. doi:10.1038/srep41241. PMID 28120890.
62. Moscardó García M et al. Importance of the biomass formulation for cancer metabolic modeling and drug prediction. *iScience* 2021;24:103110. doi:10.1016/j.isci.2021.103110. PMID 34622163.
63. Lewis NE et al. Omic data from evolved E. coli are consistent with computed optimal growth from genome-scale models (pFBA). *Mol Syst Biol* 2010;6:390. doi:10.1038/msb.2010.47. PMID 20664636.
64. Segrè D, Vitkup D, Church GM. Analysis of optimality in natural and perturbed metabolic networks (MOMA). *PNAS* 2002;99:15112. doi:10.1073/pnas.232349399. PMID 12415116.
65. Dai Z, Yang S, Xu L et al. Identification of cancer-associated metabolic vulnerabilities by modeling multi-objective optimality in metabolism. *Cell Commun Signal* 2019;17:124. doi:10.1186/s12964-019-0439-y. PMID 31601242.
66. Thompson NA et al. Combinatorial CRISPR screen identifies fitness effects of gene paralogues. *Nat Commun* 2021;12:1302. doi:10.1038/s41467-021-21478-9. PMID 33637726.
67. Parrish PCR et al. Discovery of synthetic lethal and tumor suppressor paralog pairs in the human genome. *Cell Rep* 2021;36:109597. doi:10.1016/j.celrep.2021.109597. PMID 34469736.
68. Ito T et al. Paralog knockout profiling identifies DUSP4 and DUSP6 as a digenic dependence in MAPK pathway-driven cancers. *Nat Genet* 2021;53:1664. doi:10.1038/s41588-021-00967-z. PMID 34857952.
69. Flister et al. Genetic mapping and predictive modeling of paralog synthetic lethality. *Cell Rep* 2025. doi:10.1016/j.celrep.2025.116512. PMID 41205176.
70. De Kegel B et al. Comprehensive prediction of robust synthetic lethality between paralog pairs in cancer cell lines. *Cell Syst* 2021;12:1144. doi:10.1016/j.cels.2021.08.006. PMID 34529928.
71. Tsitsiridis G et al. CORUM: the comprehensive resource of mammalian protein complexes—2022. *Nucleic Acids Res* 2023. doi:10.1093/nar/gkac1015. PMID 36382402.
72. Giurgiu M et al. CORUM: the comprehensive resource of mammalian protein complexes—2019. *Nucleic Acids Res* 2019. doi:10.1093/nar/gky973. PMID 30357367.
73. Meldal BHM et al. Complex Portal 2022: new curation frontiers. *Nucleic Acids Res* 2022;50:D578. doi:10.1093/nar/gkab991. PMID 34718729.
74. Nusinow DP et al. Quantitative proteomics of the Cancer Cell Line Encyclopedia. *Cell* 2020. doi:10.1016/j.cell.2019.12.023. PMID 31978347.
75. Gonçalves E et al. Widespread post-transcriptional attenuation of genomic copy-number variation in cancer. *Cell Syst* 2017;5:386. doi:10.1016/j.cels.2017.08.013. PMID 29032074.
76. Wainberg M et al. A genome-wide atlas of co-essential modules assigns function to uncharacterized genes. *Nat Genet* 2021. doi:10.1038/s41588-021-00840-z. PMID 33859415.
77. Di Filippo M, Damiani C, Pescini D. GPRuler: metabolic gene-protein-reaction rules automatic reconstruction. *PLoS Comput Biol* 2021;17:e1009550. doi:10.1371/journal.pcbi.1009550. PMID 34748537.
78. Ryu JY, Kim HU, Lee SY. Framework and resource for more than 11,000 gene-transcript-protein-reaction associations in human metabolism. *PNAS* 2017. doi:10.1073/pnas.1713050114. PMID 29078384.
79. Lim et al. Genome-scale knockout simulation and clustering analysis of drug-resistant breast cancer cells reveal drug sensitization targets. *PNAS* 2025;122:e2425384122. doi:10.1073/pnas.2425384122. PMID 40560621.
80. Corsello SM et al. Discovering the anticancer potential of non-oncology drugs by systematic viability profiling (PRISM). *Nat Cancer* 2020;1:235. doi:10.1038/s43018-019-0018-6. PMID 32613204.
81. Iorio F et al. A landscape of pharmacogenomic interactions in cancer (GDSC). *Cell* 2016;166:740. doi:10.1016/j.cell.2016.06.017. PMID 27397505.
82. Rees MG et al. Correlating chemical sensitivity and basal gene expression reveals mechanism of action (CTRP). *Nat Chem Biol* 2016;12:109. doi:10.1038/nchembio.1986. PMID 26656090.
83. Seashore-Ludlow B et al. Harnessing connectivity in a large-scale small-molecule sensitivity dataset (CTRP v2). *Cancer Discov* 2015;5:1210. doi:10.1158/2159-8290.CD-15-0235. PMID 26482930.
84. Dempster JM et al. Gene expression has more power for predicting in vitro cancer cell vulnerabilities than genomics. bioRxiv 2020. doi:10.1101/2020.02.21.959627. [preprint]
85. Li H et al. The landscape of cancer cell line metabolism. *Nat Med* 2019;25:850. doi:10.1038/s41591-019-0404-8. PMID 31068703.
86. Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. *PLoS One* 2015;10:e0118432. doi:10.1371/journal.pone.0118432. PMID 25738806.
87. McDermott MBA et al. A closer look at AUROC and AUPRC under class imbalance. NeurIPS 2024; arXiv:2401.06091. [abstract only]
88. Kapoor S, Narayanan A. Leakage and the reproducibility crisis in machine-learning-based science. *Patterns* 2023;4:100804. doi:10.1016/j.patter.2023.100804. PMID 37720327.
89. Guil F, García R, García JM. Adding metabolic tasks to human GEM models to improve the study of gene targets and their associated toxicities. *Sci Rep* 2024;14:17265. doi:10.1038/s41598-024-68073-8. PMID 39068208. (Added 2026-10-26; abstract read.)
90. Barrena N, Valcárcel LV, Olaverri-Mendizabal D, Apaolaza I, Planes FJ. Synthetic lethality in large-scale integrated metabolic and regulatory network models of human cells. *NPJ Syst Biol Appl* 2023;9:32. doi:10.1038/s41540-023-00296-3. PMID 37454223. (Added 2026-10-26; abstract read; Human1 + regulatory networks, gMCS essentiality vs in-vitro silencing data.)
91. Hart T et al. High-resolution CRISPR screens reveal fitness genes and genotype-specific cancer liabilities. *Cell* 2015;163:1515. doi:10.1016/j.cell.2015.11.015. PMID 26627737. (Split out of ref 24, 2026-10-26.)
92. McSherry F, Najork M. Computing information retrieval performance measures efficiently in the presence of tied scores. *ECIR 2008*, LNCS 4956. (Tie-aware expected AP used in DESIGN.md §6; computer-science venue, not in PubMed; title/venue as given by the statistics reviewer — verify against the Springer record before M1.)

Other DepMap resources cited: Tsherniak A et al. Defining a cancer dependency map. *Cell* 2017;170:564. doi:10.1016/j.cell.2017.06.010. PMID 28753430. Meyers RM et al. (CERES). *Nat Genet* 2017;49:1779. doi:10.1038/ng.3984. PMID 29083409. Ghandi M et al. (CCLE 2019). *Nature* 2019;569:503. doi:10.1038/s41586-019-1186-3. PMID 31068700. Thiele I et al. (Recon 2). *Nat Biotechnol* 2013;31:419. doi:10.1038/nbt.2488. PMID 23455439. Brunk E et al. (Recon3D). *Nat Biotechnol* 2018;36:272. doi:10.1038/nbt.4072. PMID 29457794.

Non-peer-reviewed sources used (and flagged in text): DepMap forum (probability-of-dependency usage; common-essential list definitions); Human-GEM GitHub release notes and `data/testResults/` (v1.x/v2.x counts, errata, graded essentiality); FluxGAT (arXiv 2403.18666, not read).
