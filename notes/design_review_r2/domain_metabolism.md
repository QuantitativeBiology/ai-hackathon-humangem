# Design review r2: cancer metabolism / CRISPR-screen domain
Reviewed DESIGN.md v1 (W2 §ladder, W3 B-add, W4 genomic constraint and backgrounds, Appendix A, §16). No label values were used.

## 1. Round-1 must-fixes
| R1 item | Status | Comment |
|---|---|---|
| M1 serum units, no nucleosides, fraction scheme, calibration | **Mostly addressed, one defect (N1 below)** | Nucleosides/hypoxanthine removed from serum. Log-fraction f ∈ {0.001..1} of a biomass requirement. M3 calibration gate 0.01–0.06 /h. Micronutrient floor. All good. The definition of req_c breaks for substitutable nutrients (N1). |
| M2 transporter gate, cystine | Addressed | M1g (Asp/Glu on SLC1A1/2/3, TPM>=1), cystine with cysteine closed unless 2-ME. Fine. |
| M3 sterol control | Addressed in the controls | Sterol branch separated from mevalonate core, B-add added. The ordering consequence is N2. |
| M4 gene-loss drop rule | Addressed in principle | One-at-a-time drop and log. Two defects: the definition of "lost" (N3) and order-dependence (see N3). |
| M5 backgrounds | Addressed | Metabolism-linked table with counts, a-priori min n, pooled within-lineage contrast, KRAS as negative control, MTAP as expected failure, IDH dropped, ENO1 as in-silico only. Good. |

## 2. New must-fix

**N1. req_c is zero for exactly the components the serum model needs (§W2 "Requirement").**
req_c is defined as "the minimal uptake of c at WT growth on the base medium", via FVA-min. That is the right definition for true auxotrophies (the micronutrient floor: linoleate, tocopherol, B12, lipoate). It is zero for every substitutable component: cholesterol, palmitate, stearate, oleate, carnitine, ethanolamine (the cell can synthesise them). Then f × req = 0 and the serum model supplies no lipids or cholesterol. The sterol-branch control and the S1/S2 contrast would silently be no-ops, and the log-fraction sweep would not vary anything for these components.
Fix: define req_c for substitutable components as the uptake carried in a pFBA solution at >= 99% of WT growth with c open at 1000 and all other inputs as in the base medium (i.e. the amount the biomass would take up if it were free). For true auxotrophs keep FVA-min. State which rule applies per component in `serum_model.tsv` and unit-test that no serum component has req = 0 unless flagged.

**N2. Heterozygous damaging mutations must not trigger a gene knockout (§W4 genomic constraint).**
"Carries a damaging mutation (matrix value >= 1)" probably includes heterozygous calls (DepMap damaging matrix: 1 = heterozygous, 2 = homozygous; the design says the codes are verified at M1). A single-allele truncation in a metabolic enzyme is usually tolerated and often a passenger or has a wild-type allele; knocking it out in the model would create many false losses (and many contradictions to log). The rule should be: homozygous call (value 2), or heterozygous with LOH/relCN < 0.5 on that segment, or heterozygous with TPM < 1. Pre-register this, and keep ">= 1" only as a sensitivity. Also make the drop procedure order-independent: test each loss alone, then jointly, and when the joint set kills growth drop the minimal subset by a fixed order (e.g. descending expression), logged.

## 3. Should-fix (not blocking)
- **S-a. Interaction of W2 and B-add.** The mevalonate core (HMGCR, MVK, FDPS...) will turn non-essential as soon as cholesterol is supplied (after N1 is fixed), and it is restored only under W3 B-add. W2's TG choice of serum variant is made without B-add, so a serum variant could be penalised or rewarded for this artefact. Report the mevalonate-core panel as a diagnostic for every W2 step, and re-check the W2 choice once B-add is fixed (the M6 factorial already crosses W2 and W3; state this explicitly).
- **S-b. Purine genes in Appendix A.** De novo purine genes (PPAT, GART, PFAS, PAICS, ADSL, ATIC, IMPDH2) are listed only in the medium-dependent group. Say in the panel that the expectation is essential in RPMI/DMEM + serum, and relaxed in F-12 / HPLM. Same for CAD/DHODH/UMPS under F-12 thymidine (the pyrimidine ring is not rescued by thymidine alone; only uridine/cytidine rescues).
- **S-c.** The Asn uptake gate is not defined, but ASNS is a main medium positive control. Not needed for freeze, but note that SLC1A5/SLC38A5-type uptake is not gated, so the ASNS control tests medium only.

## 4. Open questions in my area
- **R2-5 (micronutrient floor in serum-free media).** Keep the floor for all media, flagged. The floor represents carrier-bound or trace inputs (lipoate, B12, tocopherol, retinoate) that exist in every real medium, including B27/BSA supplements and the cell's own stores; without it the model cannot grow and the line is lost, not an informative "non-growth". Do not create a no-floor stratum. Instead: (i) the floor is fixed at a low f (same f as the serum sweep, minimum 0.01), (ii) serum-free lines are a reporting stratum anyway (already planned), and (iii) the floor is excluded from the "serum-free vs serum" lipid contrast by being identical in both. Report the list of floor components per parent.
- **R2-6 (transporter gate: medium or context).** Attribute it to **context** (an expression-informed constraint), not to medium. Reason: it uses line expression, so it changes between lines on the same medium, and it would inflate "medium alone" if counted as medium. Keep M1 (availability) as the medium lever; report M1g as a separate step "M1g − M1" in the Shapley shares with the context constraints. Do not use the gate in the "medium alone on generic parent" estimate (the generic parent has no line expression).
- **Rejected hypoxia (O4) and TIFM (O3).** Agree with both rejections: no labels or composition on disk, and hypoxia is a culture-condition change DepMap does not document per line. HPLM/Plasmax cover the physiological counterfactual. Keep one sentence in the report that these were not tested.

## 5. Verdict
Revise minimally, then freeze. Two small must-fixes (N1 requirement definition; N2 definition of "lost" gene and order-independent drop). Everything else from round 1 is adequately addressed.
