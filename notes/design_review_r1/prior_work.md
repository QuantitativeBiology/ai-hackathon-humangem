# Design review r1: prior work and reuse (reviewer: prior-work angle)

Scope: DESIGN.md v0.1 against LITERATURE.md, PROMPT.md, local Human-GEM v2.1.1 repo and venv. No essentiality performance was computed.

## 1. Citation spot-check (LITERATURE.md refs that DESIGN.md relies on)

Checked via PubMed (metadata/citation lookup) and web fetches. "OK" = DOI/PMID resolves and title/journal/year match what LITERATURE says.

| Ref | Used in DESIGN for | Result |
|---|---|---|
| [1] Robinson 2020, PMID 32209698 | H1 baseline, 57-task essentiality | OK (PMID resolves to Sci Signal 13:eaaz1482) |
| [2] Luo 2026, PNAS 123(15) e2516511123, PMID 41950094 | H1.2, W5 positive control | OK. PubMed abstract says nothing about essentiality; the MCC/OR-isozyme claims rest on the main text read via PMC, so they cannot be re-verified from the abstract. Acceptable, but DESIGN H1.2 should cite "main text, Fig 2D" not the paper generally |
| [5] Vieira 2022 (Troppo), PMID 35749559 | context for W4 thresholds, best-of-grid MCC 0.44 | DOI OK. DISCREPANCY: the PubMed abstract says "over 6000 models ... for **733** CCLE cell lines" and that **MCF7 was used as the reference to choose the best parameter combinations**. LITERATURE says "739 Achilles lines" and "best-of-grid, optimistic". The 739 vs 733 is unverified, and the "maximum over the grid" description may be wrong (the grid was tuned on one line and then applied). Fix wording and re-read the full text |
| [20] Gustafsson 2023 ftINIT, PMID 36719923 | F1, Hart reproduction | OK. Web summary confirms 891 DepMap lines, MCC, ftINIT ~ tINIT |
| [10] Lagziel 2019, PMID 31039782 | B6, W2 positive control (ASNS) | OK. The abstract mentions "162 metabolic genes" with predictive models; LITERATURE quotes "32% of variable metabolic genes" and "26% medium-associated". Plausibly consistent but state which denominator |
| [11] Rossiter 2021, PMID 33651980 | W2 HPLM check | OK |
| [12] Neiswender 2026 Nature, PMID 42557315 | W2 serum, NextGen | OK (published 2026-08-05; volume/pages 657:789 not verified) |
| [13] Bernstein 2023, PMID 37888487 | AUPRC for GEMs, GPR error | OK (resolves) |
| [16] Merzbacher 2025 Flux Cone Learning, PMID 41006282 | W7 H-c | OK (resolves) |
| [17] Kim 2026 IJMS, PMID 42278581 | hybrid ML-FBA precedent | OK; web summary confirms 50 breast lines, Recon 2M.2, tINIT, MOMA, sensitivity 0.37 to 0.55, DepMap labels |
| [9] Gopalakrishnan 2023, PMID 36566974; [50] CellFie PMID 34761247; [69] Flister PMID 41205176 | W4 stability, T5, W5 | resolve to the stated first author/journal/year |
| [79] Lim 2025 PNAS, PMID 40560621 | W6 knockout-flux clustering | OK; abstract confirms the knockout-flux UMAP/k-means design, MCF7 only |
| [33] FlowGAT, PMID 38448436 | W7 | OK |
| [14],[15],[29],[76] and CORUM [71],[72] | W5 | titles/DOIs resolve (CORUM 2022 is NAR 2023 Jan 6, as cited) |
| [6] DeepMeta Cell Rep 2025, PMID 40616842 | B6 framing | Lookup by journal/year ambiguous; PMID/DOI not independently confirmed here. LITERATURE already flags that only the preprint was read |

Hygiene issues in LITERATURE.md (not scientific errors):
- Reference list skips numbers 3-4 ("numbers not used") and lists Hart 2015 inside ref 24. DESIGN §11 requires every reference to be in LITERATURE with a checked DOI; renumber, or delete the gap.
- LITERATURE §3 and §0.2 state "No published study reports per-gene or per-lineage AUPRC ... at full scale with lineage hold-out." My searches (PubMed, web) did not find one, but the statement is absence-of-evidence; phrase it as "we found none (searched 2026-10-26, terms listed)".
- No citation errors found among the ~16 checked beyond the Vieira numbers.

## 2. Missed or under-used relevant work (2022-2026) and positioning gaps

Found by search; none is in LITERATURE.md. Titles were read from the fetched page or PubMed record. Cite only after the literature agent reads them.

1. Huijer, Jiao, Chen, Yu. Improved flux profiling in genome-scale modeling of human cell metabolism. Cell Rep Methods 2026;6(1):101275. doi:10.1016/j.crmeth.2025.101275 (PMC12853175). Human-GEM + tINIT cell-line models, medium-component constraints, regression-based exchange fluxes (REGP vs CORE). No essentiality, but it is the most recent Human-GEM precedent for quantitative exchange bounds. It directly bears on W2 M3 (concentration-proportional bounds) and on the claim in LITERATURE §4.3 that few uptake-rate sources exist. Cite it in W2 and in §15 Q5.
2. Cheng et al. Identifying essential genes in GEMs of consensus molecular subtypes of colorectal cancer. PLoS ONE 2023;18:e0286032. doi:10.1371/journal.pone.0286032. Recon3D + CORDA, DMEM vs Ham's vs with/without cholesterol uptake, checked against DepMap. Medium-aware essentiality and the cholesterol-uptake control that W2 uses; it is a prior instance of that control.
3. Guil, García, García. Adding metabolic tasks to human GEM models to improve the study of gene targets and toxicities. Sci Rep 2024;14:17265. doi:10.1038/s41598-024-68073-8 (PMID 39068208). Task-based vulnerability (gMCS) in Human GEMs; relevant to W3's "tasks as criteria" and to the separate note that tasks change which targets are found.
4. npj Syst Biol Appl 2023, "Synthetic lethality in large-scale integrated metabolic and regulatory network models of human cells" (doi:10.1038/s41540-023-00296-3; search result says it used Human1 v1.14.0). Prior Human1 KO/SL validation; authors, details not read.
5. bioRxiv 10.1101/593277, "On the inconsistent treatment of gene-protein-reaction rules in context-specific metabolic models" (authors and publication status not verified). Directly relevant to W4 F0/F1 (how GPRs map to reaction presence: AND=min, OR=max/sum) and to the claim in DESIGN §W4 that F4 uses AND=min, OR=sum. State the GPR-to-score convention as a pre-registered choice and cite this.
6. Not found (searched, state as such in LITERATURE): a published DepMap-scale resource of cell-line-specific Human1/Human2 models with released SBML besides the Human-GEM repo's own `code/DepMapGeneEss` pipeline (see §3). Metabolic Atlas cell-line models: I could not verify what is hosted there; the literature agent should check metabolicatlas.org/downloads for HPA-based cell-line/tissue Human-GEM models and say whether they are DepMap-keyed. DESIGN §8 W8 claims no novelty framing; add a sentence on what is new relative to that resource.

Positioning statements in DESIGN that are accurate: "no published generic H1/H2 FBA vs DepMap AUPRC baseline"; Kim 2026 random-CV caveat; Richelle-lineage evidence. The statement in §0 of LITERATURE ("not tested with ftINIT at DepMap scale") is questionable: Gustafsson 2023 built ftINIT/tINIT models for 891 DepMap lines and scored MCC; what is untested is variance partitioning by method vs tissue. Reword so it does not overclaim.

## 3. Reuse of existing code and resources

On disk (checked):
- `data/Human-GEM-main/code/DepMapGeneEss/` (MATLAB): `generate_DepMap_models_ftINIT.m`, `getTaskEssentialGenes.m`, `evaluateDepMapEssentialityPredictions.m`, `PrepDepMapData.m`. This is the Robinson/Gustafsson DepMap pipeline (tINIT/ftINIT model per line, task-essential genes, comparison to Achilles gene effect). DESIGN never mentions it. Use it (a) as the specification of the published protocol that B7 below must reproduce, (b) to copy gene-ID handling and task-essential logic. MATLAB RAVEN is not available, so read it, do not run it.
- `code/test/gradedEssentiality.py` plus `estimateEssentialGenes.py`, `evaluateHart2015Essentiality.py`, `diagnoseTaskEssential.py` (Python, cobra + `raven_toolbox`): already implements viability-vs-capability task scoring and growth-ratio on Ham's medium. DESIGN §10 uses it only as a positive control; it should also be the reference implementation for T1/T3 (task-augmented essentiality) so W3 does not re-implement the same logic. Add a unit test that the project's `nextgem.knockouts` growth ratio equals `growth_ratios()` of this file on the Hart lines.
- `raven_toolbox` 3.0.0b1 is installed with `init/ftinit.py`, `init/taskfill.py`, `tasks/check.py`, `tasks/tasklist.py`, `omics`, `curation`, `gapfilling`, `conditions`. The DESIGN names only `init.ftinit`; use `tasks.check` for the T-modes and `omics`/`init.score` for expression-to-reaction scoring, instead of writing own versions. Check `raven_toolbox.curation`/`manipulation` before writing own GPR-edit code in W5.
- Not installed: `troppo`, `corda`, `memote`, `cellfie`. DESIGN says "own or troppo if it installs". Decide at M3 and pin; for GIMME/FASTCORE prefer `troppo` (published, benchmarked in Vieira 2022 against Human-GEM) because reimplementing makes the F2/F3 results less comparable to the literature. CellFie is MATLAB (and a Python re-implementation exists in Troppo/elsewhere; unverified): T5 says "CellFie defaults" but owns the implementation; add a toy test against the published CellFie output or drop the claim "CellFie defaults".
- COBRApy `single_gene_deletion` is used by the repo scripts; ensure nextgem's own sweep is cross-checked against it on a sample (DESIGN §10 has a verify step; name this one explicitly).
- Human-GEM `code/GPRs`, `code/modelCuration`: read before W5 for the existing GPR-fix tooling and the v1-to-v2 diff generator.

## 4. Baselines: is the strongest feasible published comparator included?

Not fully. Findings:

1. The strongest published human-cancer GEM comparator that is reproducible is the Human-GEM protocol: per-line ftINIT ("1+0", 1 TPM) model, then 57-task essentiality (Robinson 2020, Gustafsson 2023, repo `DepMapGeneEss`). DESIGN B0-B5 are generic models only; F1 appears in W4 as one candidate and P1 is compared only to G\* (C1). A reader would ask "does P1 beat what Human-GEM already does?" This is the key positioning contrast.
2. B6 (Lagziel-style RF) omits culture medium as a feature, although Lagziel [10] found medium to be the top feature for 11% of models, and the GEM arms get medium information in W2. C5 (P2 vs B6) is then unfair to B6.
3. DeepMeta [6] is a published per-gene predictor on DepMap labels; its code is probably not reproducible at low cost, so it is not a required baseline. Use its reported per-gene AUROC 0.61 only as a sanity anchor for B6/B2 on training data (not a comparison).
4. Troppo's reference pipeline [5] is the other published human DepMap-scale GEM benchmark; F2/F3 via `troppo` would let W4 report a number in the same family.

## 5. Concrete changes, ranked

### Must-fix
1. **Add baseline B7 "published Human-GEM protocol"** (DESIGN §7, §9.3 T1, §6.4): ftINIT 1+0, global 1 TPM, RNA only, medium M0, with both (i) growth-ratio score and (ii) 57-task essentiality score (Robinson/Gustafsson/repo `DepMapGeneEss`). Report it in every table and add contrast C6 "mAP_line(P1) - mAP_line(B7)" at T1 (Holm family becomes C1-C6). Do not let it be the same object as F1 under whatever medium/tasks P1 picks; fix its configuration now, label it as published-protocol, and tune nothing.
2. **State the existing resource and code in DESIGN** (§0/§2 or a new §2.4 "Prior art and reuse"): the Human-GEM `DepMapGeneEss` pipeline, `gradedEssentiality.py`, `raven_toolbox` modules, and what is new relative to them (medium ladder, a-priori task set with controls, evidence-based GPR curation with a curation-vs-constraint factorial, lineage hold-out, AUPRC grid, SBML for 26Q3). Without this the novelty claim of the deliverable and W8 is unanchored.
3. **Give B6 the same side information** (§7): add medium (parsed base medium + serum flag, from the same `lineage_split.tsv` `medium_status` source) and keep "no lineage features". Otherwise C5 over-credits the GEM.
4. **Correct the Vieira 2022 description** in LITERATURE (733 CCLE lines; parameters chosen on the MCF7 reference; "best-of-grid" label to be checked) before DESIGN relies on "0.44" as optimistic; this affects the framing in §W4 thresholds.

### Should-fix
5. W2 (§8): cite Huijer 2026 and Cheng 2023; in the M3 step consider a one-line justification of why availability-only (1000/0) bounds are used when a regression-based flux-bound method for Human-GEM exists, and treat it as out of scope for lack of exchange data for DepMap lines.
6. W4 (§8): add the GPR-to-reaction-score convention as a pre-registered choice (AND=min, OR=max versus sum) and cite the GPR-treatment preprint (10.1101/593277, verify authors/status). Currently F4 uses OR=sum and F1/F2/F3 use whatever the tool does; differences confound the method-vs-tissue test (H4.1).
7. W3 (§8): implement the task-augmented score by calling/refactoring `gradedEssentiality.py` and `raven_toolbox.tasks.check`, with a test of equality to the repo logic on Hart lines. Also add Guil 2024 as the precedent for tasks beyond biomass.
8. T5/CellFie (§8 W3): either validate the own CellFie-like task calling against the published CellFie output on a few lines, or rename "CellFie-inspired".
9. W4 F2/F3: pin `troppo` at M3 if installable (benchmarked in the Human-GEM context); only reimplement as fallback. Same for CORDA.
10. LITERATURE §0/§3: reword "not tested with ftINIT at DepMap scale" (Gustafsson tested ftINIT on 891 DepMap lines for MCC); the untested element is method-vs-tissue variance partitioning.
11. LITERATURE: add a short "searched, not found" list (databases, date, query strings) to support the absence claims (no per-gene AUPRC on DepMap for GEMs; no curation-vs-constraint ablation; no flux-vs-omics head-to-head).

### Optional
12. Verify and, if confirmed, mention any Metabolic Atlas cell-line Human-GEM models in the positioning paragraph.
13. Use Lagziel and DeepMeta reported numbers as order-of-magnitude anchors for B6 on training data only; do not tune B6 to match them.
14. Renumber the LITERATURE reference list (gap at 3-4).
15. Gonçalves 2022 [18] and Gopalakrishnan 2023 [9] are only abstract-read; DESIGN W4 stability and W5 ProCan evidence rely on them. Have the literature agent read the full text before M1.
