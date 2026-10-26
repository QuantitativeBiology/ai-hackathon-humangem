# Design review r2: prior work and reuse

Checked DESIGN.md v1 (commit d84a593) and the corrected LITERATURE.md against my r1 must-fixes and PI_decisions.md.

## 1. Round-1 must-fixes

| r1 must-fix | Status | Evidence |
|---|---|---|
| MF1 B7 + C6 | Addressed (one wording caveat, below) | DESIGN §8 B7 row: ftINIT config pinned to `estimateEssentialGenes.py` (1+1, removeGenes, resolve_ties, gap/time limits), B7g growth ratio and B7t task essentiality, nothing tuned. C6 = P1 vs B7g and B7t listed as an estimate (§6.3). One-shot Hart test with Jaccard >= 0.95 (§12/§13). I confirmed in the repo that `estimateEssentialGenes.py` uses `1+1` while the MATLAB `DepMapGeneEss` default is `1+0`. |
| MF2 prior-art section | Addressed | New DESIGN §2 names `DepMapGeneEss`, `gradedEssentiality.py`, `estimateEssentialGenes.py`, `taskEssentialGenes.py`, `raven_toolbox` modules, cobra, with use mode ("read, not run" for MATLAB) and an equality test; lists 7 items that are new. |
| MF3 medium in B6 | Addressed | B6 row: "medium class + serum flag"; B2m medium-aware prior also added. |
| MF4 LITERATURE corrections | Addressed | Vieira row: 733 CCLE lines, parameters chosen on MCF7 reference, "selected configuration, not lineage-held-out". §0.4 and §3.2 reworded: Gustafsson ran ftINIT on 891 DepMap lines; only method-vs-tissue is untested. Numbering gap 3-4 filled (Huijer, Cheng). Added Guil 2024 (89), Barrena 2023 (90). §8.1 "searched, not found". |

Re-verified in PubMed this round: DeepMeta = Wu T et al., Cell Rep 44(7):115945, PMID 40616842, doi 10.1016/j.celrep.2025.115945 (first author "Wu", correction is right). Huijer PMID 41529686 and doi OK. Cheng PMID 37205704 OK. Barrena PMID 37454223, doi 10.1038/s41540-023-00296-3, npj Syst Biol Appl 9:32, OK (Human1 + regulatory networks, gMCS, in-vitro silencing data). Guil PMID 39068208 OK. The unverified GPR preprint (10.1101/593277) is correctly kept out of the citations and flagged.

Residual wording issue (should-fix, not blocking): DESIGN §2 and §8 call B7 "the published Human-GEM protocol". The Robinson/Gustafsson DepMap protocol (MATLAB `DepMapGeneEss`) uses series 1+0; B7 is the repository Python protocol (1+1). Rename B7 to "Human-GEM repository protocol (v2.1.1 Python, ftINIT 1+1)" in §2/§8/§6.3 and in the abstract sentence, and note the 1+0 vs 1+1 difference in the F1 row, so the claim "beats the published protocol" is not overstated. Optionally one panel-only run at 1+0 to bound the difference.

## 2. Should C6 join the gatekeeping sequence (R2-1)?

Yes, but at the END and for B7g only.
- C6 is the contrast that answers "better than what Human-GEM users already do". It should be confirmatory, otherwise the headline "optimised GEMs beat existing practice" rests on an estimate.
- Do not place it before C3, C4 or C5: a failure of C6 (P1 is built on F1-like extraction and may sit close to B7) would block the context-specificity and best-predictor claims that are separate questions. Sequence: C1 -> C1b -> C3 -> C4 -> C5 -> C6g.
- Declare B7g (growth ratio, the same score type as P1) as the single confirmatory comparator. B7t stays an estimate; this avoids choosing the weaker comparator post hoc. State that C6 is evaluated against the T1-sealed value after C5.
- Alternatively keep the main chain and give C6g its own one-sided alpha 0.025 as a second chain conditional on C1 passing. Either is acceptable; the end-of-chain option needs no extra alpha.

R2-3 (B7 on all lines vs only validation + test): B7 must run on every test and validation line (needed for C6 and G* style selection). Training lines other than the panel are only needed for the F1 training gate; running B7 there is optional from a prior-work view. No objection to the cheaper option provided the F1 training-gate gap is stated.

## 3. New must-fix issues

None from the prior-work angle. Should-fix only:
1. B7 naming (above).
2. If the H1 arm is full, say whether B7 is also run on Human1 (H1 parent at repo-equivalent settings) or only on H2; the "within-parent" contrasts imply a B7(H1) is needed for symmetry. State the choice in §8.
3. B6 neighbour-enzyme features are described as "Lagziel-style"; state in §8 that they are a re-implementation, not Lagziel's code, and that the medium feature uses the same tiered medium annotation as W2 (so B6 and GEM arms see identical medium information, including imputed-medium lines).
4. Verification of the Barrena/Guil/Huijer/Cheng entries is abstract-level only (as LITERATURE says); keep claims in DESIGN to what the abstracts state.

## 4. Verdict

Freeze from the prior-work angle, after the B7 renaming and the C6 placement decision are written into §6.3/§8.
