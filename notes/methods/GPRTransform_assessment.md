# GPRTransform (Machado, Herrgård & Rocha 2016): assessment

Label-free (hook probe blocked; no `labels.*`; locked tier untouched). Code: `results/gprt/proto/`; outputs: `results/gprt/`.

## 1. The transformation
Source: PLoS Comput Biol 12(10):e1005140.
- **Enzymes.** Each GPR is converted to disjunctive normal form (DNF). Each term is one enzyme, that is, a set of subunit genes.
- **Reactions.** Gene-associated reactions become irreversible. A reaction with k enzymes becomes k copies.
- **Genes.** Each gene is a pseudo-metabolite. Every copy that contains it consumes it at −1 per unit flux; subunit copy numbers are ignored. A usage reaction u_g (→ gene, bounds [0, ∞)) supplies it.
- u_g is the gene's total enzyme flux over all its reactions and complexes. KO: u_g = 0. Unbounded u leaves FBA unchanged.

**Caveat (ours).** In the paper, and in `reframed.gpr_transform`, each isozyme copy inherits the reaction's bound. A k-isozyme reaction can therefore carry k × its bound. Our default "arm" encoding is exact:
- the copies keep the original bounds and consume a catalysis token;
- each enzyme is an unbounded arm (subunits → token).

## 2. Analyses it enables
All are LPs (D2-compatible): gene-pFBA (min Σu); gene-(l)MOMA; essentiality as min u_g > 0 at 10% growth (= KO ratio < 0.1, plus a usage range); gene-E-Flux (u_g ≤ e_g/max e); gene-GIMME (min Σ max(0, cutoff − e_g)·u_g); gene-level MCS/EFMs.

## 3. Difference from our methods
- **Boolean KO.** u_g = 0 blocks exactly the DNF terms that contain g. Generic KO is therefore **identical by construction** (§5).
- **F4 (AND = min, OR = sum)** already grades isozyme KOs. Gene-level bounds differ only in that one capacity is *shared*:
  - across all reactions of a promiscuous gene (Human2: median 2 reactions per gene, max 257);
  - across all complexes containing the gene.

  Scores differ only where these shared budgets bind.
- **F2.** Gene-GIMME charges low-expression genes per subunit and per reaction, instead of one min-aggregated cost per reaction. It is still an extraction method.
- **Ties.** A KO that leaves the WT optimum feasible leaves the score unchanged, so no growth-ratio variant can rank the ~1,750 genes closing no reaction or the ~2,440 unused at a line optimum (§6). Only a *cost* score (Δ min Σw·u) is continuous, for the ~430 used genes.
- **A4.** "Off isozymes do not buffer" is u_g = 0, but "a re-opened reaction keeps its full GPR" has no gene-bound equivalent: presence families would change semantics; bound-type use does not.

## 4. Regulation, repository and validation
**Regulation.** **No regulatory network**: no TF rules, no rFBA/SR-FBA/PROM coupling (SR-FBA, OptORF, TIGER are cited only as MILP alternatives). "Regulation" means over/under-expression bounds on u_g relative to WT gene-pFBA usage, and cRegMCS up/down-regulation designs: gene-level handles, nothing more.

**Repository** (commit 4952ea6, 2016-11-21; cloned to `results/gprt/tmp/GPRTransform`, not executed): Python 2 notebooks plus MATLAB; pins **framed v0.2**, Gurobi 6.5, CPLEX 12.6.3. **No license** (no LICENSE file; GitHub `null`, so all rights reserved; the article is CC-BY). Unmaintained since 2016. `reframed` (Apache-2.0, active 2026-06) has `gpr_transform`; we reimplemented from the paper.

**Validation in the paper.** *E. coli* only: reduced iAF1260 (1,532 → 3,853 reactions) and the core model.
- **No yeast work and no experimental essentiality benchmark**; essentiality appears only as FVA of u_g (S1 Fig).
- **Fluxes of 24 KO mutants (Ishii 2007):** gene-pFBA/MOMA/lMOMA beat their reaction-level versions.
- **Transcriptomics (Ishii; Gerosa):** nothing beat pFBA; gene-GIMME > GIMME; gene-E-Flux **no** gain.

## 5. Prototype and equivalence
Ham's medium (M0), GLPK. Outputs: `sizes_*.json`, `equivalence_*.parquet`.

| | Human1 | Human2 |
|---|---|---|
| Reactions before → after | 12,971 → 38,604 | 12,748 → 36,229 |
| Metabolites before → after | 8,455 → 15,329 | 8,374 → 15,045 |
| Isozyme reactions | 3,987 | 3,801 |
| Complex reactions | 698 | 713 |
| Mixed AND/OR (non-DNF) | 121 (5) | 129 (17) |
| Max DNF terms (V-ATPase) / subunits | 192 / 45 | 192 / 45 |
| Transform CPU | 8.8 s | 6.3 s |
| WT solve, cold original / transformed / warm | 6.8 / 16.1 / 0.20 s | 6.6 / 18.6 / 0.16 s |
| KO sweep, every gene solved | 23 CPU-min | 16 CPU-min |
| Mismatches (ratio > 1e-6 or lethal call) | **0 / 2,887** | **0 / 2,870** |
| Max \|Δratio\| | 7e-11 | 8e-11 |
| Lethal genes (baseline = transformed) | 149 = 149 | 154 = 154 |

The baseline Gurobi sweep takes 3–3.5 CPU-min.

**Paper encoding on Human2:** WT 93.45 vs 74.76 (+25%); **67 ratios differ** (max 0.20; OXPHOS, ATP synthase, glycolysis, TYMP); same lethal set. Cause: isozyme capacity relaxation, not a GPR edge case (A4 does not apply to generic models).

## 6. Pilot: gene-level E-Flux (variant A)
Lines ACH-000957, ACH-002024, ACH-003383; F4's R1-RNA lin_g = 2^s_g; u_g ≤ max(10, 1000·lin_g/max lin). Reference: the F4 KO of the same lines (R1/real/none).

| Per line | Gene-E-Flux | F4 | Generic |
|---|---|---|---|
| Distinct ratio_r | 29 / 32 / 31 | 31 / 33 / 33 | 32 |
| Lethal | 154 | 154 | 154 |
| Partial | 89 / 94 / 93 | 87 / 87 / 89 | 91 |
| Tied at 1 | 2,627 / 2,622 / 2,623 | 2,629 / 2,629 / 2,627 | 2,625 |

- WT 0.99–1.08 (≈ F4). 2,773–2,800 of 2,854 genes sit at the floor (max lin ≈ 2×10⁴), so expression barely enters.
- 88–91 genes differ from F4 (Spearman 0.96–0.97); no lethal call changes. Isozyme genes (2,206): 36–39 change; **7–10 graded here but 1 in F4** (PKM, PKLR, PGK2, PGM1, PGAM2/4, ENO2/3, SLC29A1, DECR1, MECR, CRISP3).
- 5.3–5.6 CPU-min per line (~425 KO LPs after the exact zero-usage shortcut) + 0.4 CPU-min parent load/transform; an F4 sweep ≈ 2.7.

## 7. Candidate variants
All parameters are label-free defaults.

**A. F6-GE, gene-level E-Flux** (§6; KO and accepted losses u_g = 0; no repair; failed if WT ≤ 1e-6). Ties: none gained (measured). Row A8 = {R1, R3} × {real, perm} × 30 lines = 120 builds × 5.8 CPU-min + 60 genomic re-sweeps ≈ **17 CPU-h**.

**B. Gene-pFBA KO-cost tie-break.**
- **Formulation:**
  - w_g = 1 + max(0, −s_g);
  - C* = min Σw·u at growth ≥ 0.9μ*;
  - C_g = the same, with u_g = 0, at growth ≥ 0.9μ_g;
  - score = round(ratio_r − 10⁻³·min(1, (C_g/μ_g)/(C*/μ*) − 1), 6).
- **Ordering:** lethal genes still rank below non-lethal ones, because the smallest non-lethal ratio is 0.505.
- **Ties:** only the ~430 used genes can move, so at most ~190 of ~2,625 ties break.
- **Cost:** ≈ 6 CPU-min per line. 30 panel builds plus 2 generic models ≈ **3.5 CPU-h**.

**C. Gene-GIMME extraction** (F2 variant): min Σ max(0, −s_g)·u_g at ≥ 0.9μ*; keep reactions with merged flux > 1e-6, s_r ≥ 0 or no GPR; then the F2 pipeline (repair, A4, Boolean KO on the original). Ties: as F2. 120 builds × 3.5 CPU-min + 60 re-sweeps ≈ **8.5 CPU-h**.

## 8. Recommendation
**Do not adopt as F6 or as a W5 lever.** Generic KO is provably identical; gene-level E-Flux under F4's ratified scaling adds no distinct values at ~2× F4 cost; the transform does not curate GPRs. Lesson: avoid `reframed.gpr_transform`'s isozyme bounds.

If the PI wants a tie analysis, **B** is the only candidate with a mechanism (≤ 3.5 CPU-h, exploratory). Its §17 amendment would name B a **secondary, exploratory** score never replacing ratio_r in §6.1/§6.3, fix its formulation, parameters, GLPK and the arm encoding, add a Stage-A row, budget and FROZEN.tsv entries, and be dated before any B score meets labels. Adopting A as F6 would also need a W4 family row, a `DECLARED_SOLVER` entry, a KO/loss clause extending A4(a), a failed-build rule, variance-partition/P1 eligibility and a §13 update.

**Runtime used:** ≈ 1.7 CPU-h on GLPK; no Gurobi.
