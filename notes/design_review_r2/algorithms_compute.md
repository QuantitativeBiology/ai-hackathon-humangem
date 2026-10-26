# Design review round 2: algorithms and compute

Reviewer: CBM algorithms / HPC reviewer (also wrote round 1). This round reviewed DESIGN.md v1 (commit d84a593) and `notes/design_review_r1/PI_decisions.md`.

All round-1 probes had finished; no background jobs were running. One new label-free probe was run (`probe3.py`, 1 Gurobi session, about 9 CPU-min): F0 pruning plus the v1 repair LP on 3 lines, with sweep timing. Scripts and JSON are in `/tmp/claude-10105/-data-tmpA-nextgem/c5c0f45b-d54c-4e57-bdef-3d0657510b12/scratchpad/review_r1_compute/`.

**Verdict: REVISE (minor).** All six round-1 must-fixes are addressed in substance. There are 5 new, narrow must-fixes (§3), each a one-paragraph edit. None changes the architecture. Once they are applied I would accept a freeze without another algorithms round.

## 1. Round-1 must-fixes: status

| R1 | v1 resolution | Status |
|---|---|---|
| MF1: medium infeasible | micronutrient floor, per-(medium × parent) growth precheck, M1 redefined (W2) | **Addressed**, but the bound scheme for the floor creates a new defect (NMF1) |
| MF2: pruning repair | L1 repair LP, per-loss viability for the genomic constraint (W4) | **Addressed.** Probe on the 3 lines F0 had killed: the repair LP solves in 0.2–0.45 s, re-opens 25–37 reactions (1–2% of the 1,522–2,257 closed), and restores WT growth to 41–65% of generic. The "viability tasks feasible" clause needs one sentence of specification (SF-a). |
| MF3: ftINIT spec, cost, acceptance | B7 pinned to the repository's Python configuration; Jaccard ≥ 0.95 / AUPRC ± 0.01 one-shot test; measured cost | **Addressed**, but the §11.4 MILP failure rule contradicts the B7 configuration (NMF2), and F1 under a non-M0 medium has no growth guarantee (NMF3) |
| MF4: run matrix | Stages A/B/C; ≤ 4 configurations per line, ≤ 1 ftINIT | **Partly.** Stage B/C are fine. Stage A is still not enumerated, and read literally it is about 5× the budgeted builds (NMF4) |
| MF5: MOMA | QP `add_moma` once; drop above 1 CPU-h per model | **Addressed.** Note: QP needs Gurobi, so O-MOMA takes Gurobi slots, not GLPK |
| MF6: task modes | D-BS/D-ER product demands; NGAM from max MAR03964 | **Addressed** |

## 2. Probe results (new)

Setup: Human2, Ham's, Gurobi with 1 thread, first 3 default-expression rows, TPM < 1 → off.

| Line | off genes | closed rxns | repair LP | re-opened | WT / generic | genes to solve | Gurobi sweep | lethal calls |
|---|---|---|---|---|---|---|---|---|
| ACH-001113 | 836 | 1,522 | 0.45 s | 27 | 0.65 | 291 | 32 s | 159 |
| ACH-001289 | 903 | 2,257 | 0.37 s | 25 | 0.41 | 287 | 43 s | 157 |
| ACH-001339 | 873 | 1,861 | 0.21 s | 37 | 0.64 | 296 | 35 s | 163 |

- **Real LP unit cost.** A per-line build, repair and sweep costs about **1 CPU-min** of LP. My round-1 figure of 3.5 min was the generic sweep, which has slower KOs.
- **The Python overhead was 70–220 s per line, larger than the LPs.** It came from the zero-flux scan done with cobra `knock_out()` loops.
- **Fix:** precompute, for each line's off-set, a gene → closed-reaction map by vectorised GPR evaluation. The per-line cost is then about 1–1.5 CPU-min (SF-b).
- **Memory:** RSS was 1.45 GB per process (model + expression). ftINIT RSS is still unmeasured.

## 3. New must-fix

**NMF1. The floor and serum f-scaling caps growth and erases the continuous score (W2 M1/M2(f), "Biomass micronutrient floor", "Requirement").**
- The floor components (tocopherols, retinoate, lipoate, linoleate, linolenate, cobalamin) **cannot be synthesised**. Bounding them at f × req with f ∈ {0.001, 0.01, 0.1} caps WT growth at f × the unconstrained growth.
- Every non-lethal KO then returns ratio ≈ 1, and the partial ratios are lost.
- At f = 0.001, WT growth ≈ 1e-3 × the unconstrained growth, which sits on the 1e-3 × M0 precheck threshold, so the precheck passes or fails by rounding.
- In addition, req_c ("FVA-min at WT growth on the base medium") is undefined: the base medium alone does not grow (round-1 probe).

Fix:
- (i) Bound the floor components at a fixed, non-limiting level (e.g. 10 × req, or 1000), independent of f.
- (ii) Apply f only to **synthesisable** serum components (fatty acids beyond the essential pair, cholesterol, carnitine, ethanolamine). For these, partial import is meaningful: it makes FASN/ACACA and the sterol branch partially required.
- (iii) Compute req on base medium + floor open.
- (iv) Add a precheck assertion that no floor component is binding (its shadow price is 0 at WT).

**NMF2. The MILP failure rule contradicts B7 (§11.4 vs §8 B7).** §11.4 says "an MILP that stops at its time limit with gap > 1e-4 counts as a failed build". The B7 configuration (`time_limit = 900`, `prove_abs_gap = 1`, `resolve_ties`) is documented in `estimateEssentialGenes.py` to hit the limit:
- RPE1, step 1;
- the GBM tie-break phase.

That is 2 of 5 Hart lines. With an absolute gap of about 1 on objectives of order 10³, the relative gap can exceed 1e-4. Under v1, B7 would therefore record about 20–40% failed lines, be flagged "failing", get the no-information score, and C6 would be meaningless. The Hart one-shot test would also contradict itself, since it compares against repository values built by accepting those incumbents.

Fix: for ftINIT, define failure as "no incumbent, or infeasible". A time-limited incumbent is **accepted, as in the repository**, with the achieved gap, the step and `status = time_limit` recorded. Report the share of time-limited builds per family; the 5% flag counts only real failures. The 1e-4 rule can stay for non-ftINIT MILPs, if any.

**NMF3. F1 models evaluated off-M0 have no growth guarantee (W4 F1, §10.1).** ftINIT gap-fills the 57 tasks with their own (Ham's) inputs. The per-line F1 build is B7, and the "growth on line medium" variant is panel- and group-only. So any combination of F1 with a W2 medium (the factorial, P1) can leave the line non-growing, giving no-information scores and a >5% flag driven by the medium, not by extraction.

Fix:
- Apply the **same repair LP** after extraction whenever an F1 model is evaluated on a medium other than M0. Re-open template reactions; cost < 1 s per line, measured.
- Log the repairs.
- B7 itself stays unrepaired, on M0.

**NMF4. Stage A is not enumerated, and as written it is about 5× its budget (W4 "Staged run matrix", §13).** "Builds: families × thresholds × modalities" (± genomic), plus the variance-partition builds, plus the F5 4×4×3 grid, gives:
- about 9,000 LP builds;
- **up to 1,800 ftINIT builds** (F1 × 3 thresholds × 3 modalities × ± genomic × 100 lines, ≈ 660 MILP CPU-h);
- 4,800 F5 sweeps.

§13 budgets about 2,100 builds and about 400 ftINIT. Enumerate Stage A in a frozen table, for example:

| Family | Builds |
|---|---|
| F1 | R1/R3 × real/permuted (the 400 variance builds), plus R1 protein on the panel lines with ProCan (about 60). No R2. ± genomic is applied **post-extraction** (gene KO + repair), not as a rebuild. |
| F0, F2, F3, F4 | R1–R3 × {RNA, protein, RNA+protein} on lines with data; ± genomic post hoc |
| F5 | 24-point grid (coarse 12 + local 12); nested LOLO reuses the label-free builds |

That is about 560 ftINIT and about 8,000 LP builds, ≈ 410 CPU-h with about 205 of it MILP.

**NMF5. The resource allocation exceeds the host's RAM (§13 "Allocation").** The plan runs 14 MILP + 24 LP workers.
- LP workers measure 1.45 GB each.
- ftINIT RSS is unmeasured; deepcopy-heavy gap-filling makes 3–4 GB plausible.
- That totals about 80–90 GB, against 69 GB currently available on a shared 125 GB host (56 GB used by others) and the design's own 50 GB cap.

Fix: allocate by memory, not by sessions. A start of 12 MILP + 16 LP workers (≈ 60 GB) is plausible. Measure ftINIT peak RSS in the M3 pilot and set the final worker counts from it, as a frozen config rather than at runtime. Abort the run if available memory falls below 10 GB. The wall-clock figures below use 12 + 16.

## 4. Should-fix (new)

- **SF-a, repair with tasks.** "Viability tasks (GR, ER) feasible" cannot be part of the same LP as growth on the line's medium, because tasks are closed-system tests with their own inputs. Specify: one L1 LP per task (10 LPs) plus the growth LP, and re-open the union.
- **SF-b, per-line overhead.** Use the vectorised closure map (§2). Without it, the Python overhead triples per-line LP cost.
- **SF-c, H1 resource cost.** "+250 CPU-h" understates it if P1 includes F1: (2,030 − ≈545 already built) × 22 min ≈ **545 MILP CPU-h**. Either budget that, or state now that the H1 resource releases only the LP-family part of P1.
- **SF-d, W6 assay media.** X_KOflux in three assay media for about 1,292 lines is about 50 extra CPU-h. Raise W6 to about 90.
- **SF-e, speculative pre-builds.** P1(H1) on validation and test (≈ 200 MILP CPU-h, about 17 h on 12–14 slots) currently sits between M6 and T1, on the critical path. Its builds are label-free. Start H1 ftINIT for validation and test at M5, with the B7-equivalent configuration, and discard them if P1 has no F1.
- **SF-f, transporter gate.** The gate has 2 states, not 2³: "any of SLC1A1/2/3 ≥ 1" is one boolean.

## 5. Open questions in my area

**R2-3 (2,940 vs 2,650 CPU-h).** 2,940 is realistic for v1's matrix. My re-estimate is ≈ 2,950 (table below), but the composition differs from §13:
- Stage A is about 410, not 360, after NMF4 trimming. It would be 660+ untrimmed.
- Stage B and other LP work are cheaper (≈ 100, not 250), at the measured 1–1.5 CPU-min per line.

On B7 coverage: the training gate (TG) for MILP families is panel-based by design, so B7 on non-panel training lines feeds no decision. Recommendation:
- Build B7 now on panel + validation + test, about 650 lines and ≈ 240 CPU-h.
- **Defer** the other ~640 training lines until P1 is known. If P1 contains F1, W8 needs them anyway; if not, they are skipped, saving ≈ 235 MILP CPU-h and ≈ 17 h on the MILP lane before M5.
- C6 only needs validation and test.

**R2-5 (floor for serum-free media).** Keep the floor for every line; with NMF1 it is non-limiting. Otherwise all 160 serum-free lines fail and are scored as no-information, which tests nothing.
- Source the floor from parsed supplements where documented: B27 and N2 carry tocopherol, retinyl acetate, linoleic and linolenic acid; DMEM/F12 carries lipoate and B12. Flag `floor_source = modelling` only where nothing is documented.
- Sensitivity: instead of "no floor, expected non-growth", use a **floor-free biomass** (B-loo of the floor components, which W3 already defines). It removes the requirement rather than supplying it, and stays informative.
- Never put cholesterol in the floor. The serum-free sterol control depends on that.

**R2-6 (transporter gate attribution).** Attribute the gate to **context**, not medium. It uses line expression, which is the defining input of the context lever. PROMPT's "error attributable to medium misspecification alone" must not include an expression-derived term.
- "Medium alone" = M0 → documented medium (M1/M2), ungated, on the generic parent.
- In the Shapley factorial the gate is either part of the context factor or its own fifth factor (2^5 = 32 cells, still label-free and cheap at ~1.5 CPU-min per line; about 55 extra CPU-h on validation and test).

**R2-2 (failed models).** Not my area, but one algorithmic point. With NMF2 and NMF3, genuine failures should be rare (repair costs < 1 s). The no-information score is acceptable only if the failure share is reported per method. I would also add a "line excluded for all methods" sensitivity row.

## 6. Compute re-estimate for v1's run matrix (with NMF1–5)

Unit costs:
- LP build + repair + sweep: 1.5 CPU-min per line (probe, with SF-b);
- generic sweep: 3.5 CPU-min;
- ftINIT: 22 CPU-min per line (repository) + 0.5 h prep;
- RSS: LP 1.45 GB measured; MILP 3.5 GB assumed.

| Block | CPU-h | MILP |
|---|---|---|
| W1 + W2 (≈ 1,450 per-medium sweeps across f, gate and parents; prechecks) | 40 | 0 |
| W3 | 25 | 5 (QP MOMA) |
| W5 (enumeration; 2×2 on extracted sets) | 60 | 0 |
| Prep, pilots, Hart test | 20 | 10 |
| Stage A (enumerated as in NMF4) | 410 | 205 |
| B7 per line, all 1,292 lines (or ≈ 650 now, the rest deferred) | 475 (240) | 475 (240) |
| Stage B (3 LP configurations × 1,292) | 100 | 0 |
| Stage C (150 group models + stability) | 95 | 85 |
| H1 arm, P1(H1) on validation + test | 200 | 190 |
| Factorial (panel 24 × 100 → top 3 on training → P1 2^4 on validation + test) | 335 | 0 |
| W6 (incl. assay media) | 90 | 0 |
| W7 cap | 300 | 50 |
| W8 (750 new F1 builds + validation of 2,030) | 330 | 275 |
| Subtotal | 2,480 | 1,295 |
| Re-runs and fixes (+20%) | 500 | 260 |
| **Total** | **≈ 2,950** (≈ 2,720 with B7 deferral; +545 MILP for a full H1 F1 resource) | **≈ 1,550** (≈ 1,270 deferred) |

Wall-clock at 12 MILP slots and 16 LP workers:

| Lane | Load | Wall |
|---|---|---|
| MILP before M6 (B7 + Stage A + Stage C + prep, ×1.2) | ≈ 930 CPU-h (≈ 650 with deferral) | ≈ 78 h (≈ 54 h) |
| LP before M6 (W2/W3/W5 + Stage A/B LP, ×1.2) | ≈ 600 CPU-h | ≈ 38 h, hidden under the MILP lane |
| Between M6 and T1 | factorial 335 CPU-h on LP ≈ 21 h; P1(H1) ≈ 17 h on MILP (0 h if pre-built per SF-e) | ≈ 21 h |
| After T1: W8 | 330 CPU-h, mostly MILP | ≈ 28 h |
| After T1: W7 | ≤ 300 CPU-h | ≈ 19 h of compute, but agent iteration dominates (≈ 2 days) |

**Critical path:**

M3 Hart test (1.5 h) → MILP lane, starting at M2 with prep, B7 and Stage A F1, which are label-free (54–78 h) → M6 factorial (21 h) → T1 → W7 ‖ W8 (≈ 48 h, agent-bound) → T2 → REPORT.

That is about **125–150 h of compute-bound wall time (5–6 days)**, inside v1's D+11 schedule with its 1.5× buffer. The tight point is M5 at D+4–6. It holds only if the MILP lane starts at M2 and B7 training lines are deferred.

**Disk:** unchanged, ≤ 10 GB.

**Licence:** 12 of the 16 sessions go to MILP, set by memory, not by the licence. The other 4 cover QP MOMA, verification and interactive work.
