# Label-exposure register

Every viewing of DepMap label *values* (or label-derived summaries) by any agent or the PI. Required by DESIGN.md §11.5. Agents listed here must not author or approve a priori lists, candidate sets or decision thresholds.

## Entries (as of DESIGN v2, 2026-10-26)

- the baseline review read the CRISPR header only;
- the literature agent saw the all-line DepMap common-essential count;
- the **round-1 statistics review**:
  - evaluated the retired seed 20261026 under the v0 rules;
  - named the old test tier (Pancreas, Soft Tissue, PNS, Prostate, Eye) and validation tier (CNS/Brain, Ovary/Fallopian Tube, Bladder/Urinary Tract, Cervix);
  - computed label summaries on the old tier: positives per line, eligible-gene counts and a prevalence/SD remark;
  - also computed all-line prevalences and the common-essential share of positives;
- the **round-1 data reviewer** computed ASNS prevalence by medium and per-library positive rates on all lines;
- the **design author**, while inspecting a file header with `head`, printed about 10 L-prob values for one line (ACH-002730);
- the **all-line label copy** `ge_U.parquet` (CRISPRGeneEffect restricted to U, all 1,361 lines) and other derived scratch files (`testtier.py`, `splits.pkl`). These were written by a round-1 reviewer in the shared scratchpad and **deleted by the PI on 2026-10-26**;
- the **PI's header-only read** of a label file while testing the hook.

**Consequences.** The retired seed's tiers are known, which is why the seed is replaced by a user-supplied one (§4.1). All other exposures are aggregate (no gene × line information on the new test tier) and method-symmetric. Agents in the register **must not author or approve** a priori lists, candidate sets or decision thresholds; the gene panels and the drug list are produced by unexposed agents at M1 (open question F-4 on the lists already in this document). Whoever writes an a priori list attests in the commit message that they are not in the register.

## Later entries

- 2026-10-26, PI: two hook-liveness probes after restart were blocked (no values seen). Earlier pre-activation probes printed the gene-name header row only.
- 2026-10-26, PI: read reviewer summaries containing all-line aggregate label statistics (common-essential share of positives 86.5%; ASNS essentiality 50% DMEM vs 6% RPMI; per-library positive rates). No gene × line values. The PI therefore does not author a priori lists (F-4); an unexposed agent re-approves them.
