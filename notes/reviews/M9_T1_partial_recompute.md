# T1 partial recompute (PI, 2026-10-28)

The test labels can be read only by the frozen scripts/evaluate_test.py (tripwire), so no non-author agent can recompute C1 from labels. This check recomputes what the unsealed tables allow, from results/m6/test/unsealed/T1/per_lineage.tsv (P1[H1] rows, delta_vs_ref = lineage mean of per-line ΔtAP vs G*):

- lineages: 7; all positive: True
- lineage-equal mean of the 7 lineage deltas = 0.0497345; harness C1 estimate = 0.0497345; |diff| = 4.4e-09
- t-interval (df = 6) from the 7 lineage means = [0.04737, 0.05209]; harness transfer_C1 = [0.04737, 0.05209]
- exact sign-flip p over 2^7 with all 7 positive = 1/128 = 0.0078125; harness = 0.0078125

Not recomputed: the line-bootstrap CI, the 10,000-flip p, C1b and C3. These need per-line test values, which the unsealed bundle does not contain. Further evidence: the rehearsal on validation (same frozen code) reproduced the independently reviewed M6 validation look exactly (C1 = +0.0517).
