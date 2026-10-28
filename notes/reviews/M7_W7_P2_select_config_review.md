# M7 W7_P2_select look config review (DESIGN 17 A6, non-author)

Verdict: APPROVE, conditional on adding the one FROZEN row below before the look.

- Label-read probe (`head -c 50` on the CRISPR gene-effect file): BLOCKED by the label guard. I did not open the W7_shortlist_1 result/decision files, data/locked/ or any raw label file.
- Spec diff vs look_spec_W7_shortlist_1.json (formatted JSON diff): only look_id, candidates (adds "W7:B6" first), decision_rule and notes differ. Methods, labels, genesets, B, seed, reference and scores file are identical.
- Candidates are [W7:B6, W7:W7-05], simplest first, as required (A11(g): order = simplicity).
- Rule text matches DESIGN W7 step 7 (P2 = best validated idea; within one paired SE -> the simpler), section 9 (one-SE) and A11(g) (gate vs G* before one-SE). RESEARCH_TREE D13 states the same rule and the same candidate order.
- Harness (scripts/evaluate_val.py decide_shortlist): gate = CI lower bound > 0 vs G*. best = max delta_lin among eligible. Each eligible r gets diff = paired mean (best - r) and SE = SD of the bootstrap replicates; within_1se = diff <= SE. selected = first eligible in spec order with within_1se.
  - Both eligible, W7-05 best: B6 diff = D(W7-05 - B6). If it is <= SE, B6 is selected; otherwise B6 fails and W7-05 is selected. This is exactly step 7.
  - Ties go to B6. If B6 is best, it is selected. If only one passes the gate, that one is selected.
  - Caveat, by design (A11(g)): a candidate that fails the gate vs G* cannot be P2. If neither passes, nothing is selected, and the result reports best_point_estimate instead.
- `--print-config`: config_sha = 6b19d5315d864f2ccc1d7fd9a9ae3cef83589d915fa66ec7eb69e46eedfa19ba (matches the expected 6b19d531...). The look_spec input hash is f30ffbbc...
- `--dry-run`: runs (10.7 s) with synthetic labels. It selects nothing, as expected, since synthetic deltas are about -0.0006 and not eligible. The harness runs for the 2-candidate spec, and the winner's-curse term is reported (m = 2).
- `--print-frozen-candidates`: exactly one row is needed. The tab-separated columns are kind=file, path=results/w7/val/look_spec_W7_P2_select.json, sha256=f30ffbbc0e88b651eec7541a83182f998c6e89d8753d2eb160f00b46777485f0, date=2026-10-28, note="W7 validation look: input".
- `--preflight-only` currently FAILS only for that missing FROZEN row ("no file row in committed HEAD:FROZEN.tsv"). The run stops at the first failure, so I cannot confirm the other frozen inputs and code hashes pass; re-run preflight after the row is committed.
- Nothing was committed or edited besides this report.
