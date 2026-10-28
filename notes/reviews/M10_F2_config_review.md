# M10_F2 look config review (A6, non-author, post hoc)

Verdict: PASS on content; the look cannot start until the 13 FROZEN.tsv rows below are committed.

- Label-guard probe (head -c 50 on the raw CRISPR gene-effect file): BLOCKED by the hook (it also blocked the compound command it sat in). No label or locked files were read.
- Spec vs D18/A15: matches. Single candidate F2, reference H1, P1 descriptive and never selected, m=1 (term 0), L_prob / all / U / none, B=2000, seed=11, workstream and milestone both M10 (A13). Spec keys follow the M6 format. The rule text gives F2-P1 as an exact point estimate with no CI and says no test evaluation follows.
- Notes state "Post hoc / exploratory (DESIGN 17 A15; D18, fixed before any F2 validation scoring)" and carry the A15 blindness statement.
- Scores parquet: 567,468 rows x 6 columns (ModelID, entrez, H1, H2, P1, F2). It has no label columns and no NaNs or duplicate keys. It has 198 lines. Its keys are identical to m6_val_scores.parquet. H1, H2 and P1 are bitwise equal to the M6 columns (np.array_equal). The F2 values are in [0,1]. The harness load_scores check passed in --dry-run.
- Provenance: all 25 pins match the committed files' SHA-256, as do the spec's scores sha (933f7e92...) and look_spec sha (ccda3232...).
- --print-config: config_sha c665f0c5f0899e92769384319d0bfa5b3aaa149b78e131d0d1a3ca89f5a7acd0, as expected.
- --dry-run (synthetic labels, 0 ledger rows, validation_labels_read=false): runs clean. The synthetic F2-H1 delta is -5.8e-5, so the gate is not passed on synthetic labels. This is only a plumbing check.
- --preflight-only FAILS, with the only error "no file row in the committed HEAD:FROZEN.tsv" for these 13 paths:
  results/m10/val/look_spec_M10_F2.json; results/m10/val/f2_val_scores.parquet; results/w7/sweeps/sweeps.tsv; results/w7/sweeps/mev.tsv; results/f2/mev_check.tsv; results/f2/scores_checks.json; scripts/w7_artefact_sweeps.py; scripts/f2_common.py; scripts/f2_scores.py; scripts/f2_score_train.py; scripts/f2_val_spec.py; results/f2/train/summary.json; results/f2/train/contrasts.tsv
  The rows can be taken verbatim from `evaluate_val.py --look M10_F2 --print-frozen-candidates` (type file, sha, date 2026-10-28, note "M10 validation look: input"). They must be committed before the look.
- Caution: the FROZEN commit must not change any of these files, or the sha pins and config_sha break. Re-run --print-config and --preflight-only afterwards.
