# Guard changelog (user decision D7, 2026-10-27)

Guard-internal changes to frozen code are logged here instead of as DESIGN.md §17 amendments. Guard-internal means test renames, scanner and check-tool fixes, and test-harness plumbing. This applies only while every frozen label, metric, split and list SHA is unchanged and no reportable number changes. Each entry still needs a passing full test suite, `check_frozen.py`, and a FROZEN.tsv refresh with `--amendment A6`. Retired test IDs go in results/test_retirements.tsv with amendment id A6.

| Date | Files | Change | Why | Frozen label/metric/split/list SHAs unchanged? |
|---|---|---|---|---|
| 2026-10-27 | tests/test_gurobi_pool.py | Two pre-review WIP tests renamed and strengthened: test_env_override → test_size_is_global_not_per_caller; test_slot_blocks_other_process → test_slot_blocks_other_process_and_crash_releases. Retired in results/test_retirements.tsv under A6. | M3 W4 review fixes (global pool size, crash release) | yes |
| 2026-10-27 | src/nextgem/scan.py, tests/test_scanner.py | The bare-integer ModelID heuristic (IDs that lost their zero padding) now needs ≥ 10 rows (MIN_BARE_INT_ROWS). New test test_small_integer_columns_are_not_bare_modelids. | False positive: a 2-row W5 diagnostic table (results/w5/blocked_W5a_pid30_human1.tsv) with integer group/rank columns 76, 612 was read as test ModelIDs ACH-000076/ACH-000612; the file has no ModelIDs. | yes |
