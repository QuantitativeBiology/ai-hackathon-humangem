# Guard changelog (user decision D7, 2026-10-27)

Guard-internal changes to frozen code are logged here instead of as DESIGN.md §17 amendments. Guard-internal means test renames, scanner and check-tool fixes, and test-harness plumbing. This applies only while every frozen label, metric, split and list SHA is unchanged and no reportable number changes. Each entry still needs a passing full test suite, `check_frozen.py`, and a FROZEN.tsv refresh with `--amendment A6`. Retired test IDs go in results/test_retirements.tsv with amendment id A6.

| Date | Files | Change | Why | Frozen label/metric/split/list SHAs unchanged? |
|---|---|---|---|---|
