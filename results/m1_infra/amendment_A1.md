# Amendment A1: draft for the PI

This is the text for DESIGN.md §17. It is drafted by the M1 infrastructure agent, as tightened by the independent review (`notes/reviews/M1_infra_review.md` §3).

## 1. Row to paste into the §17 amendment table

Replace the placeholder row `| — | — | — | — | — |` with the line below. Keep it on **one line**: `make_splits.py` accepts post-freeze DESIGN.md changes only as rows of this table, and the placeholder row is the only line that may be removed.

```
| 2026-10-26 | 4.1 (steps 1-3), 12 | A1: The split code did not exist at the design freeze 61b0beb, so V-1's diff check runs against the reviewed M1 code commit C (`make_splits.py --code-commit C`), with freeze ⊑ C ⊑ HEAD, a clean tree, and C..HEAD touching no src/scripts/tests/DESIGN.md/FROZEN.tsv. In addition: (i) DESIGN.md may change after 61b0beb only by rows added to this table; (ii) before drawing, splits.draw/checks and the accepted attempt must equal the functions extracted verbatim from scripts/design_numbers.py at 61b0beb (seeds 1-1000 and S..S+99), together with the 20 reference draws and the 1-1000 acceptance statistics; (iii) the draw inputs are pinned: input files against data/checksums.tsv at 61b0beb; harmonized_MS_CCLE_Gygi.csv and CTRPAUCMatrix.csv, which that file lacks, through the split.inputs_sha256 row of FROZEN.tsv; per-lineage expression/proteomics counts equal to coverage_by_lineage of the frozen JSON; and FROZEN.tsv rows split.* (inputs, line→lineage map, per-lineage counts, draw parameters) plus the SHA-256 values of make_splits.py, splits.py and cohort.py; (iv) production mode refuses --reference-commit ≠ 61b0beb, --repo and NEXTGEM_* root overrides, and sandbox seeds must be ≥ 500000. Deviation recorded: the user stated S = 38 before the freeze commit (RESEARCH_TREE D4), whereas §4.1 step 2 requires after. The design authors did not know S, and S was kept out of all files until 5a8c218. Every post-reveal author (M1 infrastructure, metrics, D5 list re-derivation, this amendment) attests that they did not evaluate draw() on seeds 38-137 and were not told the tiers. Test IDs retired by the M1 review are listed in results/test_retirements.tsv. | The literal V-1 rule cannot pass, because 61b0beb contains no make_splits.py. The pins make the production tiers a function of pre-reveal artefacts only (review M1-M3, §3 a-e). | Before the split draw; no label or draw has been computed |
```

## 2. `results/splits/SEED_SOURCE.md` (commit after C, before the draw)

```
# Split seed source (DESIGN.md §4.1 step 2, amendment A1)

- Seed S = 38, supplied by the user by message (timestamp: <paste the message timestamp>).
- Design freeze commit: 61b0beba6f098c55b62c2881b8bacec96de780db (61b0beb).
- Revealed in the record after the freeze in commit 5a8c218 (RESEARCH_TREE D4). The user stated it before the freeze; see A1.
- Code commit C (reviewed M1 infrastructure): <paste C>.
```

`make_splits.py` checks two things in this file: that it is committed and unmodified, and that it contains `38` and `61b0beb`.

## 3. Order of operations

1. Paste the A1 row into DESIGN.md §17. Commit it **together with** all M1 code and `FROZEN.tsv`: this is commit C.
   - Stage everything under src/, scripts/ and tests/ (`make_splits.py` refuses untracked files there), plus `FROZEN.tsv`, `results/m1_infra/`, `results/test_runs*` and `results/test_retirements.tsv`.
   - If any frozen file is edited before committing, refresh its row first: `.venv/bin/python scripts/check_frozen.py --add <path>`.
2. Commit `results/splits/SEED_SOURCE.md`.
3. Draw the split: `.venv/bin/python scripts/make_splits.py --seed 38 --freeze-commit 61b0beb --code-commit <C>`.
   - The `split_started` ledger row is committed automatically.
   - Then commit `results/splits/`.
4. Build the labels: `.venv/bin/python scripts/build_labels.py`.
   - The kappa ledger row is committed automatically.
   - Then commit `FROZEN.tsv` and `results/labels/`.
