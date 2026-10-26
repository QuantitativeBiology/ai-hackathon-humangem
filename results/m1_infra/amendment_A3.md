# Amendment A3: draft for the PI

Add this as one new row of the DESIGN.md §17 amendment table. Commit it together with, or before, the five changed frozen files. Their FROZEN.tsv rows already name `A3`.

```
| 2026-10-26 | 11.1 (content scanner, pre-commit test-keyed rule) | A3: The canonical split file results/splits/lineage_split.tsv (ModelID, tier, label-free covariates; written once by make_splits.py) is exempt from the scanner, at that path only and only while its SHA-256 equals SPLIT_SHA256 and split_meta.json. Every copy, every modified version and every other table keyed by test ModelIDs (e.g. per-line results) is still flagged. Changed files: src/nextgem/scan.py, scripts/scan_label_copies.py, tests/test_scanner.py; and tests/test_frozen.py, tests/test_splits.py (two assertions that assumed no amendment rows and no drawn split). | After the draw, the pre-commit test-keyed rule blocked committing the split itself: an integer covariate column (n_screens) parses as a gene column, so the file looked like a long ModelID x gene table. | After the split draw (639f3a9 → 62000d3); before any label build or scoring; the split is unchanged |
```
