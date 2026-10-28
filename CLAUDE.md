# nextgem

Constraint-based modelling with the Human-GEM genome-scale metabolic model, integrated with DepMap (CRISPR, expression, CN, mutations).

## Environment
- Python: always use `/data/tmpA/nextgem/.venv/bin/python` (conda prefix env, Python 3.11). Do not use the base or `py` conda envs (no cobra).
- Pinned versions in `requirements.txt` (regenerate with `.venv/bin/python -m pip freeze --exclude-editable > requirements.txt` after installing packages).
- Solver: Gurobi 13 (WLS academic license at `~/gurobi.lic`) is cobra's default; GLPK is the fallback. The WLS license needs internet access at solve time.
- Jupyter kernel: "Python (nextgem)".

## Layout
- `src/nextgem/` — project package, installed editable (`import nextgem`). `model.py`: `load_model()`, `set_hams_medium()`.
- `scripts/` — runnable scripts (e.g. `scripts/fba_sanity_check.py`).
- `notebooks/` — exploratory notebooks.
- `results/` — analysis outputs, committed to git (bulky intermediates go under `results/**/tmp/` or `*.large.*`, which are ignored and must be regenerable).
- Research record: `PROMPT.md` (goal), `RESEARCH_TREE.md` (status, decisions — resume from here), `DESIGN.md` (pre-registered design), `LESSONS.md`, `DATA.md`.
- `data/` — git-ignored inputs; provenance in `data/README.md`. Human-GEM 2.1.1 in `data/Human-GEM-main/`, DepMap Public 26Q3 in `data/depmap/`. Treat `data/` as read-only.

## Modelling notes
- Load `data/Human-GEM-main/model/Human-GEM.xml` once (~15 s); use `with model:` or `model.copy()` for perturbations.
- The raw model has all 1,654 exchanges open (cells can import ATP etc.) — always set a medium first. `set_hams_medium()` opens the 44 inputs of the `GR` task in `metabolicTasks_Essential.txt`.
- Biomass reaction: `MAR13082`. With uptake bounds of 1000, growth values are in arbitrary units.
- Expected on Ham's medium: growth ≈ 74.76 for Human2 (Human1 ≈ 65.57), 0 anaerobically, CAD/UMPS KOs lethal (`scripts/fba_sanity_check.py`).
- Human-GEM's own essentiality/tINIT pipelines live in `data/Human-GEM-main/code/` (e.g. `code/test/gradedEssentiality.py`).

## DepMap notes
- Gene columns are `SYMBOL (EntrezID)`; join to Human-GEM via `geneEntrezID` in `model/genes.tsv`.
- `CRISPRGeneEffect.csv` is indexed by `ModelID`; expression has several rows per model — filter `IsDefaultEntryForModel == "Yes"`.
- Files are large (up to ~630 MB); read only needed columns (`usecols`) and write intermediates as Parquet.

## Labels and leakage (see DESIGN.md)
- Never read raw DepMap label files (CRISPR gene effect/dependency, screen gene effect, common-essential/non-essential lists) directly. Labels are accessed only through `src/nextgem/labels.py` inside project scripts, fold-scoped.
- A PreToolUse hook (`.claude/hooks/guard_labels.py`) blocks direct reads; never work around it. New agent sessions start with a probe that a label read is blocked.
- Never run the real split seed or write label summaries into committed files outside the pre-registered pipeline. Test-tier results stay sealed in `data/locked/` until the design says otherwise.
- Don't leave data or label copies in the scratchpad; write intermediates under `results/` via scripts.

## Housekeeping
- `/data/tmpA` is a shared NFS volume that was ~98% full (~100 GB free); keep outputs compact.
- The label hook splits shell commands on newlines before stripping quotes: use one `-m` per paragraph in commit messages and avoid protected path names in them.
- Pin BLAS/OpenMP to one thread in every worker and analysis script (`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`, or `threadpoolctl`): an unpinned numpy pilot spun ~40 threads on this 48-core shared host and burned ~7 CPU-h (W6, 2026-10-28).
- Background agents stopped by an API limit resume intact via SendMessage; have them check their own background processes first so they don't start duplicate runs.
