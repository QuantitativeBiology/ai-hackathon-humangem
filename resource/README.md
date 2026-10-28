# nextgem P1 resource v1.0.0 (release candidate (P1 per D11; built before T1))

P1 (RESEARCH_TREE D11): generic Human-GEM parent + W3 B-loo:cofactor_pool_biomass (the cofactor/vitamin pool species removed from the biomass reaction MAR13082) on the W2 documented medium M2(f=1) (base medium + supplements + serum model at f = 1), gate-free; no context family, no pruning, no repaired reactions; objective MAR13082.

P1 models differ across lines only by the line's medium signature, so this resource holds **one SBML model per
(parent, medium signature)** and a manifest that maps every line to its model:

- `v1.0.0/medium-<signature>_<parent>-<Human-GEM version>_P1.xml.gz`: SBML L3V1 FBC v2, gzipped. Notes
  record the pipeline, Human-GEM release, medium signature/keys/tiers, W3 variant, git commit, input SHA-256s and
  flags. Closed reactions are kept with bounds (0, 0).
- `models.tsv`: one row per model file (sha256, size, counts, flags, validation gates).
- `manifest.tsv`: one row per (line, parent): the line's own medium key, text, source and stratum (tier),
  fallback flag (imputed / sheet missing / unmappable lines get classic RPMI + 10% FBS), lineage, and its model file.

| Parent | Models | Gates passed | Lines |
|---|---|---|---|
| human2 | 20 | 19 | 1361 |
| human1 | 20 | 20 | 1361 |

A per-line file (`<ModelID>_<parent>-<version>_P1.xml.gz`, the shared model plus the line's notes) is written on
demand: `python scripts/w8_line_sbml.py --lines ACH-000001 --parent human2`; in Python,
`nextgem.resource.load_line("resource", "ACH-000001", "human2")`.

Flags (RESEARCH_TREE D11): W3-D1, W2-SHUF, W2-MOPEN on every model; W3-MEV on every serum medium (mevalonate core
non-lethal). See results/w8/README.md.

Licence: CC BY 4.0, inherited from Human-GEM. Medium bounds come from the W2 documented-media library (vendor sheets; no HPLM/Plasmax composition is used). The resource contains no per-protein values and no DepMap screen data.
