"""Frozen project constants (DESIGN.md §5, §7, §11.2)."""

from __future__ import annotations

# §7: fixed seed list. SEEDS[0] drives bootstraps/permutations (spawned), SEEDS[0:5] random baselines,
# SEEDS[0:3] ML/CV shuffles. The split seed is NOT here: it is user-supplied (§4.1) and passed to
# scripts/make_splits.py on the command line.
SEEDS = [11, 23, 37, 41, 53, 67, 79, 83, 97, 101]

# §4.1: the v0/v1 split seed, exposed in review round 1. Any S whose attempt seeds S..S+99 touch
# 20261026..20261125 is refused.
RETIRED_SPLIT_SEED = 20261026
RETIRED_SEED_RANGE = (20261026, 20261125)  # inclusive
MAX_SPLIT_ATTEMPTS = 100

# §5.1 / §6.8: label variants. name -> (source column, comparison, threshold)
#   L-prob : CRISPRGeneDependency > 0.5 (primary unless the kappa rule switches to L-ge05)
#   L-ge05 : Chronos gene effect < -0.5 (sensitivity; primary if kappa(L-prob, L-ge05) < 0.5 on training)
#   L-ge10 : Chronos gene effect < -1   (sensitivity)
#   L-ge06 : Chronos gene effect < -0.6 (secondary only; Human1/Human2-paper protocols, §6.8)
LABELS = {
    "L_prob": ("dep_prob", ">", 0.5),
    "L_ge05": ("gene_effect", "<", -0.5),
    "L_ge06": ("gene_effect", "<", -0.6),
    "L_ge10": ("gene_effect", "<", -1.0),
}
KAPPA_RULE_THRESHOLD = 0.5  # §5.1: L-ge05 becomes primary if kappa(L-prob, L-ge05) < 0.5 on training lines
CE_FIT_FRACTION = 0.90  # §5.3: CE_fit = positive in >= 90% of fitting lines
# §6.8: Human1-paper robustness grid for gene-effect thresholds: 21 steps from -1.0 to 0.0
GE_THRESHOLD_GRID = [round(-1.0 + 0.05 * i, 2) for i in range(21)]

# §11.2: label-derived sources; accessors require purpose in LABEL_DERIVED_PURPOSES.
LABEL_DERIVED_SOURCES = [
    "DeKegel2021_paralogs",
    "Ito2021_paralogs",
    "ProjectScore_Behan2019",
    "CoRe",
    "Hart_CEG_NEG",
    "DEMETER2",
    "DepMap_older_releases",
    "DepMap_common_essential_and_control_lists",
    "DepMap_coessentiality_outside_training_folds",
]
LABEL_DERIVED_PURPOSES = ("reporting", "annotation")

# §4.1 / §4.2 tiers
TIERS = ("train", "validation", "test")

# §11.3 look caps
VALIDATION_LOOK_CAP = 40
VALIDATION_LOOK_CAP_W7 = 20
