"""W2 media library: formulation sheets -> Human-GEM exchange bounds (DESIGN.md §9 W2, PI decisions 23, 27).

Label-free. Nothing here opens a label-bearing file.

Layers
------
1. **Sheets** (`parse_vendor_sheet`, `parse_hplm`, `parse_plasmax`): each formulation on disk (four Gibco/ATCC
   sheets, HPLM = Cantor 2017 Table S1 sheet "HPLM formulation", Plasmax = Vande Voorde 2019 Table S1) is read
   into rows of (raw component, concentration).
2. **Component map** (`RAW_COMPONENTS`): every raw component string -> canonical metabolites with stoichiometry
   (salts split into ions), or a non-nutrient class (buffer, dye, protein, antibiotic). Explicit and versioned
   (`MAPPING_VERSION`).
3. **Exchange map** (`EXCHANGES`): canonical metabolite -> exchange reaction id, the same id in Human1 v1.19.0 and
   Human2 v2.1.1; `verify_exchange_map(model)` checks each id and its metabolite name per parent. Canonical
   metabolites with no exchange in either parent (selenite, nitrate, 4-aminobenzoate, ...) are listed as
   ``no_exchange`` and reported per medium.
4. **Formulation strings** (`parse_formulation`): DepMap medium text -> base medium (sheet, mixture or
   unmappable class), serum %, and parsed supplements (DESIGN parser rules).
5. **Line assignment** (`assign_line_media`): screen > onboarded > imputed tiers, attribution set, fallback and
   unmappable strata (same tier logic as scripts/design_numbers.py, whose counts are asserted in tests).
6. **Bounds** (`uptake_bounds`): the variant ladder M-open, M0, M1, M2(f), M3, the +g transporter gate and the
   HPLM/Plasmax counterfactuals, from a `Medium` plus a per-parent `ParentConstraints` (floor and req_c).

Uptake-bound scheme (availability; DESIGN W2 variant ladder)
-------------------------------------------------------------
- M-open: every exchange of the parent at 1000 (negative control).
- M0: Ham's GR-task inputs at 1000 (``nextgem.model.hams_medium``).
- M1: every mapped component of the documented medium at 1000, every other uptake 0; background O2 and H2O at
  1000; cystine instead of cysteine (A9: cysteine opens only when 2-mercaptoethanol and cystine are both present;
  the sheet's documented free cysteine is used only by the sensitivity rule, see `CYSTEINE_RULE`); plus the biomass micronutrient floor at ``FLOOR_MULT x req_c`` (1000 when the
  medium itself supplies the component).
- M2(f): M1 plus, if the medium contains serum, each synthesisable serum component at ``f x req_c`` (1000 if the
  medium itself supplies it) and the serum trace metals at 1000. Serum-free media get no serum.
- M3: the best M2 plus glucose, glutamine and O2 at calibrated rates (``calibration`` mapping); secretion free.
- +g (gate closed): aspartate and glutamate uptake 0. Gate open = no change.
- M-HPLM, M-Plasmax: the HPLM/Plasmax composition treated as a serum-containing medium under M2(f).

Secretion is never constrained here (cobra ``model.medium`` only sets uptake bounds).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from nextgem import paths

MAPPING_VERSION = "w2-media-map-v2 (2026-10-27; A9 fixes)"

# ----------------------------------------------------------------------------------------------- fixed parameters
UPTAKE_OPEN = 1000.0
FLOOR_MULT = 10.0  # DESIGN W2: floor bound = 10 x req_c (raised to 100 x req_c if the non-binding precheck fails)
FLOOR_MULT_FALLBACK = 100.0
SERUM_F_LEVELS = (0.001, 0.01, 0.1, 1.0)
PRECHECK_REL = 1e-3  # WT growth > 1e-3 x M0 growth (DESIGN W2 prechecks)
FLOOR_NONBINDING_RTOL = 1e-6
PFBA_FRACTION = 0.99  # req_c for synthesisable serum components: pFBA at >= 99% of WT growth
CALIBRATION_RANGE = (0.01, 0.06)  # h^-1, DESIGN W2 M3 calibration window
GATE_GENES = ("SLC1A1", "SLC1A2", "SLC1A3")
GATE_ENTREZ = {"SLC1A1": 6505, "SLC1A2": 6506, "SLC1A3": 6507}
GATE_THRESHOLD = 1.0  # log2(TPM+1) >= 1 in any of the three opens the gate
GATED = ("aspartate", "glutamate")
BACKGROUND = ("O2", "H2O")  # environment, not formulation: open in every documented variant (as in M0)
IMPUTED_TEXT = "RPMI + 10% FBS"  # DESIGN W2: imputed / fallback medium
# DESIGN M1 / amendment A9(b): 'design_literal' (default) = cysteine uptake closed unless 2-mercaptoethanol is
# present, in every medium (free cysteine autoxidises to cystine in medium; Ishii & Bannai 1985, 1987).
# 'documented' (sensitivity) = cysteine also open where the formulation lists free cysteine (DMEM/F12, HPLM, Plasmax).
CYSTEINE_RULE = "design_literal"
CYSTEINE_RULES = ("design_literal", "documented")

# Serum model (DESIGN W2). f-scaled synthesisable components; trace metals at availability 1000.
SERUM_F_SCALED = ("palmitate", "stearate", "oleate", "cholesterol", "carnitine", "ethanolamine")
SERUM_TRACE_METALS = ("Fe2+", "Fe3+", "Zn2+", "Cu2+")  # selenium: no selenite exchange in either parent
SERUM_UNREPRESENTABLE = ("selenite",)

# ------------------------------------------------------------------------------------------ canonical -> exchange
# Same reaction id in both parents (checked by verify_exchange_map); value = (exchange id, expected metabolite name
# in the parent, lower-case).
EXCHANGES: dict[str, tuple[str, str]] = {
    # amino acids
    "glycine": ("MAR09067", "glycine"), "alanine": ("MAR09061", "alanine"), "arginine": ("MAR09066", "arginine"),
    "asparagine": ("MAR09062", "asparagine"), "aspartate": ("MAR09070", "aspartate"),
    "cystine": ("MAR09363", "cystine"), "cysteine": ("MAR09065", "cysteine"), "glutamate": ("MAR09071", "glutamate"),
    "glutamine": ("MAR09063", "glutamine"), "histidine": ("MAR09038", "histidine"),
    "hydroxyproline": ("MAR08386", "trans-4-hydroxy-l-proline"), "isoleucine": ("MAR09039", "isoleucine"),
    "leucine": ("MAR09040", "leucine"), "lysine": ("MAR09041", "lysine"), "methionine": ("MAR09042", "methionine"),
    "phenylalanine": ("MAR09043", "phenylalanine"), "proline": ("MAR09068", "proline"),
    "serine": ("MAR09069", "serine"), "threonine": ("MAR09044", "threonine"),
    "tryptophan": ("MAR09045", "tryptophan"), "tyrosine": ("MAR09064", "tyrosine"), "valine": ("MAR09046", "valine"),
    # vitamins and cofactors
    "biotin": ("MAR09109", "biotin"), "choline": ("MAR09083", "choline"),
    "pantothenate": ("MAR09145", "pantothenate"), "folate": ("MAR09146", "folate"),
    "nicotinamide": ("MAR09378", "nicotinamide"), "pyridoxine": ("MAR09144", "pyridoxine"),
    "pyridoxal": ("MAR09400", "pyridoxal"), "riboflavin": ("MAR09143", "riboflavin"),
    "thiamin": ("MAR09159", "thiamin"), "cobalamin": ("MAR09269", "aquacob(iii)alamin"),
    "inositol": ("MAR09361", "inositol"), "ascorbate": ("MAR09158", "ascorbate"),
    # ions
    "Ca2+": ("MAR09082", "ca2+"), "Mg2+": ("MAR13072", "mg2+"), "K+": ("MAR09081", "k+"), "Na+": ("MAR09077", "na+"),
    "chloride": ("MAR09150", "chloride"), "HCO3-": ("MAR09078", "hco3-"), "Pi": ("MAR09072", "pi"),
    "sulfate": ("MAR09074", "sulfate"), "Fe2+": ("MAR09076", "fe2+"), "Fe3+": ("MAR09096", "fe3+"),
    "Cu2+": ("MAR13073", "cu2+"), "Zn2+": ("MAR09080", "zinc"), "NH4+": ("MAR09073", "nh3"),
    # carbon sources and other small molecules
    "glucose": ("MAR09034", "glucose"), "pyruvate": ("MAR09133", "pyruvate"), "GSH": ("MAR09351", "gsh"),
    "hypoxanthine": ("MAR09358", "hypoxanthine"), "thymidine": ("MAR09423", "thymidine"),
    "uridine": ("MAR09439", "uridine"), "uracil": ("MAR09437", "uracil"), "urate": ("MAR09075", "urate"),
    "adenine": ("MAR09253", "adenine"), "adenosine": ("MAR09254", "adenosine"),
    "linoleate": ("MAR09035", "linoleate"), "linolenate": ("MAR09036", "linolenate"),
    "lipoate": ("MAR09167", "lipoic acid"), "putrescine": ("MAR11426", "putrescine"),
    "2-hydroxybutyrate": ("MAR09216", "2-hydroxybutyrate"), "3-hydroxybutyrate": ("MAR09134", "(r)-3-hydroxybutanoate"),
    "3-hydroxyisobutyrate": ("MAR10185", "3-hydroxyisobutyrate"), "acetate": ("MAR09086", "acetate"),
    "acetone": ("MAR09243", "acetone"), "acetylglycine": ("MAR10196", "acetyl-glycine"),
    "2-aminobutyrate": ("MAR10206", "(s)-2-aminobutanoate"), "betaine": ("MAR09341", "betaine"),
    "carnitine": ("MAR09292", "l-carnitine"), "acetylcarnitine": ("MAR09920", "o-acetylcarnitine"),
    "citrate": ("MAR09286", "citrate"), "citrulline": ("MAR09201", "citrulline"), "creatine": ("MAR09290", "creatine"),
    "creatinine": ("MAR09460", "creatinine"), "formate": ("MAR09318", "formate"), "fructose": ("MAR09139", "fructose"),
    "galactose": ("MAR09140", "galactose"), "glycerol": ("MAR09085", "glycerol"), "lactate": ("MAR09135", "l-lactate"),
    "malonate": ("MAR09164", "malonate"), "ornithine": ("MAR09087", "ornithine"), "succinate": ("MAR09415", "succinate"),
    "taurine": ("MAR09418", "taurine"), "urea": ("MAR09438", "urea"), "homocysteine": ("MAR09853", "homocysteine"),
    "5-oxoproline": ("MAR09025", "5-oxoproline"), "carnosine": ("MAR08644", "carnosine"),
    # serum lipids and floor components
    "palmitate": ("MAR00611", "palmitate"), "stearate": ("MAR00639", "stearate"), "oleate": ("MAR00650", "oleate"),
    "cholesterol": ("MAR09285", "cholesterol"), "ethanolamine": ("MAR09084", "ethanolamine"),
    "phosphoethanolamine": ("MAR09849", "ethanolamine-phosphate"),
    "alpha-tocopherol": ("MAR09151", "alpha-tocopherol"), "gamma-tocopherol": ("MAR09153", "gamma-tocopherol"),
    "retinoate": ("MAR09404", "retinoate"), "retinol": ("MAR09147", "retinol"),
    # hormones in supplements
    "cortisol": ("MAR09293", "cortisol"), "estradiol": ("MAR09314", "estradiol-17beta"),
    "triiodothyronine": ("MAR09427", "triiodothyronine"), "testosterone": ("MAR09429", "testosterone"),
    "progesterone": ("MAR09393", "progesterone"), "corticosterone": ("MAR09294", "corticosterone"),
    # background
    "O2": ("MAR09048", "o2"), "H2O": ("MAR09047", "h2o"),
}
# canonical metabolites that are real nutrients but have no exchange in either parent (representability gaps)
NO_EXCHANGE = {"selenite", "nitrate", "4-aminobenzoate", "Mn2+", "vanadate", "methyl acetoacetate", "BrdU",
               "2-mercaptoethanol"}
# non-nutrient classes (listed, never mapped)
NON_NUTRIENT = {"buffer", "dye", "protein", "antibiotic", "proprietary"}

# ----------------------------------------------------------------------------------------------- raw components
# raw component string (sheet text, stripped) -> list of (canonical, stoichiometry) or a non-nutrient class name.
RAW_COMPONENTS: dict[str, list[tuple[str, float]] | str] = {
    # ---- Gibco/ATCC sheets
    "Glycine": [("glycine", 1)], "L-Alanine": [("alanine", 1)],
    "L-Arginine": [("arginine", 1)], "L-Arginine hydrochloride": [("arginine", 1), ("chloride", 1)],
    "L-Asparagine": [("asparagine", 1)], "L-Asparagine-H2O": [("asparagine", 1)],
    "L-Asparagine (freebase)": [("asparagine", 1)], "L-Aspartic acid": [("aspartate", 1)],
    "L-Cystine 2HCl": [("cystine", 1), ("chloride", 2)],
    "L-Cysteine hydrochloride-H2O": [("cysteine", 1), ("chloride", 1)],
    "L-Glutamic Acid": [("glutamate", 1)], "L-Glutamine": [("glutamine", 1)],
    "L-Histidine": [("histidine", 1)], "L-Histidine hydrochloride-H2O": [("histidine", 1), ("chloride", 1)],
    "L-Hydroxyproline": [("hydroxyproline", 1)], "L-Isoleucine": [("isoleucine", 1)], "L-Leucine": [("leucine", 1)],
    "L-Lysine hydrochloride": [("lysine", 1), ("chloride", 1)], "L-Methionine": [("methionine", 1)],
    "L-Phenylalanine": [("phenylalanine", 1)], "L-Proline": [("proline", 1)], "L-Serine": [("serine", 1)],
    "L-Threonine": [("threonine", 1)], "L-Tryptophan": [("tryptophan", 1)],
    "L-Tyrosine disodium salt dihydrate": [("tyrosine", 1), ("Na+", 2)],
    "L-Tyrosine disodium salt": [("tyrosine", 1), ("Na+", 2)], "L-Valine": [("valine", 1)],
    "Biotin": [("biotin", 1)], "Choline chloride": [("choline", 1), ("chloride", 1)],
    "D-Calcium pantothenate": [("pantothenate", 2), ("Ca2+", 1)], "Folic Acid": [("folate", 1)],
    "Niacinamide": [("nicotinamide", 1)], "Para-Aminobenzoic Acid": [("4-aminobenzoate", 1)],
    "Pyridoxine hydrochloride": [("pyridoxine", 1), ("chloride", 1)],
    "Pyridoxal hydrochloride": [("pyridoxal", 1), ("chloride", 1)], "Riboflavin": [("riboflavin", 1)],
    "Thiamine hydrochloride": [("thiamin", 1), ("chloride", 2)],  # thiamine chloride hydrochloride (MW 337) "Vitamin B12": [("cobalamin", 1)],
    "i-Inositol": [("inositol", 1)],
    "Calcium nitrate (Ca(NO3)2 4H2O)": [("Ca2+", 1), ("nitrate", 2)],
    "Calcium Chloride (CaCl2) (anhyd.)": [("Ca2+", 1), ("chloride", 2)],
    "Magnesium Sulfate (MgSO4) (anhyd.)": [("Mg2+", 1), ("sulfate", 1)],
    "Magnesium Chloride (anhydrous)": [("Mg2+", 1), ("chloride", 2)],
    "Potassium Chloride (KCl)": [("K+", 1), ("chloride", 1)], "Potassium Nitrate (KNO3)": [("K+", 1), ("nitrate", 1)],
    "Sodium Bicarbonate (NaHCO3)": [("Na+", 1), ("HCO3-", 1)], "Sodium Chloride (NaCl)": [("Na+", 1), ("chloride", 1)],
    "Sodium Phosphate dibasic (Na2HPO4) anhydrous": [("Na+", 2), ("Pi", 1)],
    "Sodium Phosphate monobasic (NaH2PO4-H2O)": [("Na+", 1), ("Pi", 1)],
    "Ferric Nitrate (Fe(NO3)3\"9H2O)": [("Fe3+", 1), ("nitrate", 3)],
    "Ferrous sulfate (FeSO4-7H2O)": [("Fe2+", 1), ("sulfate", 1)],
    "Cupric sulfate (CuSO4-5H2O)": [("Cu2+", 1), ("sulfate", 1)],
    "Zinc sulfate (ZnSO4-7H2O)": [("Zn2+", 1), ("sulfate", 1)],
    "Sodium Selenite (Na2SeO3-5H20)": [("Na+", 2), ("selenite", 1)],
    "D-Glucose (Dextrose)": [("glucose", 1)], "Glutathione (reduced)": [("GSH", 1)], "HEPES": "buffer",
    "Phenol Red": "dye", "Sodium Pyruvate": [("pyruvate", 1), ("Na+", 1)],
    "Hypoxanthine Na": [("hypoxanthine", 1), ("Na+", 1)], "Linoleic Acid": [("linoleate", 1)],
    "Lipoic Acid": [("lipoate", 1)], "Putrescine 2HCl": [("putrescine", 1), ("chloride", 2)],
    "Thymidine": [("thymidine", 1)],
    # ---- HPLM (Cantor 2017 Table S1, sheet "HPLM formulation"; concentrations in uM)
    "Glucose": [("glucose", 1)], "Alanine": [("alanine", 1)], "Arginine": [("arginine", 1)],
    "Asparagine": [("asparagine", 1)], "Aspartate": [("aspartate", 1)], "Cysteine": [("cysteine", 1)],
    "Cystine": [("cystine", 1)], "Glutamate": [("glutamate", 1)], "Glutamine": [("glutamine", 1)],
    "Histidine": [("histidine", 1)], "Isoleucine": [("isoleucine", 1)], "Leucine": [("leucine", 1)],
    "Lysine": [("lysine", 1)], "Methionine": [("methionine", 1)], "Phenylalanine": [("phenylalanine", 1)],
    "Proline": [("proline", 1)], "Serine": [("serine", 1)], "Threonine": [("threonine", 1)],
    "Tryptophan": [("tryptophan", 1)], "Tyrosine": [("tyrosine", 1)], "Valine": [("valine", 1)],
    "Choline": [("choline", 1)], "Folate": [("folate", 1)], "myo-Inositol": [("inositol", 1)],
    "p-Aminobenzoate": [("4-aminobenzoate", 1)], "Pantothenate": [("pantothenate", 1)],
    "Pyridoxine": [("pyridoxine", 1)], "Thiamine": [("thiamin", 1)], "Vitamin B-12": [("cobalamin", 1)],
    "2-hydroxybutyrate": [("2-hydroxybutyrate", 1)], "3-hydroxybutyrate": [("3-hydroxybutyrate", 1)],
    "4-hydroxyproline": [("hydroxyproline", 1)], "Acetate": [("acetate", 1)], "Acetone": [("acetone", 1)],
    "Acetylglycine": [("acetylglycine", 1)], "Alpha-aminobutyrate": [("2-aminobutyrate", 1)],
    "Betaine": [("betaine", 1)], "Carnitine": [("carnitine", 1)], "Citrate": [("citrate", 1)],
    "Citrulline": [("citrulline", 1)], "Creatine": [("creatine", 1)], "Creatinine": [("creatinine", 1)],
    "Formate": [("formate", 1)], "Fructose": [("fructose", 1)], "Galactose": [("galactose", 1)],
    "Glutathione": [("GSH", 1)], "Glycerol": [("glycerol", 1)], "Hypoxanthine": [("hypoxanthine", 1)],
    "Lactate": [("lactate", 1)], "Malonate": [("malonate", 1)], "Ornithine": [("ornithine", 1)],
    "Pyruvate": [("pyruvate", 1)], "Succinate": [("succinate", 1)], "Taurine": [("taurine", 1)],
    "Urea": [("urea", 1)], "Uric acid": [("urate", 1)], "CaCl2": [("Ca2+", 1), ("chloride", 2)],
    "KCl": [("K+", 1), ("chloride", 1)], "MgCl2": [("Mg2+", 1), ("chloride", 2)], "MgSO4": [("Mg2+", 1), ("sulfate", 1)],
    "NaCl": [("Na+", 1), ("chloride", 1)], "NaHCO3": [("Na+", 1), ("HCO3-", 1)], "Na2HPO4": [("Na+", 2), ("Pi", 1)],
    "Ca(NO3)2.4H2O": [("Ca2+", 1), ("nitrate", 2)], "NH4Cl": [("NH4+", 1), ("chloride", 1)], "Phenol red": "dye",
    "Niacinamide ": [("nicotinamide", 1)], "Riboflavin ": [("riboflavin", 1)],
    # ---- Plasmax (Vande Voorde 2019 Table S1; uM; footnote digits stripped)
    "L-Glutamate": [("glutamate", 1)], "L-Cysteine": [("cysteine", 1)], "alpha-Aminobutyrate": [("2-aminobutyrate", 1)],
    "L-Citrulline": [("citrulline", 1)], "L-Cystine": [("cystine", 1)], "L-Homocysteine": [("homocysteine", 1)],
    "4-Hydroxy-L-proline": [("hydroxyproline", 1)], "L-Ornitine": [("ornithine", 1)],
    "L-Pyroglutamate": [("5-oxoproline", 1)], "L-Acetyl glycine": [("acetylglycine", 1)],
    "L-Carnosine": [("carnosine", 1)], "N-Trimethylglycine (betaine)": [("betaine", 1)],
    "Acetyl carnitine": [("acetylcarnitine", 1)], "D-Glucose": [("glucose", 1)],
    "3-Hydroxyisobutyrate": [("3-hydroxyisobutyrate", 1)], "Methyl acetoacetate": [("methyl acetoacetate", 1)],
    "Uracil": [("uracil", 1)], "Urate": [("urate", 1)], "Uridine": [("uridine", 1)],
    "Ammonium Chloride": [("NH4+", 1), ("chloride", 1)], "Calcium Chloride": [("Ca2+", 1), ("chloride", 2)],
    "Magnesium Sulfate": [("Mg2+", 1), ("sulfate", 1)], "Potassium Chloride": [("K+", 1), ("chloride", 1)],
    "Sodium Bicarbonate": [("Na+", 1), ("HCO3-", 1)], "Sodium Chloride": [("Na+", 1), ("chloride", 1)],
    "Sodium Phosphate monobasic": [("Na+", 1), ("Pi", 1)],
    "Ammonium Metavanadate": [("NH4+", 1), ("vanadate", 1)], "Cupric Sulfate": [("Cu2+", 1), ("sulfate", 1)],
    "Ferric Nitrate": [("Fe3+", 1), ("nitrate", 3)], "Ferric Sulfate": [("Fe3+", 2), ("sulfate", 3)],
    "Manganous Chloride": [("Mn2+", 1), ("chloride", 2)], "Sodium Selenite": [("Na+", 2), ("selenite", 1)],
    "Zinc Sulfate": [("Zn2+", 1), ("sulfate", 1)], "Ascorbate": [("ascorbate", 1)], "D-Biotin": [("biotin", 1)],
    "D-Pantothenic acid hemicalcium": [("pantothenate", 1), ("Ca2+", 0.5)], "Vitamin B12": [("cobalamin", 1)],
    "L-Lysine": [("lysine", 1)], "L-Tyrosine": [("tyrosine", 1)], "2-Hydroxybutyrate": [("2-hydroxybutyrate", 1)],
    "3-Hydroxybutyrate": [("3-hydroxybutyrate", 1)], "Calcium Nitrate": [("Ca2+", 1), ("nitrate", 2)],
    "Magnesium Chloride": [("Mg2+", 1), ("chloride", 2)],
}
# Notes recorded in exchange_map.tsv for mappings that are not literal.
MAPPING_NOTES = {
    "cobalamin": "cyanocobalamin (vitamin B12) represented by the parent's only cobalamin exchange, "
                 "aquacob(III)alamin (as in the GR task)",
    "Pi": "phosphate salts -> Pi",
    "NH4+": "ammonium -> NH3 exchange (the parents carry ammonia as NH3)",
    "hydroxyproline": "L-hydroxyproline / 4-hydroxy-L-proline -> trans-4-hydroxy-L-proline",
    "2-aminobutyrate": "alpha-aminobutyrate -> (S)-2-aminobutanoate",
    "3-hydroxybutyrate": "3-hydroxybutyrate -> (R)-3-hydroxybutanoate",
    "5-oxoproline": "pyroglutamate = 5-oxoproline",
    "lactate": "lactate -> L-lactate",
    "selenite": "no selenite exchange in either parent (selenate MAR09169 is a different species): representability gap",
    "nitrate": "no nitrate exchange in either parent",
    "4-aminobenzoate": "no 4-aminobenzoate exchange in either parent",
    "Mn2+": "no manganese exchange in either parent",
    "vanadate": "no vanadate exchange in either parent",
    "methyl acetoacetate": "no exchange in either parent",
}

SHEET_FILES = {  # base medium -> (file in data/media, sheet name or None)
    "RPMI": ("rpmi1640_ATCC.xlsx", None),
    "DMEM": ("DMEM_HighGlucose_NoGlutamine_NoPhenolRed.xlsx", None),
    "DMEM:F12": ("DMEM_F12_NoGlutamine_NoHEPES.xlsx", None),
    "IMDM": ("IMDM_NoPhenolRed.xlsx", None),
}
HPLM_FILE, HPLM_SHEET = "mmc1.xlsx", "HPLM formulation"
CLASSIC_SHEET = "Classic synthetic media"  # Cantor 2017 Table S1: BME, MEM, DMEM 11965, RPMI 1640 11875 (uM)
RPMI_ATCC = "RPMI-ATCC"  # the on-disk ATCC 30-2001 sheet (4.5 g/L glucose, 10 mM HEPES, 1 mM pyruvate)
# A9(d): the RPMI base is classic RPMI 1640 (Gibco 11875; Moore 1967): the ATCC sheet without pyruvate and HEPES,
# with glucose and bicarbonate at the 11875 values taken from Cantor 2017 Table S1 "Classic synthetic media".
RPMI_CLASSIC_DROP = ("Sodium Pyruvate", "HEPES")
RPMI_CLASSIC_FROM_CANTOR = {"D-Glucose (Dextrose)": "Glucose", "Sodium Bicarbonate (NaHCO3)": "NaHCO3"}
# A9(a): sheet errors corrected at parse time (recorded in the `correction` column of compositions.tsv).
SHEET_CORRECTIONS = {
    ("IMDM", "L-Cystine 2HCl"): dict(mw=313.0, note="IMDM sheet uses MW 240 (free cystine) for cystine 2HCl; corrected "
                                               "to the salt MW 313 -> 292 uM (Gibco IMDM; A9a)"),
}
PLASMAX_FILE = "aau7314_sm.pdf"

# Standard-formulation defaults for the on-disk "NoGlutamine" sheets (DESIGN ambiguity, flagged
# `gln_standard_default`): the DepMap text names plain DMEM / DMEM:F12, whose standard Gibco formulations
# (11965, 11320) contain 4 mM and 2.5 mM L-glutamine. Sensitivity: `add_standard_glutamine=False`.
STANDARD_GLN_UM = {"DMEM": 4000.0, "DMEM:F12": 2500.0}

# ------------------------------------------------------------------------------------------------- media classes
SHEET_BASES = {"RPMI", "DMEM", "DMEM:F12", "MEM", "EMEM", "IMDM", "F12", "McCoy's 5A", "L-15"}  # = design_numbers


def media_dir() -> Path:
    return paths.data_root() / "media"


# =========================================================================================== 1. sheet parsing
def _canon_rows(raw: str) -> list[tuple[str, float]] | str:
    key = raw if raw in RAW_COMPONENTS else raw.strip()
    if key not in RAW_COMPONENTS:
        raise KeyError(f"raw component not in RAW_COMPONENTS: {raw!r}")
    return RAW_COMPONENTS[key]


def parse_vendor_sheet(path: Path, medium: str) -> pd.DataFrame:
    """Gibco/ATCC formulation sheet -> rows (medium, section, raw_component, mw, mg_per_L, uM, source_row)."""
    df = pd.read_excel(path, header=None)
    assert str(df.iloc[0, 0]).strip() == "Components" and "mg/L" in str(df.iloc[0, 2]), f"unexpected header in {path}"
    rows, section = [], None
    for i in range(1, len(df)):
        name, mw, mgl, mm = df.iloc[i, 0], df.iloc[i, 1], df.iloc[i, 2], df.iloc[i, 3]
        if pd.isna(name):
            continue
        name = str(name).strip()
        if pd.isna(mw) and pd.isna(mgl):
            section = name
            continue
        mw, mgl = float(mw), float(mgl)
        um = float(mm) * 1000 if pd.notna(mm) else mgl / mw * 1000
        corr = SHEET_CORRECTIONS.get((medium, name))
        if corr:
            mw = corr["mw"]
            um = mgl / mw * 1000
        rows.append(dict(medium=medium, section=section, raw_component=name, mw=mw, mg_per_L=mgl, uM=um,
                         source=f"{path.name}:row{i + 1}", correction=corr["note"] if corr else ""))
    return pd.DataFrame(rows)


def parse_hplm(path: Path) -> pd.DataFrame:
    """HPLM (Cantor 2017 Table S1, sheet 'HPLM formulation'), concentrations in uM."""
    df = pd.read_excel(path, sheet_name=HPLM_SHEET, header=None)
    assert "All concentrations in uM" in str(df.iloc[4, 0])
    rows, section = [], "Glucose"
    for i in range(6, len(df)):
        name, conc = df.iloc[i, 0], df.iloc[i, 1]
        if pd.isna(name):
            if rows:  # first blank row after the table ends it (ion summary follows)
                break
            continue
        name = str(name).strip()
        if name == "Small ion":
            break
        if pd.isna(conc):
            section = name
            continue
        rows.append(dict(medium="HPLM", section=section, raw_component=name, mw=math.nan, mg_per_L=math.nan,
                         uM=float(conc), source=f"{path.name}:{HPLM_SHEET}:row{i + 1}"))
    return pd.DataFrame(rows)


_PLASMAX_LINE = re.compile(r"^\s*(?:(?:Proteinogenic|Non-proteinogenic) Amino Acids|Amino Acids Derivatives|Other "
                           r"Components|Inorganic Salts|Trace Elements|Vitamins)?\s*(?P<name>.+?)\s+(?P<p>NA|[\d.]+)"
                           r"\s+(?P<h>NA|[\d.]+\*?)\s*$")
_PLASMAX_SECTIONS = ["Proteinogenic Amino Acids", "Non-proteinogenic Amino Acids", "Amino Acids Derivatives",
                     "Other Components", "Inorganic Salts", "Trace Elements", "Vitamins"]


def _strip_footnote(name: str) -> str:
    """'L-Alanine1' -> 'L-Alanine'; 'Vitamin B126' -> 'Vitamin B12'; '(betaine) 2' -> '(betaine)'."""
    # Symbol-font glyphs in the PDF text: U+F061 U+F02D = "alpha-" (A9a: keep the alpha in the name)
    name = name.replace("\uf061\uf02d", "alpha-")
    name = re.sub(r"[\uf000-\uf0ff]", "", name).strip()
    for cand in (name, re.sub(r"\s*\d$", "", name)):
        if cand in RAW_COMPONENTS:
            return cand
    raise KeyError(f"Plasmax component not in RAW_COMPONENTS: {name!r}")


def plasmax_text(path: Path) -> str:
    import pypdf  # pinned in requirements.txt; only needed for the Plasmax supplement

    reader = pypdf.PdfReader(str(path))
    return "\n".join(reader.pages[i].extract_text() for i in (1, 2))  # Table S1 spans PDF pages 2-3


def parse_plasmax(path: Path) -> pd.DataFrame:
    """Plasmax column of Vande Voorde 2019 Table S1 (PDF pages 2-3), uM. Rows with Plasmax 'NA' are absent."""
    text = plasmax_text(path)
    start = text.index("Plasmax  HPLM")
    end = text.index("* present as dibasic salt")
    rows, section = [], None
    for k, line in enumerate(text[start:end].splitlines()[1:]):
        for s in _PLASMAX_SECTIONS:
            if s in line:
                section = s
        m = _PLASMAX_LINE.match(line)
        if not m:
            continue
        name = _strip_footnote(m.group("name"))
        if m.group("p") == "NA":
            continue
        rows.append(dict(medium="Plasmax", section=section, raw_component=name, mw=math.nan, mg_per_L=math.nan,
                         uM=float(m.group("p")), source=f"{path.name}:TableS1:line{k + 1}"))
    return pd.DataFrame(rows)


def parse_cantor_classic(path: Path, column: str = "RPMI 1640") -> dict[str, float]:
    """One column of Cantor 2017 Table S1 'Classic synthetic media' (uM), keyed by metabolite name."""
    df = pd.read_excel(path, sheet_name=CLASSIC_SHEET, header=None)
    assert "All concentrations in uM" in str(df.iloc[1, 0])
    col = next(j for j in range(df.shape[1]) if str(df.iloc[3, j]).strip() == column)
    assert str(df.iloc[4, col]).strip() == "11875" or column != "RPMI 1640"
    out = {}
    for i in range(6, len(df)):
        name, v = df.iloc[i, 0], df.iloc[i, col]
        if pd.isna(name) or str(name).strip() == "Small ion":
            if out:
                break
            continue
        out[str(name).strip()] = float(str(v).rstrip("*"))
    return out


def classic_rpmi(atcc_rows: pd.DataFrame, cantor: dict[str, float], cantor_file: str) -> pd.DataFrame:
    """Classic RPMI 1640 (A9d) derived from the ATCC sheet rows (see RPMI_CLASSIC_*)."""
    rows = atcc_rows[~atcc_rows.raw_component.isin(RPMI_CLASSIC_DROP)].copy()
    rows["medium"] = "RPMI"
    rows["correction"] = ""
    for raw, cname in RPMI_CLASSIC_FROM_CANTOR.items():
        i = rows.index[rows.raw_component == raw]
        assert len(i) == 1, raw
        um = cantor[cname]
        rows.loc[i, "uM"] = um
        rows.loc[i, "mg_per_L"] = um * rows.loc[i, "mw"] / 1000
        rows.loc[i, "correction"] = f"classic RPMI 1640 (Gibco 11875) value from {cantor_file}:{CLASSIC_SHEET}"
    rows["source"] = rows.source + " (derived: classic RPMI 1640, A9d)"
    return rows


def read_compositions(root: Path | None = None) -> pd.DataFrame:
    """All formulation rows: four vendor sheets (RPMI as the on-disk ATCC sheet `RPMI-ATCC` and the derived classic
    RPMI 1640 `RPMI`), HPLM, Plasmax."""
    root = root or media_dir()
    parts = [parse_vendor_sheet(root / f, RPMI_ATCC if base == "RPMI" else base) for base, (f, _s) in SHEET_FILES.items()]
    parts.insert(1, classic_rpmi(parts[0], parse_cantor_classic(root / HPLM_FILE), HPLM_FILE))
    parts += [parse_hplm(root / HPLM_FILE), parse_plasmax(root / PLASMAX_FILE)]
    out = pd.concat(parts, ignore_index=True)
    out["correction"] = out.correction.fillna("")
    return out


def canonical_composition(rows: pd.DataFrame) -> pd.DataFrame:
    """Per-medium canonical metabolite concentrations (uM): raw components x stoichiometry, summed."""
    out = []
    for r in rows.itertuples():
        mapped = _canon_rows(r.raw_component)
        if isinstance(mapped, str):
            out.append(dict(medium=r.medium, canonical=f"<{mapped}>", uM=r.uM, raw_component=r.raw_component))
            continue
        for canon, n in mapped:
            out.append(dict(medium=r.medium, canonical=canon, uM=r.uM * n, raw_component=r.raw_component))
    df = pd.DataFrame(out)
    return (df.groupby(["medium", "canonical"], as_index=False)
              .agg(uM=("uM", "sum"), raw_components=("raw_component", lambda s: "; ".join(sorted(set(s))))))


def mapping_status(canonical: str) -> str:
    if canonical.startswith("<"):
        return "non_nutrient"
    if canonical in EXCHANGES:
        return "mapped"
    if canonical in NO_EXCHANGE:
        return "no_exchange"
    raise KeyError(f"canonical metabolite with no status: {canonical!r}")


def verify_exchange_map(model) -> pd.DataFrame:
    """Check every EXCHANGES id in a parent: exists, is an exchange, metabolite name matches. One row each."""
    ex_ids = {r.id for r in model.exchanges}
    rows = []
    for canon, (rid, name) in EXCHANGES.items():
        ok_rxn = rid in ex_ids
        met_name = list(model.reactions.get_by_id(rid).metabolites)[0].name if rid in model.reactions else None
        rows.append(dict(canonical=canon, exchange=rid, is_exchange=ok_rxn, met_name=met_name,
                         name_ok=bool(met_name) and met_name.lower() == name))
    return pd.DataFrame(rows)


# ===================================================================================== 2. formulation strings
@dataclass
class Supplement:
    token: str
    kind: str  # nutrient | serum | non_nutrient | cysteine_rule | unrepresentable
    components: dict[str, float | None] = field(default_factory=dict)  # canonical -> uM (None = unstated amount)
    note: str = ""


@dataclass
class ParsedText:
    text: str
    base: str  # text before the first ' + ' (design_numbers.base_of)
    base_class: str  # sheet:<X> | mixture_of_sheets | organoid | stem | other_defined
    parts: list[tuple[str, float]]  # [(sheet base, fraction)] for sheet bases and two-sheet mixtures
    serum_pct: float
    serum_type: str  # FBS | horse | none
    supplements: list[Supplement]
    flags: list[str]

    @property
    def serum_free(self) -> bool:
        return self.serum_pct == 0


def base_of(text: str) -> str:
    return text.split(" + ")[0].strip()


def base_class(text: str) -> str:
    """Identical to scripts/design_numbers.py::base_class (asserted in tests)."""
    b = base_of(text)
    if re.search(r"organoid|OPAC|OCMI", text, re.I):
        return "organoid"
    if re.search(r"Tumor Stem", text, re.I):
        return "stem"
    if b in SHEET_BASES:
        return "sheet:" + ("MEM/EMEM" if b in ("MEM", "EMEM") else b)
    parts = re.sub(r"\s*\(.*\)", "", b).split(":")
    if len(parts) == 2 and all(p in SHEET_BASES for p in parts):
        return "mixture_of_sheets"
    return "other_defined"


B27 = ("biotin", "alpha-tocopherol", "retinol", "corticosterone", "galactose", "ethanolamine", "GSH", "carnitine",
       "linoleate", "linolenate", "progesterone", "putrescine", "selenite", "triiodothyronine")
N2 = ("progesterone", "putrescine", "selenite")
NEAA = ("glycine", "alanine", "asparagine", "aspartate", "glutamate", "proline", "serine")  # MEM NEAA, 100 uM each
MEM_EAA = ("arginine", "cystine", "histidine", "isoleucine", "leucine", "lysine", "methionine", "phenylalanine",
           "threonine", "tryptophan", "tyrosine", "valine")
_NUM = r"(\d+(?:\.\d+)?)"
_UNIT_UM = {"mm": 1000.0, "um": 1.0, "nm": 1e-3, "pm": 1e-6}


def _conc_um(token: str, mw: float | None = None) -> float | None:
    """Concentration in uM from '2mM', '100uM', '500uM', '15 mg/L', '1g/l' (mg/L needs mw); None if unstated."""
    t = token.lower().replace(" ", "")
    m = re.match(_NUM + r"(mm|um|nm|pm)", t)
    if m:
        return float(m.group(1)) * _UNIT_UM[m.group(2)]
    m = re.match(_NUM + r"(mg|g|ug|ng)/(l|ml)", t)
    if m and mw:
        mass = float(m.group(1)) * {"g": 1e3, "mg": 1.0, "ug": 1e-3, "ng": 1e-6}[m.group(2)]  # mg
        per_l = mass * (1000.0 if m.group(3) == "ml" else 1.0)  # mg/L
        return per_l / mw * 1000.0
    return None


# token regex -> (kind, canonical list, mw for mass units, note)
_SUPPLEMENTS: list[tuple[str, str, tuple[str, ...], float | None, str]] = [
    (r"glutamax", "nutrient", ("glutamine",), None, "GlutaMAX (alanyl-glutamine) -> glutamine (DESIGN parser rule)"),
    (r"(l-)?glutamine", "nutrient", ("glutamine",), 146.0, ""),
    (r"sodium pyruvate|pyruvate", "nutrient", ("pyruvate", "Na+"), 110.0, ""),
    (r"\bmem neaa\b|\bneaa\b", "nutrient", NEAA, None, "MEM NEAA: Gly, Ala, Asn, Asp, Glu, Pro, Ser (100 uM each at 1x)"),
    (r"\beaa\b", "nutrient", MEM_EAA, None, "MEM essential amino acids"),
    (r"insulin-transferrin-(sodium )?selenite?|insulin-transferrin-selenium|\bits\b", "nutrient", ("selenite",), None,
     "ITS: insulin and transferrin are proteins (not mapped); selenite has no exchange"),
    (r"hepes", "non_nutrient", (), None, "buffer"),
    (r"sodium bicarbonate", "nutrient", ("HCO3-", "Na+"), 84.0, ""),
    (r"sodium hydrocortisone succinate|hydrocortisone", "nutrient", ("cortisol",), 362.5, "hydrocortisone = cortisol"),
    (r"beta-estradiol", "nutrient", ("estradiol",), 272.4, ""),
    (r"triiodo-l-thyronine", "nutrient", ("triiodothyronine",), 651.0, ""),
    (r"testosterone", "nutrient", ("testosterone",), 288.4, ""),
    (r"phosphorylethanolamine", "nutrient", ("phosphoethanolamine",), 141.1, ""),
    (r"ethanolamine", "nutrient", ("ethanolamine",), 61.1, ""),
    (r"sodium selenite|selenium", "unrepresentable", ("selenite",), None, "selenite has no exchange"),
    (r"mercaptoethanol", "cysteine_rule", (), None, "2-ME: cysteine uptake allowed (DESIGN parser rule)"),
    (r"\bb-?27\b", "nutrient", B27, None, "B27 (Brewer 1993 composition; amounts unstated): retinyl acetate -> "
     "retinol, DL-alpha-tocopherol -> alpha-tocopherol; proteins not mapped; selenite has no exchange"),
    (r"\bn-?2\b", "nutrient", N2, None, "N2: insulin and transferrin not mapped; selenite has no exchange"),
    (r"linoleic acid", "nutrient", ("linoleate",), 280.4, ""),
    (r"naco3|sodium carbonate", "nutrient", ("HCO3-", "Na+"), 106.0, "sodium carbonate -> HCO3- (bicarbonate buffer)"),
    (r"glutathione", "nutrient", ("GSH",), 307.3, ""),
    (r"hypoxanthine", "nutrient", ("hypoxanthine",), 136.1, ""),
    (r"\badenine\b", "nutrient", ("adenine",), 135.1, ""),
    (r"\batp\b", "nutrient", ("adenosine", "Pi"), 507.2,
     "ATP -> adenosine + phosphate (A9h: extracellular ATP is hydrolysed by ecto-nucleotidases, not taken up)"),
    (r"thymidine", "nutrient", ("thymidine",), 242.2, ""),
    (r"proline", "nutrient", ("proline",), 115.1, ""),
    (r"cystine", "nutrient", ("cystine",), 240.3, ""),
    (r"glucose", "nutrient", ("glucose",), 180.2, ""),
    (r"cyanocobalamin|vitamin b12", "nutrient", ("cobalamin",), 1355.4, ""),
    (r"5-bromo-2'-deoxyuridine", "unrepresentable", ("BrdU",), None, "BrdU (thymidine analogue) has no exchange"),
    (r"insulin|transferrin|catalase|bsa|egf|epidermal growth factor|fgf|gm-csf|g-csf|il-6|il-2|\blif\b|heparin|"
     r"thyroid stimulating hormone|cholera toxin|human recombinant", "non_nutrient", (), None, "protein/peptide"),
    (r"pen-strep|gentam[iy]cin|\bg4[01]8\b", "non_nutrient", (), None, "antibiotic / selection agent"),
    (r"singlequot|keratinocyte growth supplement|\bbpe\b", "non_nutrient", (), None, "proprietary kit / extract"),
    (r"^dmem:f12", "non_nutrient", (), None, "base-medium component of an unmappable (stem) formulation"),
]


def _parse_supplement(token: str) -> Supplement:
    low = token.lower().strip().replace("µ", "u").replace("μ", "u")
    m = re.match(r"^" + _NUM + r"\s*%\s*(fbs|fcs|horse serum|nu-serum)$", low)
    if m:
        return Supplement(token, "serum", {}, m.group(2))
    for pat, kind, canon, mw, note in _SUPPLEMENTS:
        if re.search(pat, low):
            conc = _conc_um(low, mw)
            mult = 1.0
            mx = re.match(r"^" + _NUM + r"\s*x\b", low)
            if canon in (NEAA, MEM_EAA):
                if mx:
                    mult = float(mx.group(1))
                conc = conc if conc is not None else (100.0 * mult if canon is NEAA else None)
            comps = {c: conc for c in canon}
            return Supplement(token, kind, comps, note)
    raise ValueError(f"unparsed supplement token {token!r}")


_MIX = re.compile(r"^(?P<a>[^:()]+):(?P<b>[^:()]+?)\s*(?:\((?P<ra>\d+):(?P<rb>\d+)\))?$")


def parse_formulation(text: str) -> ParsedText:
    """DepMap medium text -> base, mixture fractions, serum %, supplements. Raises on any unparsed token."""
    text = text.strip()
    base = base_of(text)
    cls = base_class(text)
    rest = text[len(text.split(" + ")[0]):]
    tokens = [t.strip() for t in re.split(r"\+", rest) if t.strip()]
    flags: list[str] = []
    parts: list[tuple[str, float]] = []
    if cls.startswith("sheet:"):
        parts = [(base, 1.0)]
    elif cls == "mixture_of_sheets":
        m = _MIX.match(base)
        a, b = m.group("a").strip(), m.group("b").strip()
        ra, rb = (float(m.group("ra")), float(m.group("rb"))) if m.group("ra") else (1.0, 1.0)
        if not m.group("ra"):
            flags.append("mixture_ratio_unstated_assumed_1:1")
        parts = [(a, ra / (ra + rb)), (b, rb / (ra + rb))]
    sups = [_parse_supplement(t) for t in tokens]
    serum = [s for s in sups if s.kind == "serum"]
    serum_pct = sum(float(re.match(_NUM, s.token.strip()).group(1)) for s in serum)
    serum_type = "none" if not serum else ("FBS" if all(s.note in ("fbs", "fcs") for s in serum) else
                                           "+".join(sorted({s.note for s in serum})))
    if serum and serum_type != "FBS":
        flags.append("non_FBS_serum_treated_as_serum")
    if base == "DMEM" and not any("pyruvate" in s.components for s in sups):
        flags.append("dmem_pyruvate_unstated_default_none")
    if base == "DMEM" and not any("glucose" in s.components for s in sups):
        flags.append("dmem_glucose_unstated_default_high")
    return ParsedText(text, base, cls, parts, serum_pct, serum_type, [s for s in sups if s.kind != "serum"], flags)


# ===================================================================================== 3. media as compositions
@dataclass
class Medium:
    """A medium for bound construction: canonical metabolite -> uM (None = present, amount unstated)."""

    key: str
    components: dict[str, float | None]
    serum: bool
    allow_cysteine: bool = False
    flags: tuple[str, ...] = ()
    serum_pct: float = 10.0  # parsed serum % (used only by the A9(f) sensitivity f x pct/10)

    def present(self) -> set[str]:
        return {c for c, v in self.components.items() if v is None or v > 0}


def sheet_media(comp: pd.DataFrame) -> dict[str, dict[str, float]]:
    """canonical_composition() -> {medium: {canonical: uM}} (non-nutrients dropped)."""
    out: dict[str, dict[str, float]] = {}
    for r in comp.itertuples():
        if r.canonical.startswith("<"):
            continue
        out.setdefault(r.medium, {})[r.canonical] = float(r.uM)
    return out


def build_medium(parsed: ParsedText, sheets: dict[str, dict[str, float]], add_standard_glutamine: bool = True,
                 dmem_pyruvate: bool = False, key: str | None = None, cysteine_rule: str | None = None,
                 rpmi_base: str = "RPMI") -> Medium:
    """Documented medium for a parsed text whose sheet(s) are on disk (linear blend for mixtures)."""
    comps: dict[str, float | None] = {}
    flags = list(parsed.flags)
    sheet_cys = False
    for base, frac in parsed.parts:
        if base not in sheets:
            raise KeyError(f"sheet for {base!r} not on disk")
        base_comps = dict(sheets[rpmi_base if base == "RPMI" else base])
        if base == "RPMI" and rpmi_base != "RPMI":
            flags.append(f"rpmi_base:{rpmi_base}")
        if add_standard_glutamine and base in STANDARD_GLN_UM and "glutamine" not in base_comps:
            base_comps["glutamine"] = STANDARD_GLN_UM[base]
            flags.append(f"gln_standard_default:{base}")
        if dmem_pyruvate and base == "DMEM":
            base_comps["pyruvate"] = 1000.0
            flags.append("dmem_with_pyruvate_sensitivity")
        sheet_cys |= base_comps.get("cysteine", 0) > 0
        for c, v in base_comps.items():
            comps[c] = (comps.get(c) or 0.0) + v * frac
    allow_cys = False
    for s in parsed.supplements:
        if s.kind == "cysteine_rule":
            allow_cys = True
        for c, v in s.components.items():
            if v is None:
                comps.setdefault(c, None)
            else:
                comps[c] = (comps.get(c) or 0.0) + v
    rule = cysteine_rule or CYSTEINE_RULE
    assert rule in CYSTEINE_RULES, rule
    if rule == "documented":
        allow_cys |= sheet_cys
        flags.append("cysteine_documented_sensitivity")
    return Medium(key or parsed.text, comps, serum=not parsed.serum_free, allow_cysteine=allow_cys,
                  flags=tuple(sorted(set(flags))), serum_pct=parsed.serum_pct)


def physiological_medium(name: str, sheets: dict[str, dict[str, float]], cysteine_rule: str | None = None) -> Medium:
    """M-HPLM / M-Plasmax: composition as published, treated as serum-containing (DESIGN: + the M2 serum model).

    The cysteine rule applies as in every medium (A9b): the listed free cysteine is closed by default.
    """
    comps = dict(sheets[name])
    rule = cysteine_rule or CYSTEINE_RULE
    documented = rule == "documented"
    return Medium(f"M-{name}" + ("|cys-documented" if documented else ""), comps, serum=True,
                  allow_cysteine=documented and comps.get("cysteine", 0) > 0,
                  flags=("physiological_counterfactual",) + (("cysteine_documented_sensitivity",) if documented else ()))


# ===================================================================================== 4. bounds
@dataclass
class ParentConstraints:
    """Per-parent floor and requirement table (from scripts/w2_medium_floor.py)."""

    floor: tuple[str, ...]  # canonical metabolites
    req: dict[str, float]  # canonical -> req_c (floor: FVA-min; serum: pFBA uptake)
    floor_mult: float = FLOOR_MULT

    def floor_bound(self, c: str) -> float:
        return min(self.floor_mult * self.req[c], UPTAKE_OPEN)  # never above the availability bound


VARIANTS = ("M-open", "M0", "M1", "M2", "M3")


def uptake_bounds(model, variant: str, medium: Medium | None = None, pc: ParentConstraints | None = None,
                  f: float | None = None, gate_open: bool = True, calibration: dict[str, float] | None = None,
                  floor_override: float | None = None, serum_pct_scaling: bool = False) -> dict[str, float]:
    """cobra medium dict {exchange id: max uptake} for one ladder variant (see module docstring).

    ``floor_override`` sets every floor component that the medium does not supply to that bound (used by the
    non-binding precheck: 1000 vs FLOOR_MULT x req_c).
    ``serum_pct_scaling`` (A9f sensitivity, run only at the TG-chosen f): the serum components use
    f_eff = f x serum_pct / 10 instead of f (10% serum = the reference).
    """
    from nextgem.model import hams_medium

    if variant == "M-open":
        return {r.id: UPTAKE_OPEN for r in model.exchanges}
    if variant == "M0":
        return hams_medium(model)
    if variant not in ("M1", "M2", "M3"):
        raise ValueError(f"unknown variant {variant!r}")
    assert medium is not None and pc is not None
    present = medium.present()
    b: dict[str, float] = {}
    for c in sorted(present):
        if c in EXCHANGES and (c != "cysteine" or medium.allow_cysteine):
            b[EXCHANGES[c][0]] = UPTAKE_OPEN
    if medium.allow_cysteine and "cystine" in present:
        # DESIGN M1 / A9(b): with 2-mercaptoethanol, extracellular cystine is reduced to cysteine -> cysteine uptake
        b[EXCHANGES["cysteine"][0]] = UPTAKE_OPEN
    for c in BACKGROUND:
        b[EXCHANGES[c][0]] = UPTAKE_OPEN
    for c in pc.floor:
        if c not in present:
            b[EXCHANGES[c][0]] = floor_override if floor_override is not None else pc.floor_bound(c)
    if variant in ("M2", "M3") and medium.serum:
        if f is None:
            raise ValueError("M2/M3 need f")
        f_eff = f * medium.serum_pct / 10.0 if serum_pct_scaling else f
        for c in SERUM_F_SCALED:
            if c not in present:
                b[EXCHANGES[c][0]] = f_eff * pc.req[c]
        for c in SERUM_TRACE_METALS:
            b[EXCHANGES[c][0]] = UPTAKE_OPEN
    if variant == "M3":
        if not calibration:
            raise ValueError("M3 needs calibration rates {glucose, glutamine, O2}")
        for c, v in calibration.items():
            b[EXCHANGES[c][0]] = v
    if not gate_open:
        for c in GATED:
            b.pop(EXCHANGES[c][0], None)
    return {k: v for k, v in b.items() if v > 0}


def bounds_signature(bounds: dict[str, float]) -> str:
    """Stable hash of a bounds dict (distinct sweeps are run per signature)."""
    import hashlib

    s = ";".join(f"{k}={v:.12g}" for k, v in sorted(bounds.items()))
    return hashlib.sha256(s.encode()).hexdigest()[:16]


# ===================================================================================== 5. line assignment
def screened_lines() -> list[str]:
    dm = paths.depmap_dir()
    return sorted(pd.read_csv(dm / "CRISPRScreenMap.csv").ModelID.unique())


def assign_line_media(on_disk_sheets: set[str]) -> pd.DataFrame:
    """Per screened line: tier, formulation text, base class, stratum, attribution flag, medium used.

    Tier logic is scripts/design_numbers.py's (QC-passing combined screens -> ModelCondition text in the column
    named FormulationID, modal text; else Model.csv text in FormulationID; else imputed). Content assertions check
    the documented column swap.
    """
    dm = paths.depmap_dir()
    screened = screened_lines()
    model = pd.read_csv(dm / "Model.csv").set_index("ModelID")
    smap = pd.read_csv(dm / "CRISPRScreenMap.csv")
    ssm = pd.read_csv(dm / "ScreenSequenceMap.csv")
    ok = ssm[(ssm.PassesQC == True) & (ssm.ExcludeFromCRISPRCombined == False) & (ssm.ScreenType != "pDNA")]  # noqa: E712
    ok = ok[ok.ScreenID.isin(smap.ScreenID)]
    mc = pd.read_csv(dm / "ModelCondition.csv", usecols=["ModelConditionID", "FormulationID", "GrowthMedia",
                                                        "SerumFreeMedia"])
    assert mc.FormulationID.dropna().str.contains(r"FBS|RPMI|DMEM", regex=True).mean() > 0.5
    assert mc.GrowthMedia.dropna().str.startswith("MF-").mean() > 0.9
    assert model.FormulationID.dropna().str.contains(r"FBS|RPMI|DMEM", regex=True).mean() > 0.5
    assert model.OnboardedMedia.dropna().str.startswith("MF-").mean() > 0.9
    scr = ok.merge(mc, on="ModelConditionID", how="left").dropna(subset=["FormulationID"])
    sets = scr.groupby("ModelID").FormulationID.agg(lambda s: list(s.value_counts().index))
    mf_of = scr.groupby("ModelID").GrowthMedia.agg(lambda s: ";".join(sorted(set(s.dropna()))))
    rows = []
    for mid in screened:
        mf = ""
        if mid in sets.index:
            texts = sets[mid]
            if len({base_of(t) for t in texts}) > 1:
                tier, text = "screen_mixed", None
            else:
                tier, text = "screen", texts[0]
                mf = mf_of.get(mid, "")
            n_texts = len(texts)
        elif pd.notna(model.loc[mid, "FormulationID"]):
            tier, text, n_texts = "onboarded", model.loc[mid, "FormulationID"], 1
            mf = model.loc[mid, "OnboardedMedia"] if pd.notna(model.loc[mid, "OnboardedMedia"]) else ""
        else:
            tier, text, n_texts = "imputed", None, 0
        bcls = base_class(text) if text else tier
        onboarded = model.loc[mid, "FormulationID"] if pd.notna(model.loc[mid, "FormulationID"]) else None
        rows.append(dict(ModelID=mid, lineage=model.loc[mid, "OncotreeLineage"], tier=tier, text=text,
                         n_screen_texts=n_texts, mf_id=mf, base_class=bcls, onboarded_text=onboarded))
    df = pd.DataFrame(rows)
    on_disk_cls = {"sheet:" + b for b in on_disk_sheets}
    df["attribution_set"] = (df.tier == "screen") & df.base_class.isin(on_disk_cls)

    def stratum(r) -> str:
        if r.tier in ("imputed",):
            return "imputed"
        if r.tier == "screen_mixed":
            return "unmappable"
        if r.base_class in ("organoid", "stem", "other_defined"):
            return "unmappable"
        if r.base_class.startswith("sheet:"):
            if r.base_class in on_disk_cls:
                return "attribution" if r.tier == "screen" else "onboarded_sensitivity"
            return "fallback_sheet_missing"
        if r.base_class == "mixture_of_sheets":
            p = parse_formulation(r.text)
            return "mixture_on_disk" if all(b in on_disk_sheets for b, _ in p.parts) else "fallback_sheet_missing"
        raise AssertionError(r)

    df["stratum"] = df.apply(stratum, axis=1)
    documented = df.stratum.isin(["attribution", "onboarded_sensitivity", "mixture_on_disk"])
    df["medium_text_used"] = df.text.where(documented, IMPUTED_TEXT)
    df["medium_source"] = df.stratum.map({"attribution": "documented", "onboarded_sensitivity": "documented",
                                          "mixture_on_disk": "documented_mixture", "imputed": "imputed_flagged",
                                          "fallback_sheet_missing": "imputed_flagged_fallback",
                                          "unmappable": "imputed_flagged_unmappable"})
    parsed = {t: parse_formulation(t) for t in df.text.dropna().unique()}
    df["serum_free"] = df.text.map(lambda t: parsed[t].serum_free if isinstance(t, str) else pd.NA)
    df["serum_pct"] = df.text.map(lambda t: parsed[t].serum_pct if isinstance(t, str) else pd.NA)
    return df


def mf_crosscheck(df: pd.DataFrame) -> pd.DataFrame:
    """MF-ID family (MF-XXX) x base medium of the text: the modal base per family and the disagreeing lines."""
    d = df[df.mf_id.astype(bool) & df.text.notna()].copy()
    d["mf_family"] = d.mf_id.str.split(";").str[0].str.extract(r"^(MF-\d{3})")[0]
    d["base"] = d.text.map(base_of)
    tab = d.groupby(["mf_family", "base"]).size().rename("n").reset_index()
    modal = tab.sort_values("n", ascending=False).drop_duplicates("mf_family").set_index("mf_family").base
    tab["family_modal_base"] = tab.mf_family.map(modal)
    tab["agrees"] = tab.base == tab.family_modal_base
    return tab.sort_values(["mf_family", "n"], ascending=[True, False])


def gate_states(lines: list[str]) -> pd.DataFrame:
    """+g gate per line from default expression rows: open if any SLC1A1/2/3 log2(TPM+1) >= 1; NA if no expression."""
    f = paths.depmap_dir() / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    header = pd.read_csv(f, nrows=0).columns
    cols = {g: next(c for c in header if c.startswith(f"{g} (") and c.endswith(f"({GATE_ENTREZ[g]})"))
            for g in GATE_GENES}
    idcols = [c for c in ("ModelID", "IsDefaultEntryForModel") if c in header]
    e = pd.read_csv(f, usecols=idcols + list(cols.values()))
    e = e[e.IsDefaultEntryForModel == "Yes"].set_index("ModelID")
    assert e.index.is_unique
    e = e.rename(columns={v: k for k, v in cols.items()})[list(GATE_GENES)]
    out = pd.DataFrame(index=pd.Index(lines, name="ModelID"))
    out = out.join(e)
    has = out[list(GATE_GENES)].notna().all(axis=1)
    out["gate_open"] = (out[list(GATE_GENES)] >= GATE_THRESHOLD).any(axis=1).where(has, pd.NA)
    return out.reset_index()


def shuffled_assignment(df: pd.DataFrame, seeds, key_col: str = "medium_key") -> pd.DataFrame:
    """Negative control: permute medium keys among attribution-set lines within lineage x serum-free status."""
    import numpy as np

    a = df[df.attribution_set].copy()
    out = a[["ModelID", "lineage", "serum_free", key_col]].copy()
    for s in seeds:
        rng = np.random.default_rng(s)
        sh = a[key_col].copy()
        for _, idx in a.groupby(["lineage", "serum_free"]).groups.items():
            idx = list(idx)
            sh.loc[idx] = rng.permutation(a.loc[idx, key_col].to_numpy())
        out[f"shuffled_seed{s}"] = sh.values
    return out


# ===================================================================================== 6. library helpers
def medium_key(parsed: ParsedText, med: Medium, sheets: dict[str, dict[str, float]]) -> str:
    """Readable key from bound-relevant content: base (or mixture) + mapped additions vs the base sheet(s)."""
    base_set: set[str] = set()
    for b, _ in parsed.parts:
        base_set |= set(sheets[b])
        if b in STANDARD_GLN_UM:
            base_set.add("glutamine")
    adds = sorted(c for c in med.present() - base_set if c in EXCHANGES)
    label = parsed.base if parsed.base_class.startswith("sheet:") else \
        ":".join(f"{b}{round(f, 3)}" for b, f in parsed.parts)
    key = label + "".join(f"+{a}" for a in adds)
    if med.allow_cysteine and "cysteine" not in base_set:
        key += "+cys(2-ME)"
    return key + ("" if med.serum else "|serum-free")


def load_sheets() -> dict[str, dict[str, float]]:
    return sheet_media(canonical_composition(read_compositions()))


def documented_media(texts, sheets) -> dict[str, Medium]:
    """medium_key -> Medium for every text whose sheet(s) are on disk (first text per key is representative)."""
    out: dict[str, Medium] = {}
    for t in sorted(texts):
        p = parse_formulation(t)
        if p.parts and all(b in sheets for b, _ in p.parts):
            med = build_medium(p, sheets)
            med.key = medium_key(p, med, sheets)
            out.setdefault(med.key, med)
    return out


def sensitivity_media(sheets) -> dict[str, Medium]:
    """Pre-declared sensitivity variants (A9; flagged; not part of the ladder). M1 and M2(f) are swept for each."""
    p_rpmi, p_dmem, p_df12 = (parse_formulation(f"{b} + 10% FBS") for b in ("RPMI", "DMEM", "DMEM:F12"))
    out = {
        # A9(d): the on-disk ATCC-modified RPMI (with 1 mM pyruvate) instead of classic RPMI 1640
        "SENS:RPMI-ATCC": build_medium(p_rpmi, sheets, key="SENS:RPMI-ATCC", rpmi_base=RPMI_ATCC),
        # DESIGN W2 parser rule: DMEM with unstated pyruvate -> none; with-pyruvate sensitivity
        "SENS:DMEM+pyruvate": build_medium(p_dmem, sheets, dmem_pyruvate=True, key="SENS:DMEM+pyruvate"),
        # A9(c): no standard glutamine added to plain DMEM / DMEM:F12
        "SENS:DMEM-no-standard-Gln": build_medium(p_dmem, sheets, add_standard_glutamine=False,
                                                  key="SENS:DMEM-no-standard-Gln"),
        "SENS:DMEM:F12-no-standard-Gln": build_medium(p_df12, sheets, add_standard_glutamine=False,
                                                      key="SENS:DMEM:F12-no-standard-Gln"),
        # A9(b): documented free cysteine open
        "SENS:DMEM:F12-cysteine-documented": build_medium(p_df12, sheets, key="SENS:DMEM:F12-cysteine-documented",
                                                          cysteine_rule="documented"),
    }
    for m in out.values():
        m.flags = m.flags + ("sensitivity",)
    return out


def physiological_sensitivity_media(sheets) -> dict[str, Medium]:
    """A9(b) sensitivity for the physiological media: HPLM / Plasmax with their listed free cysteine open."""
    return {m.key: m for m in (physiological_medium(n, sheets, cysteine_rule="documented") for n in ("HPLM", "Plasmax"))}


def floor_constraints(floor_tsv: Path, req_tsv: Path, mult_tsv: Path | None = None) -> dict[str, ParentConstraints]:
    """Per-parent ParentConstraints; ``mult_tsv`` (parent, floor_mult) is written by the non-binding precheck."""
    fl = pd.read_csv(floor_tsv, sep="\t")
    rq = pd.read_csv(req_tsv, sep="\t")
    mult = {}
    if mult_tsv is not None and Path(mult_tsv).exists():
        mult = pd.read_csv(mult_tsv, sep="\t").set_index("parent").floor_mult.to_dict()
    out = {}
    for p, g in rq.groupby("parent"):
        out[p] = ParentConstraints(floor=tuple(sorted(fl[fl.parent == p].component)),
                                   req=dict(zip(g.component, g.req_c)), floor_mult=float(mult.get(p, FLOOR_MULT)))
    return out
