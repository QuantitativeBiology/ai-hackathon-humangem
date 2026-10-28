"""W6 — information in flux states beyond the inputs (DESIGN.md §9 W6, §6.5, §13; amendments A2, A9, A11, A12).

Label-free library: nothing here opens a CRISPR label file or a label-derived list (DESIGN §11.2). The response
data (PRISM, GDSC2, CTRP drug screens; CCLE metabolomics; OncoTree / inferred molecular subtypes) are W6 targets
only (§11.2 "Drug screens are W6 targets only").

Contents
--------
1. The a priori compound table (``notes/apriori_lists_final.md`` §4 under A2, used verbatim): a parser of the §4 table
   and a transcription of the targets / activating genes per listed compound, cross-checked against the parsed tokens.
2. Screen-name matching (PRISM 20Q2 secondary AUC → 24Q2 primary LFC fallback, GDSC2 AUC, CTRP AUC).
3. Assay media (DESIGN W6 "Assay medium"): PRISM from its SOP, GDSC from "Screen Medium", metabolomics (and, DESIGN
   being silent, CTRP and subtypes) from the line's culture medium in the frozen W2 library (A9).
4. Flux-state simulation per *model state* (arm × medium): WT pFBA, the single-gene KO e-vector, and the KO pFBA
   flux of the curated drug-target genes (compute SF6). A generic model has one state per medium, so its features
   are line-invariant except through the medium. That is a property of the features, not a bug.
5. Comparator blocks: X_RNA, X_prot, and GPR-mapped reaction expression X_rxnRNA / X_rxnProt (E-Flux capacity:
   AND = min, OR = sum, as ``nextgem.context.capacity``).
6. Matched-capacity models: per-block PCA (20; 50 as sensitivity) computed in *state space* from the block Gram
   matrix (exact; handles repeated states), ridge regression with lineage one-hot covariate, LOLO folds; within-
   lineage stratified 5-fold kNN for subtypes (SEEDS[0:3]).
7. Controls (permuted within lineage; synthetic nonlinear injection; flux as a linear function of expression) and the
   pre-registered decision statistics (one-sided Wilcoxon of per-unit Δ, Holm across the 3 analyses, above the
   permuted control).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from nextgem.config import SEEDS

PROVISIONAL_TAG = "provisional (pre-P1)"
APRIORI_FILE = Path("notes") / "apriori_lists_final.md"  # A2: replaces notes/apriori_lists_approved.md
APRIORI_SECTION = "# 4. Final metabolic-drug list"

MIN_LINES = 100  # DESIGN W6: a compound needs >= 100 lines with features
PCA_K = 20
PCA_K_SENS = 50
EIG_REL_TOL = 1e-9  # PCA components with eigenvalue <= tol * max are numerically zero (constant blocks)
RIDGE_ALPHAS = tuple(float(a) for a in np.logspace(-3, 5, 17))
KNN_K = 5
SUBTYPE_MIN_CLASS = 5  # DESIGN W6: lineages with >= 2 subtypes of >= 5 lines
SUBTYPE_FOLDS = 5
CV_SEEDS = tuple(SEEDS[0:3])  # §7: ML and CV shuffles
PERM_SEEDS = tuple(SEEDS[0:5])  # §7: permuted controls
FIT_FREE_MIN_LINES = 10  # lines per lineage for a within-lineage Spearman rho
WITHIN_LINEAGE_MIN = 10
ALPHA_ONE_SIDED = 0.025  # §6.3 convention (one-sided 0.025)
LETHAL_REL = 1e-6  # KO growth <= LETHAL_REL * WT growth: lethal KO
ZERO_FLUX_TOL = 1e-9
N_SYNTH = 20  # synthetic positive-control units per run
NEG_DIM = 200  # columns of the synthetic "flux = linear function of expression" block


# ===================================================================================== 1. compound table
def _split_outside_parens(s: str, seps=(",", ";", " and ")) -> list[str]:
    out, cur, depth, i = [], "", 0, 0
    while i < len(s):
        c = s[i]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        if depth == 0:
            hit = next((sep for sep in seps if s.startswith(sep, i)), None)
            if hit is not None:
                out.append(cur)
                cur = ""
                i += len(hit)
                continue
        cur += c
        i += 1
    out.append(cur)
    return out


def _cut_sentence(s: str) -> str:
    """Drop trailing prose: the first '. ' outside parentheses ends the compound list of a cell."""
    depth = 0
    for i, c in enumerate(s):
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif c == "." and depth == 0 and s[i + 1:i + 2] == " ":
            return s[:i]
    return s


def _strip_groups(s: str) -> str:
    s = re.sub(r"\[[^\]]*\]", "", s)  # citations [LIT ...]
    prev = None
    while prev != s:  # nested-safe removal of parenthetical groups
        prev, s = s, re.sub(r"\([^()]*\)", "", s)
    return s


def parse_apriori_compounds(text: str) -> pd.DataFrame:
    """Compound tokens of the §4 table, verbatim, in file order: columns class, token, row, pos."""
    if APRIORI_SECTION not in text:
        raise ValueError("§4 section not found in the a priori list file")
    sec = text.split(APRIORI_SECTION, 1)[1].split("\n# 5.", 1)[0]
    rows = []
    in_table = False
    for line in sec.splitlines():
        if line.startswith("| Class"):
            in_table = True
            continue
        if not in_table:
            continue
        if not line.startswith("|"):
            if rows:
                break
            continue
        if set(line.replace("|", "").strip()) <= {"-"}:
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 2:
            raise ValueError(f"unexpected §4 row: {line[:80]}")
        cls, cell = cells
        cell = _cut_sentence(re.sub(r"\[[^\]]*\]", "", cell))
        for pos, tok in enumerate(_split_outside_parens(cell)):
            t = re.sub(r"\s+", " ", _strip_groups(tok)).strip().rstrip(".").strip()
            if t:
                rows.append({"class": cls, "token": t, "row": len({r["class"] for r in rows}), "pos": pos})
    df = pd.DataFrame(rows)
    if df.empty or df.token.duplicated().any():
        raise ValueError("compound parse failed or duplicated tokens")
    return df


def positive_control_text(text: str) -> str:
    sec = text.split(APRIORI_SECTION, 1)[1].split("\n# 5.", 1)[0]
    m = re.search(r"\*\*Positive-control pairs\*\*:(.*)", sec)
    if not m:
        raise ValueError("positive-control paragraph not found")
    return m.group(1)


def _c(targets=(), src="table", act=(), mod=(), flags=(), aliases=()):
    return {"targets": list(targets), "target_source": src if targets else "none", "activators": list(act),
            "modifiers": list(mod), "flags": list(flags), "aliases": list(aliases)}


# Transcription of §4 (key = parsed token, verbatim). target_source: "table" = stated in the token's own
# parenthesis or a parenthesis shared with the compounds listed before it in the same group (the table's
# convention, e.g. "simvastatin, lovastatin (...), ..., rosuvastatin (HMGCR)"); "class" = the row's class names one
# enzyme family; "literature" = not stated in the table (assigned here, flagged in the output). Aliases are screen
# names (normalised by ``norm_name`` before matching); the token itself and its '/'-separated parts are always tried.
COMPOUNDS: dict[str, dict] = {
    # Antifolates
    "methotrexate": _c(["DHFR"], mod=["SLC19A1", "FPGS"]),
    "pemetrexed": _c(["TYMS", "DHFR", "GART"]),
    "pralatrexate": _c(["DHFR"], "class"),
    "raltitrexed": _c(["TYMS"], "literature"),
    "pyrimethamine": _c(["DHFR"], "class"),
    "trimethoprim": _c(["DHFR"], "class", flags=["bacterial-selective DHFR inhibitor"]),
    # TYMS
    "5-FU": _c(["TYMS"], act=["UMPS", "UPP1", "UCK2", "TYMP", "TK1"], mod=["DPYD"], flags=["prodrug"],
               aliases=["fluorouracil", "5-fluorouracil"]),
    "floxuridine": _c(["TYMS"], act=["TK1"], flags=["prodrug"]),
    "capecitabine": _c(["TYMS"], act=["CES1", "CES2", "CDA", "TYMP"], flags=["prodrug"]),
    "trifluridine+tipiracil": _c(["TYMS"], mod=["TYMP"], act=["TK1"], flags=["screened as trifluridine"],
                                 aliases=["trifluridine", "trifluridine/tipiracil", "TAS-102"]),
    # RNR
    "hydroxyurea": _c(["RRM2"], aliases=["hydroxycarbamide"]),
    "gemcitabine": _c(["RRM1"], act=["DCK"], flags=["prodrug"]),
    "triapine": _c(["RRM2"]),
    # DHODH
    "brequinar": _c(["DHODH"], "class"),
    "teriflunomide": _c(["DHODH"], "class"),
    "leflunomide": _c(["DHODH"], flags=["non-enzymatic conversion to teriflunomide"]),
    "BAY-2402234": _c(["DHODH"], "class"),
    # IMPDH
    "mycophenolic acid / mofetil": _c(["IMPDH1", "IMPDH2"], "class", act=["CES1"],
                                      aliases=["mycophenolic acid", "mycophenolate mofetil"]),
    "ribavirin": _c(["IMPDH1", "IMPDH2"], "class", act=["ADK"], flags=["prodrug"]),
    "mizoribine": _c(["IMPDH1", "IMPDH2"], "class", act=["ADK"], flags=["prodrug"]),
    # Purine analogues
    "6-mercaptopurine": _c(["PPAT"], "literature", act=["HPRT1"], flags=["prodrug", "'de novo purine' target -> PPAT"],
                           aliases=["mercaptopurine", "6-MP"]),
    "6-thioguanine": _c(["PPAT"], "literature", act=["HPRT1"], flags=["prodrug", "'de novo purine' target -> PPAT"],
                        aliases=["thioguanine", "tioguanine"]),
    "fludarabine": _c(act=["DCK"], flags=["prodrug"]),
    "cladribine": _c(act=["DCK"], flags=["prodrug"]),
    "clofarabine": _c(act=["DCK"], flags=["prodrug"]),
    "azaserine": _c(["PPAT", "PFAS"], "literature", flags=["glutamine amidotransferases"]),
    "allopurinol": _c(["XDH"]),
    # Pyrimidine analogues
    "cytarabine": _c(act=["DCK"], flags=["prodrug"], aliases=["ara-c"]),
    "decitabine": _c(flags=["no model target (DNMT1)"]),
    "azacitidine": _c(flags=["no model target (DNMT1)"], aliases=["5-azacytidine"]),
    # Glutaminase / amidotransferase
    "CB-839/telaglenastat": _c(["GLS"]),
    "BPTES": _c(["GLS"]),
    "acivicin": _c(["CTPS1"]),
    # Statins and FDPS
    "simvastatin": _c(["HMGCR"], flags=["lactone prodrug"]),
    "lovastatin": _c(["HMGCR"], flags=["lactone prodrug"]),
    "atorvastatin": _c(["HMGCR"]),
    "fluvastatin": _c(["HMGCR"]),
    "pitavastatin": _c(["HMGCR"]),
    "rosuvastatin": _c(["HMGCR"]),
    "zoledronate": _c(["FDPS"], aliases=["zoledronic acid"]),
    "alendronate": _c(["FDPS"], aliases=["alendronic acid", "alendronate sodium"]),
    "ibandronate": _c(["FDPS"], aliases=["ibandronic acid"]),
    "risedronate": _c(["FDPS"], aliases=["risedronic acid"]),
    # SQLE
    "terbinafine": _c(["SQLE"]),
    # FASN / ACLY / SCD / ACC
    "TVB-3166": _c(["FASN"]),
    "TVB-2640/denifanstat": _c(["FASN"]),
    "GSK2194069": _c(["FASN"]),
    "orlistat": _c(["FASN"]),
    "cerulenin": _c(["FASN"]),
    "BMS-303141": _c(["ACLY"]),
    "SB-204990": _c(["ACLY"]),
    "NDI-091143": _c(["ACLY"]),
    "bempedoic acid": _c(["ACLY"], act=["SLC27A2"], flags=["prodrug", "liver-restricted activation"]),
    "A939572": _c(["SCD"]),
    "CAY10566": _c(["SCD"]),
    "MF-438": _c(["SCD"]),
    "ND-646": _c(["ACACA", "ACACB"]),
    "firsocostat": _c(["ACACA", "ACACB"]),
    "BAY ACCi": _c(["ACACA", "ACACB"], aliases=["BAY ACCi"]),
    # CPT1 / FAO
    "etomoxir": _c(["CPT1A"], flags=["off-target at high dose"]),
    "perhexiline": _c(["CPT1A"], "class"),
    # Glycolysis and transporters
    "2-deoxyglucose": _c(["HK1", "HK2"], aliases=["2-deoxy-D-glucose"]),
    "lonidamine": _c(["HK2"], "literature"),
    "3-bromopyruvate": _c(["HK2", "GAPDH"], "literature"),
    "TT-232/TLN-232": _c(["PKM"]),
    "shikonin": _c(["PKM"]),
    "AZD3965": _c(["SLC16A1"]),
    "BAY-876": _c(["SLC2A1"]),
    "GSK2837808A": _c(["LDHA", "LDHB"]),
    "GNE-140": _c(["LDHA", "LDHB"]),
    # OXPHOS ("complex I" resolved to MT-ND1: every complex-I subunit KO closes the same reaction set)
    "IACS-010759": _c(["MT-ND1"]),
    "rotenone": _c(["MT-ND1"]),
    "phenformin": _c(["MT-ND1"]),
    "metformin": _c(["MT-ND1"], flags=["weak, mechanism debated"]),
    "antimycin A": _c(["MT-CYB"]),
    "atovaquone": _c(["MT-CYB"]),
    "oligomycin-A": _c(["ATP5F1A"], aliases=["oligomycin A"]),
    "CPI-613/devimistat": _c(["PDHA1", "DLD", "DLST"]),
    # Redox / ferroptosis
    "buthionine sulfoximine": _c(["GCLC", "GCLM"], aliases=["buthionine sulphoximine", "BSO",
                                                            "L-buthionine-(S,R)-sulfoximine"]),
    "sulfasalazine": _c(["SLC7A11"]),
    "erastin": _c(["SLC7A11"], flags=["also hits VDAC"]),
    "disulfiram": _c(["ALDH2"]),
    "elesclomol": _c(),
    "N-acetyl cysteine": _c(aliases=["acetylcysteine", "N-acetylcysteine"]),
    "α-lipoic acid": _c(aliases=["lipoic acid", "alpha-lipoic acid", "thioctic acid"]),
    # NAD
    "FK866/daporinad": _c(["NAMPT"], mod=["NAPRT"], aliases=["APO866"]),
    "GNE-617": _c(["NAMPT"], mod=["NAPRT"]),
    "GMX1778/CHS-828": _c(["NAMPT"], mod=["NAPRT"]),
    "STF-118804": _c(["NAMPT"], mod=["NAPRT"]),
    # 1C / SHMT
    "SHIN1": _c(["SHMT1", "SHMT2"]),
    "SHIN2": _c(["SHMT1", "SHMT2"]),
    # Polyamines
    "eflornithine/DFMO": _c(["ODC1"]),
    # Mutant IDH (negative-expectation class)
    "ivosidenib": _c(["IDH1"], flags=["negative-expectation class"], aliases=["AG-120"]),
    "AGI-5198": _c(["IDH1"], flags=["negative-expectation class"]),
    "enasidenib": _c(["IDH2"], flags=["negative-expectation class"], aliases=["AG-221"]),
    "AGI-6780": _c(["IDH2"], flags=["negative-expectation class"]),
    "vorasidenib": _c(["IDH1", "IDH2"], "literature", flags=["negative-expectation class"], aliases=["AG-881"]),
    # Off-model controls (excluded from the primary test; reported separately)
    "tigecycline": _c(flags=["off-model control (mitoribosome)"]),
    "doxycycline": _c(flags=["off-model control (mitoribosome)"]),
    "chloramphenicol": _c(flags=["off-model control (mitoribosome)"]),
    "linezolid": _c(flags=["off-model control (mitoribosome)"]),
    "bortezomib": _c(flags=["off-model control (no model target)"]),
    "belinostat": _c(flags=["off-model control (no model target)"]),
    "alisertib": _c(flags=["off-model control (no model target)"]),
    "mitoxantrone": _c(flags=["off-model control (no model target)"]),
}
OFF_MODEL_CLASS = "Off-model controls"

# Positive-control pairs of §4, transcribed (genes, compound tokens, note).
POSITIVE_CONTROL_PAIRS: list[tuple[tuple[str, ...], tuple[str, ...], str]] = [
    (("DHODH",), ("brequinar",), ""),
    (("DHFR",), ("methotrexate",), ""),
    (("TYMS",), ("5-FU", "floxuridine", "pemetrexed"), ""),
    (("IMPDH2",), ("mycophenolic acid / mofetil",), ""),
    (("RRM1", "RRM2"), ("hydroxyurea", "gemcitabine"), ""),
    (("HMGCR",), ("simvastatin", "lovastatin", "atorvastatin", "fluvastatin", "pitavastatin", "rosuvastatin"), ""),
    (("SQLE",), ("terbinafine",), ""),
    (("FASN",), ("GSK2194069", "TVB-2640/denifanstat"), ""),
    (("ACACA",), ("ND-646",), ""),
    (("GLS",), ("CB-839/telaglenastat", "BPTES"), ""),
    (("SLC16A1",), ("AZD3965",), ""),
    (("SLC2A1",), ("BAY-876",), ""),
    (("CPT1A",), ("etomoxir",), ""),
    (("NAMPT",), ("FK866/daporinad",), "NAPRT modifier"),
    (("GCLC", "GCLM"), ("buthionine sulfoximine",), ""),
    (("SLC7A11",), ("sulfasalazine", "erastin"), ""),
    (("ATP5F1A",), ("oligomycin-A",), ""),
    (("PDHA1", "DLD", "DLST"), ("CPI-613/devimistat",), ""),
    (("HK1", "HK2"), ("2-deoxyglucose",), ""),
    (("ODC1",), ("eflornithine/DFMO",), ""),
    (("SHMT1", "SHMT2"), ("SHIN1", "SHIN2"), ""),
    (("IDH1",), ("ivosidenib", "AGI-5198"), "restricted to IDH_mut lines"),
]


def compound_table(text: str) -> pd.DataFrame:
    """Parsed §4 tokens joined to the transcription; raises if the two disagree (verbatim check)."""
    parsed = parse_apriori_compounds(text)
    missing = set(parsed.token) - set(COMPOUNDS)
    extra = set(COMPOUNDS) - set(parsed.token)
    if missing or extra:
        raise ValueError(f"compound transcription differs from §4: missing {sorted(missing)}, extra {sorted(extra)}")
    rows = []
    for _, r in parsed.iterrows():
        c = COMPOUNDS[r.token]
        rows.append({"compound": r.token, "class": r["class"], "targets": ";".join(c["targets"]),
                     "target_source": c["target_source"], "activators": ";".join(c["activators"]),
                     "modifiers": ";".join(c["modifiers"]), "flags": ";".join(c["flags"]),
                     "aliases": ";".join(c["aliases"]), "off_model_control": r["class"] == OFF_MODEL_CLASS})
    return pd.DataFrame(rows)


def positive_control_pairs_long() -> pd.DataFrame:
    rows = [{"gene": g, "compound": c, "note": note} for gs, cs, note in POSITIVE_CONTROL_PAIRS for g in gs for c in cs]
    return pd.DataFrame(rows)


def target_genes(table: pd.DataFrame) -> list[str]:
    """Sorted curated drug-target symbols (the X_KOflux gene set); off-model controls have none."""
    s = set()
    for t in table.targets:
        s |= {x for x in str(t).split(";") if x}
    return sorted(s)


# ===================================================================================== 2. screen names
def norm_name(s) -> str:
    s = str(s).lower().replace("α", "alpha").replace("β", "beta")
    return re.sub(r"[^a-z0-9]", "", s)


def candidate_names(token: str) -> list[str]:
    """Names tried for a compound, in priority order: aliases, the token, its '/'-separated parts."""
    c = COMPOUNDS[token]
    names = list(c["aliases"]) + [token] + [p.strip() for p in re.split(r"/", token) if p.strip()]
    out, seen = [], set()
    for n in names:
        k = norm_name(n)
        if k and k not in seen:
            seen.add(k)
            out.append(k)
    return out


def match_compounds(tokens, index: dict[str, set], n_lines: dict) -> pd.DataFrame:
    """For each token, the screen id with the most lines among the ids matched by its first matching name.

    ``index``: normalised screen name -> set of screen ids; ``n_lines``: screen id -> lines with data.
    """
    rows = []
    for t in tokens:
        hit_name, ids = None, set()
        for k in candidate_names(t):
            if k in index:
                hit_name, ids = k, index[k]
                break
        if not ids:
            rows.append({"compound": t, "screen_id": None, "matched_name": None, "n_ids": 0, "n_lines_screen": 0})
            continue
        best = sorted(ids, key=lambda i: (-n_lines.get(i, 0), str(i)))[0]
        rows.append({"compound": t, "screen_id": best, "matched_name": hit_name, "n_ids": len(ids),
                     "n_lines_screen": int(n_lines.get(best, 0))})
    return pd.DataFrame(rows)


def name_index(pairs) -> dict[str, set]:
    """pairs: iterable of (screen id, name string possibly listing synonyms separated by , ; |)."""
    idx: dict[str, set] = {}
    for sid, names in pairs:
        if names is None or (isinstance(names, float) and math.isnan(names)):
            continue
        for part in re.split(r"[,;|]", str(names)):
            k = norm_name(part)
            if k:
                idx.setdefault(k, set()).add(sid)
    return idx


# ===================================================================================== 3. assay media
PRISM_MEDIUM_KEY = "RPMI"
PRISM_MEDIUM_SOURCE = ("data/prism/PRISM_SOP_Supplemental_Materials.pdf, 'PRISM Assay' materials and 'Reagent setup': "
                       "RPMI 1640, no phenol red (Gibco 11835055) + 10% heat-inactivated FBS -> W2 key 'RPMI' (classic "
                       "RPMI 1640 + 10% FBS; phenol red is a non-nutrient)")
GDSC_MEDIUM_KEYS = {"R": "RPMI", "D/F12": "DMEM:F12"}
ASSAY_MAPS = ("prism", "gdsc", "culture")
DATASET_MAP = {"PRISM": "prism", "GDSC2": "gdsc", "CTRP": "culture", "metabolomics": "culture", "subtypes": "culture"}


def assay_medium_map(kind: str, lines, culture: pd.Series, gdsc_screen_medium: pd.Series | None = None) -> pd.Series:
    """line -> W2 medium key for one assay map; lines without a medium are absent (NaN dropped)."""
    lines = list(lines)
    if kind == "prism":
        return pd.Series(PRISM_MEDIUM_KEY, index=pd.Index(lines, name="ModelID"))
    if kind == "culture":
        return culture.reindex(lines).dropna()
    if kind == "gdsc":
        if gdsc_screen_medium is None:
            raise ValueError("gdsc map needs the GDSC Screen Medium per line")
        return gdsc_screen_medium.reindex(lines).map(GDSC_MEDIUM_KEYS).dropna()
    raise ValueError(kind)


# ===================================================================================== 4. simulation
class SolverInfraError(RuntimeError):
    """Solver / licence / environment error: aborts the job with no output (DESIGN §11.4)."""


@dataclass(frozen=True)
class Arm:
    """A feature-producing model configuration. ``variant`` is a W3 Variant JSON (or None)."""

    name: str
    parent: str = "human1"
    variant: dict | None = None
    bounds_variant: str = "M2"
    f: float | None = 1.0
    gate_open: bool = True

    def variant_obj(self):
        from nextgem import w3

        return w3.Variant.from_json(self.variant) if self.variant else w3.Variant("T0")


def state_id(arm: str, medium_key: str) -> str:
    return f"{arm}|{medium_key}"


def reaction_var_names(model) -> list[tuple[str, str, str]]:
    return [(r.id, r.forward_variable.name, r.reverse_variable.name) for r in model.reactions]


def _flux_vector(model, names) -> np.ndarray:
    pv = model.solver.primal_values
    return np.array([pv[f] - pv[r] for _, f, r in names], dtype=np.float64)


def pfba(model, names, biomass: str) -> tuple[float, np.ndarray | None, str]:
    """FBA then minimal-L1 flux at the FBA optimum (biomass >= mu (1 - 1e-9)). The current medium/variant/KO apply.

    Every LP gets the 300 s limit with one retry (``nextgem.w3.solve_status``). Returns (mu, flux or None, status);
    infeasible -> (0, None, 'infeasible').
    """
    from nextgem import w3

    model.objective = biomass
    st = w3.solve_status(model)
    if st != "optimal":
        return 0.0, None, st
    mu = float(model.solver.objective.value)
    with model:
        bio = model.reactions.get_by_id(biomass)
        bio.lower_bound = max(0.0, mu * (1 - 1e-9))
        from optlang.symbolics import Zero

        obj = model.problem.Objective(Zero, direction="min", sloppy=True)
        model.objective = obj
        coefs = {}
        for r in model.reactions:
            coefs[r.forward_variable] = 1.0
            coefs[r.reverse_variable] = 1.0
        model.objective.set_linear_coefficients(coefs)
        st2 = w3.solve_status(model)
        if st2 != "optimal":
            raise SolverInfraError(f"pFBA L1 step not optimal ({st2}) after a feasible FBA")
        v = _flux_vector(model, names)
    return mu, v, "optimal"


def ko_groups_for(symbols: list[str], sym2gene: dict, closures: dict) -> dict[frozenset, list[str]]:
    """Target symbols grouped by identical KO closure sets (one KO pFBA per group)."""
    out: dict[frozenset, list[str]] = {}
    for s in symbols:
        g = sym2gene.get(s)
        if g is None or g not in closures:
            continue
        out.setdefault(frozenset(closures[g]), []).append(s)
    return out


def simulate_state(model, names, genes, closures, bounds: dict, variant, target_syms: list[str], sym2gene: dict,
                   biomass: str, do_ko: bool = True, do_koflux: bool = True) -> dict:
    """One model state: WT pFBA flux, KO ratio_r for every gene, KO pFBA flux per target group.

    KO pFBA shortcut (exact): if no reaction closed by the KO carries flux in the WT pFBA solution, that solution stays
    feasible at the same growth and remains L1-minimal over the smaller feasible set, so it is the KO pFBA solution.
    """
    import time

    from nextgem import w3

    t0, c0 = time.time(), time.process_time()
    out: dict = {"n_lp_koflux": 0, "n_shortcut_koflux": 0}
    with model:
        model.medium = bounds
        variant.apply(model)
        mu, v, st = pfba(model, names, biomass)
        out.update(wt_growth=mu, wt_status=st, flux=v)
        koflux, ko_growth = {}, {}
        if do_koflux and v is not None and mu > 0:
            carrying = {rid for (rid, _, _), x in zip(names, v) if abs(x) > ZERO_FLUX_TOL}
            for clos, syms in ko_groups_for(target_syms, sym2gene, closures).items():
                if not clos or not (set(clos) & carrying):
                    vk, mk = v, mu
                    out["n_shortcut_koflux"] += 1
                else:
                    with model:
                        model.genes.get_by_id(sym2gene[syms[0]]).knock_out()
                        mk, vk, stk = pfba(model, names, biomass)
                    out["n_lp_koflux"] += 1
                    if vk is None:
                        vk, mk = np.zeros_like(v), 0.0
                for s in syms:
                    koflux[s], ko_growth[s] = vk, mk
        out["koflux"], out["ko_growth"] = koflux, ko_growth
    if do_ko:
        r = w3.ko_sweep(model, genes, closures, bounds, variant)
        out.update(ratio_r=r["ratio_r"], ko_wt_growth=r["wt_growth"], ko_n_solved=r["n_solved"],
                   shortcut_verify_max_dev=r["shortcut_verify_max_dev"])
    out["wall_s"], out["cpu_s"] = time.time() - t0, time.process_time() - c0
    return out


def flux_features(v: np.ndarray, mu: float) -> np.ndarray:
    """X_flux transform: asinh(v / mu) (biomass-normalised; asinh is a fixed, fold-free variance stabiliser)."""
    if v is None or not mu > 0:
        return np.full(0 if v is None else len(v), np.nan)
    return np.arcsinh(v / mu)


def koflux_features(vk: np.ndarray, mu_wt: float) -> np.ndarray:
    """X_KOflux transform: KO pFBA flux normalised by the *WT* biomass flux (finite for lethal KOs), asinh."""
    return np.arcsinh(vk / mu_wt)


# ===================================================================================== 5. comparator blocks
def _compile_vec(node, cols: dict, X: np.ndarray):
    """E-Flux capacity on vectors (lines): AND = min over known children, OR = sum (unknown if any child unknown)."""
    if isinstance(node, str):
        j = cols.get(node)
        return None if j is None else X[:, j]
    op, ch = node
    vals = [_compile_vec(c, cols, X) for c in ch]
    if op == "and":
        known = [v for v in vals if v is not None]
        return np.minimum.reduce(known) if known else None
    if any(v is None for v in vals):
        return None
    return np.add.reduce(vals)


def reaction_capacity(gprs: dict, lin: pd.DataFrame) -> pd.DataFrame:
    """Lines x reactions E-Flux capacity from linear gene values (columns = model gene ids); reactions whose GPR
    is unknown (an OR child without data, or no gene with data) are dropped."""
    cols = {g: j for j, g in enumerate(lin.columns)}
    X = lin.to_numpy(dtype=np.float64)
    out = {}
    for rid, node in gprs.items():
        v = _compile_vec(node, cols, X)
        if v is not None and np.all(np.isfinite(v)):
            out[rid] = v
    return pd.DataFrame(out, index=lin.index)


def log_capacity(cap: pd.DataFrame) -> pd.DataFrame:
    return np.log2(cap.clip(lower=0) + 1.0)


def impute_line_floor(X: pd.DataFrame, q: float = 0.01) -> pd.DataFrame:
    """Missing = not detected: each line's NaNs get that line's own q-quantile (no cross-line statistic)."""
    fl = X.quantile(q, axis=1)
    return X.T.fillna(fl).T


# ===================================================================================== 6. state-space PCA
@dataclass
class Block:
    """Features of one block as state rows; lines map to states through ``state_of`` (one map per assay).

    ``G`` = S S^T over all states (float64). Omics blocks have one state per line.
    """

    name: str
    states: list
    G: np.ndarray
    p: int
    constant: bool = False

    @staticmethod
    def from_matrix(name: str, S: pd.DataFrame) -> "Block":
        A = S.to_numpy(dtype=np.float64)
        if not np.all(np.isfinite(A)):
            raise ValueError(f"{name}: non-finite features")
        A = A[:, A.std(axis=0) > 0] if len(A) > 1 else A
        G = A @ A.T
        return Block(name, list(S.index), G, A.shape[1], constant=A.shape[1] == 0)

    def index_of(self, keys) -> np.ndarray:
        pos = {s: i for i, s in enumerate(self.states)}
        return np.array([pos[k] for k in keys], dtype=np.int64)


def pca_fit_transform(G: np.ndarray, train_idx: np.ndarray, all_idx: np.ndarray, k: int) -> np.ndarray:
    """PCA scores (up to k components) of the rows ``all_idx`` (state indices, repeats allowed) fitted on the
    training rows ``train_idx`` (repeats = weights), from the Gram matrix only. Exactly equals sklearn PCA on the
    explicit row matrix up to component signs. Components with eigenvalue <= EIG_REL_TOL * max are dropped."""
    u, w = np.unique(train_idx, return_counts=True)
    n = float(w.sum())
    if len(u) < 2:
        return np.zeros((len(all_idx), 0))
    Guu = G[np.ix_(u, u)]
    a = Guu @ w / n  # a_j = <S_j, mu>
    b = float(w @ Guu @ w) / n ** 2  # <mu, mu>
    C = Guu - a[:, None] - a[None, :] + b  # <S_i - mu, S_j - mu>
    sw = np.sqrt(w)
    K = (sw[:, None] * C * sw[None, :]) / n  # covariance-normalised (eigenvalues = PCA variances)
    lam, V = np.linalg.eigh(K)
    order = np.argsort(lam)[::-1]
    lam, V = lam[order], V[:, order]
    if lam[0] <= 0:
        return np.zeros((len(all_idx), 0))
    keep = lam > EIG_REL_TOL * lam[0]
    keep[k:] = False
    lam, V = lam[keep], V[:, keep]
    # sign convention: largest-|loading| state weight positive (deterministic)
    sgn = np.sign(V[np.abs(V).argmax(axis=0), np.arange(V.shape[1])])
    V = V * sgn
    t_rows = np.unique(all_idx)
    Gtu = G[np.ix_(t_rows, u)]
    at = Gtu @ w / n
    Ct = Gtu - at[:, None] - a[None, :] + b  # <S_t - mu, S_j - mu>
    Z = (Ct * sw[None, :]) @ V / np.sqrt(lam * n)
    pos = {s: i for i, s in enumerate(t_rows)}
    return Z[[pos[s] for s in all_idx]]


# ===================================================================================== 7. models and metrics
def spearman(x, y) -> float:
    from scipy.stats import spearmanr

    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return np.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def within_lineage_rho(x: pd.Series, y: pd.Series, lineage: pd.Series, min_lines: int = WITHIN_LINEAGE_MIN) -> tuple[float, int]:
    """n-weighted mean (Fisher z) of per-lineage Spearman rho over lineages with >= min_lines and non-constant x, y."""
    df = pd.DataFrame({"x": x, "y": y, "lin": lineage}).dropna()
    zs, ns = [], []
    for _, d in df.groupby("lin"):
        if len(d) < min_lines:
            continue
        r = spearman(d.x, d.y)
        if np.isfinite(r):
            zs.append(np.arctanh(np.clip(r, -0.999999, 0.999999)))
            ns.append(len(d))
    if not zs:
        return np.nan, 0
    return float(np.tanh(np.average(zs, weights=ns))), len(zs)


@dataclass
class View:
    """A block seen through a line -> state map (optionally permuted within lineage)."""

    block: Block
    state_of: pd.Series  # line -> state key (for omics: the line itself)
    tag: str  # cache key: block name | map | perm

    def idx(self, lines) -> np.ndarray:
        return self.block.index_of(self.state_of.loc[list(lines)].to_numpy())


def permute_within_lineage(state_of: pd.Series, lineage: pd.Series, seed: int) -> pd.Series:
    """Shuffle the line -> state assignment among lines of the same lineage (multiset per lineage preserved)."""
    rng = np.random.default_rng(seed)
    out = state_of.copy()
    lin = lineage.reindex(state_of.index)
    for _, ix in state_of.groupby(lin.fillna("__none__"), sort=True).groups.items():
        ix = sorted(ix)
        out.loc[ix] = state_of.loc[ix].to_numpy()[rng.permutation(len(ix))]
    return out


class PCACache:
    """PCA scores per (view tag, fold key, k), fitted on the view's universe minus the held-out lines."""

    def __init__(self):
        self._c: dict = {}
        self._p: dict = {}

    def scores(self, view: View, universe: list, heldout: set, fold_key, k: int) -> pd.DataFrame:
        """Fitted once at max(k, PCA_K_SENS) components; the leading k columns are the k-component PCA."""
        kmax = max(k, PCA_K_SENS)
        key = (view.tag, fold_key)
        if key not in self._c or self._c[key][1] < kmax:
            train = [l for l in universe if l not in heldout]
            Z = pca_fit_transform(view.block.G, view.idx(train), view.idx(universe), kmax)
            self._c[key] = (pd.DataFrame(Z, index=universe), kmax)
        Z = self._c[key][0]
        return Z.iloc[:, :k]

    def positions(self, tag, fold_key, Z) -> dict:
        key = ("pos", tag, fold_key)
        if key not in self._p:
            self._p[key] = {l: i for i, l in enumerate(Z.index)}
        return self._p[key]

    def heldout(self, lineage: pd.Series, L) -> set:
        key = ("held", L)
        if key not in self._p:
            self._p[key] = set(lineage.index[lineage == L])
        return self._p[key]


def _design(views: list[View], universes: dict, cache: PCACache, heldout: set, fold_key, k: int, lines) -> np.ndarray:
    parts = []
    for v in views:
        Z = cache.scores(v, universes[v.tag], heldout, fold_key, k)
        pos = cache.positions(v.tag, fold_key, Z)
        parts.append(Z.to_numpy()[[pos[l] for l in lines]])
    return np.hstack(parts) if parts else np.zeros((len(lines), 0))


def ridge_loo_fit_predict(Xtr: np.ndarray, ytr: np.ndarray, Xte: np.ndarray, alphas=RIDGE_ALPHAS) -> tuple[np.ndarray, float]:
    """Ridge with intercept; alpha chosen by exact leave-one-out squared error (as sklearn RidgeCV's default efficient
    LOO; first alpha on ties). Uses one SVD of the centred training matrix."""
    n = len(ytr)
    xm, ym = Xtr.mean(0), ytr.mean()
    Xc, yc = Xtr - xm, ytr - ym
    if Xc.shape[1] == 0:
        return np.full(len(Xte), ym), float("nan")
    U, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    Uty = U.T @ yc
    U2 = U ** 2
    best, best_a, best_w = np.inf, None, None
    for a in alphas:
        f = sv ** 2 / (sv ** 2 + a)
        fitted = U @ (f * Uty)
        h = U2 @ f + 1.0 / n
        loo = np.mean(((yc - fitted) / (1.0 - h)) ** 2)
        if loo < best - 1e-15 * abs(best if np.isfinite(best) else 1.0):
            best, best_a = loo, a
            best_w = Vt.T @ ((sv / (sv ** 2 + a)) * Uty)
    return (Xte - xm) @ best_w + ym, float(best_a)


def lolo_ridge(y: pd.Series, views: list[View], universes: dict, lineage: pd.Series, cache: PCACache,
               k: int = PCA_K) -> pd.Series:
    """Held-out predictions under leave-one-lineage-out: per fold, block PCA fitted without the held-out lineage,
    scores standardised on training lines, ridge (alpha by exact LOO over RIDGE_ALPHAS) with the lineage one-hot
    as covariate."""
    y = y.dropna()
    lin = lineage.reindex(y.index)
    lin_arr = lin.to_numpy()
    yv = y.to_numpy(float)
    pred = np.full(len(y), np.nan)
    idx = np.asarray(y.index)
    for L in sorted(lin.dropna().unique()):
        te_m = lin_arr == L
        tr_m = ~te_m & pd.notna(lin_arr)
        if tr_m.sum() < 10:
            continue
        held = cache.heldout(lineage, L)
        Xtr = _design(views, universes, cache, held, L, k, idx[tr_m])
        Xte = _design(views, universes, cache, held, L, k, idx[te_m])
        if Xtr.shape[1]:
            mu, sd = Xtr.mean(0), Xtr.std(0)
            ok = sd > 1e-12
            Xtr, Xte = (Xtr[:, ok] - mu[ok]) / sd[ok], (Xte[:, ok] - mu[ok]) / sd[ok]
        levels = np.array(sorted(set(lin_arr[tr_m])))
        D_tr = (lin_arr[tr_m][:, None] == levels[None, :]).astype(float)
        D_te = np.zeros((te_m.sum(), len(levels)))
        p, _ = ridge_loo_fit_predict(np.hstack([Xtr, D_tr]), yv[tr_m], np.hstack([Xte, D_te]))
        pred[te_m] = p
    return pd.Series(pred, index=y.index)


def heldout_scores(y: pd.Series, pred: pd.Series, lineage: pd.Series) -> dict:
    ok = pred.dropna().index
    rho = spearman(pred.loc[ok], y.loc[ok])
    wl, nl = within_lineage_rho(pred.loc[ok], y.loc[ok], lineage)
    return {"rho": rho, "rho_within": wl, "n_lineages_within": nl, "n_pred": len(ok)}


def subtype_scores(labels: pd.Series, views: list[View], k: int = PCA_K, seeds=CV_SEEDS, cache: dict | None = None) -> dict:
    """Within-lineage subtype task: stratified 5-fold kNN macro-F1 (mean over seeds; PCA per training fold) and the
    silhouette of the labels on the standardised PCA scores of all task lines. ``cache`` (per task) reuses a view's
    fold PCA across feature sets that share it (identical folds for identical labels)."""
    from sklearn.metrics import f1_score, silhouette_score
    from sklearn.model_selection import StratifiedKFold
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.preprocessing import StandardScaler

    labels = labels.dropna()
    lines = list(labels.index)
    yv = labels.to_numpy()
    cache = {} if cache is None else cache
    lk = hash(tuple(lines))

    def pca(v, ltr, lall, key):
        ck = (v.tag, lk, key, k)
        if ck not in cache:
            cache[ck] = pca_fit_transform(v.block.G, v.idx(ltr), v.idx(lall), k)
        return cache[ck]

    f1s = []
    for seed in seeds:
        pred = np.empty(len(lines), dtype=object)
        skf = StratifiedKFold(n_splits=SUBTYPE_FOLDS, shuffle=True, random_state=seed)
        for fi, (tr, te) in enumerate(skf.split(np.zeros(len(lines)), yv)):
            ltr, lte = [lines[i] for i in tr], [lines[i] for i in te]
            Ztr, Zte = [], []
            for v in views:
                Z = pca(v, ltr, ltr + lte, (seed, fi))
                Ztr.append(Z[:len(ltr)])
                Zte.append(Z[len(ltr):])
            Xtr, Xte = np.hstack(Ztr), np.hstack(Zte)
            if Xtr.shape[1] == 0:
                Xtr, Xte = np.zeros((len(ltr), 1)), np.zeros((len(lte), 1))
            sc = StandardScaler().fit(Xtr)
            clf = KNeighborsClassifier(n_neighbors=KNN_K).fit(sc.transform(Xtr), yv[tr])
            pred[te] = clf.predict(sc.transform(Xte))
        f1s.append(f1_score(yv, pred.astype(yv.dtype), average="macro"))
    Zs = [pca(v, lines, lines, "all") for v in views]
    Z = np.hstack(Zs)
    if Z.shape[1] == 0:
        sil = np.nan
    else:
        Zs_ = (Z - Z.mean(0)) / np.where(Z.std(0) > 0, Z.std(0), 1.0)
        sil = float(silhouette_score(Zs_, yv)) if len(set(yv)) > 1 else np.nan
    return {"knn_macro_f1": float(np.mean(f1s)), "silhouette": sil}


# ===================================================================================== 8. controls and decisions
def wilcoxon_greater(d) -> float:
    """One-sided Wilcoxon signed-rank p for median(d) > 0 (zero differences dropped; all-zero or n < 1 -> p = 1)."""
    from scipy.stats import wilcoxon

    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    if (d != 0).sum() < 1:
        return 1.0
    return float(wilcoxon(d, alternative="greater", zero_method="wilcox").pvalue)


def holm(p: dict) -> dict:
    """Holm step-down adjusted p-values (monotone, capped at 1)."""
    keys = sorted(p, key=lambda k: (p[k], str(k)))
    m = len(keys)
    out, running = {}, 0.0
    for i, k in enumerate(keys):
        running = max(running, min(1.0, (m - i) * p[k]))
        out[k] = running
    return out


def decision_rows(deltas: dict[str, pd.Series], perm: dict[str, dict[int, pd.Series]], alpha: float = ALPHA_ONE_SIDED) -> pd.DataFrame:
    """One row per analysis: n units, median Δ, one-sided Wilcoxon p, Holm across the analyses, the permuted
    control's medians, 'above permuted' (real median > every permuted-seed median) and the verdict."""
    p = {a: wilcoxon_greater(d) for a, d in deltas.items()}
    ph = holm(p)
    rows = []
    for a, d in deltas.items():
        d = d.dropna()
        pm = {s: float(np.nanmedian(x)) if len(x.dropna()) else np.nan for s, x in perm.get(a, {}).items()}
        med = float(np.median(d)) if len(d) else np.nan
        above = bool(pm) and all(np.isfinite(v) and med > v for v in pm.values())
        rows.append({"analysis": a, "n_units": len(d), "median_delta": med, "mean_delta": float(d.mean()) if len(d) else np.nan,
                     "n_positive": int((d > 0).sum()), "n_zero": int((d == 0).sum()), "p_wilcoxon_one_sided": p[a],
                     "p_holm": ph[a], "perm_medians": ";".join(f"{s}:{v:.4g}" for s, v in pm.items()),
                     "above_permuted": above,
                     "verdict": "supported" if (ph[a] < alpha and med > 0 and above) else "not supported"})
    return pd.DataFrame(rows)


def synthetic_positive_response(z_state: pd.Series, state_of: pd.Series, seed: int, snr: float = 1.0) -> pd.Series:
    """Injected nonlinear signal: y = tanh(2 z) + noise with var(signal)/var(noise) = snr; z standardised over lines."""
    rng = np.random.default_rng(seed)
    z = z_state.reindex(state_of.to_numpy()).to_numpy(float)
    z = (z - z.mean()) / (z.std() if z.std() > 0 else 1.0)
    s = np.tanh(2 * z)
    sd = s.std() if s.std() > 0 else 1.0
    y = s + rng.normal(0, sd / math.sqrt(snr), len(s))
    return pd.Series(y, index=state_of.index)


def random_state_direction(block_matrix: pd.DataFrame, seed: int) -> pd.Series:
    """Projection of each state on a random unit direction of the block's (non-constant) feature space."""
    rng = np.random.default_rng(seed)
    A = block_matrix.to_numpy(float)
    A = A[:, A.std(0) > 0]
    if A.shape[1] == 0:
        return pd.Series(0.0, index=block_matrix.index)
    w = rng.normal(size=A.shape[1])
    return pd.Series((A - A.mean(0)) @ (w / np.linalg.norm(w)), index=block_matrix.index)


def linear_of_expression(rna: pd.DataFrame, seed: int, dim: int = NEG_DIM) -> pd.DataFrame:
    """Synthetic negative control: a 'flux' block that is an exact linear function of expression (random map)."""
    rng = np.random.default_rng(seed)
    A = rna.to_numpy(float)
    W = rng.normal(size=(A.shape[1], dim)) / math.sqrt(A.shape[1])
    return pd.DataFrame((A - A.mean(0)) @ W, index=rna.index, columns=[f"neg{j}" for j in range(dim)])


def fit_free_rho(x: pd.Series, y_sens: pd.Series, lineage: pd.Series) -> tuple[float, int]:
    """Within-lineage rho between a line-specific target score and sensitivity (fit-free); NaN when the score is
    constant within every lineage (e.g. generic-model e under one assay medium)."""
    return within_lineage_rho(x, y_sens, lineage, min_lines=FIT_FREE_MIN_LINES)


@dataclass
class RunConfig:
    tag: str
    label: str
    arms: list[Arm] = field(default_factory=list)
    note: str = ""

    @staticmethod
    def from_json(d: dict, project: Path) -> "RunConfig":
        import json

        arms = []
        for a in d["arms"]:
            if a.get("context_family"):
                raise NotImplementedError("line-specific (context-model) states need a per-line state builder; "
                                          "W6 currently simulates one state per (arm, medium)")
            var = None
            if a.get("variant"):
                spec = json.loads((project / a["variants_file"]).read_text())
                pool = [x for k in spec for x in (spec[k] if isinstance(spec[k], list) else [])]
                hits = [x for x in pool if x["name"] == a["variant"]]
                if len(hits) != 1:
                    raise ValueError(f"variant {a['variant']!r} not found once in {a['variants_file']}")
                var = hits[0]
            arms.append(Arm(a["name"], a.get("parent", "human1"), var, a.get("bounds_variant", "M2"), a.get("f", 1.0),
                            a.get("gate_open", True)))
        return RunConfig(d["tag"], d["label"], arms, d.get("note", ""))


# ===================================================================================== 9. worker side (spawn pool)
_W: dict = {}


def symbol_map(genes_tsv: Path, model_genes) -> dict[str, str]:
    """Gene symbol -> model gene id (symbols that map to one model gene only)."""
    t = pd.read_csv(genes_tsv, sep="\t", dtype=str)
    t = t[t.genes.isin(set(model_genes))]
    vc = t.geneSymbols.value_counts()
    return {s: g for s, g in zip(t.geneSymbols, t.genes) if isinstance(s, str) and vc[s] == 1}


def init_worker(parent: str, genes_tsv: str) -> None:
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem import w3
    from nextgem.model import load_model

    m = load_model(parent)
    _W.update(model=m, names=reaction_var_names(m), genes=sorted(g.id for g in m.genes), closures=w3.gene_closures(m),
              sym2gene=symbol_map(Path(genes_tsv), [g.id for g in m.genes]))


def worker_meta() -> dict:
    return {"rxn_ids": [n[0] for n in _W["names"]], "genes": _W["genes"]}


def state_task(args) -> dict:
    """(state_key, bounds, variant JSON or None, target symbols, do_ko, do_koflux) -> results with sparse fluxes."""
    from nextgem import w3
    from nextgem.knockouts import WT_GROWTH_TOL

    key, bounds, var_json, targets, do_ko, do_koflux = args
    variant = w3.Variant.from_json(var_json) if var_json else w3.Variant("T0")
    r = simulate_state(_W["model"], _W["names"], _W["genes"], _W["closures"], bounds, variant, targets, _W["sym2gene"],
                       w3.BIOMASS, do_ko=do_ko, do_koflux=do_koflux)
    if r["flux"] is None or not r["wt_growth"] > WT_GROWTH_TOL:
        raise RuntimeError(f"{key}: WT does not grow ({r['wt_status']}, {r['wt_growth']})")

    def sparse(v):
        nz = np.flatnonzero(np.abs(v) > ZERO_FLUX_TOL)
        return nz.astype(np.int32), v[nz]

    r["flux"] = sparse(r["flux"])
    r["koflux"] = {s: sparse(v) for s, v in r["koflux"].items()}
    r["key"] = key
    return r
