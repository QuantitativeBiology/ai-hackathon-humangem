"""Label-free per-line covariates for the split file (DESIGN.md §4.1 step 6, §8 W2 medium tiers).

Re-implements, for the screened lines, the label-free logic of ``scripts/design_numbers.py`` (cohort,
omics coverage, library class, medium tiers). ``tests/test_splits.py`` checks every count against
``results/design/design_numbers.json``. Opens no label-bearing file: the screened set comes from
``CRISPRScreenMap.csv``.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from nextgem import paths

SHEET_BASES = {"RPMI", "DMEM", "DMEM:F12", "MEM", "EMEM", "IMDM", "F12", "McCoy's 5A", "L-15"}


@dataclass
class Cohort:
    screened: list[str]
    lineage: pd.Series  # ModelID -> OncotreeLineage (screened)
    model_type: pd.Series
    expr: set[str]
    cn: set[str]
    nusinow: set[str]
    harmonized_gygi: set[str]
    procan: set[str]
    metabolomics: set[str]
    prism: set[str]
    ctrp: set[str]
    library_class: pd.Series  # screened -> class
    n_screens: pd.Series
    screen_type: pd.Series
    medium_source: pd.Series  # screen / screen_mixed / onboarded / imputed
    base_medium_class: pd.Series
    extra: dict = field(default_factory=dict)

    @property
    def proteomics_any(self) -> set[str]:
        return self.nusinow | self.harmonized_gygi | self.procan

    def lineage_counts(self) -> dict[str, int]:
        return {k: int(v) for k, v in self.lineage.value_counts().items()}

    def by_lineage(self, ids: set[str]) -> dict[str, int]:
        return {L: int(len(set(g.index) & ids)) for L, g in self.lineage.groupby(self.lineage)}


def _default_ids(path: Path) -> set[str]:
    df = pd.read_csv(path, usecols=["ModelID", "IsDefaultEntryForModel"])
    df = df[df.IsDefaultEntryForModel == "Yes"]
    assert df.ModelID.is_unique, path
    return set(df.ModelID)


def _lib_class(s: frozenset) -> str:
    if "Avana" in s:
        return "Avana-containing"
    if s == {"KY"}:
        return "KY-only"
    if s == {"Humagne-CD"}:
        return "Humagne-only"
    return "other"


def _base_of(text: str) -> str:
    return text.split(" + ")[0].strip()


def _base_class(text: str) -> str:
    b = _base_of(text)
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


INPUT_FILES = [  # every file load() reads, relative to the data root (pinned by make_splits, review M3)
    "depmap/CRISPRScreenMap.csv",
    "depmap/Model.csv",
    "depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
    "depmap/OmicsCNGeneWGS.csv",
    "depmap/harmonized_MS_CCLE_Gygi.csv",
    "depmap/ScreenSequenceMap.csv",
    "depmap/ModelCondition.csv",
    "proteomics_nusinow/protein_quant_current_normalized.csv.gz",
    "proteomics_goncalves/Protein_matrix_averaged_20250211.tsv",
    "metabolomics_ccle/CCLE_metabolomics_20190502.csv",
    "prism/repurposing_24Q2/Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv",
    "prism/repurposing_20Q2_secondary/prism-repurposing-20q2-secondary-screen-cell-line-info.csv",
    "ctrp/CTRPAUCMatrix.csv",
]


def load(data_root: Path | None = None) -> Cohort:
    """Every input in INPUT_FILES is mandatory (no silent 'if it exists' fallbacks; review M3)."""
    D = Path(data_root) if data_root else paths.data_root()
    missing = [f for f in INPUT_FILES if not (D / f).exists()]
    if missing:
        raise FileNotFoundError(f"split inputs missing: {missing}")
    DM = D / "depmap"
    screened = sorted(pd.read_csv(DM / "CRISPRScreenMap.csv").ModelID.unique())
    model = pd.read_csv(DM / "Model.csv").set_index("ModelID")
    scr = model.loc[screened]

    expr = _default_ids(DM / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
    cn = _default_ids(DM / "OmicsCNGeneWGS.csv")
    with gzip.open(D / "proteomics_nusinow/protein_quant_current_normalized.csv.gz", "rt") as f:
        nus_cols = f.readline().strip().split(",")
    nus_names = {re.sub(r"_TenPx\d+$", "", c.strip('"')) for c in nus_cols if "_TenPx" in c and not c.startswith('"TenPx')}
    ccle2ach = model.reset_index().dropna(subset=["CCLEName"]).groupby("CCLEName").ModelID.agg(set)
    nus = set().union(*[ccle2ach[n] for n in nus_names if n in ccle2ach.index])
    harm_path = DM / "harmonized_MS_CCLE_Gygi.csv"
    harm = set(pd.read_csv(harm_path, usecols=[0]).iloc[:, 0])
    pc_path = D / "proteomics_goncalves/Protein_matrix_averaged_20250211.tsv"
    pc = pd.read_csv(pc_path, sep="\t", usecols=[1], skiprows=3, header=None).iloc[:, 0]
    sid2ach = model.reset_index().dropna(subset=["SangerModelID"]).groupby("SangerModelID").ModelID.agg(set)
    procan: set[str] = set().union(*[sid2ach[s] for s in pc if s in sid2ach.index])
    metab = set(pd.read_csv(D / "metabolomics_ccle/CCLE_metabolomics_20190502.csv", usecols=["DepMap_ID"]).DepMap_ID.dropna())
    prism24 = set(pd.read_csv(D / "prism/repurposing_24Q2/Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv", usecols=["depmap_id"]).depmap_id.dropna())
    prism20 = set(pd.read_csv(D / "prism/repurposing_20Q2_secondary/prism-repurposing-20q2-secondary-screen-cell-line-info.csv", usecols=["depmap_id"]).depmap_id.dropna())
    ctrp_path = D / "ctrp/CTRPAUCMatrix.csv"
    ctrp = set(pd.read_csv(ctrp_path, usecols=[0]).iloc[:, 0])

    smap = pd.read_csv(DM / "CRISPRScreenMap.csv")
    ssm = pd.read_csv(DM / "ScreenSequenceMap.csv")
    ssm_ok = ssm[(ssm.PassesQC == True) & (ssm.ExcludeFromCRISPRCombined == False) & (ssm.ScreenType != "pDNA")]  # noqa: E712
    ssm_ok = ssm_ok[ssm_ok.ScreenID.isin(smap.ScreenID)]
    lib = ssm_ok.groupby("ModelID").Library.agg(lambda s: frozenset(s)).map(_lib_class)
    stype = ssm_ok.groupby("ModelID").ScreenType.agg(lambda s: "+".join(sorted(set(s))))
    nscreens = ssm_ok.groupby("ModelID").ScreenID.nunique()

    mc = pd.read_csv(DM / "ModelCondition.csv", usecols=["ModelConditionID", "FormulationID", "GrowthMedia", "SerumFreeMedia"])
    # content assertions: the text formulation lives in the column named FormulationID (names swapped)
    assert mc.FormulationID.dropna().str.contains(r"FBS|RPMI|DMEM", regex=True).mean() > 0.5
    assert mc.GrowthMedia.dropna().str.startswith("MF-").mean() > 0.9
    assert model.FormulationID.dropna().str.contains(r"FBS|RPMI|DMEM", regex=True).mean() > 0.5
    scrmed = ssm_ok.merge(mc, on="ModelConditionID", how="left").dropna(subset=["FormulationID"])
    scrmed_sets = scrmed.groupby("ModelID").FormulationID.agg(lambda s: list(s.value_counts().index))
    tier, base = {}, {}
    for mid in screened:
        if mid in scrmed_sets.index:
            texts = scrmed_sets[mid]
            if len({_base_of(t) for t in texts}) > 1:
                tier[mid], text = "screen_mixed", None
            else:
                tier[mid], text = "screen", texts[0]
        elif pd.notna(model.loc[mid, "FormulationID"]):
            tier[mid], text = "onboarded", model.loc[mid, "FormulationID"]
        else:
            tier[mid], text = "imputed", None
        base[mid] = _base_class(text) if text else tier[mid]

    idx = pd.Index(screened, name="ModelID")
    return Cohort(
        screened=screened, lineage=scr.OncotreeLineage, model_type=scr.ModelType,
        expr=expr, cn=cn, nusinow=nus, harmonized_gygi=harm, procan=procan, metabolomics=metab,
        prism=prism24 | prism20, ctrp=ctrp,
        library_class=lib.reindex(idx).fillna("no_mapped_screen"),
        n_screens=nscreens.reindex(idx).fillna(0).astype(int),
        screen_type=stype.reindex(idx).fillna("none"),
        medium_source=pd.Series(tier).reindex(idx), base_medium_class=pd.Series(base).reindex(idx),
    )


def line_table(c: Cohort) -> pd.DataFrame:
    """One row per screened line with the §4.1 step-6 covariates (tier/stratum added by make_splits)."""
    df = pd.DataFrame(index=pd.Index(c.screened, name="ModelID"))
    df["lineage"] = c.lineage
    df["library_class"] = c.library_class
    df["n_screens"] = c.n_screens
    df["screen_type"] = c.screen_type
    df["model_type"] = c.model_type
    df["medium_source"] = c.medium_source
    df["base_medium_class"] = c.base_medium_class
    for name, s in [("has_expression", c.expr), ("has_cn", c.cn), ("has_nusinow", c.nusinow), ("has_procan", c.procan),
                    ("has_proteomics_any", c.proteomics_any), ("has_metabolomics", c.metabolomics),
                    ("has_prism", c.prism), ("has_ctrp", c.ctrp)]:
        df[name] = df.index.isin(s)
    return df.reset_index()
