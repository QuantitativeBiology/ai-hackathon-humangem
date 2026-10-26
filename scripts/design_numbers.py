"""Compute every numeric data fact quoted in DESIGN.md -> results/design/design_numbers.json.

Label policy (DESIGN.md §11.1, PI decisions r1-14/26, r2-14/15): this script opens NO label-bearing
file (CRISPRGeneEffect, CRISPRGeneDependency, ScreenGeneEffect, common-essential/control lists).
The screened-line set comes from CRISPRScreenMap.csv (screen -> model map, label-free); the gene
universe comes from results/baseline/*_gene_mapping.tsv. It computes no performance and no prevalence.
The split section simulates the draw function on seeds 1..1000 only; the real split seed is supplied
by the user after the M1 freeze commit (DESIGN.md §4.1) and is never evaluated here.

Run: .venv/bin/python scripts/design_numbers.py   (~1-2 min; reads only a few columns of the large omics files)
"""

from __future__ import annotations

import gzip
import json
import re
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data"
DM = D / "depmap"
OUT = ROOT / "results" / "design" / "design_numbers.json"

RETIRED_SPLIT_SEED = 20261026  # exposed in review round 1; make_splits.py refuses it
SIM_SEEDS = range(1, 1001)
MIN_LINES_ELIGIBLE = 15
STRATA = {"S1": (0, 8), "S2": (8, 16), "S3": (16, 23)}  # rank slices of eligible lineages
ALLOC_TEST = {"S1": 2, "S2": 2, "S3": 3}  # PI decision 3: 7 test lineages
ALLOC_VAL = {"S1": 1, "S2": 1, "S3": 2}  # 4 validation lineages

# Pre-registered thresholds (DESIGN.md §W4); log2(TPM+1) scale for expression, relative CN for CN.
CN_DEL = 0.3
CN_DEL_SENS = (0.1, 0.5)
CN_AMP = 2.0
EXPR_SILENT = 1.0  # log2(TPM+1) < 1  <=>  TPM < 1

SHEET_BASES = {"RPMI", "DMEM", "DMEM:F12", "MEM", "EMEM", "IMDM", "F12", "McCoy's 5A", "L-15"}
# formulation-sheet file patterns in data/media/ (case-insensitive) -> base medium class
SHEET_FILE_PATTERNS = {"RPMI": r"^rpmi", "DMEM": r"^dmem_highglucose", "DMEM:F12": r"^dmem_f12",
                       "IMDM": r"^imdm", "MEM/EMEM": r"^e?mem_", "F12": r"^(f12|ham)", "McCoy's 5A": r"mccoy",
                       "L-15": r"^(l-?15|leibovitz)"}
CN_LOH_PROXY = 0.5  # no LOH file on disk: heterozygous damaging + relCN < 0.5 is the LOH proxy

PANEL_GENES = ["PPAT", "GART", "PFAS", "PAICS", "ADSL", "ATIC", "IMPDH2", "GMPS", "TYMS", "DHFR",
               "DHODH", "UMPS", "CAD", "SQLE", "LSS", "DHCR7", "HMGCR", "MVK", "HMGCS1", "FDPS",
               "FASN", "ACACA", "SCD", "ASNS", "GOT1", "GOT2", "PHGDH", "PSAT1", "NAMPT", "SLC7A11",
               "ACADL", "HADHA", "ECHS1", "CYP4F2", "BBOX1", "TMLHE"]

out: dict = {"_meta": {"script": "scripts/design_numbers.py", "generated": time.strftime("%Y-%m-%d %H:%M"),
                       "label_values_read": False, "performance_computed": False}}


def header(path: Path) -> list[str]:
    return list(pd.read_csv(path, nrows=0).columns)


def gene_cols(cols: list[str], symbols: list[str]) -> dict[str, str]:
    """Map symbol -> 'SYMBOL (Entrez)' column (exact symbol match)."""
    res = {}
    for c in cols:
        m = re.match(r"^(\S+) \((\d+)\)$", c)
        if m and m.group(1) in symbols:
            res[m.group(1)] = c
    return res


def default_rows(path: Path, symbols: list[str]) -> pd.DataFrame:
    cols = header(path)
    gc = gene_cols(cols, symbols)
    df = pd.read_csv(path, usecols=["ModelID", "IsDefaultEntryForModel", *gc.values()])
    df = df[df.IsDefaultEntryForModel == "Yes"].drop(columns="IsDefaultEntryForModel")
    assert df.ModelID.is_unique
    return df.set_index("ModelID").rename(columns={v: k for k, v in gc.items()})


# ---------------------------------------------------------------- cohort / lineages
screened = pd.Series(sorted(pd.read_csv(DM / "CRISPRScreenMap.csv").ModelID.unique()))
model = pd.read_csv(DM / "Model.csv")
model = model.set_index("ModelID")
scr = model.loc[screened.values]
lin_counts = scr.OncotreeLineage.value_counts()
eligible = lin_counts[lin_counts >= MIN_LINES_ELIGIBLE]
rare = lin_counts[lin_counts < MIN_LINES_ELIGIBLE]
out["cohort"] = {
    "n_screened_lines": int(len(screened)),
    "n_lineages": int(len(lin_counts)),
    "lineage_counts": {k: int(v) for k, v in lin_counts.items()},
    "n_eligible_lineages": int(len(eligible)),
    "n_lines_eligible_lineages": int(eligible.sum()),
    "rare_lineages": {k: int(v) for k, v in rare.items()},
    "n_lines_rare": int(rare.sum()),
    "model_type_screened": {k: int(v) for k, v in scr.ModelType.value_counts().items()},
}

# ---------------------------------------------------------------- gene universe
maps = {m: pd.read_csv(ROOT / f"results/baseline/{m}_gene_mapping.tsv", sep="\t", dtype=str) for m in ("human1", "human2")}
scr_genes = {m: set(df.depmap_column.dropna()) for m, df in maps.items()}
U = scr_genes["human1"] | scr_genes["human2"]
Ui = scr_genes["human1"] & scr_genes["human2"]
out["genes"] = {
    "human1_model_genes": int(len(maps["human1"])),
    "human2_model_genes": int(len(maps["human2"])),
    "human1_screened": int(len(scr_genes["human1"])),
    "human2_screened": int(len(scr_genes["human2"])),
    "U_union": int(len(U)),
    "U_intersection": int(len(Ui)),
    "U_only_human1": int(len(scr_genes["human1"] - scr_genes["human2"])),
    "U_only_human2": int(len(scr_genes["human2"] - scr_genes["human1"])),
}

# ---------------------------------------------------------------- baseline KO facts (label-free)
ko = {}
for m in ("human1", "human2"):
    k = pd.read_parquet(ROOT / f"results/baseline/{m}_ko.parquet")
    r = k.ratio.clip(0, 1).round(6)
    bysym = dict(zip(k.symbol, r))
    ko[m] = {
        "genes": int(len(k)), "lethal_lt_0.01": int((r < 0.01).sum()), "ratio_eq_1": int((r == 1).sum()),
        "distinct_ratio_r": int(r.nunique()), "solver_values": sorted(k.solver.unique().tolist()),
        "panel_ratio_r": {g: (None if g not in bysym else float(bysym[g])) for g in PANEL_GENES},
    }
summ = pd.read_csv(ROOT / "results/baseline/summary.tsv", sep="\t").set_index("model")
for m in ("human1", "human2"):
    ko[m]["genes_closing_no_reaction"] = int(summ.loc[m, "genes_closing_no_reaction"])
    ko[m]["knockouts_wall_s_15workers"] = float(summ.loc[m, "knockouts_wall_s"])
    ko[m]["addressable_genes_ko_closes_reaction"] = ko[m]["genes"] - ko[m]["genes_closing_no_reaction"]
out["baseline_ko"] = ko

# Ham's (GR task) inputs and task lists
tasks_dir = D / "Human-GEM-main/data/metabolicTasks"
ess = pd.read_csv(tasks_dir / "metabolicTasks_Essential.txt", sep="\t", header=None, names=range(30), dtype=str, keep_default_na=False)
gr_start = ess.index[ess[1] == "GR"][0]
gr_inputs = []
for i in range(gr_start, len(ess)):
    if i > gr_start and ess.loc[i, 1] != "":
        break
    if ess.loc[i, 4]:  # columns: 0 '', 1 ID, 2 DESCRIPTION, 3 SHOULD FAIL, 4 IN
        gr_inputs.append(ess.loc[i, 4])
task_ids = ess[(ess[1] != "") & (ess[0] != "#") & (ess[1] != "ID")][1]
groups = ess[(ess[0] == "#") & (ess[1] != "") & (~ess[1].str.contains(":"))][1]
cf = pd.read_csv(tasks_dir / "metabolicTasks_CellfieConsensus.txt", sep="\t", header=None, names=range(30), dtype=str, keep_default_na=False)
out["tasks"] = {
    "essential_n_tasks": int(len(task_ids)),
    "essential_by_category": {k: int(v) for k, v in task_ids.value_counts().items()},
    "essential_n_groups": int(len(groups)),
    "essential_group_names": groups.tolist(),
    "cellfie_n_tasks": int(((cf[1] != "") & (cf[1] != "ID")).sum()),
    "hams_gr_inputs_n": len(gr_inputs),
    "hams_gr_inputs": gr_inputs,
}

# Human-GEM repository Hart-2015 reference values (read, not typed)
summary_md = (D / "Human-GEM-main/data/testResults/gene-essential_summary.md").read_text()
hart = {}
for line in summary_md.splitlines():
    cells = [c.strip() for c in line.strip().strip("|").split("|")]
    if len(cells) == 9 and cells[0] not in ("cellLine", "---"):
        hart[cells[0]] = {"growthAUPRC": float(cells[6]), "baseRate": float(cells[7]),
                          "viabilityMCC": float(cells[3]) if cells[3] else None}
out["hart_reference"] = hart

# ---------------------------------------------------------------- omics coverage of screened lines
S = set(screened)
LOSS_GENES = ["KEAP1", "STK11", "ME2", "FH", "SDHA", "SDHB", "SDHC", "SDHD", "ENO1"]
expr = default_rows(DM / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", LOSS_GENES +
                    ["NAPRT", "ASS1", "MTAP", "ENO1", "ME2", "PHGDH", "SLC1A3", "SLC1A1", "SLC1A2"])
cn = default_rows(DM / "OmicsCNGeneWGS.csv", ["MTAP", "ENO1", "ME2", "PHGDH", "KEAP1", "STK11", "FH",
                                               "SDHA", "SDHB", "SDHC", "SDHD", "CDKN2A", "SMAD4"])
dmg = default_rows(DM / "OmicsSomaticMutationsMatrixDamaging.csv", ["KEAP1", "STK11", "ME2", "FH", "SDHA", "SDHB", "SDHC", "SDHD", "ENO1"])
hot = default_rows(DM / "OmicsSomaticMutationsMatrixHotspot.csv", ["IDH1", "IDH2", "NFE2L2", "KRAS", "BRAF", "PIK3CA"])

# proteomics
with gzip.open(D / "proteomics_nusinow/protein_quant_current_normalized.csv.gz", "rt") as f:
    nus_cols = f.readline().strip().split(",")
nus_names = {re.sub(r"_TenPx\d+$", "", c.strip('"')) for c in nus_cols if "_TenPx" in c and not c.startswith('"TenPx')}
ccle2ach = model.reset_index().dropna(subset=["CCLEName"]).groupby("CCLEName").ModelID.agg(set)
nus = set().union(*[ccle2ach[n] for n in nus_names if n in ccle2ach.index])
harm = set(pd.read_csv(DM / "harmonized_MS_CCLE_Gygi.csv", usecols=[0]).iloc[:, 0]) if (DM / "harmonized_MS_CCLE_Gygi.csv").exists() else set()
pc_path = D / "proteomics_goncalves/Protein_matrix_averaged_20250211.tsv"
procan = set()
procan_unmapped = None
if pc_path.exists():
    pc = pd.read_csv(pc_path, sep="\t", usecols=[1], skiprows=3, header=None).iloc[:, 0]
    sid2ach = model.reset_index().dropna(subset=["SangerModelID"]).groupby("SangerModelID").ModelID.agg(set)
    procan = set().union(*[sid2ach[s] for s in pc if s in sid2ach.index])
    procan_unmapped = int(sum(s not in sid2ach.index for s in pc))
metab = set(pd.read_csv(D / "metabolomics_ccle/CCLE_metabolomics_20190502.csv", usecols=["DepMap_ID"]).DepMap_ID.dropna())
prism24 = set(pd.read_csv(D / "prism/repurposing_24Q2/Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv", usecols=["depmap_id"]).depmap_id.dropna())
prism20 = set(pd.read_csv(D / "prism/repurposing_20Q2_secondary/prism-repurposing-20q2-secondary-screen-cell-line-info.csv", usecols=["depmap_id"]).depmap_id.dropna())
ctrp_path = D / "ctrp/CTRPAUCMatrix.csv"
ctrp = set(pd.read_csv(ctrp_path, usecols=[0]).iloc[:, 0]) if ctrp_path.exists() else set()
E, C = set(expr.index), set(cn.index)
# GDSC (xlsx; openpyxl pinned per PI decision r2-17): join on SangerModelID, cross-check on COSMICID
gdsc_cov = {}
sid2ach_all = model.reset_index().dropna(subset=["SangerModelID"]).groupby("SangerModelID").ModelID.agg(set)
for tag, fn in (("gdsc1", "GDSC1_fitted_dose_response_27Oct23.xlsx"), ("gdsc2", "GDSC2_fitted_dose_response_27Oct23.xlsx")):
    g = pd.read_excel(D / "gdsc" / fn, usecols=["SANGER_MODEL_ID", "COSMIC_ID"])
    ids = set().union(*[sid2ach_all[x] for x in g.SANGER_MODEL_ID.dropna().unique() if x in sid2ach_all.index])
    gdsc_cov[f"{tag}_lines_total"] = int(g.SANGER_MODEL_ID.nunique())
    gdsc_cov[f"{tag}_unmapped_sanger_ids"] = int(sum(x not in sid2ach_all.index for x in g.SANGER_MODEL_ID.dropna().unique()))
    gdsc_cov[f"{tag}_screened"] = len(set(screened) & ids)
    gdsc_cov[f"_{tag}_set"] = ids
gdsc_any = gdsc_cov.pop("_gdsc1_set") | gdsc_cov.pop("_gdsc2_set")
gdsc_cov["gdsc_any_screened"] = len(set(screened) & gdsc_any)
gdsc_cov["gdsc_any_screened_with_expr"] = len(set(screened) & gdsc_any & E)
cld = pd.read_excel(D / "gdsc" / "Cell_Lines_Details.xlsx", sheet_name="Cell line details", usecols=["Screen Medium"])
gdsc_cov["gdsc_screen_medium_values"] = {str(k): int(v) for k, v in cld["Screen Medium"].value_counts(dropna=False).items()}
out["coverage_screened"] = {
    "expression_default": len(S & E), "no_expression": len(S - E),
    "cn_wgs_default": len(S & C), "no_cn": len(S - C),
    "mutation_matrix": len(S & set(dmg.index)),
    "expr_and_cn": len(S & E & C),
    "nusinow_raw_mapped_total": len(nus), "nusinow_screened": len(S & nus), "nusinow_screened_with_expr": len(S & nus & E),
    "harmonized_gygi_total": len(harm), "harmonized_gygi_screened": len(S & harm),
    "procan_mapped_total": len(procan), "procan_screened": len(S & procan), "procan_unmapped_ids": procan_unmapped,
    "any_proteomics_screened": len(S & (nus | harm | procan)),
    "metabolomics_total": len(metab), "metabolomics_screened": len(S & metab),
    "prism_24q2_or_20q2_screened": len(S & (prism24 | prism20)),
    "ctrp_screened": len(S & ctrp),
    **gdsc_cov,
}
lin_of = scr.OncotreeLineage
cov_by_lin = {}
for L, n in lin_counts.items():
    ids = set(lin_of[lin_of == L].index)
    cov_by_lin[L] = {"screened": int(n), "expr": len(ids & E), "cn": len(ids & C),
                     "proteomics_any": len(ids & (nus | harm | procan)), "nusinow": len(ids & nus),
                     "metabolomics": len(ids & metab), "prism": len(ids & (prism24 | prism20))}
out["coverage_by_lineage"] = cov_by_lin

# ---------------------------------------------------------------- resource scope (W8)
expr_all = set(expr.index)
ma = model.loc[sorted(expr_all & set(model.index))]
noncancer = ma.OncotreeLineage.isin(["Fibroblast", "Normal"]) | ma.OncotreePrimaryDisease.astype(str).str.contains("Non-Cancerous", case=False)
out["resource_scope"] = {"default_expression_models": len(expr_all),
                         "excluded_noncancer": int(noncancer.sum()),
                         "resource_models": int((~noncancer).sum()),
                         "resource_models_screened": int(len(set(ma[~noncancer].index) & S)),
                         "resource_models_not_benchmarkable": int(len(set(ma[~noncancer].index) - S))}

# ---------------------------------------------------------------- screens: library, type, medium tiers
smap = pd.read_csv(DM / "CRISPRScreenMap.csv")
ssm = pd.read_csv(DM / "ScreenSequenceMap.csv")
ssm_ok = ssm[(ssm.PassesQC == True) & (ssm.ExcludeFromCRISPRCombined == False) & (ssm.ScreenType != "pDNA")]  # noqa: E712
ssm_ok = ssm_ok[ssm_ok.ScreenID.isin(smap.ScreenID)]
lib = ssm_ok.groupby("ModelID").Library.agg(lambda s: frozenset(s))


def lib_class(s: frozenset) -> str:
    if "Avana" in s:
        return "Avana-containing"
    if s == {"KY"}:
        return "KY-only"
    if s == {"Humagne-CD"}:
        return "Humagne-only"
    return "other"


libc = lib.map(lib_class)
stype = ssm_ok.groupby("ModelID").ScreenType.agg(lambda s: "+".join(sorted(set(s))))
nscreens = ssm_ok.groupby("ModelID").ScreenID.nunique()
out["screens"] = {
    "library_class_counts": {k: int(v) for k, v in libc.reindex(screened).fillna("no_mapped_screen").value_counts().items()},
    "screen_type_counts": {k: int(v) for k, v in stype.reindex(screened).fillna("none").value_counts().items()},
    "models_with_2plus_screens": int((nscreens.reindex(screened).fillna(0) >= 2).sum()),
}

mc = pd.read_csv(DM / "ModelCondition.csv", usecols=["ModelConditionID", "FormulationID", "GrowthMedia", "SerumFreeMedia"])
# content assertion: the text formulation lives in the column named FormulationID (names swapped)
assert mc.FormulationID.dropna().str.contains(r"FBS|RPMI|DMEM", regex=True).mean() > 0.5
assert mc.GrowthMedia.dropna().str.startswith("MF-").mean() > 0.9
assert model.FormulationID.dropna().str.contains(r"FBS|RPMI|DMEM", regex=True).mean() > 0.5
scrmed = ssm_ok.merge(mc, on="ModelConditionID", how="left").dropna(subset=["FormulationID"])
# modal text per model (data r2: several formulations sharing one base -> modal text, flagged)
scrmed_sets = scrmed.groupby("ModelID").FormulationID.agg(lambda s: list(s.value_counts().index))


def base_of(text: str) -> str:
    return text.split(" + ")[0].strip()


def base_class(text: str) -> str:
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


tier, text_of, base_cls, multi_supp = {}, {}, {}, {}
for mid in screened:
    if mid in scrmed_sets.index:
        texts = scrmed_sets[mid]
        bases = {base_of(t) for t in texts}
        if len(bases) > 1:
            tier[mid] = "screen_mixed"
            text_of[mid] = None
        else:
            tier[mid] = "screen"
            text_of[mid] = texts[0]  # modal formulation text
            multi_supp[mid] = len(texts) > 1
    elif pd.notna(model.loc[mid, "FormulationID"]):
        tier[mid] = "onboarded"
        text_of[mid] = model.loc[mid, "FormulationID"]
    else:
        tier[mid] = "imputed"
        text_of[mid] = None
    base_cls[mid] = base_class(text_of[mid]) if text_of[mid] else tier[mid]
tier_s, base_s = pd.Series(tier), pd.Series(base_cls)
disagree = sum(1 for mid in screened if tier[mid] == "screen" and pd.notna(model.loc[mid, "FormulationID"])
               and model.loc[mid, "FormulationID"] != text_of[mid])
both = sum(1 for mid in screened if tier[mid] == "screen" and pd.notna(model.loc[mid, "FormulationID"]))
serum_free = sum(1 for mid in screened if text_of[mid] and "FBS" not in text_of[mid])
kyonly = set(libc[libc == "KY-only"].index)
out["medium"] = {
    "tier_counts": {k: int(v) for k, v in tier_s.value_counts().items()},
    "screen_vs_onboarded_text_disagree": int(disagree), "screen_and_onboarded_both_present": int(both),
    "base_class_counts": {k: int(v) for k, v in base_s.value_counts().items()},
    "mappable_to_requested_sheets": int(base_s.str.startswith("sheet:").sum()),
    "mappable_incl_mixtures": int((base_s.str.startswith("sheet:") | (base_s == "mixture_of_sheets")).sum()),
    "no_fbs_in_text": int(serum_free),
    "screen_multi_supplement": int(sum(multi_supp.values())),
    "screen_tier_sheet_mappable": int(((tier_s == "screen") & base_s.str.startswith("sheet:")).sum()),
    "screen_tier_sheet_or_mixture": int(((tier_s == "screen") & (base_s.str.startswith("sheet:") | (base_s == "mixture_of_sheets"))).sum()),
    "onboarded_tier_sheet_or_mixture": int(((tier_s == "onboarded") & (base_s.str.startswith("sheet:") | (base_s == "mixture_of_sheets"))).sum()),
    "unmappable_organoid_stem_other": int(base_s.isin(["organoid", "stem", "other_defined"]).sum()),
    "kyonly_lines": len(kyonly), "kyonly_without_screen_medium": int(sum(tier[m] != "screen" for m in kyonly)),
    "media_files_present": sorted(p.name for p in (D / "media").glob("*")),
}
# within-lineage medium contrast power: lineages with >=8 lines in each of two sheet base media
lin_med = pd.DataFrame({"lin": lin_of, "base": base_s.reindex(lin_of.index)})
pow_lins = []
for L, g in lin_med.groupby("lin"):
    vc = g.base[g.base.str.startswith("sheet:")].value_counts()
    if (vc >= 8).sum() >= 2:
        pow_lins.append(L)
out["medium"]["lineages_with_2_media_ge8"] = pow_lins
out["medium"]["n_lineages_with_2_media_ge8"] = len(pow_lins)
# PI decision r2-9: attribution set = screen-tier lines whose base-medium sheet is on disk
media_files = [p.name.lower() for p in (D / "media").glob("*")]
on_disk = sorted(b for b, pat in SHEET_FILE_PATTERNS.items() if any(re.search(pat, f) for f in media_files))
on_disk_cls = {"sheet:" + b for b in on_disk}
attr = (tier_s == "screen") & base_s.isin(on_disk_cls)
out["medium"]["sheet_bases_on_disk"] = on_disk
out["medium"]["sheet_bases_missing"] = sorted(set(SHEET_FILE_PATTERNS) - set(on_disk))
out["medium"]["attribution_set_screen_tier_on_disk"] = int(attr.sum())
out["medium"]["attribution_set_by_base"] = {k: int(v) for k, v in base_s[attr].value_counts().items()}
out["medium"]["fallback_lines_sheet_missing_any_tier"] = int((base_s.str.startswith("sheet:") & ~base_s.isin(on_disk_cls)).sum())
pow2 = []
for L, g in lin_med[attr.reindex(lin_med.index).fillna(False).astype(bool)].groupby("lin"):
    vc = g.base.value_counts()
    if (vc >= 8).sum() >= 2:
        pow2.append(L)
out["medium"]["lineages_with_2_media_ge8_on_disk_sheets"] = pow2
out["medium"]["n_lineages_with_2_media_ge8_on_disk_sheets"] = len(pow2)

# ---------------------------------------------------------------- genomic backgrounds (label-free)
def cnv(g):
    return cn[g] if g in cn.columns else pd.Series(dtype=float)


def dm(g):
    return dmg[g] if g in dmg.columns else pd.Series(dtype=float)


def defined(*series_):
    idx = set()
    for s in series_:
        idx |= set(s.dropna().index)
    return idx


def lost(g):
    """r2-11: homozygous damaging (code 2), or heterozygous damaging (code 1) with LOH proxy
    (relCN < 0.5) or TPM < 1, or homozygous deletion (relCN < 0.3)."""
    d, c = dm(g), cnv(g)
    e = expr[g] if g in expr.columns else pd.Series(dtype=float)
    idx = d.index.union(c.index)
    d, c, e = d.reindex(idx), c.reindex(idx), e.reindex(idx)
    return (d >= 2) | ((d == 1) & ((c < CN_LOH_PROXY) | (e < EXPR_SILENT))) | (c < CN_DEL)


def lost_sens(g):
    """v1 rule, kept as sensitivity: any damaging allele or relCN < 0.3."""
    return (dm(g) >= 1) | (cnv(g) < CN_DEL)


groups_def = {
    "MTAP_del": (lambda: cnv("MTAP") < CN_DEL, "relCN(MTAP)<0.3", "expected FBA failure (MTA/PRMT5)"),
    "ENO1_del": (lambda: cnv("ENO1") < CN_DEL, "relCN(ENO1)<0.3", "ENO2"),
    "ME2_loss": (lambda: lost("ME2"), "loss rule r2-11 (ME2)", "ME3"),
    "NAPRT_silenced": (lambda: expr["NAPRT"] < EXPR_SILENT, "log2(TPM+1)(NAPRT)<1", "NAMPT"),
    "ASS1_silenced": (lambda: expr["ASS1"] < EXPR_SILENT, "log2(TPM+1)(ASS1)<1", "arginine import"),
    "NRF2_active": (lambda: lost("KEAP1").reindex(hot.index.union(lost("KEAP1").index)).fillna(False) | (hot["NFE2L2"] >= 1).reindex(hot.index.union(lost("KEAP1").index)).fillna(False), "KEAP1 loss (r2-11) or NFE2L2 hotspot", "SLC7A11/cystine, NADPH"),
    "STK11_loss": (lambda: lost("STK11"), "loss rule r2-11 (STK11)", "CPS1/pyrimidine"),
    "FH_SDH_loss": (lambda: pd.concat([lost(g) for g in ("FH", "SDHA", "SDHB", "SDHC", "SDHD")], axis=1).fillna(False).any(axis=1), "loss rule r2-11 (FH/SDHx)", "TCA truncation"),
    "IDH_mut": (lambda: (hot["IDH1"] >= 1) | (hot["IDH2"] >= 1), "IDH1/2 hotspot", "2-HG (needs added reaction)"),
    "PHGDH_amp": (lambda: cnv("PHGDH") >= CN_AMP, "relCN(PHGDH)>=2", "serine synthesis"),
    "KRAS_hotspot_negctrl": (lambda: hot["KRAS"] >= 1, "KRAS hotspot", "signalling negative control"),
}
bg = {}
group_members = {}
for name, (fn, rule, mech) in groups_def.items():
    s = fn().astype(bool)
    s = s[s.index.isin(S)]
    members = set(s[s].index)
    group_members[name] = members
    known = set(s.index)
    elig_lins, n_in_elig = [], 0
    for L in lin_counts.index:
        ids = set(lin_of[lin_of == L].index) & known
        k = len(ids & members)
        if k >= 3 and len(ids) - k >= 3:
            elig_lins.append(L)
            n_in_elig += k
    by_lin = Counter(lin_of[list(members)])
    bg[name] = {"rule": rule, "mechanism": mech, "n_screened_members": len(members),
                "n_screened_status_known": len(known),
                "top_lineages": dict(sorted(by_lin.items(), key=lambda kv: (-kv[1], kv[0]))[:4]),  # ties broken by name (deterministic)
                "within_lineage_eligible_lineages": len(elig_lins), "members_in_eligible_lineages": n_in_elig}
for g, (lo, hi) in {"MTAP": CN_DEL_SENS}.items():
    bg["MTAP_del"]["sens_relCN_0.1"] = int(((cnv(g) < lo) & cnv(g).index.isin(S)).sum())
    bg["MTAP_del"]["sens_relCN_0.5"] = int(((cnv(g) < hi) & cnv(g).index.isin(S)).sum())
out["backgrounds"] = bg
out["loss_rule_sensitivity_v1_rule_counts"] = {
    g: int(lost_sens(g).fillna(False)[lambda s: s].index.isin(S).sum()) for g in ("ME2", "KEAP1", "STK11", "FH")}
out["loss_rule_primary_counts"] = {
    g: int(lost(g).fillna(False)[lambda s: s].index.isin(S).sum()) for g in ("ME2", "KEAP1", "STK11", "FH")}
out["mutation_matrix_values"] = {"damaging_STK11": {str(k): int(v) for k, v in dmg["STK11"].value_counts().items()}}

# ---------------------------------------------------------------- split simulation (label-free)
elig_sorted = sorted(eligible.items(), key=lambda kv: (-kv[1], kv[0]))
strata = {s: [l for l, _ in elig_sorted[a:b]] for s, (a, b) in STRATA.items()}
expr_by_lin = {L: len(set(lin_of[lin_of == L].index) & E) for L in lin_counts.index}
prot_any = nus | harm | procan
prot_by_lin = {L: len(set(lin_of[lin_of == L].index) & prot_any) for L in lin_counts.index}


def draw(seed: int):
    rng = np.random.default_rng(seed)
    test, val = [], []
    for s in ("S1", "S2", "S3"):
        perm = list(rng.permutation(strata[s]))
        test += perm[:ALLOC_TEST[s]]
        val += perm[ALLOC_TEST[s]:ALLOC_TEST[s] + ALLOC_VAL[s]]
    return test, val


def checks(test, val):
    ntest_expr = sum(expr_by_lin[L] for L in test)
    return {
        "a_min10_expr_each": all(expr_by_lin[L] >= 10 for L in test + val),
        "b_test_ge250_expr_lines": ntest_expr >= 250,
        "c_test_ge40_proteomic_lines": sum(prot_by_lin[L] for L in test) >= 40,
        "d_not_both_heme_heldout": not ({"Lymphoid", "Myeloid"} <= set(test + val)),
    }


sims = []
for sd in SIM_SEEDS:
    t, v = draw(sd)
    c = checks(t, v)
    ntest = sum(expr_by_lin[L] for L in t)
    sims.append({"ok": all(c.values()), **c, "test_expr": ntest,
                 "val_expr": sum(expr_by_lin[L] for L in v),
                 "largest_share": max(expr_by_lin[L] for L in t) / ntest,
                 "test_prot": sum(prot_by_lin[L] for L in t),
                 "test": t})
sd = pd.DataFrame(sims)
acc = sd[sd.ok]
feas = {}
for name, members in group_members.items():
    ok = 0
    for t in acc.test:
        tl = set(lin_of[lin_of.isin(t)].index)
        k = 0
        for L in t:
            ids = set(lin_of[lin_of == L].index)
            m = len(ids & members)
            if m >= 3 and len(ids) - m >= 3:
                k += m
        ok += k >= 5
    feas[name] = round(ok / len(acc), 3)
q = lambda s, p: float(np.quantile(s, p))  # noqa: E731
out["split_simulation"] = {
    "note": f"{len(sd)} simulated draws on seeds 1..1000; the real seed is user-supplied after the M1 freeze commit and is never evaluated here; {RETIRED_SPLIT_SEED} is retired",
    "strata": strata, "alloc_test": ALLOC_TEST, "alloc_val": ALLOC_VAL,
    "acceptance_rate": round(float(sd.ok.mean()), 3),
    "constraint_fail_rate": {k: round(float((~sd[k]).mean()), 3) for k in ("a_min10_expr_each", "b_test_ge250_expr_lines", "c_test_ge40_proteomic_lines", "d_not_both_heme_heldout")},
    "accepted_test_expr_lines_median_p05_p95": [q(acc.test_expr, .5), q(acc.test_expr, .05), q(acc.test_expr, .95)],
    "accepted_val_expr_lines_median_p05_p95": [q(acc.val_expr, .5), q(acc.val_expr, .05), q(acc.val_expr, .95)],
    "accepted_largest_lineage_share_median_p95": [round(q(acc.largest_share, .5), 3), round(q(acc.largest_share, .95), 3)],
    "accepted_test_proteomic_lines_median_p05": [q(acc.test_prot, .5), q(acc.test_prot, .05)],
    "background_test_feasibility_share": feas,
    "training_eligible_lineages_per_draw": 23 - 7 - 4,
    "reference_draws_seeds_1_to_20": {int(sd_): {"test": list(draw(sd_)[0]), "val": list(draw(sd_)[1])} for sd_ in range(1, 21)},
}

OUT.parent.mkdir(parents=True, exist_ok=True)
tmp = OUT.with_suffix(".tmp")
tmp.write_text(json.dumps(out, indent=1, default=lambda o: int(o) if isinstance(o, np.integer) else str(o)))
tmp.replace(OUT)
print(f"wrote {OUT}")
