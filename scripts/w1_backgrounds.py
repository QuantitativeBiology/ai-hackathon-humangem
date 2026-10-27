"""M1 leftover: evaluable genomic-background groups per tier (DESIGN.md §4.1 step 7, W4 "Genomic backgrounds";
group definitions from notes/apriori_lists_final.md §2 per amendment A2). Label-free.

Thresholds (fixed a priori in the list file): deep deletion relative CN < 0.3; amplification CN >= 2.0;
expression-null log2(TPM+1) < 1.0; "TPM < x" / "TPM >= x" are read as TPM, i.e. log2(TPM+1) vs log2(x+1).
Damaging = damaging-allele count >= 1 (OmicsSomaticMutationsMatrixDamaging); hotspot = hotspot matrix > 0;
IDH_mut from the long mutation table (IDH1 R132*, IDH2 R140*/R172*). Default entry per model.
Within a tier, a lineage is *eligible* for a group if it has >= 3 members and >= 3 non-members (DESIGN W4);
a group is *evaluable* in a tier if >= 2 lineages are eligible and those lineages hold >= 10 members.
Member status is evaluated on screened lines; a line with the defining omics missing is "unknown" (not a
member, not a non-member).

Usage: .venv/bin/python scripts/w1_backgrounds.py
Outputs: results/w1/background_groups.tsv (counts per group x tier), results/w1/background_members.tsv
(group x screened line membership: 1 / 0 / NA; label-free).
"""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w1_common import DEPMAP, OUT_DIR, PROJECT, atomic_write_tsv, current, run_main  # noqa: E402

DM = PROJECT / "data" / "depmap"
GENE_COL = re.compile(r"^(?P<symbol>.+) \((?P<entrez>\d+)\)$")
NULL = 1.0  # log2(TPM+1) < 1


def tpm(x):
    return np.log2(x + 1.0)


def read_genes(fname: str, symbols: set[str], all_cols: bool = False) -> pd.DataFrame:
    head = list(pd.read_csv(DM / fname, nrows=0).columns)
    keep = [c for c in head if (m := GENE_COL.match(c)) and m["symbol"] in symbols]
    use = ["ModelID", "IsDefaultEntryForModel"] + (([c for c in head if GENE_COL.match(c)]) if all_cols else keep)
    df = pd.read_csv(DM / fname, usecols=use)
    df = df[df.IsDefaultEntryForModel.astype(str).str.lower().isin(["yes", "true"])].set_index("ModelID")
    df = df.drop(columns="IsDefaultEntryForModel")
    df.columns = [GENE_COL.match(c)["symbol"] for c in df.columns]
    return df


def main() -> int:
    t0, c0 = time.time(), time.process_time()
    from nextgem import labels

    split = labels.load_split().set_index("ModelID")
    ids = split.index
    E = read_genes("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                   {"MTAP", "ASS1", "ASNS", "NAPRT", "ME2", "SQLE", "HK1", "HK2", "LDHB", "SLC7A11", "CDS1"}).reindex(ids)
    C = read_genes("OmicsCNGeneWGS.csv", {"MTAP", "CDKN2A", "KEAP1", "STK11", "PHGDH", "ME2", "VHL", "ENO1", "SMARCA4"}).reindex(ids)
    Dall = read_genes("OmicsSomaticMutationsMatrixDamaging.csv", set(), all_cols=True)
    burden = Dall.sum(axis=1)  # damaging-allele burden over all default models (NEG-3 rule)
    b90 = float(np.percentile(burden, 90))
    D = Dall[[g for g in ("KEAP1", "STK11", "VHL", "SMARCA4", "TTN") if g in Dall.columns]].reindex(ids)
    del Dall
    H = read_genes("OmicsSomaticMutationsMatrixHotspot.csv", {"NFE2L2", "KRAS"}).reindex(ids)
    mut = pd.read_csv(DM / "OmicsSomaticMutations.csv", usecols=["ModelID", "IsDefaultEntryForModel", "HugoSymbol", "ProteinChange"])
    mut = mut[mut.IsDefaultEntryForModel.astype(str).str.lower().isin(["yes", "true"]) & mut.HugoSymbol.isin(["IDH1", "IDH2"])]
    idh_ids = set(mut[((mut.HugoSymbol == "IDH1") & mut.ProteinChange.astype(str).str.match(r"p\.R132"))
                      | ((mut.HugoSymbol == "IDH2") & mut.ProteinChange.astype(str).str.match(r"p\.R(140|172)"))].ModelID)
    has_mutcalls = set(mut.ModelID) | set(H.dropna(how="all").index)

    def known(*series):
        ok = pd.Series(True, index=ids)
        for s in series:
            ok &= s.notna()
        return ok

    def any_known(*series):
        ok = pd.Series(False, index=ids)
        for s in series:
            ok |= s.notna()
        return ok

    def union(conds_and_data):
        """OR of conditions: member if any known condition is true; non-member only if all are known false."""
        mem = pd.Series(False, index=ids)
        allknown = pd.Series(True, index=ids)
        for cond, data in conds_and_data:
            mem |= cond.fillna(False) & data.notna()
            allknown &= data.notna()
        return pd.Series(np.where(mem, 1.0, np.where(allknown, 0.0, np.nan)), index=ids)

    def single(cond, data):
        return pd.Series(np.where(data.notna(), cond.astype(float), np.nan), index=ids)

    g = {}
    g["MTAP_loss"] = union([(C.MTAP < 0.3, C.MTAP), (E.MTAP < NULL, E.MTAP)])
    neg1 = (C.CDKN2A < 0.3) & (C.MTAP >= 0.5) & (E.MTAP >= tpm(2))
    g["NEG-1_CDKN2A_del_MTAP_intact"] = pd.Series(np.where(known(C.CDKN2A, C.MTAP, E.MTAP), neg1.astype(float), np.nan), index=ids)
    g["ASS1_silenced"] = single(E.ASS1 < NULL, E.ASS1)
    g["ASNS_low"] = single(E.ASNS < tpm(4), E.ASNS)
    g["NRF2_active"] = union([(D.KEAP1 >= 1, D.KEAP1), (C.KEAP1 < 0.3, C.KEAP1), (H.NFE2L2 > 0, H.NFE2L2)])
    g["STK11_loss"] = union([(D.STK11 >= 1, D.STK11), (C.STK11 < 0.3, C.STK11)])
    g["NAPRT_silenced"] = single(E.NAPRT < NULL, E.NAPRT)
    g["PHGDH_amp"] = single(C.PHGDH >= 2.0, C.PHGDH)
    g["ME2_loss"] = union([(C.ME2 < 0.3, C.ME2), (E.ME2 < NULL, E.ME2)])
    g["SQLE_null"] = single(E.SQLE < NULL, E.SQLE)
    hk = (E.HK1 < tpm(2)) & (E.HK2 >= tpm(4))
    g["HK1low_HK2high"] = pd.Series(np.where(known(E.HK1, E.HK2), hk.astype(float), np.nan), index=ids)
    g["LDHB_null"] = single(E.LDHB < NULL, E.LDHB)
    g["SLC7A11_null"] = single(E.SLC7A11 < NULL, E.SLC7A11)
    g["CDS1_null"] = single(E.CDS1 < NULL, E.CDS1)
    g["VHL_loss"] = union([(D.VHL >= 1, D.VHL), (C.VHL < 0.3, C.VHL)])
    g["ENO1_del (in-silico only)"] = single(C.ENO1 < 0.3, C.ENO1)
    idh = pd.Series([1.0 if m in idh_ids else (0.0 if m in has_mutcalls else np.nan) for m in ids], index=ids)
    g["IDH_mut (secondary, in-silico only)"] = idh
    ttn = (D.TTN >= 1) & (pd.Series(burden).reindex(ids) <= b90)
    g["NEG-3_TTN_damaging_nonhyper"] = pd.Series(np.where(D.TTN.notna(), ttn.astype(float), np.nan), index=ids)
    g["NEG-4_SMARCA4_loss"] = union([(D.SMARCA4 >= 1, D.SMARCA4), (C.SMARCA4 < 0.3, C.SMARCA4)])
    g["NEG-5_KRAS_hotspot (secondary)"] = single(H.KRAS > 0, H.KRAS)
    members = pd.DataFrame(g)

    rows = []
    for name in members:
        for tier in ("train", "validation", "test", "all"):
            sel = split.index if tier == "all" else split.index[split.tier == tier]
            m = members.loc[sel, name]
            lin = split.loc[sel, "lineage"]
            elig = []
            for L, mm in m.groupby(lin):
                nm, nn = int((mm == 1).sum()), int((mm == 0).sum())
                if nm >= 3 and nn >= 3:
                    elig.append((L, nm))
            n_in_elig = sum(n for _, n in elig)
            rows.append({"group": name, "tier": tier, "n_lines": len(sel), "n_members": int((m == 1).sum()),
                         "n_nonmembers": int((m == 0).sum()), "n_unknown": int(m.isna().sum()),
                         "n_eligible_lineages": len(elig), "members_in_eligible_lineages": n_in_elig,
                         "eligible_lineages": ";".join(f"{L}({n})" for L, n in elig),
                         "evaluable": bool(len(elig) >= 2 and n_in_elig >= 10)})
    tab = pd.DataFrame(rows)
    tab.loc[tab.group.str.contains("in-silico"), "evaluable"] = False
    tab["note"] = np.where(tab.group.str.contains("in-silico"), "in-silico only per apriori_lists_final §2",
                           np.where(tab.group.str.startswith("NEG"), "negative control", ""))
    atomic_write_tsv(tab, OUT_DIR / "background_groups.tsv")
    atomic_write_tsv(members.reset_index().rename(columns={"index": "ModelID"}), OUT_DIR / "background_members.tsv")
    current().output(OUT_DIR / "background_groups.tsv", OUT_DIR / "background_members.tsv")
    wide = tab.pivot(index="group", columns="tier", values="evaluable")
    print(tab[tab.tier != "all"].pivot(index="group", columns="tier", values="n_members").join(wide, rsuffix="_eval").to_string())
    print(f"damaging burden 90th percentile = {b90}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run_main(main, inputs=[DEPMAP / f for f in (
        "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", "OmicsCNGeneWGS.csv", "OmicsSomaticMutationsMatrixDamaging.csv",
        "OmicsSomaticMutationsMatrixHotspot.csv", "OmicsSomaticMutations.csv")] + [PROJECT / "results" / "splits" / "lineage_split.tsv"]))
