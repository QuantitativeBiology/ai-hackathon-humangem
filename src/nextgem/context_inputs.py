"""Label-free inputs for W4 context-specific models (DESIGN.md §9 W4).

Everything here reads omics *inputs* only (expression, WGS CN, damaging mutations, ProCan
proteomics, the split file's label-free covariates). No label-bearing file is opened.

Conventions (DESIGN §9 W4, "Thresholds" and "Modalities"):

* RNA values are DepMap log2(TPM+1) on the default row of each model.
* A gene's *signed score* is ``s_g = value_g - threshold_g`` (log2 units); ``s_g < 0`` is "off".
  NaN means "no evidence" (the gene has no column in the source, or, for protein, it is a
  transporter/membrane protein whose absence is not evidence, see below).
* Threshold rules (RNA): R1 global TPM 1 (``log2(1+1) = 1``); R2 the global 25th percentile of all
  (model gene, training line) values; R3 the gene's mean over training lines, clipped to the
  global [P25, P75] (a label-free statistic, so no LOLO). The protein rules mirror R2/R3 on
  log2 ProCan intensities; protein R1 is "detected" (implementation choice, README §Open items).
* Protein absence (not detected in a line, while quantified in the dataset) counts as evidence only
  for *soluble enzymes*: genes with no transport reaction and no "Cell membrane"/"Inner
  mitochondria" localisation in Human-GEM ``genes.tsv`` (a proxy, flagged). Transporters and
  membrane proteins use RNA only.
* RNA+protein: protein overrides RNA wherever the protein score is defined.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import pyarrow.csv as pacsv

from nextgem import paths
from nextgem.config import SEEDS
from nextgem.genes import parse_depmap_columns

EXPR_FILE = "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
CN_FILE = "OmicsCNGeneWGS.csv"
DAMAGING_FILE = "OmicsSomaticMutationsMatrixDamaging.csv"
PROCAN_FILE = "proteomics_goncalves/Protein_matrix_averaged_20250211.tsv"

R1_LOG2 = 1.0  # TPM 1 -> log2(1 + 1)
THRESHOLD_RULES = ("R1", "R2", "R3")
MODALITIES = ("rna", "protein", "rna_protein")

# Loss rule (PI decision r2-11; same constants as scripts/design_numbers.py)
CN_DEL = 0.3
CN_DEL_SENS = (0.1, 0.5)
CN_LOH_PROXY = 0.5
EXPR_SILENT = 1.0

PANEL_LINEAGES = 10
PANEL_LINES_PER_LINEAGE = 3
TRANSPORT_SUBSYSTEM = "Transport reactions"
MEMBRANE_COMPARTMENTS = {"Cell membrane", "Inner mitochondria"}


# ----------------------------------------------------------------------------- gene table
def gene_table(model, genes_tsv=None) -> pd.DataFrame:
    """One row per model gene: entrez (primary), symbol, uniprot, transporter, soluble."""
    genes_tsv = genes_tsv or (paths.data_root() / "Human-GEM-main" / "model" / "genes.tsv")
    t = pd.read_csv(genes_tsv, sep="\t", dtype=str).set_index("genes")
    ids = sorted(g.id for g in model.genes)
    df = pd.DataFrame(index=pd.Index(ids, name="gene"))
    df["entrez"] = t.geneEntrezID.reindex(df.index).str.split(";").str[0]
    df["symbol"] = t.geneSymbols.reindex(df.index)
    df["uniprot"] = t.geneUniProtID.reindex(df.index)
    comps = t.compartments.reindex(df.index).fillna("")
    transporter = {g.id: any(r.subsystem == TRANSPORT_SUBSYSTEM for r in g.reactions) for g in model.genes}
    df["transporter"] = df.index.map(transporter).astype(bool)
    df["membrane"] = comps.map(lambda s: bool(set(s.split(";")) & MEMBRANE_COMPARTMENTS))
    df["soluble"] = ~df.transporter & ~df.membrane
    return df


# ----------------------------------------------------------------------------- readers
def _read_depmap_matrix(fname: str, genes: pd.DataFrame, lines=None) -> pd.DataFrame:
    """Default rows of a DepMap ``SYMBOL (Entrez)`` matrix, restricted to model genes.

    Returns lines x model gene (Ensembl ids); genes without a column are absent from the result.
    """
    path = paths.depmap_dir() / fname
    hdr = pd.read_csv(path, nrows=0).columns
    dm = parse_depmap_columns(hdr)
    dm = dm[dm.entrez.isin(set(genes.entrez.dropna()))]
    meta = ["ModelID", "IsDefaultEntryForModel"]
    tab = pacsv.read_csv(path, convert_options=pacsv.ConvertOptions(include_columns=meta + list(dm.column)))
    df = tab.to_pandas()
    df = df[df.IsDefaultEntryForModel == "Yes"].drop(columns="IsDefaultEntryForModel")
    assert df.ModelID.is_unique, fname
    df = df.set_index("ModelID")
    if lines is not None:
        df = df.reindex([x for x in lines if x in df.index])
    by_entrez = dict(zip(dm.entrez, dm.column))
    out = {g: df[by_entrez[e]] for g, e in genes.entrez.items() if e in by_entrez}
    return pd.DataFrame(out, index=df.index).astype("float64")


def read_expression(genes: pd.DataFrame, lines=None) -> pd.DataFrame:
    return _read_depmap_matrix(EXPR_FILE, genes, lines)


def read_cn(genes: pd.DataFrame, lines=None) -> pd.DataFrame:
    return _read_depmap_matrix(CN_FILE, genes, lines)


def read_damaging(genes: pd.DataFrame, lines=None) -> pd.DataFrame:
    return _read_depmap_matrix(DAMAGING_FILE, genes, lines)


def read_procan(genes: pd.DataFrame, lines=None) -> pd.DataFrame:
    """ProCan-DepMapSanger log2 intensities (lines x model gene), mapped Sanger ID -> ModelID.

    Proteins map to model genes by UniProt accession (``genes.tsv`` geneUniProtID). A protein
    quantified in the dataset but NaN in a line is "not detected" in that line; genes whose protein
    is not in the dataset are absent from the result. A Sanger model mapping to several ModelIDs
    is assigned to each of them.
    """
    path = paths.data_root() / PROCAN_FILE
    head = pd.read_csv(path, sep="\t", header=None, nrows=2, dtype=str)
    uniprot = head.iloc[0, 2:].tolist()
    body = pd.read_csv(path, sep="\t", header=None, skiprows=3, dtype={0: str, 1: str})
    body = body.set_index(1).drop(columns=0)
    body.columns = uniprot
    model_csv = pd.read_csv(paths.depmap_dir() / "Model.csv", usecols=["ModelID", "SangerModelID"]).dropna()
    rows = []
    for _, r in model_csv.iterrows():
        if r.SangerModelID in body.index and (lines is None or r.ModelID in set(lines)):
            rows.append((r.ModelID, r.SangerModelID))
    idx = pd.Index([a for a, _ in rows], name="ModelID")
    vals = body.loc[[s for _, s in rows]].to_numpy(dtype="float64")
    prot = pd.DataFrame(vals, index=idx, columns=uniprot)
    prot = prot.loc[~prot.index.duplicated()]
    up2gene = {u: g for g, u in genes.uniprot.dropna().items()}
    keep = [u for u in prot.columns if u in up2gene]
    prot = prot[keep]
    prot.columns = [up2gene[u] for u in keep]
    return prot.T.groupby(level=0).max().T  # two accessions -> one gene: keep the larger value


# ----------------------------------------------------------------------------- thresholds
@dataclass
class ThresholdStats:
    """Label-free threshold statistics from training lines (fitted on fitting lines only, §4.2)."""

    modality: str  # "rna" or "protein"
    p25: float
    p75: float
    gene_mean: pd.Series
    floor_value: float  # value used for "not detected" proteins (below every threshold)
    n_lines: int

    def thresholds(self, rule: str, genes) -> pd.Series:
        genes = pd.Index(genes)
        if rule == "R1":
            return pd.Series(R1_LOG2 if self.modality == "rna" else np.nan, index=genes)
        if rule == "R2":
            return pd.Series(self.p25, index=genes)
        if rule == "R3":
            return self.gene_mean.reindex(genes).clip(self.p25, self.p75).fillna(self.p25)
        raise ValueError(rule)


def threshold_stats(values: pd.DataFrame, modality: str) -> ThresholdStats:
    """``values``: training lines x genes (NaN = not detected for protein)."""
    flat = values.to_numpy().ravel()
    flat = flat[~np.isnan(flat)]
    p25, p75 = np.percentile(flat, [25, 75])
    return ThresholdStats(modality, float(p25), float(p75), values.mean(axis=0, skipna=True),
                          float(flat.min()) - 1.0, int(values.shape[0]))


def signed_scores(
    line: str,
    rule: str,
    modality: str,
    genes: pd.DataFrame,
    rna: pd.DataFrame | None,
    rna_stats: ThresholdStats | None,
    prot: pd.DataFrame | None = None,
    prot_stats: ThresholdStats | None = None,
) -> pd.Series:
    """Signed gene scores ``s_g = value - threshold`` for one line (index: all model genes).

    NaN = no evidence. Protein: not-detected soluble enzymes get ``floor_value - threshold`` (< 0);
    non-soluble genes get NaN; protein R1 maps detected -> +1, not detected (soluble) -> -1.
    """
    out = pd.Series(np.nan, index=genes.index, dtype="float64")
    s_rna = s_prot = None
    if modality in ("rna", "rna_protein"):
        v = rna.loc[line].reindex(genes.index)
        s_rna = v - rna_stats.thresholds(rule, genes.index)
    if modality in ("protein", "rna_protein"):
        if prot is None or line not in prot.index:
            raise KeyError(f"{line}: no proteomics")
        v = prot.loc[line].reindex(genes.index)
        quantified = genes.index.isin(prot.columns)
        detected = v.notna()
        absent_ok = quantified & ~detected.to_numpy() & genes.soluble.to_numpy()
        if rule == "R1":
            s_prot = pd.Series(np.where(detected, 1.0, np.where(absent_ok, -1.0, np.nan)), index=genes.index)
        else:
            thr = prot_stats.thresholds(rule, genes.index)
            vv = v.where(detected, np.where(absent_ok, prot_stats.floor_value, np.nan))
            s_prot = vv - thr
        # transporters / membrane proteins: protein evidence is not used at all (RNA only)
        s_prot[~genes.soluble.to_numpy()] = np.nan
    if modality == "rna":
        out = s_rna
    elif modality == "protein":
        out = s_prot
    else:
        out = s_prot.where(s_prot.notna(), s_rna)
    return out.astype("float64")


# ----------------------------------------------------------------------------- genomic loss rule
def lost_genes(
    line: str,
    dmg: pd.DataFrame | None,
    cn: pd.DataFrame | None,
    rna: pd.DataFrame | None,
    cn_del: float = CN_DEL,
    rule: str = "r2-11",
) -> tuple[list[str], dict]:
    """Model genes lost in ``line`` under the loss rule (DESIGN W4, PI decision r2-11).

    r2-11: damaging code 2; or code 1 and (relCN < 0.5 [LOH proxy] or log2(TPM+1) < 1); or relCN < cn_del.
    v1 (sensitivity): any damaging allele (code >= 1) or relCN < cn_del.
    Lines without CN use mutations only (flag ``cn_missing``).
    """
    flags = {"cn_missing": cn is None or line not in cn.index, "mut_missing": dmg is None or line not in dmg.index}
    d = dmg.loc[line] if not flags["mut_missing"] else pd.Series(dtype=float)
    c = cn.loc[line] if not flags["cn_missing"] else pd.Series(dtype=float)
    e = rna.loc[line] if (rna is not None and line in rna.index) else pd.Series(dtype=float)
    idx = d.index.union(c.index)
    d, c, e = d.reindex(idx), c.reindex(idx), e.reindex(idx)
    if rule == "r2-11":
        lost = (d >= 2) | ((d == 1) & ((c < CN_LOH_PROXY) | (e < EXPR_SILENT))) | (c < cn_del)
    elif rule == "v1":
        lost = (d >= 1) | (c < cn_del)
    else:
        raise ValueError(rule)
    return sorted(lost[lost.fillna(False)].index), flags


# ----------------------------------------------------------------------------- panel and permutation
def stage_a_panel(split: pd.DataFrame, seed: int = SEEDS[0]) -> pd.DataFrame:
    """Stage-A panel: 10 eligible training lineages x 3 training lines with expression (DESIGN W4).

    Implementation of the frozen prose ("10 training lineages x 3 lines, sampled with SEEDS[0]; if a
    lineage has fewer lines, the gap is filled from the largest training lineages"):
    ``rng = default_rng(seed)``; the eligible training lineages (>= 15 screened lines, i.e. not the
    rare-line pool) sorted by name are permuted and the first 10 taken; within each, 3 lines are drawn
    without replacement from its expression-profiled training lines sorted by ModelID. A shortfall is
    filled from the largest training lineages (size desc, name asc), drawing among their unused lines.
    Only training lines are ever drawn (asserted).
    """
    tr = split[(split.tier == "train") & split.has_expression.astype(bool)]
    sizes = split[split.tier == "train"].lineage.value_counts()
    eligible = sorted(L for L, n in sizes.items() if n >= 15)
    rng = np.random.default_rng(seed)
    chosen = [eligible[i] for i in rng.permutation(len(eligible))[:PANEL_LINEAGES]]
    rows, short = [], 0
    for L in chosen:
        pool = sorted(tr[tr.lineage == L].ModelID)
        k = min(PANEL_LINES_PER_LINEAGE, len(pool))
        pick = [pool[i] for i in sorted(rng.choice(len(pool), k, replace=False))] if k else []
        short += PANEL_LINES_PER_LINEAGE - k
        rows += [(m, L, "drawn") for m in pick]
    if short:
        used = {r[0] for r in rows}
        for L in sorted(eligible, key=lambda x: (-sizes[x], x)):
            pool = sorted(set(tr[tr.lineage == L].ModelID) - used)
            take = min(short, len(pool))
            if take:
                pick = [pool[i] for i in sorted(rng.choice(len(pool), take, replace=False))]
                rows += [(m, L, "fill") for m in pick]
                short -= take
            if not short:
                break
    panel = pd.DataFrame(rows, columns=["ModelID", "lineage", "how"])
    meta = split.set_index("ModelID")
    assert (meta.loc[panel.ModelID, "tier"] == "train").all()
    for c in ("has_cn", "has_procan", "base_medium_class"):
        panel[c] = meta.loc[panel.ModelID, c].to_numpy()
    return panel


def permute_across_lines(values: pd.DataFrame, seed: int = SEEDS[0]) -> pd.DataFrame:
    """Permuted-expression control: each gene's values permuted independently across the lines."""
    rng = np.random.default_rng(seed)
    arr = values.to_numpy().copy()
    for j in range(arr.shape[1]):
        arr[:, j] = arr[rng.permutation(arr.shape[0]), j]
    return pd.DataFrame(arr, index=values.index, columns=values.columns)


@dataclass
class LineInputs:
    """Everything one line build needs (picklable; sent to workers)."""

    line: str
    scores: pd.Series  # signed scores s_g (all model genes; NaN = no evidence)
    lost: list = field(default_factory=list)  # genes lost under the loss rule (pre-drop)
    expr_for_drop: pd.Series | None = None  # log2(TPM+1) used to order joint-drop (desc)
    flags: dict = field(default_factory=dict)


# ----------------------------------------------------------------------------- A5: F5 LCP scores (amendment A4 d)
LCP_GRID = [(1.0, w_p, w_c, tau) for w_p in (0.0, 0.5, 1.0) for w_c in (0.0, 1.0) for tau in (-1.5, -1.0, -0.5, 0.0)]


@dataclass
class ZStats:
    """Per-gene training-line statistics for the A5 z-scores (fitting lines only, like R3)."""

    rna_mean: pd.Series
    rna_sd: pd.Series
    prot_mean: pd.Series  # over training lines with ProCan, "not detected" set to the protein floor value
    prot_sd: pd.Series
    prot_floor: float


def _z(v: pd.Series, mean: pd.Series, sd: pd.Series) -> pd.Series:
    sd = sd.reindex(v.index)
    z = (v - mean.reindex(v.index)) / sd.where(sd > 0)
    return z.fillna(0.0)  # no column, zero training variance, or no value: no evidence (z = 0)


def lcp_scores(line: str, w_r: float, w_p: float, w_c: float, tau: float, genes: pd.DataFrame,
               rna: pd.DataFrame, prot: pd.DataFrame | None, zs: ZStats, accepted_losses) -> pd.Series:
    """Signed F5 score ``s_g = w_r z_RNA + w_p z_prot + w_c loss - tau``; the gene is off iff s_g < 0.

    * z_RNA = (log2(TPM+1) - training mean) / training SD (per gene).
    * z_prot = (log2 ProCan - training mean) / training SD per gene, with "not detected" set to the
      protein floor value in both the line and the training distribution; only soluble enzymes
      quantified in ProCan (README §5 item 7); 0 for other genes and for lines without proteomics.
    * loss = -1 for genes in the line's accepted-loss set (loss rule after the order-independent
      drop), else 0, so the term pushes lost genes towards "off".
    """
    idx = genes.index
    z_rna = _z(rna.loc[line].reindex(idx), zs.rna_mean, zs.rna_sd) if line in rna.index else pd.Series(0.0, index=idx)
    z_prot = pd.Series(0.0, index=idx)
    if prot is not None and line in prot.index:
        quant = idx.isin(prot.columns) & genes.soluble.to_numpy()
        v = prot.loc[line].reindex(idx).fillna(zs.prot_floor)
        z_prot = _z(v, zs.prot_mean, zs.prot_sd).where(quant, 0.0)
    loss = pd.Series(np.where(idx.isin(list(accepted_losses)), -1.0, 0.0), index=idx)
    return (w_r * z_rna + w_p * z_prot + w_c * loss - tau).astype("float64")
