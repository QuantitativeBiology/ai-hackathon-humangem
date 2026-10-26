"""Map Human-GEM genes (Ensembl) to Entrez IDs and symbols using the model's own ``genes.tsv``.

Nothing is dropped: every model gene gets one row with a ``map_status``:
  one_to_one       one Entrez ID listed
  one_to_many      several Entrez IDs listed; ``entrez`` is the first (primary) one, all in ``entrez_all``
  many_to_one      the Entrez ID is shared with another model gene
  unmapped         no Entrez ID in genes.tsv
  not_in_table     gene is in the SBML but missing from genes.tsv
"""

import re
from pathlib import Path

import cobra
import pandas as pd

DEPMAP_COLUMN = re.compile(r"^(?P<symbol>.+) \((?P<entrez>\d+)\)$")  # DepMap gene columns: "SYMBOL (Entrez)"


def read_genes_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).set_index("genes")


def parse_depmap_columns(columns) -> pd.DataFrame:
    """DepMap ``SYMBOL (Entrez)`` column labels -> DataFrame(symbol, entrez, column)."""
    rows = [(m["symbol"], m["entrez"], c) for c in columns if (m := DEPMAP_COLUMN.match(c))]
    return pd.DataFrame(rows, columns=["symbol", "entrez", "column"])


def gene_mapping(model: cobra.Model, genes_tsv: Path, depmap_columns=None) -> pd.DataFrame:
    """One row per model gene: gene, symbol, entrez, entrez_all, n_entrez, map_status, n_reactions.

    If ``depmap_columns`` (gene column labels only, no values) is given, adds
    ``depmap_column`` (label matched on the primary Entrez ID), ``depmap_symbol_differs``
    and, for unmapped genes, ``depmap_symbol_candidate`` (label matched on symbol; not
    used as a mapping, reported for review only).
    """
    table = read_genes_tsv(genes_tsv)
    df = pd.DataFrame({"gene": sorted(g.id for g in model.genes)})
    df["n_reactions"] = df.gene.map({g.id: len(g.reactions) for g in model.genes})
    df["in_table"] = df.gene.isin(table.index)
    df["symbol"] = df.gene.map(table["geneSymbols"])
    df["entrez_all"] = df.gene.map(table["geneEntrezID"])
    ids = df.entrez_all.str.split(";")
    df["n_entrez"] = ids.str.len().fillna(0).astype(int)
    df["entrez"] = ids.str[0]

    shared = df.entrez.notna() & df.entrez.duplicated(keep=False)
    df["map_status"] = "one_to_one"
    df.loc[df.n_entrez > 1, "map_status"] = "one_to_many"
    df.loc[shared, "map_status"] = "many_to_one"
    df.loc[df.n_entrez == 0, "map_status"] = "unmapped"
    df.loc[~df.in_table, "map_status"] = "not_in_table"

    if depmap_columns is not None:
        dm = parse_depmap_columns(depmap_columns)
        by_entrez = dm.set_index("entrez")
        df["depmap_column"] = df.entrez.map(by_entrez["column"])
        df["depmap_symbol_differs"] = df.depmap_column.notna() & (
            df.entrez.map(by_entrez["symbol"]) != df.symbol
        )
        by_symbol = dm.drop_duplicates("symbol", keep=False).set_index("symbol")
        df["depmap_symbol_candidate"] = df.symbol.map(by_symbol["column"]).where(df.entrez.isna())
    return df.drop(columns="in_table")


def mapping_summary(mapping: pd.DataFrame) -> dict[str, int]:
    summary = {"genes": len(mapping), **mapping.map_status.value_counts().to_dict()}
    if "depmap_column" in mapping:
        summary["entrez_in_depmap_crispr"] = int(mapping.depmap_column.notna().sum())
        summary["depmap_symbol_differs"] = int(mapping.depmap_symbol_differs.sum())
    return summary
