"""Sanity-check FBA on Human-GEM under Ham's medium: growth, O2/glucose dependence, known KOs.

Usage: .venv/bin/python scripts/fba_sanity_check.py
"""

import pandas as pd

from nextgem import BIOMASS_REACTION, load_model, set_hams_medium
from nextgem.model import HUMAN_GEM_DIR

# Expected under Ham's medium (thymidine but no uridine supplied): CAD/UMPS lethal, rest viable.
KO_GENES = ["CAD", "UMPS", "TYMS", "DHFR", "GAPDH", "PGK1", "HPRT1", "LDHA"]


def exchange_by_name(model, name):
    return next(r for r in model.exchanges if next(iter(r.metabolites)).name == name and r.lower_bound < 0)


def main():
    model = load_model()
    set_hams_medium(model)
    model.objective = BIOMASS_REACTION
    print(f"Solver: {model.solver.interface.__name__} | medium exchanges: {len(model.medium)}")

    wt = model.slim_optimize(error_value=0)
    print(f"Wild-type growth: {wt:.4f}")

    for label, rxn in [("no glucose", exchange_by_name(model, "glucose")), ("anaerobic", exchange_by_name(model, "O2"))]:
        with model:
            rxn.lower_bound = 0
            print(f"Growth, {label}: {model.slim_optimize(error_value=0):.4f}")

    symbols = pd.read_csv(HUMAN_GEM_DIR / "model" / "genes.tsv", sep="\t", index_col=0)["geneSymbols"]
    for symbol in KO_GENES:
        with model:
            model.genes.get_by_id(symbols[symbols == symbol].index[0]).knock_out()
            print(f"KO {symbol:6s}: growth ratio = {model.slim_optimize(error_value=0) / wt:.3f}")


if __name__ == "__main__":
    main()
