"""nextgem package. Model/knockout helpers are imported lazily (PEP 562), so light modules such as
``nextgem.labels``, ``nextgem.ledger`` or ``nextgem.paths`` (used by git hooks) do not import cobra."""

_LAZY = {
    "NonOptimalError": "nextgem.knockouts",
    "configure": "nextgem.knockouts",
    "essential_genes": "nextgem.knockouts",
    "single_gene_knockouts": "nextgem.knockouts",
    "BIOMASS_REACTION": "nextgem.model",
    "GEMS": "nextgem.model",
    "gem_of": "nextgem.model",
    "hams_medium": "nextgem.model",
    "hams_medium_exchanges": "nextgem.model",
    "load_model": "nextgem.model",
    "set_hams_medium": "nextgem.model",
}

__all__ = sorted(_LAZY)


def __getattr__(name):
    if name in _LAZY:
        import importlib

        value = getattr(importlib.import_module(_LAZY[name]), name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'nextgem' has no attribute {name!r}")
