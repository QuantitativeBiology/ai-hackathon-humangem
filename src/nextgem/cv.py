"""LOLO cross-validation harness (DESIGN.md §4.2, §11.1). All CV goes through :func:`lolo`.

Fitting code receives only ``labels.fold(k)`` (training minus the held-out lineage; for nested LOLO,
minus an inner lineage too) and the IDs of the held-out lines. The held-out labels never reach it:
the harness scores the predictions itself with a caller-supplied pure ``score`` function.
"""

from __future__ import annotations

from typing import Callable

import pandas as pd

from nextgem import labels


def lolo(fit_predict: Callable[[pd.DataFrame, list[str]], pd.DataFrame],
         score: Callable[[pd.DataFrame, pd.DataFrame], object],
         *, columns: list[str] | None = None, outer=None, folds: list[str] | None = None) -> dict:
    """Leave-one-lineage-out over the training folds.

    ``fit_predict(fit_labels, heldout_model_ids)`` -> predictions (ModelID, entrez, score).
    ``score(predictions, heldout_labels)`` -> any result; it runs inside the harness.
    ``outer``: for nested (inner) LOLO within training-minus-``outer``; folds then exclude ``outer``.
    Returns {fold lineage: score result}.
    """
    all_folds = labels.folds()
    if outer is not None:
        outer_name = all_folds[outer] if isinstance(outer, int) else outer
        fl = [f for f in (folds or all_folds) if f != outer_name]
    else:
        outer_name, fl = None, list(folds or all_folds)
    results = {}
    for k in fl:
        if outer_name is None:
            fit = labels.fold(k, columns=columns)
            held = labels._heldout(k, columns=columns)
        else:
            fit = labels.fold(outer_name, inner=k, columns=columns)
            held = labels._heldout(outer_name, inner=k, columns=columns)
        ids = sorted(held.ModelID.unique())
        if set(ids) & set(fit.ModelID):
            raise AssertionError(f"fold {k}: held-out lines leak into the fitting labels")
        preds = fit_predict(fit, ids)
        results[k] = score(preds, held)
    return results
