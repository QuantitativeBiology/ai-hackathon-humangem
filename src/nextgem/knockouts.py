"""Wild-type growth and GPR-aware single-gene knockouts for a given medium and objective.

A knockout marks the gene non-functional and closes every reaction whose GPR then
evaluates to False (``cobra.Gene.knock_out``): isozymes (OR) buffer, complex
subunits (AND) do not. Growth is the optimal objective value.

Status policy (nothing non-optimal silently becomes "non-essential"):
- ``infeasible`` KO LP: recorded as growth 0 and ratio 0, status kept as ``infeasible``.
- any other non-optimal KO status (time limit, numeric, unbounded, ...): ``NonOptimalError``
  (``strict=False`` keeps such rows with NaN growth and their status instead).
- wild type infeasible or not growing: ratios NaN, status ``wt_no_growth`` (with a warning);
  any other non-optimal wild-type status: ``NonOptimalError``.

``ratio = ko_growth / wt_growth`` is the raw value; ``ratio_r`` is ratio rounded to 6 decimals
and clipped to [0, 1]. Raw ratios carry solver noise (~1e-13) that varies between runs and
creates spurious ranks among ties at 1, so all downstream scoring uses ``ratio_r``.

Parallel runs use ``spawn`` worker processes that each unpickle the configured model
once, so each worker holds one solver environment (one Gurobi WLS session) and solves
single-threaded.
"""

import math
import pickle
import time
import warnings
from collections.abc import Iterable, Mapping
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

import cobra
import pandas as pd
from optlang.interface import INFEASIBLE, OPTIMAL

from nextgem.model import BIOMASS_REACTION, hams_medium

COLUMNS = ["gene", "wt_growth", "ko_growth", "ratio", "ratio_r", "status", "solver", "runtime"]
WT_GROWTH_TOL = 1e-6  # wild-type growth at or below this is "no growth"
RATIO_DECIMALS = 6
ACCEPTED_STATUSES = {OPTIMAL, INFEASIBLE}


class NonOptimalError(RuntimeError):
    """A solve ended with a status other than optimal or infeasible."""


def rounded_ratio(ratio: pd.Series) -> pd.Series:
    """Ratio rounded to ``RATIO_DECIMALS`` and clipped to [0, 1] (``+ 0.0`` turns -0.0 into 0.0)."""
    return ratio.round(RATIO_DECIMALS).clip(0.0, 1.0) + 0.0

_WORKER_MODEL: cobra.Model | None = None


def configure(
    model: cobra.Model,
    medium: Mapping[str, float] | None = None,
    objective: str | dict = BIOMASS_REACTION,
    tasks_file: Path | None = None,
) -> None:
    """Set medium (default: Ham's medium from the ``GR`` task) and objective in place.

    ``medium`` is a cobra medium dict {exchange id: max uptake}; every other uptake is closed.
    ``tasks_file`` defaults to the task file of the model's own release (``model.tasks_file_for``).
    """
    if medium is None:
        medium = hams_medium(model, tasks_file=tasks_file)
    unknown = sorted(set(medium) - {r.id for r in model.reactions})
    if unknown:
        raise KeyError(f"Medium reactions not in model: {unknown[:10]} ({len(unknown)} total)")
    model.medium = dict(medium)
    model.objective = objective


def solver_name(model: cobra.Model) -> str:
    return model.solver.interface.__name__.split(".")[-1].removesuffix("_interface")


def _solve(model: cobra.Model) -> tuple[float, str]:
    model.slim_optimize(error_value=math.nan)
    status = model.solver.status
    if status == OPTIMAL:
        return model.solver.objective.value, status
    return (0.0 if status == INFEASIBLE else math.nan), status


def _knockout(model: cobra.Model, gene_id: str) -> tuple[str, float, str, float]:
    start = time.perf_counter()
    with model:
        model.genes.get_by_id(gene_id).knock_out()
        growth, status = _solve(model)
    return gene_id, growth, status, time.perf_counter() - start


def limit_threads(model: cobra.Model, threads: int = 1) -> None:
    """Cap solver threads (Gurobi's default concurrent LP uses many cores per solve)."""
    if solver_name(model) == "gurobi":
        model.solver.problem.Params.Threads = threads


def _init_worker(model_bytes: bytes) -> None:
    global _WORKER_MODEL
    _WORKER_MODEL = pickle.loads(model_bytes)
    limit_threads(_WORKER_MODEL)  # solver parameters do not survive pickling; one thread per worker


def _worker_chunk(gene_ids: list[str]) -> list[tuple[str, float, str, float]]:
    return [_knockout(_WORKER_MODEL, g) for g in gene_ids]


def wild_type_growth(model: cobra.Model) -> tuple[float, str]:
    """Optimal objective value and solver status of the (configured) model."""
    return _solve(model)


def check_genes(model: cobra.Model, genes: Iterable[str]) -> list[str]:
    """Return ``genes`` as a list, raising KeyError if any is not a model gene id."""
    genes = list(genes)
    unknown = [g for g in genes if g not in model.genes]
    if unknown:
        raise KeyError(f"Not model genes: {unknown[:10]} ({len(unknown)} total)")
    return genes


def single_gene_knockouts(
    model: cobra.Model,
    genes: Iterable[str] | None = None,
    medium: Mapping[str, float] | None = None,
    objective: str | dict = BIOMASS_REACTION,
    tasks_file: Path | None = None,
    solver: str | None = None,
    processes: int = 1,
    chunksize: int = 25,
    strict: bool = True,
) -> pd.DataFrame:
    """Wild-type growth and KO/WT growth ratio for each gene (default: all model genes).

    ``model`` is left unchanged. ``medium``/``objective``/``tasks_file`` as in ``configure``;
    ``solver`` (e.g. "gurobi", "glpk") overrides the model's solver for this run.
    Returns one row per gene with columns ``COLUMNS``; see the module docstring for the status policy.
    """
    genes = check_genes(model, sorted(g.id for g in model.genes) if genes is None else genes)
    if solver is not None and solver_name(model) != solver:
        model = model.copy()  # solver switches are not reliably undone by ``with model:``
        model.solver = solver
    with model:
        configure(model, medium, objective, tasks_file)
        wt, wt_status = wild_type_growth(model)
        if wt_status not in ACCEPTED_STATUSES:
            raise NonOptimalError(f"Wild-type solve status {wt_status!r}")
        name = solver_name(model)
        if processes > 1 and len(genes) > chunksize:
            chunks = [genes[i : i + chunksize] for i in range(0, len(genes), chunksize)]
            with ProcessPoolExecutor(
                min(processes, len(chunks)),
                mp_context=get_context("spawn"),
                initializer=_init_worker,
                initargs=(pickle.dumps(model),),
            ) as pool:
                rows = [row for chunk in pool.map(_worker_chunk, chunks) for row in chunk]
        else:
            rows = [_knockout(model, g) for g in genes]

    df = pd.DataFrame(rows, columns=["gene", "ko_growth", "status", "runtime"])
    bad = df[~df.status.isin(ACCEPTED_STATUSES)]
    if strict and not bad.empty:
        raise NonOptimalError(
            f"{len(bad)} knockouts ended non-optimal: {bad.status.value_counts().to_dict()}; "
            f"genes {bad.gene.tolist()[:10]}"
        )
    df["wt_growth"] = wt
    df["solver"] = name
    if wt_status == OPTIMAL and wt > WT_GROWTH_TOL:
        df["ratio"] = df.ko_growth / wt  # infeasible KOs have ko_growth 0 -> ratio 0
    else:
        warnings.warn(f"Wild type does not grow (growth={wt}, status={wt_status}); ratios set to NaN")
        df["ratio"] = math.nan
        df["status"] = "wt_no_growth"
    df["ratio_r"] = rounded_ratio(df.ratio)
    return df[COLUMNS]


def essential_genes(result: pd.DataFrame, threshold: float = 0.01, nan: str = "raise") -> set[str]:
    """Genes with rounded KO/WT growth ratio (``ratio_r``) below ``threshold``.

    NaN ratios (wild type not growing, or non-optimal rows kept with ``strict=False``) are never
    silently called non-essential: ``nan="raise"`` (default) raises, ``nan="exclude"`` drops them.
    """
    if nan not in ("raise", "exclude"):
        raise ValueError(f"nan must be 'raise' or 'exclude', got {nan!r}")
    n_nan = int(result.ratio_r.isna().sum())
    if n_nan and nan == "raise":
        raise ValueError(f"{n_nan} NaN ratios (statuses {result.status[result.ratio_r.isna()].unique().tolist()})")
    return set(result.gene[result.ratio_r < threshold])
