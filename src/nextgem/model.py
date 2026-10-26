"""Load Human-GEM and constrain it to a defined growth medium."""

from pathlib import Path

import cobra
from raven_toolbox.tasks.tasklist import parse_task_list

PROJECT_ROOT = Path(__file__).resolve().parents[2]
HUMAN_GEM_DIR = PROJECT_ROOT / "data" / "Human-GEM-main"
MODEL_FILE = HUMAN_GEM_DIR / "model" / "Human-GEM.xml"
ESSENTIAL_TASKS_FILE = HUMAN_GEM_DIR / "data" / "metabolicTasks" / "metabolicTasks_Essential.txt"

BIOMASS_REACTION = "MAR13082"  # Generic human cell biomass reaction
GROWTH_TASK_ID = "GR"  # "Growth on Ham's media (biomass production)"


def load_model(path: Path = MODEL_FILE) -> cobra.Model:
    """Read Human-GEM from SBML (~15 s); load once and use ``model.copy()`` or ``with model:``."""
    return cobra.io.read_sbml_model(str(path))


def hams_medium_exchanges(model: cobra.Model) -> set[str]:
    """Exchange reaction ids for the inputs of Human-GEM's ``GR`` task (Ham's medium).

    Task inputs are ``name[compartment]`` tokens resolved by metabolite name and
    compartment, as in ``code/test/gradedEssentiality.py::media_exchange_ids``.
    """
    task = next(t for t in parse_task_list(str(ESSENTIAL_TASKS_FILE)) if t.id == GROWTH_TASK_ID)
    exchange_ids: set[str] = set()
    unresolved: list[str] = []
    for token, _lb, _ub in task.inputs:
        name, _, comp = token.rpartition("[") if "[" in token else (token, "", "")
        name, comp = name.strip().lower(), comp.rstrip("]") or None
        matched = {
            rxn.id
            for met in model.metabolites
            if met.name and met.name.lower() == name and (comp is None or met.compartment == comp)
            for rxn in met.reactions
            if rxn.boundary
        }
        if not matched:
            unresolved.append(token)
        exchange_ids |= matched
    if unresolved:
        raise ValueError(f"Unresolved Ham's medium inputs: {unresolved}")
    return exchange_ids & {rxn.id for rxn in model.exchanges}


def set_hams_medium(model: cobra.Model, uptake_bound: float = 1000.0) -> None:
    """Close all uptakes except Ham's medium components (in place).

    With the default bound of 1000, uptake is effectively unconstrained, so growth is
    in arbitrary units; pass measured rates per exchange for quantitative predictions.
    """
    model.medium = {rxn_id: uptake_bound for rxn_id in hams_medium_exchanges(model)}
