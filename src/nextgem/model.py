"""Load Human-GEM (Human1 = v1.19.0, Human2 = v2.1.1) and constrain it to a defined growth medium."""

from dataclasses import dataclass
from pathlib import Path

import cobra
from raven_toolbox.tasks.tasklist import parse_task_list

PROJECT_ROOT = Path(__file__).resolve().parents[2]
HUMAN_GEM_DIR = PROJECT_ROOT / "data" / "Human-GEM-main"  # Human2, v2.1.1
HUMAN1_DIR = PROJECT_ROOT / "data" / "Human1"  # last v1.x release, see scripts/download_human1.py
MODEL_FILE = HUMAN_GEM_DIR / "model" / "Human-GEM.xml"
# Human2 task file. Functions never default to it: they derive the file from the model release (tasks_file_for).
ESSENTIAL_TASKS_FILE = HUMAN_GEM_DIR / "data" / "metabolicTasks" / "metabolicTasks_Essential.txt"

BIOMASS_REACTION = "MAR13082"  # Generic human cell biomass reaction (same id in Human1 and Human2)
GROWTH_TASK_ID = "GR"  # "Growth on Ham's media (biomass production)"


@dataclass(frozen=True)
class GEM:
    """Files of one Human-GEM release; both releases share the repo layout."""

    name: str
    version: str
    root: Path

    @property
    def sbml(self) -> Path:
        return self.root / "model" / "Human-GEM.xml"

    @property
    def genes_tsv(self) -> Path:
        return self.root / "model" / "genes.tsv"

    @property
    def essential_tasks(self) -> Path:
        return self.root / "data" / "metabolicTasks" / "metabolicTasks_Essential.txt"


GEMS = {
    "human1": GEM("human1", "v1.19.0", HUMAN1_DIR),
    "human2": GEM("human2", "v2.1.1", HUMAN_GEM_DIR),
}


RELEASE_NOTE = "nextgem_release"  # model.notes key naming the GEMS entry a model was loaded from


def load_model(path: Path | str = MODEL_FILE) -> cobra.Model:
    """Read Human-GEM from SBML (~15 s); load once and use ``model.copy()`` or ``with model:``.

    ``path`` is an SBML file or a key of ``GEMS`` ("human1", "human2"). Models loaded from a
    ``GEMS`` release are tagged (``model.notes[RELEASE_NOTE]``, kept by copy and pickle) so that
    release-specific files such as the task list are derived from the model, not assumed.
    """
    if isinstance(path, str) and path in GEMS:
        path = GEMS[path].sbml
    model = cobra.io.read_sbml_model(str(path))
    release = next((g.name for g in GEMS.values() if g.sbml.resolve() == Path(path).resolve()), None)
    if release is not None:
        model.notes[RELEASE_NOTE] = release
    return model


def gem_of(model: cobra.Model) -> GEM:
    """The ``GEMS`` release a model was loaded from; raises if the model is untagged."""
    release = model.notes.get(RELEASE_NOTE)
    if release not in GEMS:
        raise ValueError(
            f"Model {model.id!r} has no release tag (notes[{RELEASE_NOTE!r}]); load it with load_model() "
            "from a GEMS release or pass tasks_file explicitly"
        )
    return GEMS[release]


def tasks_file_for(model: cobra.Model, tasks_file: Path | None = None) -> Path:
    """``tasks_file`` if given, else the essential-task file of the model's own release."""
    return Path(tasks_file) if tasks_file is not None else gem_of(model).essential_tasks


def task_input_exchanges(
    model: cobra.Model, task_id: str = GROWTH_TASK_ID, tasks_file: Path | None = None
) -> set[str]:
    """Exchange reaction ids for the inputs of a metabolic task (default ``GR``, Ham's medium).

    Task inputs are ``name[compartment]`` tokens resolved by metabolite name and
    compartment, as in ``code/test/gradedEssentiality.py::media_exchange_ids``.
    ``tasks_file`` defaults to the task file of the model's own release (``tasks_file_for``).
    """
    tasks_file = tasks_file_for(model, tasks_file)
    task = next(t for t in parse_task_list(str(tasks_file)) if t.id == task_id)
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
        raise ValueError(f"Unresolved inputs of task {task_id}: {unresolved}")
    return exchange_ids & {rxn.id for rxn in model.exchanges}


def hams_medium_exchanges(model: cobra.Model, tasks_file: Path | None = None) -> set[str]:
    """Exchange reaction ids for the inputs of Human-GEM's ``GR`` task (Ham's medium)."""
    return task_input_exchanges(model, GROWTH_TASK_ID, tasks_file)


def hams_medium(
    model: cobra.Model, uptake_bound: float = 1000.0, tasks_file: Path | None = None
) -> dict[str, float]:
    """Ham's medium as a cobra medium dict (exchange id -> maximal uptake)."""
    return {rxn_id: uptake_bound for rxn_id in sorted(hams_medium_exchanges(model, tasks_file))}


def set_hams_medium(
    model: cobra.Model, uptake_bound: float = 1000.0, tasks_file: Path | None = None
) -> None:
    """Close all uptakes except Ham's medium components (in place).

    With the default bound of 1000, uptake is effectively unconstrained, so growth is
    in arbitrary units; pass measured rates per exchange for quantitative predictions.
    """
    model.medium = hams_medium(model, uptake_bound, tasks_file)
