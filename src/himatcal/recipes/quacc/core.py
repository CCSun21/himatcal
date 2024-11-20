from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

from quacc import job
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from quacc.types import OptParams, OptSchema


@job
def relax_job(
    atoms: Atoms,
    calc,
    relax_cell: bool = False,
    opt_params: OptParams | None = None,
    additional_fields: dict[str, Any] | None = None,
) -> OptSchema:
    opt_defaults = {"fmax": 0.05, "max_steps": 1000}

    # Ensure opt_params is converted to a dictionary if it's not already
    opt_flags = recursive_dict_merge(
        opt_defaults, dict(opt_params) if opt_params else {}
    )

    # Make sure that the 'run_opt' method returns an appropriate type for 'opt'
    dyn = Runner(atoms, calc).run_opt(relax_cell=relax_cell, **opt_flags)

    # Ensure dyn is of the correct type for Summarize.opt
    return Summarize(
        additional_fields={"name": "MLP Relax"} | (additional_fields or {})
    ).opt(dyn)


def load_quacc_result(file_path: str) -> dict:
    with Path.open(Path(file_path)) as file:
        return json.load(file)
