from __future__ import annotations

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
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    dyn = Runner(atoms, calc).run_opt(relax_cell=relax_cell, **opt_flags)

    return Summarize(
        additional_fields={"name": "MLP Relax"} | (additional_fields or {})
    ).opt(dyn)
