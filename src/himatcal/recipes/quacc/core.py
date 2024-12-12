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
    """Perform a relaxation job on a given set of atoms using a specified calculator.

    This function executes an optimization process to relax the atomic structure,
    optionally allowing for cell relaxation. It merges default optimization parameters
    with any user-provided parameters and returns the results of the optimization.

    Args:
        atoms (Atoms): The atomic structure to be relaxed.
        calc: The calculator to be used for the relaxation process.
        relax_cell (bool, optional): Whether to relax the cell dimensions. Defaults to False.
        opt_params (OptParams | None, optional): Additional optimization parameters. Defaults to None.
        additional_fields (dict[str, Any] | None, optional): Extra fields to include in the results. Defaults to None.

    Returns:
        OptSchema: The results of the relaxation job, encapsulated in an optimization schema.

    Examples:
        result = relax_job(atoms, calculator, relax_cell=True, opt_params={"fmax": 0.01})
    """

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
