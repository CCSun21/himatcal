"""https://github.com/kumaranu/ts-workflow-examples/tree/main"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any

from quacc import strip_decorator

from himatcal.recipes.newtonnet.ts import geodesic_job, irc_job, ts_job

if TYPE_CHECKING:
    from ase.atoms import Atoms


def geodesic_ts_hess_irc_newtonnet(
    reactant: Atoms,
    product: Atoms,
    calc_kwargs1: dict[str, Any] | None = None,
    calc_kwargs2: dict[str, Any] | None = None,
    clean_up: bool = True,
) -> list[dict[str, Any]]:
    """
    Perform geodesic, transition state, and intrinsic reaction coordinate (IRC) calculations using NewtonNet.

    Parameters
    ----------
    reactant : Atoms
        The reactant structure.
    product : Atoms
        The product structure.
    calc_kwargs1 : dict
        Keyword arguments for the ASE calculator for the geodesic and IRC jobs.
    calc_kwargs2 : dict
        Keyword arguments for the ASE calculator for the TS job with custom Hessian.
    clean_up : bool, optional
        Whether to clean up raw files after completion, by default True.

    Returns
    -------
    List[Dict[str, Any]]
        List containing results from geodesic, TS, and IRC jobs.
    """
    # Set default values for calc_kwargs1 and calc_kwargs2 if they are None
    if calc_kwargs1 is None:
        calc_kwargs1 = {"hess_method": None}
    if calc_kwargs2 is None:
        calc_kwargs2 = {"hess_method": "autograd"}

    # Create NEB job
    job1 = strip_decorator(geodesic_job)(reactant, product, calc_kwargs=calc_kwargs1)
    logging.info("Created Geodesic job.")

    # Create TS job with custom Hessian
    job2 = strip_decorator(ts_job)(
        job1["highest_e_atoms"], use_custom_hessian=True, **calc_kwargs2
    )
    logging.info("Created TS job with custom Hessian.")

    # Create IRC job in forward direction
    job3 = strip_decorator(irc_job)(job2["atoms"], direction="forward", **calc_kwargs1)
    logging.info("Created IRC job in forward direction.")

    # Create IRC job in reverse direction
    job4 = strip_decorator(irc_job)(job2["atoms"], direction="reverse", **calc_kwargs1)
    logging.info("Created IRC job in reverse direction.")

    logging.info("All jobs executed successfully.")

    if clean_up:
        # Delete the raw files
        directory_patterns = ["quacc-*", "tmp*"]

        # Delete directories matching patterns
        for pattern in directory_patterns:
            for dir_path in Path.cwd().glob(pattern):
                if dir_path.is_dir():
                    shutil.rmtree(dir_path)

    return [job1, job2, job3, job4]
