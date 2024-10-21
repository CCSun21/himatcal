from __future__ import annotations

from typing import TYPE_CHECKING

import psutil
from quacc import job

from himatcal.recipes.gaussian._base import run_and_summarize

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from quacc.types import Filenames, RunSchema, SourceDirectory


@job
def relax_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    label: str = "relax",
    xc: str = "wb97xd",
    basis: str = "def2tzvp",
    freq: bool = False,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Run a Gaussian relaxation calculation.

    Args:
        atoms: Atoms - The atoms for the calculation.
        charge: int - The charge of the system.
        spin_multiplicity: int - The spin multiplicity of the system.
        label: str - (Optional) The label for the calculation (default is "relax").
        xc: str - (Optional) The exchange-correlation functional to use (default is "wb97xd").
        basis: str - (Optional) The basis set to use (default is "def2tzvp").
        freq: bool - (Optional) Whether to perform frequency calculations (default is False).
        copy_files: Union[SourceDirectory, dict[SourceDirectory, Filenames], None] - (Optional) Files to copy after the calculation.
        **calc_kwargs - Additional keyword arguments for the calculation.

    Returns:
        RunSchema - The summarized result of the Gaussian relaxation calculation.
    """

    calc_defaults = {
        "mem": "64GB",
        "chk": "Gaussian.chk",
        "nprocshared": psutil.cpu_count(logical=False),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": spin_multiplicity,
        "opt": "",
        "pop": "CM5",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "ioplist": ["2/9=2000"], # ASE issue #660
    }
    if freq:
        calc_defaults["freq"] = ""

    return run_and_summarize(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        label=label,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "Gaussian Relax"},
        copy_files=copy_files,
    )


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    label: str = "static",
    xc: str = "wb97xd",
    basis: str = "def2tzvp",
    freq: bool = True,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    xc
        Exchange-correlation functional
    basis
        Basis set
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Gaussian calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.gaussian.Gaussian][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results
    """
    calc_defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": psutil.cpu_count(logical=False),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": spin_multiplicity,
        # "force": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "pop": "CM5",
        "gfinput": "",
        "ioplist": ["6/7=3", "2/9=2000"],  # see ASE issue #660
    }
    if freq:
        calc_defaults["freq"] = ""

    return run_and_summarize(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        label=label,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "Gaussian Static"},
        copy_files=copy_files,
    )
