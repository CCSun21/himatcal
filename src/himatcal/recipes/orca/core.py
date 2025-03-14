"""Core recipes for ORCA."""

from __future__ import annotations

from typing import TYPE_CHECKING

import psutil
from quacc import job
from quacc.recipes.orca._base import run_and_summarize

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms
    from quacc.types import Filenames, RunSchema, SourceDirectory


@job
def neb_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    n_images: int = 6,
    preopt: bool = False,
    guessTS: bool = False,
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    nprocs: int | Literal["max"] | None = "max",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
) -> RunSchema:
    """
    product.xyz is necessary for NEB calculation, guessts.xyz is optional; Add a folder to copy_files to copy the product.xyz and guessts.xyz to the working directory.
    """
    additional_fields = {"name": "ORCA Static"} | (additional_fields or {})
    nprocs = psutil.cpu_count(logical=False) if nprocs in ["max", None] else nprocs
    default_inputs = [xc, basis, "neb-ts", "freq"]
    neb_block = '%NEB NEB_END_XYZFILE "product.xyz"'
    if guessTS:
        neb_block += ' NEB_TS_XYZFILE "guessts.xyz"'
    if preopt:
        neb_block += " PREOPT_ENDS TRUE"
    neb_block += f" Nimages {n_images}"
    neb_block += " END"
    default_blocks = [f"%pal nprocs {nprocs} end", neb_block]

    return run_and_summarize(
        atoms,
        charge,
        spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=orcasimpleinput,
        block_swaps=orcablocks,
        additional_fields=additional_fields,
        copy_files=copy_files,
    )


@job
def bare_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "b97-3c",
    basis: str = "def2-tzvp",
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    nprocs: int | Literal["max"] | None = "max",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
) -> RunSchema:
    additional_fields = {"name": "ORCA Static"} | (additional_fields or {})
    nprocs = psutil.cpu_count(logical=False) if nprocs in ["max", None] else nprocs
    default_inputs = [xc, basis]
    default_blocks = [f"%pal nprocs {nprocs} end"]

    return run_and_summarize(
        atoms,
        charge,
        spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=orcasimpleinput,
        block_swaps=orcablocks,
        additional_fields=additional_fields,
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "b97-3c",
    basis: str = "def2-tzvp",
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    nprocs: int | Literal["max"] | None = "max",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
) -> RunSchema:
    additional_fields = {"name": "ORCA Relax"} | (additional_fields or {})
    nprocs = psutil.cpu_count(logical=False) if nprocs in ["max", None] else nprocs
    default_inputs = [xc, basis, "opt", "freq"]
    default_blocks = [f"%pal nprocs {nprocs} end"]

    return run_and_summarize(
        atoms,
        charge,
        spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=orcasimpleinput,
        block_swaps=orcablocks,
        additional_fields=additional_fields,
        copy_files=copy_files,
    )
