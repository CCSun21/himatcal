from __future__ import annotations

from ase import Atoms
from ase.cell import Cell
from pymatgen.core import Structure


@staticmethod
def pymatgen_to_ase(structure):
    lattice = structure.lattice
    """
    Converts a pymatgen Structure to an ASE Atoms object.

    Args:
        structure: pymatgen Structure object to convert.

    Returns:
        ASE Atoms object representing the same structure.
    """

    return Atoms(
        scaled_positions=structure.frac_coords,
        numbers=structure.atomic_numbers,
        pbc=True,
        cell=Cell.fromcellpar(
            [
                lattice.a,
                lattice.b,
                lattice.c,
                lattice.alpha,
                lattice.beta,
                lattice.gamma,
            ]
        ),
    )


@staticmethod
def ase_to_pymatgen(atoms):
    """
    Converts an ASE Atoms object to a pymatgen Structure.

    Args:
        atoms: ASE Atoms object to convert.

    Returns:
        pymatgen Structure object representing the same atomic configuration.
    """

    return Structure(
        coords_are_cartesian=True,
        coords=atoms.positions,
        species=atoms.symbols,
        lattice=atoms.cell,
    )
