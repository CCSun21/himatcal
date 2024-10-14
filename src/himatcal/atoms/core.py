from __future__ import annotations

import contextlib
import logging
from pathlib import Path

from ase import Atoms
from ase.io import write

from himatcal.recipes.crest.core import iMTD_GC

PF6 = Atoms(
    symbols="PF6",
    positions=[
        [6.747, 7.453, 7.469],
        [7.944, 6.319, 7.953],
        [6.127, 6.381, 6.461],
        [5.794, 8.645, 7.001],
        [5.815, 7.032, 8.699],
        [7.617, 8.534, 8.484],
        [7.91, 7.908, 6.284],
    ],
)

Li = Atoms(symbols="Li", positions=[[0, 0, 0]])

Na = Atoms(symbols="Na", positions=[[0, 0, 0]])

elements = [
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Uub",
    "Uut",
    "Uuq",
    "Uup",
    "Uuh",
    "Uus",
    "Uuo",
]


def dock_atoms(
    ship: Atoms,
    dock: Atoms | str = "PF6",
    chg: int = 0,
    mult: int = 1,
    offset: float = 1.5,
    crest_sampling: bool = True,
    topo_change: bool = False,
):
    """
    Dock the ship🚢 atoms to the dock⚓ atoms (default is PF6).

    Parameters:
    -----------
    ship_atoms (ase.Atoms): The ship atoms.
    """
    dock_atoms_dict = {"PF6": PF6.copy(), "Li": Li.copy(), "Na": Na.copy()}

    if isinstance(dock, str):
        dock = dock_atoms_dict.get(dock, dock)

    docked_atoms = ship.copy()
    ship_atoms_center = docked_atoms.get_center_of_mass()
    ship_atoms_center[0] = max(docked_atoms.positions.T[0])
    dock_atoms_center = dock.get_center_of_mass()
    dock_atoms_center[0] = min(dock.positions.T[0])
    vector = ship_atoms_center - dock_atoms_center + [offset, 0, 0]
    dock.positions += vector
    docked_atoms.extend(dock)
    if not crest_sampling:
        return docked_atoms

    processed_atoms = None
    for _ in range(3):
        logging.info(f"Trying sampling the docked atoms using iMTD-GC the {_} time")
        with contextlib.suppress(Exception):
            processed_atoms = iMTD_GC(
                docked_atoms, chg=chg, mult=mult, topo_change=topo_change
            )
            break
    if processed_atoms is None:
        logging.info("Sampling failed, trying the docked atoms with topology change")
        for _ in range(3):
            with contextlib.suppress(Exception):
                processed_atoms = iMTD_GC(
                    docked_atoms, chg=chg, mult=mult, topo_change=True
                )
                break
        logging.info("Crest iMTD-GC Sampling failed!")
    return processed_atoms


# TODO: dock atoms using orca or xtb


def tmp_atoms(atoms, filename="tmp.xyz", create_tmp_folder=True):
    """
    Write the atoms to a temporary file in the tmp folder and return the path.

    Args:

        create_tmp_folder:
        atoms (ase.Atoms): The atoms object.
        filename (str): The filename of the temporary file.

    Returns:

        filepath (str): The path of the temporary file
    """

    _CWD = Path.cwd()
    if create_tmp_folder:
        from monty.os import makedirs_p

        tmp_path = _CWD / "tmp"
        makedirs_p(tmp_path)
        filepath = _CWD / "tmp" / filename
    else:
        filepath = _CWD / filename
    write(filepath, atoms, format="xyz")
    return filepath


def add_cell(atoms, cell_parameters):
    """
    Add cell parameters to the atoms object.

    Args:

        atoms (ase.Atoms): The atoms object.
        cell_parameters (list): The cell parameters.

    Returns:

        atoms (ase.Atoms): The atoms object with cell parameters.
    """

    atoms.set_cell(cell_parameters)
    return atoms
