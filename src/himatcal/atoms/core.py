from __future__ import annotations

import contextlib
from pathlib import Path
from typing import TYPE_CHECKING

from ase import Atoms
from ase.io import write

from himatcal.recipes.crest.core import relax, iMTD_GC

if TYPE_CHECKING:
    from typing import Literal

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

Li = Atoms(
    symbols="Li",
    positions=[[0, 0, 0]],
)

Na = Atoms(
    symbols="Na",
    positions=[[0, 0, 0]],
)

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
    ship_atoms: Atoms,
    dock_atoms: Atoms | Literal["PF6", "Li", "Na"] = "PF6",
    offset: float = 1.5,
    crest_sampling: bool = True,
    chg: int = 0,
    mult: int = 1,
):
    # TODO: this dock function is not general enough, it should be able to dock any two atoms
    """
    Dock the ship🚢 atoms to the dock⚓ atoms (default is PF6).

    Parameters:
    -----------
    ship_atoms (ase.Atoms): The ship atoms.
    """
    dock_atoms_dict = {"PF6": PF6.copy(), "Li": Li.copy(), "Na": Na.copy()}

    if isinstance(dock_atoms, str):
        dock_atoms = dock_atoms_dict.get(dock_atoms, dock_atoms)

    ship_atoms = ship_atoms.copy()
    ship_atoms_center = ship_atoms.get_center_of_mass()
    ship_atoms_center[0] = max(ship_atoms.positions.T[0])
    dock_atoms_center = dock_atoms.get_center_of_mass()
    dock_atoms_center[0] = min(dock_atoms.positions.T[0])
    vector = ship_atoms_center - dock_atoms_center
    offset = [offset, 0, 0]
    dock_atoms.positions = dock_atoms.positions + vector + offset
    ship_atoms.extend(dock_atoms)
    if crest_sampling:
        for _ in range(3):
            with contextlib.suppress(Exception):
                ship_atoms = iMTD_GC(ship_atoms, chg=chg, mult=mult)
                break
    return ship_atoms


def tmp_atoms(atoms, filename="tmp.xyz", create_tmp_folder=True):
    """
    Write the atoms to a temporary file in the tmp folder and return the path.

    Args:

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
