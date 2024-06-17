from __future__ import annotations

from ase import Atoms

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
)  # the PF6 atoms type: ase.Atoms()

Li = Atoms(
    symbols="Li",
    positions=[[0, 0, 0]],
)  # the Li atoms type: ase.Atoms()


def dock_atoms(ship_atoms, dock_atoms=None, offset=1.5):
    # TODO: this dock function is not general enough, it should be able to dock any two atoms
    """
    Dock the ship🚢 atoms to the dock⚓ atoms (default is PF6).

    Parameters:
    -----------
    ship_atoms (ase.Atoms): The ship atoms.
    """
    if dock_atoms is None:
        dock_atoms = PF6.copy()
    ship_atoms = ship_atoms.copy()
    ship_atoms_center = ship_atoms.get_center_of_mass()
    ship_atoms_center[0] = max(ship_atoms.positions.T[0])
    dock_atoms_center = dock_atoms.get_center_of_mass()
    dock_atoms_center[0] = min(dock_atoms.positions.T[0])
    vector = ship_atoms_center - dock_atoms_center
    offset = [offset, 0, 0]
    dock_atoms.positions = dock_atoms.positions + vector + offset
    ship_atoms.extend(dock_atoms)
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
    import os

    from ase.io import write

    _CWD = os.getcwd()
    if create_tmp_folder:
        from monty.os import makedirs_p

        makedirs_p(f"{_CWD}/tmp")
        filepath = os.path.join(f"{_CWD}/tmp", filename)
    else:
        filepath = os.path.join(_CWD, filename)
    write(filepath, atoms, format="xyz")
    return filepath
