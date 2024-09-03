"""core functions for GSM"""
from __future__ import annotations

import numpy as np
from ase import Atoms
from pyGSM.utilities import elements


def atoms2geom(atoms: Atoms):
    """
    Converts an Atoms object into a geometric representation suitable for further processing. This function extracts atomic symbols and positions, returning them in a structured format.

    The function takes an Atoms object as input and retrieves the atomic positions and symbols. It then constructs a list of atomic data and a geometric representation, which includes the symbols and their corresponding coordinates.

    Args:
        atoms (Atoms): The Atoms object containing atomic positions and symbols.

    Returns:
        tuple: A tuple containing:
            - List of atom objects corresponding to the atomic symbols.
            - Numpy array of atomic positions.
            - List of geometric representations combining symbols and positions.

    Examples:
        To convert an Atoms object to its geometric representation, call:
        >>> atom_list, positions, geometry = atoms2geom(atoms)
    """

    xyz = atoms.positions
    geom = np.column_stack([atoms.symbols, xyz]).tolist()
    ELEMENT_TABLE = elements.ElementData()
    atom = [ELEMENT_TABLE.from_symbol(atom) for atom in atoms.symbols]
    return atom, xyz, geom


def gsm2atoms(gsm):
    """
    convert GSM DE_GSM object to ASE Atoms objects
    """
    # string
    frames = []
    for energy, geom in zip(gsm.energies, gsm.geometries):
        at = Atoms(symbols=[x[0] for x in geom], positions=[x[1:4] for x in geom])
        at.info["energy"] = energy
        frames.append(at)

    # TS
    ts_geom = gsm.nodes[gsm.TSnode].geometry
    ts_atoms = Atoms(symbols=[x[0] for x in ts_geom], positions=[x[1:4] for x in ts_geom])

    return frames, ts_atoms
