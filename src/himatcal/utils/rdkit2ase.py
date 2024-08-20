"""https://github.com/zincware/rdkit2ase/tree/main"""

from __future__ import annotations

import functools
import io
import operator
import pathlib
import subprocess
import tempfile
from typing import Union

import ase.io
import numpy as np
import rdkit.Chem.AllChem
import rdkit.Chem.rdDetermineBonds
import rdkit.Geometry
from rdkit import Chem
from rdkit.Chem import rdDistGeom

OBJ_OR_STR = Union[str, Chem.rdchem.Mol, ase.Atoms]

OBJ_OR_STR_OR_LIST = Union[OBJ_OR_STR, list[tuple[OBJ_OR_STR, float]]]


def rdkit2ase(mol) -> ase.Atoms:
    """Convert an RDKit molecule to an ASE atoms object."""
    mol = rdkit.Chem.AddHs(mol)
    rdkit.Chem.AllChem.EmbedMolecule(mol)
    rdkit.Chem.AllChem.UFFOptimizeMolecule(mol)

    return ase.Atoms(
        positions=mol.GetConformer().GetPositions(),
        numbers=[atom.GetAtomicNum() for atom in mol.GetAtoms()],
    )


def ase2rdkit(atoms: ase.Atoms) -> rdkit.Chem.Mol:
    """Convert an ASE Atoms object to an RDKit molecule."""
    with io.StringIO() as f:
        ase.io.write(f, atoms, format="xyz")
        f.seek(0)
        xyz = f.read()
        raw_mol = rdkit.Chem.MolFromXYZBlock(xyz)

    mol = rdkit.Chem.Mol(raw_mol)
    rdkit.Chem.rdDetermineBonds.DetermineBonds(
        mol, charge=int(sum(atoms.get_initial_charges()))
    )
    return mol


def smiles2atoms(smiles: str) -> ase.Atoms:
    """
    Convert a SMILES string to an ASE Atoms object.

    Args:
        smiles (str): The SMILES string.

    Returns:
        atoms (ase.Atoms): The Atoms object.
    """
    mol = Chem.MolFromSmiles(smiles)
    return rdkit2ase(mol)


def smiles2conformers(
    smiles: str,
    numConfs: int,
    randomSeed: int = 42,
    maxAttempts: int = 1000,
) -> list[ase.Atoms]:
    """Create multiple conformers for a SMILES string.

    Args:
        smiles (str): The SMILES string.
        numConfs (int): The number of conformers to generate.
        randomSeed (int): The random seed.
        maxAttempts (int): The maximum number of attempts.

    Returns:
        images (list[ase.Atoms]): The list of conformers.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    rdDistGeom.EmbedMultipleConfs(
        mol,
        numConfs=numConfs,
        randomSeed=randomSeed,
        maxAttempts=maxAttempts,
    )

    images: list[ase.Atoms] = []

    for conf in mol.GetConformers():
        atoms = ase.Atoms(
            positions=conf.GetPositions(),
            numbers=[atom.GetAtomicNum() for atom in mol.GetAtoms()],
        )
        images.append(atoms)

    return images


def _get_cell_vectors(images: list[ase.Atoms], density: float) -> list[float]:
    """Get the box size from the molar volume.

    Attributes
    ----------
    images : list[ase.Atoms]
        All the atoms that should be packed.
    density: float
        Density of the system in kg/m^3.
    """
    molar_mass = sum(sum(atoms.get_masses()) for atoms in images)
    molar_volume = molar_mass / density / 1000  # m^3 / mol

    # convert to particles / A^3
    volume = molar_volume * ase.units.m**3 / ase.units.mol

    return [volume ** (1 / 3) for _ in range(3)]


def pack(
    data: list[list[ase.Atoms]],
    counts: list[int],
    density: float,
    seed: int = 42,
    tolerance: float = 2,
    logging: bool = False,
) -> ase.Atoms:
    """
    Pack the given molecules into a box with the specified density.

    Parameters
    ----------
    data : list[list[ase.Atoms]]
        A list of lists of ASE Atoms objects representing the molecules to be packed.
    counts : list[int]
        A list of integers representing the number of each type of molecule.
    density : float
        The target density of the packed system in kg/m^3.
    seed : int, optional
        The random seed for reproducibility, by default 42.
    tolerance : float, optional
        The tolerance for the packing algorithm, by default 2.
    logging : bool, optional
        If True, enables logging of the packing process, by default False.

    Returns
    -------
    ase.Atoms
        An ASE Atoms object representing the packed system.

    Example
    -------
    >>> from rdkit2ase import pack, smiles2conformers
    >>> water = smiles2conformers("O", 1)
    >>> ethanol = smiles2conformers("CCO", 1)
    >>> density = 1000  # kg/m^3
    >>> packed_system = pack([water, ethanol], [7, 5], density)
    >>> print(packed_system)
    Atoms(symbols='C10H44O12', pbc=True, cell=[8.4, 8.4, 8.4])
    """
    rng = np.random.default_rng(seed)
    selected_idx: list[np.ndarray] = []

    for images, count in zip(data, counts):
        selected_idx.append(
            rng.choice(range(len(images)), count, replace=len(images) < count)
        )

    images = [
        [data[category][idx] for idx in indices]
        for category, indices in enumerate(selected_idx)
    ]
    images = functools.reduce(operator.iadd, images, [])

    cell = _get_cell_vectors(images=images, density=density)

    file = f"""
tolerance {tolerance}
filetype xyz
output mixture.xyz
pbc 0 0 0 {" ".join([f"{x:.6f}" for x in cell])}
    """
    for category, indices in enumerate(selected_idx):
        for idx in indices:
            file += f"""
structure struct_{category}_{idx}.xyz
    filetype xyz

end structure
                     """

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = pathlib.Path(tmpdir)
        for category, indices in enumerate(selected_idx):
            for idx in set(indices):
                atoms = data[category][idx]
                ase.io.write(
                    tmpdir_path / f"struct_{category}_{idx}.xyz", atoms, format="xyz"
                )
        (tmpdir_path / "pack.inp").write_text(file)
        subprocess.run(
            "packmol < pack.inp",
            cwd=tmpdir_path,
            shell=True,
            check=True,
            capture_output=not logging,
        )
        atoms: ase.Atoms = ase.io.read(tmpdir_path / "mixture.xyz")

    atoms.cell = cell
    atoms.pbc = True
    return atoms
