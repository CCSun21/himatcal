"""https://github.com/zincware/rdkit2ase/tree/main"""

from __future__ import annotations

import functools
import io
import operator
import pathlib
import subprocess
import tempfile
from typing import Union

import numpy as np
from ase import Atoms, units
from ase.io import read, write
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds, rdDistGeom

OBJ_OR_STR = Union[str, Chem.rdchem.Mol, Atoms]

OBJ_OR_STR_OR_LIST = Union[OBJ_OR_STR, list[tuple[OBJ_OR_STR, float]]]


def rdkit2ase(mol) -> Atoms:
    """Convert an RDKit molecule to an ASE atoms object."""
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    return Atoms(
        positions=mol.GetConformer().GetPositions(),
        numbers=[atom.GetAtomicNum() for atom in mol.GetAtoms()],
    )


def ase2rdkit(atoms: Atoms, charge:int=0) -> Chem.Mol:
    """Convert an ASE Atoms object to an RDKit molecule."""
    with io.StringIO() as f:
        write(f, atoms, format="xyz")
        f.seek(0)
        xyz = f.read()
        raw_mol = Chem.MolFromXYZBlock(xyz)

    mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineBonds(mol, charge=charge)
    return mol


def smiles2atoms(smiles: str) -> Atoms:
    """
    Convert a SMILES string to an ASE Atoms object.

    Args:
        smiles (str): The SMILES string.

    Returns:
        atoms (Atoms): The Atoms object.
    """
    mol = Chem.MolFromSmiles(smiles)
    return rdkit2ase(mol)


def smiles2conformers(
    smiles: str,
    numConfs: int,
    randomSeed: int = 42,
    maxAttempts: int = 1000,
) -> list[Atoms]:
    """Create multiple conformers for a SMILES string.

    Args:
        smiles (str): The SMILES string.
        numConfs (int): The number of conformers to generate.
        randomSeed (int): The random seed.
        maxAttempts (int): The maximum number of attempts.

    Returns:
        images (list[Atoms]): The list of conformers.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    rdDistGeom.EmbedMultipleConfs(
        mol,
        numConfs=numConfs,
        randomSeed=randomSeed,
        maxAttempts=maxAttempts,
    )

    images: list[Atoms] = []

    for conf in mol.GetConformers():
        atoms = Atoms(
            positions=conf.GetPositions(),
            numbers=[atom.GetAtomicNum() for atom in mol.GetAtoms()],
        )
        images.append(atoms)

    return images


def _get_cell_vectors(images: list[Atoms], density: float) -> list[float]:
    """Get the box size from the molar volume.

    Attributes
    ----------
    images : list[Atoms]
        All the atoms that should be packed.
    density: float
        Density of the system in kg/m^3.
    """
    molar_mass = sum(sum(atoms.get_masses()) for atoms in images)
    molar_volume = molar_mass / density / 1000  # m^3 / mol

    # convert to particles / A^3
    volume = molar_volume * units.m**3 / units.mol

    return [volume ** (1 / 3) for _ in range(3)]


def pack(
    data: list[list[Atoms]],
    counts: list[int],
    density: float,
    seed: int = 42,
    tolerance: float = 2,
    logging: bool = False,
) -> Atoms:
    """
    Pack the given molecules into a box with the specified density.

    Parameters
    ----------
    data : list[list[Atoms]]
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
    Atoms
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
                write(tmpdir_path / f"struct_{category}_{idx}.xyz", atoms, format="xyz")
        (tmpdir_path / "pack.inp").write_text(file)
        subprocess.run(
            "packmol < pack.inp",
            cwd=tmpdir_path,
            shell=True,
            check=True,
            capture_output=not logging,
        )
        atoms: Atoms = read(tmpdir_path / "mixture.xyz")

    atoms.cell = cell
    atoms.pbc = True
    return atoms


# * Gasteger charge visualization
def plot_gasteiger_charges(mol):
    """
    Plot Gasteiger charges on a molecule.
    """
    from rdkit.Chem.Draw import SimilarityMaps

    AllChem.ComputeGasteigerCharges(mol)
    contribs = [
        float(mol.GetAtomWithIdx(i).GetProp("_GasteigerCharge"))
        for i in range(mol.GetNumAtoms())
    ]
    SimilarityMaps.GetSimilarityMapFromWeights(
        mol, contribs, colorMap="jet", contourLines=10
    )

def smiles_to_rdkit(smi, gen_3d=True, nconf=100):
    """
    Convert smiles to RDKit molecule.
    Tries to generate the lowest-energy conformer.
    """
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)

    if gen_3d:
        cids = AllChem.EmbedMultipleConfs(mol, nconf, AllChem.ETKDG())

        AllChem.MMFFSanitizeMolecule(mol)
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol)

        energies = []
        for cid in cids:
            ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=cid)
            ff.Minimize()
            energy = ff.CalcEnergy()
            energies.append(energy)

        energies = np.asarray(energies)
        min_energy_idx = np.argsort(energies)[0]

        new_mol = Chem.Mol(mol)
        new_mol.RemoveAllConformers()
        min_conf = mol.GetConformer(cids[min_energy_idx])
        new_mol.AddConformer(min_conf, assignId=True)
        mol = new_mol

    return mol
