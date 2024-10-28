from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def mol2xyz(mol, comment=None):
    """
    Convert rdkit mol to xyz file.
    from https://github.com/yanfeiguan/QM_descriptors_calculation/blob/master/lib/file_parser.py
    """
    c = mol.GetConformers()[0]
    coords = c.GetPositions()
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]

    xyz = f"{len(atoms)}\n{comment}\n"
    for a, c in zip(atoms, coords):
        xyz += "{}     {:14.9f}    {:14.9f}    {:14.9f}\n".format(a, *c)

    return xyz

def xyz2mol(xyz, smiles):
    """
    Convert xyz files to rdkit mol.
    from https://github.com/yanfeiguan/QM_descriptors_calculation/blob/master/lib/file_parser.py
    """
    lines = xyz.splitlines()
    N_atoms = int(lines[0])
    comments = lines[1]

    if N_atoms != len(lines[2:]):
        raise ValueError("Number of atoms does not match")

    mol = Chem.MolFromSmiles(smiles)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    mol = Chem.AddHs(mol, addCoords=True)
    try:
        conf = mol.GetConformers()[0]
    except Exception:
        AllChem.EmbedMultipleConfs(
            mol,
            numConfs=1,
            pruneRmsThresh=0.5,
            randomSeed=1,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
        )
        try:
            conf = mol.GetConformers()[0]
        except Exception:
            return None, None

    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    for i, coord in enumerate(lines[2:]):
        split_coord = coord.split()

        if atoms[i] != split_coord[0]:
            raise ValueError("Atom does not match")

        conf.SetAtomPosition(i, np.array(split_coord[1:]).astype("float"))

    mol.SetProp("comments", comments)
    return mol, comments
