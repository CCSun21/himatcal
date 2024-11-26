from __future__ import annotations

import json
import urllib.request
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

import requests
from ase.io import write
from chemspipy import ChemSpider
from rdkit import Chem
from rdkit.Chem import AllChem

from himatcal import SETTINGS
from himatcal.recipes.crest.core import relax as crest_relax
from himatcal.utils.rdkit.core import rdkit2ase

if TYPE_CHECKING:
    from ase import Atoms


def get_molecular_structure(
    molecular_cas: str,
    write_mol: bool = True,
    chemspider_api: str = SETTINGS.CHEMSPIDER_API_KEY,
):
    """
    Get molecular structure from CAS number, using chemspipy from RSC ChemSpider.

    This function retrieves the molecular structure corresponding to the provided CAS number from ChemSpider, processes it, and optionally writes it to a file in XYZ format.

    Args:
        molecular_cas (str): The CAS number of the molecule.
        write_mol (str | None): The file name to write the molecular structure to in XYZ format. Defaults to None.
        chemspider_api (str): The ChemSpider API key. Defaults to the value in SETTINGS.CHEMSPIDER_API_KEY.

    Returns:
        None
    """

    cs = ChemSpider(chemspider_api)
    c1 = cs.search(molecular_cas)[0]
    try:
        mol_file = StringIO(c1.mol_3d)
        mol = Chem.MolFromMolBlock(mol_file.getvalue(), removeHs=False)
        mol = Chem.AddHs(mol, addCoords=True)
        if write_mol:
            Chem.MolToXYZFile(mol, f"{molecular_cas}.xyz")
        return mol
    except Exception as e:
        return f"Unexpected error: {e}"


def consumeApi(urlPath):
    dataResponse = requests.get(urlPath)
    return None if (dataResponse.status_code != 200) else dataResponse.text


def cas2xyz(CAS_ID, relax_atoms=True):
    """
    Converts a CAS ID into an XYZ file format representation of the corresponding molecule.

    This function retrieves molecular data from the Common Chemistry API using the provided CAS ID, constructs a molecular structure, and optionally relaxes the atomic positions before saving the structure to an XYZ file.

    Args:
        CAS_ID (str): The Chemical Abstracts Service identifier for the desired molecule.
        relax_atoms (bool): A flag indicating whether to relax the atomic positions before saving. Defaults to True.

    Returns:
        None

    Raises:
        ValueError: If the CAS ID is invalid or if the API call fails.

    Examples:
        cas2xyz("50-00-0")  # Converts the CAS ID for formaldehyde to an XYZ file.
    """
    sources = [
        f"https://commonchemistry.cas.org/api/detail?cas_rn={CAS_ID}",
        f"https://www.chemicalbook.com/CAS/mol/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20210305/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20210111/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20180601/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20150408/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20200515/MOL/{CAS_ID}.mol",
    ]

    mol = None
    for source in sources:
        try:
            if "commonchemistry" in source:
                result_dict = json.loads(consumeApi(source))
                mol = Chem.MolFromInchi(result_dict["inchi"])
            else:
                request = urllib.request.Request(
                    source, headers={"User-Agent": "Mozilla/5.0"}
                )
                response = urllib.request.urlopen(request)
                if response.status == 200:
                    mol_file_content = response.read().decode("utf-8")
                    mol = Chem.MolFromMolBlock(mol_file_content, removeHs=False)
            if mol:
                break
        except Exception:
            # print(f"Failed to retrieve data from {source}: {e}")
            continue

    if not mol:
        try:
            mol = get_molecular_structure(CAS_ID)
        except Exception as e:
            raise ValueError(
                f"Could not retrieve molecular data for CAS ID {CAS_ID}"
            ) from e

    mol = Chem.AddHs(mol, addCoords=True)
    AllChem.EmbedMultipleConfs(mol, numConfs=10)
    if relax_atoms:
        atoms = rdkit2ase(mol)
        atoms_relaxed = crest_relax(atoms)
        write(f"{CAS_ID}.xyz", atoms_relaxed)
    else:
        Chem.MolToXYZFile(mol, f"{CAS_ID}.xyz")
    with Path.open(Path(f"{CAS_ID}.xyz")) as file:
        return file.read()


# TODO: include PubGrep ( https://github.com/grimme-lab/PubGrep )

# TODO: Methods or flow to relax the mols


def relax_mol(mol: Atoms, chg: int = 0, mult: int = 1, method: str = "crest", **kwargs):
    """
    Relax the molecular structure using the specified method.

    This function relaxes the atomic positions of the provided molecular structure using the specified method.

    Args:
        mol (Atoms): The molecular structure to relax.
        method (str): The relaxation method to use. Defaults to 'crest'.
        **kwargs: Additional keyword arguments for the relaxation method.

    Returns:
        Atoms: The relaxed molecular structure.

    Raises:
        ValueError: If the relaxation method is not supported.

    Examples:
        relax_mol(mol, method='crest')  # Relaxes the molecular structure using the CREST method.
    """
    if method == "crest":
        return crest_relax(mol, chg=chg, mult=mult, **kwargs)
    else:
        raise ValueError(f"Relaxation method '{method}' is not supported")
