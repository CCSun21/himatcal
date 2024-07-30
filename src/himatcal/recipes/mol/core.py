from __future__ import annotations

from io import StringIO

from chemspipy import ChemSpider
from rdkit import Chem

from himatcal import SETTINGS


def get_molecular_structure(
    molecular_cas: str,
    write_mol: str | None = None,
    chemspider_api: str = SETTINGS.CHEMSPIDER_API_KEY,
):
    """
    Get molecular structure from CAS number.

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
    except Exception as e:
        return f"Unexpected error: {e}"
