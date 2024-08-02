from __future__ import annotations

import os

from ase.io import read
from chemspipy import ChemSpider
from pydantic import Field

from himatcal import SETTINGS


def get_molecular_structure(
    molecular_cas: str = Field(None, description="CAS number"),
    chemspider_api: str = SETTINGS.CHEMSPIDER_API_KEY,
):
    cs = ChemSpider(chemspider_api)
    c1 = cs.search(molecular_cas)[0]
    try:
        ide = c1.record_id
        mol = c1.mol_3d
        name = c1.common_name
        open(f"tmp_{molecular_cas}.mol", "w").write(mol)
        molecular_atoms = read(f"./tmp_{molecular_cas}.mol")
        os.remove(f"tmp_{molecular_cas}.mol")
        return molecular_atoms
    except Exception as e:
        print(f"molecular dowload error:{e}")
