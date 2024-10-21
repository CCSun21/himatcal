from __future__ import annotations

import logging

from ase.io import read

from himatcal.recipes.electrolyte.sol_stru.genFF import genFF

atoms = read("./EC.xyz")
result = genFF(atoms=atoms, label="EC") # type: ignore
logging.info(result)

