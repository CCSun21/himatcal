"""Example of using the relax_job function to relax a structure using the mace-mp-0 recipe."""

from __future__ import annotations

import logging
from pathlib import Path

from ase.io import read, write

from himatcal.calculator.aimnet import AIMNet2ASE
from himatcal.recipes.quacc.core import relax_job

atoms_path = Path.cwd() / "input.xyz"
charge = -1
mult = 1

calc = AIMNet2ASE("aimnet2_b973c", charge=charge, mult=mult)  # * aimnet2/aimnet2_wb97m , aimnet2_b973c, aimnet2-qr

atoms = read(atoms_path)
result = relax_job(atoms, calc)
logging.info(result)

# * write the relaxed structure to a file
write("relaxed.xyz", result["atoms"])
