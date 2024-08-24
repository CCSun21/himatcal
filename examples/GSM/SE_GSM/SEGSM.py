from __future__ import annotations

import os
from pathlib import Path

from ase.io import read
from monty.os import cd
from xtb_ase import XTB

from himatcal.recipes.gsm.SE_GSM import ASE_SE_GSM
from himatcal.utils.os import labeled_dir

main_workdir = Path.cwd()
label = "IMI"
atoms = [read(f"{main_workdir}/input.xyz")]  # * list of atoms
driving_coords = [["BREAK", 1, 4]]

_CWD = labeled_dir(main_workdir, label)
# xtb test
gsm = ASE_SE_GSM(
    atom=atoms,
    driving_coords=driving_coords,
    calculator=XTB(method="gfn2-xtb", charge=-1, uhf=0, gbsa={"solvent": "acetone"}),
)

# # orca test
# from ase.calculators.orca import ORCA
# from ase.calculators.orca import OrcaProfile

# profile = OrcaProfile(command="/home/suncc/orca_6_0_0/orca")

# calculator = ORCA(
#     profile=profile,
#     charge=-1,
#     mult=1,
#     orcasimpleinput="B3LYP def2-TZVP EnGrad",
#     orcablocks="%pal nprocs 16 end \n%maxcore 1000",
# )

# gsm = ASE_SE_GSM(
#     atom=atoms,
#     driving_coords=driving_coords,
#     calculator=calculator,
# )

with cd(_CWD):
    gsm.run()
