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
atoms = read(f"{main_workdir}/input.xyz")

"""
# * set driving coordinates with format ["BREAK"/"ADD", atom1, atom2] or ["ANGLE", atom1, atom2, atom3], or ["TORSION", atom1, atom2, atom3, atom4], or ["OOP", atom1, atom2, atom3, atom4]
"""
driving_coords = [["BREAK", 1, 4]]

_CWD = labeled_dir(main_workdir, label)
#########   xTB example    ###############
gsm = ASE_SE_GSM(
    atoms=atoms,
    driving_coords=driving_coords,
    calculator=XTB(method="gfn2-xTB", charge=-1, uhf=0, gbsa={"solvent": "acetone"}),
)

#########   ORCA example   ###############
# from ase.calculators.orca import ORCA
# from ase.calculators.orca import OrcaProfile

# profile = OrcaProfile(command="/home/suncc/orca_6_0_0/orca")

# calc = ORCA(
#     profile=profile,
#     charge=-1,
#     mult=1,
#     orcasimpleinput="B3LYP g-d3 def2-TZVP EnGrad", # using EnGrad for force calculation
#     orcablocks="%pal nprocs 16 end \n%maxcore 1000",
# )

# gsm = ASE_SE_GSM(
#     atom=atoms,
#     driving_coords=driving_coords,
#     calculator=calc,
# )


######### Gaussian example ###############
# from ase.calculators.gaussian import Gaussian

# calc = Gaussian(
#     charge=-1,
#     mult=1,
#     label="IMI",
#     method="B3LYP",
#     basis="6-31G(d)",
#     scf="xqc",
#     force="",  # * remember to return force
#     nosymm="",
#     mem="64GB",
#     nprocshared=16,
# )
# gsm = ASE_SE_GSM(
#     atom=atoms,
#     driving_coords=driving_coords,
#     calculator=calc,
# )

with cd(_CWD):
    gsm.run()
