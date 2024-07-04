from __future__ import annotations

import os

from ase.io import read
from monty.os import cd
from xtb_ase import XTB
from pathlib import Path
from himatcal.recipes.gsm.SE_GSM import ASE_SE_GSM
from himatcal.utils.os import labeled_dir

main_workdir = Path.cwd()
label = "IMI"
atoms = read(f"{main_workdir}/input.xyz", ":")
driving_coords = [["BREAK", 1, 4]]

_CWD = labeled_dir(main_workdir, label)
gsm = ASE_SE_GSM(
    atom=atoms,
    driving_coords=driving_coords,
    calculator=XTB(method="gfn2-xtb", charge=-1, uhf=0, gbsa={"solvent": "acetone"}),
)

with cd(_CWD):
    gsm.run()
