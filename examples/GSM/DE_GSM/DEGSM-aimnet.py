from __future__ import annotations

import os
from pathlib import Path

from ase.io import read
from monty.os import cd

from himatcal.calculator.aimnet import AIMNet2ASE
from himatcal.recipes.gsm.DE_GSM import ASE_DE_GSM
from himatcal.utils.os import labeled_dir

main_workdir = Path.cwd()
label = "IMI"
reactant = read(f"{main_workdir}/reactant.xyz")
product = read(f"{main_workdir}/product.xyz")
charge = 1
multiplicity = 1


_CWD = labeled_dir(main_workdir, label)
#########   xTB example    ###############
gsm = ASE_DE_GSM(
    reactant=reactant,
    product=product,
    calculator=AIMNet2ASE(
        "aimnet2_b973c", charge=charge, mult=multiplicity
    ),  # * aimnet2/aimnet2_wb97m , aimnet2_b973c, aimnet2-qr
    multiplicity=multiplicity,
    fixed_product=True,
)

with cd(_CWD):
    gsm.run()
