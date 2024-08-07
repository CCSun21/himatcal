from __future__ import annotations

import uuid
from pathlib import Path

from ase.io import read
from monty.os import cd
from monty.serialization import loadfn

from himatcal.recipes.gaussian.flow import calc_free_energy
from himatcal.tools.db import save_to_db

# Calculate the free energy of EC-POF3-c1s2
workdir =Path.cwd()
atoms = read(f"{workdir}/input.xyz")
charge = 0
mult = 1
label = "EC-c0s1"

with cd(workdir):
    freeE = calc_free_energy(
        atoms, charge=charge, mult=mult, label=label, relax=False
    )
    freeE.run()
    free_energy = freeE.extract_free_energy()

# Save the result to the database
info = {
    "charge": charge,
    "multiplicity": mult,
    "free_energy": free_energy,
    "free_energy_unit": "a.u.",
    "molecule": "EC",
    "counter_ion": "POF3",
    "task_id": str(uuid.uuid4()),
    "quacc_results": loadfn(
        f"{workdir}/quacc-2024-06-14-02-13-29-893148-31118/quacc_results.json"
    ),
}

save_to_db(label, info)
