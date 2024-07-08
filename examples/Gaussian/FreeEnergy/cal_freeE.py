import uuid

from ase.io import read
from monty.os import cd
from monty.serialization import loadfn
from pathlib import Path
from himatcal.recipes.gaussian.flow import calc_free_energy
from himatcal.tools.db import save_to_db

# Calculate the free energy of EC-POF3-c1s2
workdir =Path.cwd()
atoms = read(f"{workdir}/input.xyz")
charge = 0
mult = 1
label = "EC-c0s1"

with cd(workdir):
    free_energy = calc_free_energy(
        atoms, charge=charge, mult=mult, label=label, relax=False
    )

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
        "/home/suncc/Code/himatcal/examples/Gaussian/FreeEnergy/quacc-2024-06-14-02-13-29-893148-31118/quacc_results.json"
    ),
}

save_to_db(label, info)
