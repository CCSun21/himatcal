import uuid

from ase.io import read
from monty.os import cd
from monty.serialization import loadfn

from himatcal.recipes.gaussian.flow import calc_free_energy
from himatcal.tools.db import save_to_db

# Calculate the free energy of EC-POF3-c1s2
workdir = "/home/suncc/Code/himatcal/examples/Gaussian/FreeEnergy"
atoms = read(f"{workdir}/final.xyz")
charge = 1
mult = 2
label = "EC-POF3-c1s2"

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
