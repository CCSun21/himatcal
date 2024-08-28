from __future__ import annotations

import logging
import uuid
from pathlib import Path

from ase.io import read, write
from quacc import change_settings

from himatcal.recipes.gaussian.core import relax_job
from himatcal.utils.os import get_chg_mult

CWD = Path.cwd()
label, chg, mult = get_chg_mult(Path(CWD / "EC_TS-c1s2.xyz").stem)
logging.info(f"Atomfile label: {label}, charge: {chg}, multiplicity: {mult}")
atoms = read(CWD / "EC_TS-c1s2.xyz")
job_id = uuid.uuid4()
label = str(Path(CWD / "EC_TS-c1s2.xyz").stem)
job_dir = CWD / f"himatcal-{job_id!s}"
Path.mkdir(job_dir)
logging.info(f"Creating directory for job: {job_id}")

write(job_dir / f"{label}.xyz", atoms, format="xyz")
logging.info(f"Writing xyz file: {label}")

logging.info(f"Calculating {label}")

calc_keywords = {
    "mem": "64GB",
    "label": label,
    "chk": "Gaussian.chk",
    "nprocshared": 4,
    "xc": "b3lyp",
    "basis": "6-311+G* em=GD3BJ",
    "opt": "calcall ts noeigen",
    "scf": ["maxcycle=250", "xqc"],
    "integral": "ultrafine",
    "nosymmetry": "",
    "pop": "CM5",
    "ioplist": ["2/9=2000"],
    "scrf": ["pcm", "solvent=acetone"],
}

with change_settings({"RESULTS_DIR": job_dir}):
    logging.info(f"Calculating {label}")
    result = relax_job(
        atoms, charge=chg, spin_multiplicity=mult, freq=True, **calc_keywords
    )
    logging.info(f"Completed {label}")
    write(job_dir / f"{label}_relaxed.xyz", result["atoms"], format="xyz")
