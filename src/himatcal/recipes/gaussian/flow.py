"""recipes for gaussian calculation"""

from __future__ import annotations

import contextlib
import logging
from pathlib import Path

from himatcal.recipes.gaussian.core import relax_job, static_job

Logger = logging.getLogger(__name__)


class calc_free_energy:
    def __init__(
        self,
        atoms,
        charge=0,
        mult=1,
        label="Gaussian",
        relax=False,
        solvent="acetone",
        relax_stage1_params=None,
    ):
        self.atoms = atoms
        self.charge = charge
        self.mult = mult
        self.label = label
        self.relax = relax
        self.solvent = solvent

        if relax_stage1_params:
            self.relax_stage1_params = relax_stage1_params
        else:
            # * 1. relax at 6-31G(d) level
            self.relax_stage1_params = {
                "xc": "B3LYP",
                "basis": "6-31G(d)",
                "label": f"{self.label}",
            }
            if self.solvent:
                self.relax_stage1_params["scrf"] = f"pcm, solvent={self.solvent}"

    def relax_stage1(self):
        result = relax_job(
            self.atoms,
            charge=self.charge,
            spin_multiplicity=self.mult,
            **self.relax_stage1_params,
        )
        self.atoms = result["atoms"]
        return self.atoms

    def relax_stage2(self):
        # * 2. relax at 6-311+g(d) level
        calc_params = {
            "xc": "B3LYP",
            "basis": "6-311+g(d) em=GD3BJ",
            "label": f"{self.label}",
        }
        if self.solvent:
            calc_params["scrf"] = f"pcm, solvent={self.solvent}"
        result = relax_job(
            self.atoms,
            charge=self.charge,
            spin_multiplicity=self.mult,
            freq=True,
            **calc_params,
        )
        self.atoms = result["atoms"]
        return self.atoms

    def write_relaxed(self):
        from ase.io import write

        from himatcal.utils.os import write_chg_mult_label

        relaxed_filename = write_chg_mult_label(
            f"{self.label}_relaxed", self.charge, self.mult
        )
        write(f"{relaxed_filename}.xyz", self.atoms)
        return self.atoms

    def single_point(self):
        calc_params = {
            "xc": "B3LYP",
            "basis": "6-311+g(d) em=GD3BJ",
            "scrf": "pcm, solvent=acetone",
            "label": f"{self.label}",
            "freq": "",
        }
        if self.solvent:
            calc_params["scrf"] = f"pcm, solvent={self.solvent}"
        return static_job(
            self.atoms, charge=self.charge, spin_multiplicity=self.mult, **calc_params
        )

    def run(self):
        if self.relax:
            # * 1. relax at 6-31G(d) level
            self.relax_stage1()
            # * 2. relax at 6-311+g(d) level
            self.relax_stage2()
            # * write relaxed atoms to file
            self.write_relaxed()
        # * 3. Single point calculation and get the free energy.
        result = self.single_point()
        self.result = result
        return result

    def extract_free_energy(self):
        if not hasattr(self, "result"):
            raise AttributeError("Please run the calculation first!")

        log_path = Path(self.result["dir_name"])
        with contextlib.suppress(FileNotFoundError):
            return self._extracted_from_extract_free_energy(log_path)

    def _extracted_from_extract_free_energy(self, log_path):
        import gzip

        import cclib

        if gzip_log := list(log_path.glob("*.log.gz")):
            unzip_file = gzip.decompress(Path.open(gzip_log[0], "rb").read())
            logfile = gzip_log[0].with_suffix("")
            with Path.open(logfile, "w") as f:
                f.write(unzip_file.decode())
        log_files = list(log_path.glob("*.log"))
        data = cclib.io.ccread(log_path / log_files[0])
        free_energy = data.freeenergy
        with Path.open(log_path / "free_energy.txt", "w") as f:
            f.write(f"{free_energy}")
        return free_energy
