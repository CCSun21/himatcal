from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from pydantic import BaseModel

from himatcal.recipes.gaussian.core import relax_job

if TYPE_CHECKING:
    from ase import Atoms


class RedoxPotential(BaseModel):
    molecule: Atoms | None = (None,)
    chg_mult: list[int] = ([-1, 1, 0, 2, 1, 1, 0, 2],)
    cal_type: Literal["ox", "re"] = ("ox",)
    calc_kwards: dict = {
        "opt_method": "b3lyp",
        "opt_basis": "6-311+G(d,p)",
        "sol_method": "m062x",
        "sol_basis": "6-31G*",
        "solvent": "Acetone",
    }

    # * 1. relax the molecule in low level of theory
    def relax_llot(self,chg,mult):
        calc_keywords = {
            "label": "relax_llot",
            "mem": "64GB",
            "chk": "Gaussian.chk",
            "nprocshared": 64,
            "xc": "b3lyp",
            "basis": "6-31G*",
            "opt": "",
            "scf": ["maxcycle=250", "xqc"],
            "integral": "ultrafine",
            "nosymmetry": "",
            "pop": "CM5",
            "ioplist": ["2/9=2000"]
        }
        relax_result = relax_job(
            self.molecule,
            charge=chg,
            mult=mult,
            **calc_keywords
        )
        return relax_result["atoms"]
