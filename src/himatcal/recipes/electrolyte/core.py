from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Literal

from himatcal.recipes.gaussian.core import relax_job, static_job
from himatcal.utils.os import cclib_result

if TYPE_CHECKING:
    from ase import Atoms


class RedoxPotential:
    def __init__(
        self,
        neutral_molecule: Atoms | None = None,
        charged_molecule: Atoms | None = None,
        chg_mult: list[int]
        | None = None,  # * default for oxdiation potential calculation, for reduction, set [1, 1, 0, 2]
        calc_type: Literal["ox", "re"] = "ox",
        calc_kwards: dict | None = None,
    ):
        if calc_kwards is None:
            calc_kwards = {
                "opt_xc": "b3lyp",
                "opt_basis": "6-31G* em=GD3BJ",
                "gas_xc": "b3lyp",
                "gas_basis": "6-311+G**",
                "sol_xc": "m062x",
                "sol_basis": "6-31G*",  # only for solvent gibbs free energy correction
                "solvent": "Acetone",
            }
        if chg_mult is None:
            chg_mult = [-1, 1, 0, 2]
        self.neutral_molecule = neutral_molecule
        self.charged_molecule = charged_molecule
        self.chg_mult = chg_mult
        self.calc_type = calc_type
        self.calc_kwards = calc_kwards
        # Ensure 'gas_xc' and 'gas_basis' keys exist
        self.calc_kwards.setdefault("gas_xc", "m062x")
        self.calc_kwards.setdefault("gas_basis", "6-311G**")

    # * 1. relax the molecule in low level of theory
    def relax_llot(self, chg: int, mult: int, kwargs: dict | None = None):
        # sourcery skip: class-extract-method
        """
        low level of theory calculation for neutral and charged molecule, using for structure optimization and gibbs free energy correction
        if using in solvent, please pass {"scrf": ["SMD", f"solvent={self.calc_kwards["solvent"}"]}  to kwargs
        """
        if kwargs is None:
            kwargs = {}
        calc_keywords = {
            "label": "relax_llot",
            "mem": "64GB",
            "chk": "Gaussian.chk",
            "nprocshared": 64,
            "xc": "b3lyp",
            "basis": "6-31G* em=GD3BJ",
            "opt": "",
            "scf": ["maxcycle=250", "xqc"],
            "integral": "ultrafine",
            "nosymmetry": "",
            "pop": "CM5",
            "ioplist": ["2/9=2000"],
        } | kwargs
        logging.info(f"Relaxing {chg} charge molecule in low level of theory")
        quacc_results = relax_job(
            self.molecule,
            charge=chg,
            spin_multiplicity=mult,
            freq=True,
            **calc_keywords,
        )
        cclib_results = cclib_result(Path(quacc_results["dir_name"]))
        logging.info(f"Relaxation of {chg} charge molecule in low level of theory done")
        return (quacc_results, cclib_results)

    def sp_hlot(self, chg: int, mult: int, kwargs: dict):
        """
        high level of theory single point energy calculation for neutral and charged molecule, using for gibbs the base of free energy
        if using in solvent, please pass {"scrf": ["SMD", f"solvent={self.calc_kwards["solvent"}"]}  to kwargs
        """
        calc_keywords = {
            "label": "sp_hlot",
            "mem": "64GB",
            "chk": "Gaussian.chk",
            "nprocshared": 64,
            "xc": "m062x",
            "basis": "6-311+G**",
            "scf": ["maxcycle=250", "xqc"],
            "integral": "ultrafine",
            "nosymmetry": "",
            "pop": "CM5",
            "ioplist": ["2/9=2000"],
        } | kwargs
        logging.info(
            f"Calculating single point energy for {chg} charge molecule in high level of theory"
        )
        quacc_results = static_job(
            self.molecule, charge=chg, spin_multiplicity=mult, **calc_keywords
        )
        cclib_results = cclib_result(Path(quacc_results["dir_name"]))
        logging.info(
            f"Single point energy calculation for {chg} charge molecule in high level of theory done"
        )
        return (quacc_results, cclib_results)

    def cal_energy(
        self,
        chg_status: Literal["neutral", "charged"],
        phase_status: Literal["gas", "solvent"],
    ):
        """
        Calculates the Gibbs free energy and single point energy of a molecular system.

        This function determines the energy of a molecule based on its charge and phase status, performing necessary calculations to obtain the Gibbs free energy and single point energy. It adjusts the molecule's properties according to the specified charge and phase, and logs the results.

        Args:
            chg_status: A string indicating the charge status of the molecule, either "neutral" or "charged".
            phase_status: A string indicating the phase of the molecule, either "gas" or "solvent".

        Returns:
            tuple: A tuple containing the Gibbs free energy and single point energy in eV.

        """
        if chg_status == "neutral":
            if self.neutral_molecule is None:
                raise ValueError("neutral_molecule is None")
            self.molecule = self.neutral_molecule.copy()
            self.chg = self.chg_mult[0]
            self.mult = self.chg_mult[1]
        elif chg_status == "charged":
            if self.charged_molecule is None:
                raise ValueError("charged_molecule is None")
            self.molecule = self.charged_molecule.copy()
            self.chg = self.chg_mult[2]
            self.mult = self.chg_mult[3]
        if phase_status == "solvent":
            kwargs = {
                "xc": self.calc_kwards["sol_xc"],
                "basis": self.calc_kwards["sol_basis"],
                "scrf": ["SMD", f"solvent={self.calc_kwards['solvent']}"],
            }
        elif phase_status == "gas":
            kwargs = {
                "xc": self.calc_kwards["gas_xc"],
                "basis": self.calc_kwards["gas_basis"],
            }
        # * 1. relax the molecule in low level of theory
        opt_kwargs = {
            "xc": self.calc_kwards["opt_xc"],
            "basis": self.calc_kwards["opt_basis"],
        }
        relax_results, relax_cclib_results = self.relax_llot(
            chg=self.chg, mult=self.mult, kwargs=opt_kwargs
        )
        self.molecule = relax_results["atoms"]
        # * 2.calculate the single point energy from low level of theory to get the Gibbs free energy
        sp_results, sp_cclib_results = self.sp_hlot(
            chg=self.chg, mult=self.mult, kwargs=kwargs
        )

        # * 3.calculate the single point energy at solvation level to get the spe energy

        Gibbs_energy = (
            sp_cclib_results.scfenergies[0]
            + relax_cclib_results.freeenergy * 27.211
            - relax_cclib_results.scfenergies[0]
        )
        SPE_energy = sp_cclib_results.scfenergies[0]
        logging.info(
            f"{chg_status} molecule in {phase_status} phase Gibbs free energy: {Gibbs_energy} eV, Single point energy: {SPE_energy} eV"
        )
        return Gibbs_energy, SPE_energy

    def cal_cycle(self):
        """
        calculate the potential from neutral and charged molecule, return the potential, the unit is eV
        real_potential = potential - 1.44 eV
        """
        neutral_gas_gibbs, neutral_gas_spe = self.cal_energy(
            chg_status="neutral", phase_status="gas"
        )
        charged_gas_gibbs, charged_gas_spe = self.cal_energy(
            chg_status="charged", phase_status="gas"
        )
        neutral_solvent_gibbs, neutral_solvent_spe = self.cal_energy(
            chg_status="neutral", phase_status="solvent"
        )
        charged_solvent_gibbs, charged_solvent_spe = self.cal_energy(
            chg_status="charged", phase_status="solvent"
        )

        # * \delta G_{gas} = G_{charged}^{gas} - G_{neutral}^{gas}
        delta_G_gas = charged_gas_gibbs - neutral_gas_gibbs

        # * \delta G_{solvention}(neutral) = E_{neutral}^{solvent} - E_{neutral}^{gas}
        # * \delta G_{solvention}(charged) = E_{charged}^{solvent} - E_{charged}^{gas}
        delta_G_solvention_neutral = neutral_solvent_spe - neutral_gas_spe
        delta_G_solvention_charged = charged_solvent_spe - charged_gas_spe

        if self.calc_type == "ox":
            potential = (
                delta_G_gas - delta_G_solvention_neutral + delta_G_solvention_charged
            )
        if self.calc_type == "re":
            potential = -(
                delta_G_gas - delta_G_solvention_neutral + delta_G_solvention_charged
            )

        logging.info(
            r"default potential unit is V, referring to Li/Li+, which is chosen as E_{abs} - 4.42 + 3.02 = E_{abs} - 1.4  V; To convert to SHE, please add 1.4 and then minus 4.42 V; more detail please refer to the paper Borodin, O.; Behl, W.; Jow, T. R. Oxidative Stability and Initial Decomposition Reactions of Carbonate, Sulfone, and Alkyl Phosphate-Based Electrolytes. J. Phys. Chem. C 2013, 117 (17), 8661-8682. https://doi.org/10.1021/jp400527c."
        )
        return potential - 1.4
