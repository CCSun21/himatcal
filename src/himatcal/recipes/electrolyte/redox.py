from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional

from himatcal.atoms.core import PF6, dock_atoms
from himatcal.recipes.crest.core import protonate, relax
from himatcal.recipes.electrolyte.core import RedoxPotential

if TYPE_CHECKING:
    from ase import Atoms


class RedoxCal:
    """
    RedoxCal is a class for managing and calculating redox potentials in molecular systems.

    This class defines the necessary attributes and parameters for performing redox calculations, including the molecular structure, charge multiplicities, and calculation settings. It provides a structured approach to facilitate the setup and execution of redox potential calculations.

    Attributes:
        label: A string representing the label for the redox calculation, defaulting to "redox".
        molecule: An optional Atoms object representing the molecule for the calculation.
        charge_mult: A list of integers specifying the charge multiplicities for the calculation, defaulting to [-1, 1, 0, 2, 1, 1, 0, 2].
        add_anion: A boolean indicating whether to include an anion in the system, defaulting to True.
        ions: A list containing Atoms or strings representing the ions involved in the calculation, defaulting to [PF6, 'li'].
        calc_kwards: A dictionary containing keyword arguments for the calculation methods, including optimization and solvent settings.

    """

    def __init__(
        self,
        molecule: Atoms | None = (None,),
        chg_mult: list[int] | None = None,
        add_ion: bool = (True,),
        ions: list[Atoms | str]
        | None = None,  # set [PF6, 'na'] for sodium solvent calculation
        label: str = ("redox",),
        calc_kwards: dict | None = None,
        machine_kwards: dict | None = None,
    ):
        if chg_mult is None:
            chg_mult = [-1, 1, 0, 2, 1, 1, 0, 2]
        if ions is None:
            ions = [PF6, "li"]
        if calc_kwards is None:
            calc_kwards = {
                "opt_xc": "b3lyp",
                "opt_basis": "6-311+G(d,p)",
                "sol_xc": "m062x",
                "sol_basis": "6-31G*",
                "solvent": "Acetone",
            }
        if machine_kwards is None:
            machine_kwards = {"xtb_proc": 16}
        self.molecule = molecule
        self.chg_mult = chg_mult
        self.add_ion = add_ion
        self.ions = ions
        self.label = label

        self.calc_kwards = calc_kwards
        self.machine_kwards = machine_kwards

    def get_ox(self):
        # * generate solvated molecules using anion and counter-ion
        if self.add_ion:
            neutral_molecule = dock_atoms(
                self.molecule,
                dock_atoms=self.ions[0],
                crest_sampling=True,
                chg=self.chg_mult[0],
                mult=self.chg_mult[1],
            )
            charged_molecule = dock_atoms(
                self.molecule,
                dock_atoms=self.ions[0],
                crest_sampling=True,
                chg=self.chg_mult[2],
                mult=self.chg_mult[3],
            )
        else:
            neutral_molecule = relax(
                self.molecule, chg=self.chg_mult[0], mult=self.chg_mult[1]
            )
            charged_molecule = relax(
                self.molecule, chg=self.chg_mult[2], mult=self.chg_mult[3]
            )

        # * calculate the oxidation state energies (in eV)
        redox_potential = RedoxPotential(
            neutral_molecule=neutral_molecule,
            charged_molecule=charged_molecule,
            chg_mult=self.chg_mult[:4],
            calc_type="ox",
            calc_kwards=self.calc_kwards,
        ).cal_cycle()
        return redox_potential

    def get_re(self):
        # * generate solvated molecules using anion and counter-ion
        if self.add_ion:
            neutral_molecule = protonate(
                self.molecule,
                ion=self.ions[1],
                chg=self.chg_mult[4],
                mult=self.chg_mult[5],
                threads=16,
            )
            charged_molecule = protonate(
                self.molecule,
                ion=self.ions[1],
                chg=self.chg_mult[6],
                mult=self.chg_mult[7],
                threads=16,
            )
        else:
            neutral_molecule = relax(
                self.molecule, chg=self.chg_mult[4], mult=self.chg_mult[5]
            )
            charged_molecule = relax(
                self.molecule, chg=self.chg_mult[6], mult=self.chg_mult[7]
            )

        # * calculate the oxidation state energies (in eV)
        redox_potential = RedoxPotential(
            neutral_molecule=neutral_molecule,
            charged_molecule=charged_molecule,
            chg_mult=self.chg_mult[4:8],
            calc_type="re",
            calc_kwards=self.calc_kwards,
        ).cal_cycle()
        return redox_potential

    def get_redox(self):
        oxidation_potential = self.get_ox()
        reduction_potential = self.get_re()
        return [oxidation_potential, reduction_potential]
