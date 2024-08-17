from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from himatcal.atoms.core import PF6, dock_atoms
from himatcal.recipes.crest.core import protonate, relax
from himatcal.recipes.electrolyte.core import RedoxPotential

if TYPE_CHECKING:
    from ase import Atoms


class RedoxCal:
    """
    A class to calculate the oxidation and reduction potentials of a molecular system.

    This class provides methods to compute the oxidation and reduction potentials based on the molecular structure and specified ions. It allows for the inclusion of ions and various calculation parameters to tailor the analysis.

    Attributes:
        molecule (Atoms | None): The molecular structure for the calculations.
        chg_mult (list[int] | None): Charge multiplicities for the calculations.
        add_ion (bool): Indicates whether to include an anion in the system.
        ions (list[Atoms | str] | None): Ions involved in the calculations.
        label (str): A label for the calculations.
        calc_kwards (dict | None): Keyword arguments for calculation methods.
        machine_kwards (dict | None): Machine-specific keyword arguments.

    Methods:
        get_ox(): Calculates the oxidation potential of the molecular system.
        get_re(): Calculates the reduction potential of the molecular system.
        get_redox(): Retrieves both oxidation and reduction potentials.
    """

    def __init__(
        self,
        molecule: Atoms | None = (None,),
        chg_mult: list[int] | None = None,
        add_ion: bool = (True,),
        ions: list[Atoms | str] | None = None,
        protonate_ion_string: bool = True,
        label: str = "redox",
        calc_kwards: dict | None = None,
        machine_kwards: dict | None = None,
    ):
        if chg_mult is None:
            chg_mult = [-1, 1, 0, 2, 1, 1, 0, 2]
        if ions is None:
            ions = [PF6, "Li"]  # set [PF6, 'Na'] for sodium solvent calculation
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
        self.protonate_ion_string = protonate_ion_string

        self.calc_kwards = calc_kwards
        self.machine_kwards = machine_kwards

    def get_ox(self):
        """
        Calculates the oxidation potential of a molecular system.

        This function generates the neutral and charged molecules required for calculating the oxidation potential, either by docking with ions or relaxing the molecule, depending on the presence of ions. It then computes the redox potential based on these generated molecules.

        Args:
            None

        Returns:
            float: The calculated oxidation potential in eV.

        """

        # * generate solvated molecules using anion and counter-ion
        if self.add_ion:
            logging.info("Generate and relax molecules clusters using crest")
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
            logging.info("Relaxing molecules using crest")
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
        logging.info(f"{self.label} oxidation potential: {redox_potential} eV")
        return redox_potential

    def get_re(self):
        """
        Calculates the reduction potential of a molecular system.

        This function generates the neutral and charged molecules required for calculating the reduction potential, either by protonation or relaxation, depending on the presence of ions. It utilizes a helper function to streamline the process of generating molecules and logs the resulting reduction potential.

        Args:
            add_ion: A boolean indicating whether to include an anion in the system.
            molecule: An Atoms object representing the molecule for the calculation.
            ions: A list containing Atoms or strings representing the ions involved in the calculation.
            chg_mult: A list of integers specifying the charge multiplicities for the calculation.
            calc_kwards: A dictionary containing keyword arguments for the calculation methods.

        Returns:
            float: The calculated reduction potential in eV.

        """

        def generate_molecule(
            molecule, ion, chg, mult, protonate_ion_string, threads=16
        ):
            """
            Generates a molecular structure by either protonating or docking the specified molecule.

            This function attempts to protonate the given molecule using the specified ion and charge parameters. If protonation fails, it falls back to docking the molecule with the specified ion, logging the failure of the protonation attempt.

            If the protonation will fail or results a mis-protonated molecule, please consider docking your ion into it by setting 'protonate_ion_string' to False in class RedoxCal. This will allow the docking attempt to work normally.

            Args:
                molecule: The molecular structure to be modified.
                ion: The ion used for protonation or docking.
                chg: The charge associated with the ion.
                mult: The multiplicity of the ion.
                protonate_ion_string: A boolean indicating whether to attempt protonation.
                threads (int, optional): The number of threads to use for the operation. Defaults to 16.

            Returns:
                The modified molecular structure after protonation or docking, or None if both attempts fail.
            """

            mol = (
                protonate(molecule, ion=ion, chg=chg, mult=mult, threads=threads)
                if protonate_ion_string is True and ion is str
                else None
            )
            if mol is None:
                logging.info("Protonation failed or skipped, trying docking")
                mol = dock_atoms(
                    molecule, dock_atoms=ion, crest_sampling=True, chg=chg, mult=mult
                )
            return mol

        if self.add_ion:
            logging.info("Generate and relax molecules clusters using crest")
            neutral_molecule = generate_molecule(
                self.molecule,
                self.ions[1],
                self.chg_mult[4],
                self.chg_mult[5],
                protonate_ion_string=self.protonate_ion_string,
            )
            charged_molecule = generate_molecule(
                self.molecule,
                self.ions[1],
                self.chg_mult[6],
                self.chg_mult[7],
                protonate_ion_string=self.protonate_ion_string,
            )
        else:
            logging.info("Relaxing molecules using crest")
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
        logging.info(f"{self.label} reduction potential: {redox_potential} eV")
        return redox_potential

    def get_redox(self):
        """
        Calculates the oxidation and reduction potentials of a molecular system.

        This function retrieves the oxidation and reduction potentials by calling the respective methods and returns them as a list. It provides a convenient way to access both potentials in a single call.

        Args:
            None

        Returns:
            list: A list containing the oxidation potential and reduction potential.

        """

        oxidation_potential = self.get_ox()
        reduction_potential = self.get_re()
        return [oxidation_potential, reduction_potential]
