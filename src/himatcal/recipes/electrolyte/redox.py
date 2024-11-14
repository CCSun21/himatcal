from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from pydantic import BaseModel, Field, field_validator

from himatcal.atoms.core import PF6, dock_atoms
from himatcal.recipes.crest.core import protonate, relax
from himatcal.recipes.electrolyte.core import RedoxPotential

if TYPE_CHECKING:
    from ase import Atoms

from ase import Atoms


class RedoxCal(BaseModel):
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

    molecule: Atoms | None = Field(
        None, description="The molecular structure for the calculations."
    )
    chg_mult: list[int] = Field(
        default_factory=lambda: [-1, 1, 0, 2, 1, 1, 0, 2], # * [-1, 1, 0, 2, 1, 1, 0, 2] for molecule with ion, [0, 1, 1, 2, 0, 1, -1, 2] for molecule without ion
        description="Charge multiplicities for the calculations.",
    )
    add_ion: bool = Field(
        True, description="Indicates whether to include an ion in the system."
    )
    ions: list[Atoms | str] = Field(
        default_factory=lambda: [PF6, "Li"],
        description="Ions involved in the calculations.",
    )
    protonate_ion_string: bool = Field(True, description="Indicates whether to protonate the ion.")
    label: str = Field("redox", description="A label for the calculations.")
    calc_kwards: dict = Field(
        default_factory=lambda: {
            "opt_xc": "b3lyp",
            "opt_basis": "6-311+G(d,p)",
            "sol_xc": "m062x",
            "sol_basis": "6-31G*",
            "solvent": "Acetone",
        },
        description="Keyword arguments for calculation methods.",
    )
    machine_kwards: dict = Field(
        default_factory=lambda: {"xtb_proc": 16},
        description="keyword arguments for Machine-specific.",
    )
    neutral_molecule: Atoms | None = None
    charged_molecule: Atoms | None = None

    class Config:
        arbitrary_types_allowed = True

    @field_validator("chg_mult")
    def set_chg_mult(cls, v):
        return v or [-1, 1, 0, 2, 1, 1, 0, 2]

    @field_validator("ions")
    def set_ions(cls, v):
        return v or [PF6, "Li"]

    @field_validator("calc_kwards")
    def set_calc_kwards(cls, v):
        return v or {
            "opt_xc": "b3lyp",
            "opt_basis": "6-311+G(d,p)",
            "sol_xc": "m062x",
            "sol_basis": "6-31G*",
            "solvent": "Acetone",
        }

    def prepare_ox(self):
        # * generate solvated molecules using ion and counter-ion
        if self.add_ion:
            logging.info("Generate and relax molecules clusters using crest")
            self.neutral_molecule = dock_atoms(
                self.molecule,
                dock=self.ions[0],
                crest_sampling=True,
                chg=self.chg_mult[0],
                mult=self.chg_mult[1],
            )
            self.charged_molecule = dock_atoms(
                self.molecule,
                dock=self.ions[0],
                crest_sampling=True,
                chg=self.chg_mult[2],
                mult=self.chg_mult[3],
            )
        else:
            logging.info("Relaxing molecules using crest")
            self.neutral_molecule = relax(
                self.molecule, chg=self.chg_mult[0], mult=self.chg_mult[1]
            )
            self.charged_molecule = relax(
                self.molecule, chg=self.chg_mult[2], mult=self.chg_mult[3]
            )

    def prepare_re(self):
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
                if protonate_ion_string is True and isinstance(ion, str)
                else None
            )
            if mol is None:
                logging.info("Protonation failed or skipped, trying docking")
                mol = dock_atoms(
                    molecule, dock=ion, crest_sampling=True, chg=chg, mult=mult
                )
            return mol

        if self.add_ion:
            logging.info("Generate and relax molecules clusters using crest")
            self.neutral_molecule = generate_molecule(
                self.molecule,
                self.ions[1],
                self.chg_mult[4],
                self.chg_mult[5],
                protonate_ion_string=self.protonate_ion_string,
            )
            self.charged_molecule = generate_molecule(
                self.molecule,
                self.ions[1],
                self.chg_mult[6],
                self.chg_mult[7],
                protonate_ion_string=self.protonate_ion_string,
            )
        else:
            logging.info("Relaxing molecules using crest")
            self.neutral_molecule = relax(
                self.molecule, chg=self.chg_mult[4], mult=self.chg_mult[5]
            )
            self.charged_molecule = relax(
                self.molecule, chg=self.chg_mult[6], mult=self.chg_mult[7]
            )
        logging.info("Molecule generation and relaxation complete")

    def cal_ox(self):
        """
        Calculates the oxidation potential of a molecular system.

        This function generates the neutral and charged molecules required for calculating the oxidation potential, either by docking with ions or relaxing the molecule, depending on the presence of ions. It then computes the redox potential based on these generated molecules.

        Args:
            None

        Returns:
            float: The calculated oxidation potential in eV.

        """

        # * calculate the oxidation state energies (in eV)
        logging.info("Calculating oxidation potential")
        redox_potential = RedoxPotential(
            neutral_molecule=self.neutral_molecule,
            charged_molecule=self.charged_molecule,
            chg_mult=self.chg_mult[:4],
            calc_type="ox",
            calc_kwards=self.calc_kwards,
        ).cal_cycle()
        logging.info(f"{self.label} oxidation potential: {redox_potential} eV")
        return redox_potential

    def cal_re(self):
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
        # * calculate the oxidation state energies (in eV)
        logging.info("Calculating reduction potential")
        redox_potential = RedoxPotential(
            neutral_molecule=self.neutral_molecule,
            charged_molecule=self.charged_molecule,
            chg_mult=self.chg_mult[4:8],
            calc_type="re",
            calc_kwards=self.calc_kwards,
        ).cal_cycle()
        logging.info(f"{self.label} reduction potential: {redox_potential} eV")
        return redox_potential

    def get_ox(self):
        """
        Calculates the oxidation potential of a molecular system.

        This function generates the neutral and charged molecules required for calculating the oxidation potential, either by docking with ions or relaxing the molecule, depending on the presence of ions. It then computes the redox potential based on these generated molecules.

        Args:
            None

        Returns:
            float: The calculated oxidation potential in eV.

        """
        self.prepare_ox()
        return self.cal_ox()

    def get_re(self):
        """
        Calculates the reduction potential of a molecular system.

        This function generates the neutral and charged molecules required for calculating the reduction potential, either by protonation or relaxation, depending on the presence of ions. It utilizes a helper function to streamline the process of generating molecules and logs the resulting reduction potential.

        Args:
            None

        Returns:
            float: The calculated reduction potential in eV.

        """
        self.prepare_re()
        return self.cal_re()

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


RedoxCal.model_rebuild()
