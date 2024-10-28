"""ASE pySCF calcuator, from scm"""
from __future__ import annotations

import logging
import warnings
from typing import Any, ClassVar

import ase
import ase.units
import numpy as np
import pyscf
import pyscf.dft
from ase.calculators.calculator import Calculator

"""Execute this file to see if it is working correctly."""


warnings.filterwarnings("ignore", "Since PySCF-2.3")


class PySCFCalculator(Calculator):
    """Example DFT calculator for PySCF. Only some minimal options are supported."""

    implemented_properties: ClassVar[list[str]] = ["energy", "forces", "dipole", "charges"]

    def __init__(
        self,
        xc: str = "pbe",  # https://pyscf.org/_modules/pyscf/dft/xcfun.html
        basis: str = "631g*",  # https://pyscf.org/_modules/pyscf/gto/basis.html
        symmetry: bool = True,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.xc = xc
        self.basis = basis
        self.symmetry = symmetry
        self.M = pyscf.M()
        self.dm1 = None  # density matrix
        self.ams_capabilities: Any = None

    @staticmethod
    def atoms2pyscf(atoms: ase.Atoms) -> list[tuple[str, tuple[float, float, float]]]:
        """Convert ASE Atoms to the list of symbols and coordinates required by pySCF"""
        return [
            (symbol, (xyz[0], xyz[1], xyz[2]))
            for symbol, xyz in zip(atoms.get_chemical_symbols(), atoms.get_positions())
        ]

    def calculate(self, atoms=None, properties=None, system_changes=None) -> None:
        """Perform DFT calculation, populate self.results with energy, forces, dipole, charge"""
        Calculator.calculate(
            self, atoms=atoms, properties=properties, system_changes=system_changes
        )
        self.results = {}
        if atoms is None:
            return

        pyscf_atoms = self.atoms2pyscf(atoms)
        pyscf_mol = pyscf.M(
            atom=pyscf_atoms,
            basis=self.basis,
            symmetry=self.symmetry,
        )

        self.mf = pyscf.dft.RKS(pyscf_mol)
        self.mf.xc = self.xc

        self.mf.kernel(verbose=0, dm0=self.dm1)
        self.dm1 = (
            self.mf.make_rdm1()
        )  # try to reuse density matrix for next calculation

        self.mf.analyze(verbose=0)
        _, charges = self.mf.mulliken_pop(verbose=0)

        self.results = {"energy": self.mf.energy_tot() * ase.units.Hartree}
        self.results["forces"] = (
            -np.array(self.mf.Gradients().kernel()) * ase.units.Hartree / ase.units.Bohr
        )
        self.results["dipole"] = (
            self.mf.dip_moment(verbose=0, unit="atomic") * ase.units.Bohr
        )
        self.results["charges"] = np.array(charges)


def get_calculator(**kwargs) -> PySCFCalculator:
    """Function used by the AMS-ASE interface"""
    from scm.amspipe import AMSExternalCapabilities

    calc = PySCFCalculator(**kwargs)
    calc.ams_capabilities = AMSExternalCapabilities()
    calc.ams_capabilities.apply_implemented_properties(calc.implemented_properties)
    return calc


if __name__ == "__main__":
    # Execute this file to see if the calculator is working correctly
    from ase.build import molecule as build_molecule

    atoms = build_molecule("H2O")

    atoms.calc = PySCFCalculator(
        xc="b3lyp",
        basis="631g*",
        symmetry=True,
    )
    logging.info(f"{atoms.get_chemical_symbols()}")
    logging.info(f"{atoms.get_positions()=}")
    logging.info(f"{atoms.get_potential_energy()=}")
    logging.info(f"{atoms.get_forces()=}")
