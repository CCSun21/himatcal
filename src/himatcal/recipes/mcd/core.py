from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from himatcal.utils.mcd import chem, mcd

if TYPE_CHECKING:
    from ase import Atoms


class MCD_runner:
    def __init__(
        self,
        atoms: Atoms,
        chg: int = 0,
        mult: int = 1,
        driving_coords: list | None = None,  # atoms' index, target bondlength, step, [[1, 8, 2.1, 5]]
        qcsoft: str = "gaussian",
        command: str = "g16",
        calc_kwargs: str = "#N B3LYP/6-311+g* em=GD3BJ scf(xqc) scrf(iefpcm, solvent=acetone)",
        num_relaxation: int = 5,
        step_size: float = 0.05,
        working_directory="",
        output_directory="",
        use_hessian: bool = False,
        hessian_update: str = "exact",
        unit: str = "Hartree",
    ):
        self.atoms = atoms
        self.chg = chg
        self.mult = mult
        self.driving_coords = driving_coords
        self.qcsoft = qcsoft
        self.command = command
        self.calc_kwargs = calc_kwargs
        self.num_relaxation = num_relaxation
        self.use_hessian = use_hessian
        self.hessian_update = hessian_update
        self.step_size = step_size
        if output_directory == "":
            output_directory = str(Path.cwd())
        if working_directory == "":
            working_directory = output_directory
        self.working_directory = working_directory
        self.output_directory = output_directory
        self.unit = unit

    def read_reactant(self):
        atom_list = []
        for a in self.atoms:
            atom = chem.Atom(a.symbol)
            atom.x = a.x
            atom.y = a.y
            atom.z = a.z
            atom_list.append(atom)
        reactant = chem.Molecule()
        reactant.atom_list = atom_list
        reactant.chg = self.chg
        reactant.multiplicity = self.mult
        return reactant

    def dirving_constraints(self):
        constraints = {}
        num_steps = {}
        for DC in self.driving_coords:
            constraint = (DC[0] - 1, DC[1] - 1)
            constraints[constraint] = DC[2]
            num_steps[constraint] = DC[3]
        return constraints, num_steps

    def get_calculator(self):
        if self.qcsoft == "gaussian":
            from himatcal.calculator.gaussian_mcd import Gaussian

            calculator = Gaussian(self.command)
        elif self.qcsoft == "orca":
            from himatcal.calculator.orca_mcd import Orca

            calculator = Orca(self.command)
        else:
            print(
                f"Wrong calculator (={self.qcsoft}) is given! Check the option file !!!"
            )
            calculator = None
            return calculator
        calculator.content = self.calc_kwargs
        # basis_file = os.path.join(args.input_directory,'basis') # For Effective Core Potential
        # if os.path.exists(basis_file):
        #     calculator.load_basis(basis_file)
        return calculator

    def run_MCD(self):
        reactant = self.read_reactant()
        calculator = self.get_calculator()
        constraints, num_steps = self.dirving_constraints()
        scanner = mcd.MCD(num_relaxation=self.num_relaxation, calculator=calculator)
        scanner.use_hessian = self.use_hessian
        scanner.hessian_update = self.hessian_update
        scanner.step_size = self.step_size
        scanner.log_directory = self.output_directory
        scanner.change_working_directory(self.working_directory)
        pathway = scanner.scan(
            reactant,
            constraints,
            num_steps,
            chg=self.chg,
            multiplicity=self.mult,
        )
        return pathway
