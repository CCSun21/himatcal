from __future__ import annotations

import logging
import os

from ase.io import write
from pyGSM.bin.gsm import post_processing  # type: ignore
from pyGSM.coordinate_systems import (  # type: ignore
    DelocalizedInternalCoordinates,
    Distance,
    PrimitiveInternalCoordinates,
    Topology,
)
from pyGSM.growing_string_methods import SE_GSM  # type: ignore
from pyGSM.level_of_theories.ase import ASELoT  # type: ignore
from pyGSM.molecule import Molecule  # type: ignore
from pyGSM.optimizers import eigenvector_follow  # type: ignore
from pyGSM.potential_energy_surfaces import PES  # type: ignore
from pyGSM.utilities import manage_xyz, nifty  # type: ignore
from pyGSM.utilities.cli_utils import get_driving_coord_prim  # type: ignore

from himatcal.recipes.gsm.core import atoms2geom, gsm2atoms


class ASE_SE_GSM:
    def __init__(
        self,
        atoms,
        driving_coords,
        multiplicity=1,
        calculator=None,
        cleanup_scratch=False,
    ):
        """
        Initializes the class with the specified atom and driving coordinates.

        This constructor sets up the parameters necessary for the class, including the atom to be used, the driving coordinates for simulations, and an optional calculator for performing calculations. It also allows for the option to clean up temporary files created during the process.

        Args:
            atom: The atom to be used in the calculations.
            driving_coords: A list of driving coordinates, formatted as [["BREAK", 2, 3]].
            calculator: An optional calculator instance; if not provided, a default XTB calculator is used.
            cleanup_scratch: A boolean indicating whether to clean up scratch files after computations.

        """

        self.atoms = atoms
        self.driving_coords = driving_coords  # List: driving_coords = [["BREAK", 2, 3]]
        self.multiplicity = multiplicity
        if calculator is None:
            from xtb_ase import XTB

            calculator = XTB()
        self.calculator = calculator
        self.cleanup_scratch = cleanup_scratch

    # * 1. Build the LoT
    def build_lot(self):
        """
        build ase lot from calculator
        """
        nifty.printcool(" Building the LOT")
        self.lot = ASELoT.from_options(self.calculator, geom=self.geom)

    # * 2. Build the PES
    def build_pes(self):
        nifty.printcool(" Building the PES")
        self.pes = PES.from_options(
            lot=self.lot,
            ad_idx=0,  # Adiabatic index (default: 0)
            multiplicity=self.multiplicity,
        )

    # * 3. Build the topology
    def build_topology(self):
        # * build the topology
        self.top = Topology.build_topology(
            self.xyz,
            self.atom,
        )
        # * add the driving coordinates to the topology
        driving_coord_prims = []
        for dc in self.driving_coords:
            prim = get_driving_coord_prim(dc)
            if prim is not None:
                driving_coord_prims.append(prim)

        for prim in driving_coord_prims:
            if type(prim) is Distance:
                bond = (prim.atoms[0], prim.atoms[1])
                if (
                    bond not in self.top.edges
                    and (bond[1], bond[0]) not in self.top.edges()
                ):
                    logging.info(f" Adding bond {bond} to top1")
                    self.top.add_edge(bond[0], bond[1])

    # * 4. Build the primitive internal coordinates
    def build_primitives(self):
        nifty.printcool("Building Primitive Internal Coordinates")
        self.p1 = PrimitiveInternalCoordinates.from_options(
            xyz=self.xyz,
            atoms=self.atom,
            addtr=True,  # Add TRIC
            topology=self.top,
        )

    # * 5. Build the delocalized internal coordinates
    def build_delocalized_coords(self):
        nifty.printcool("Building Delocalized Internal Coordinates")
        self.coord_obj1 = DelocalizedInternalCoordinates.from_options(
            xyz=self.xyz,
            atoms=self.atom,
            addtr=True,  # Add TRIC
            primitives=self.p1,
        )

    # * 6. Build the molecule
    def build_molecule(self):
        nifty.printcool("Building Molecule")
        self.reactant = Molecule.from_options(
            geom=self.geom,
            PES=self.pes,
            coord_obj=self.coord_obj1,
            Form_Hessian=True,
        )

    # * 7. Create the optimizer
    def create_optimizer(self):
        nifty.printcool("Creating optimizer")
        self.optimizer = eigenvector_follow.from_options(
            Linesearch="backtrack",
            OPTTHRESH=0.0005,
            DMAX=0.5,
            abs_max_step=0.5,
            conv_Ediff=0.1,
            opt_climb=True,
        )

    # * 8. Optimize the reactant
    def optimize_reactant(self):
        nifty.printcool(f"initial energy is {self.reactant.energy:5.4f} kcal/mol")
        nifty.printcool("REACTANT GEOMETRY NOT FIXED!!! OPTIMIZING")
        self.optimizer.optimize(
            molecule=self.reactant,
            refE=self.reactant.energy,
            opt_steps=50,
        )

    # * 9. Run the GSM
    def run_gsm(self):
        self.gsm = SE_GSM.from_options(
            reactant=self.reactant,
            nnodes=30,
            optimizer=self.optimizer,
            xyz_writer=manage_xyz.write_std_multixyz,
            driving_coords=self.driving_coords,
            DQMAG_MAX=1.5,  # * default value is 0.8, Maximum step size in single-ended mode
            DQMAG_MIN=0.2,  # * default value is 0.8, Maximum step size in single-ended mode
            ADD_NODE_TOL=0.1,  # * default value is 0.1, for GSM, Convergence tolerance for adding new node
            CONV_TOL=0.0005,  # * Convergence tolerance for optimizing nodes
        )
        self.gsm.go_gsm(max_iters=50, opt_steps=10, rtype=2)

    # * 10. Clean up the scratch directory
    def clean_scratch(self):
        if self.cleanup_scratch:
            cmd = f"rm scratch/growth_iters_{self.gsm.ID:03d}_*.xyz"
            os.system(cmd)
            cmd = f"rm scratch/opt_iters_{self.gsm.ID:03d}_*.xyz"
            os.system(cmd)

    def run(self):
        self.atom, self.xyz, self.geom = atoms2geom(self.atoms)
        self.build_lot()
        self.build_pes()
        self.build_topology()
        self.build_primitives()
        self.build_delocalized_coords()
        self.build_molecule()
        self.create_optimizer()
        self.optimize_reactant()
        self.run_gsm()
        post_processing(self.gsm, analyze_ICs=False, have_TS=True)

        # * 10. write the results into an extended xyz file
        string_ase, ts_ase = gsm2atoms(self.gsm)
        write(f"opt_converged_{self.gsm.ID:03d}_ase.xyz", string_ase)
        write(f"TSnode_{self.gsm.ID}.xyz", ts_ase)

        self.clean_scratch()
