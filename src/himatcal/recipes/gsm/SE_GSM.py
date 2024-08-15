from __future__ import annotations

import logging
import os
from pathlib import Path

import numpy as np
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
from pyGSM.utilities import elements, manage_xyz, nifty  # type: ignore
from pyGSM.utilities.cli_utils import get_driving_coord_prim  # type: ignore
from pyGSM.utilities.cli_utils import plot as gsm_plot  # type: ignore


class ASE_SE_GSM:
    def __init__(self, atom, driving_coords, calculator=None, cleanup_scratch=False):
        """
        Initializes the class with the specified atom and driving coordinates.

        This constructor sets up the parameters necessary for the class, including the atom to be used, the driving coordinates for simulations, and an optional calculator for performing calculations. It also allows for the option to clean up temporary files created during the process.

        Args:
            atom: The atom to be used in the calculations.
            driving_coords: A list of driving coordinates, formatted as [["BREAK", 2, 3]].
            calculator: An optional calculator instance; if not provided, a default XTB calculator is used.
            cleanup_scratch: A boolean indicating whether to clean up scratch files after computations.

        """

        self.atom = atom
        self.driving_coords = driving_coords  # List: driving_coords = [["BREAK", 2, 3]]
        if calculator is None:
            from xtb.ase.calculator import XTB

            calculator = XTB()
        self.calculator = calculator
        self.cleanup_scratch = cleanup_scratch

    def atom2geom(self):
        xyz = self.atom[0].positions
        geom = np.column_stack([self.atom[0].symbols, xyz]).tolist()
        ELEMENT_TABLE = elements.ElementData()
        atoms = [ELEMENT_TABLE.from_symbol(atom) for atom in self.atom[0].symbols]
        return atoms, xyz, geom

    # * 1. Build the LOT
    def build_lot(self):
        nifty.printcool(" Building the LOT")
        self.lot = ASELoT.from_options(self.calculator, geom=self.geom)

    # * 2. Build the PES
    def build_pes(self):
        nifty.printcool(" Building the PES")
        self.pes = PES.from_options(
            lot=self.lot,
            ad_idx=0,
            multiplicity=1,
        )

    # * 3. Build the topology
    def build_topology(self):
        # * build the topology
        self.top = Topology.build_topology(
            self.xyz,
            self.atoms,
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
            atoms=self.atoms,
            addtr=True,  # Add TRIC
            topology=self.top,
        )

    # * 5. Build the delocalized internal coordinates
    def build_delocalized_coords(self):
        nifty.printcool("Building Delocalized Internal Coordinates")
        self.coord_obj1 = DelocalizedInternalCoordinates.from_options(
            xyz=self.xyz,
            atoms=self.atoms,
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

    def run_gsm(self):
        self.gsm = SE_GSM.from_options(
            reactant=self.reactant,
            nnodes=20,
            optimizer=self.optimizer,
            xyz_writer=manage_xyz.write_std_multixyz,
            driving_coords=self.driving_coords,
            DQMAG_MAX=0.5,  # default value is 0.8
            ADD_NODE_TOL=0.01,  # default value is 0.1
            CONV_TOL=0.0005,
        )
        self.gsm.set_V0()
        self.gsm.nodes[0].gradrms = 0.0
        self.gsm.nodes[0].V0 = self.gsm.nodes[0].energy
        logging.info(f" Initial energy is {self.gsm.nodes[0].energy:1.4f}")
        self.gsm.add_GSM_nodeR()
        self.gsm.grow_string(max_iters=50, max_opt_steps=10)
        if self.gsm.tscontinue:
            self.gsm.pastts = self.gsm.past_ts()
            logging.info(f"pastts {self.gsm.pastts}")
            try:
                if self.gsm.pastts == 1:  # normal over the hill
                    self.gsm.add_GSM_nodeR(1)
                    self.gsm.add_last_node(2)
                elif self.gsm.pastts in [2, 3]:  # when cgrad is positive
                    self.gsm.add_last_node(1)
                    if (
                        self.gsm.nodes[self.gsm.nR - 1].gradrms
                        > 5.0 * self.gsm.options["CONV_TOL"]
                    ):
                        self.gsm.add_last_node(1)
            except Exception:
                logging.info("Failed to add last node, continuing.")
        self.gsm.nnodes = self.gsm.nR
        self.gsm.nodes = self.gsm.nodes[: self.gsm.nR]
        self.energies = self.gsm.energies
        if self.gsm.TSnode == self.gsm.nR - 1:
            logging.info(" The highest energy node is the last")
            logging.info(" not continuing with TS optimization.")
            self.gsm.tscontinue = False
        logging.info(f" Number of nodes is {self.gsm.nnodes}")
        logging.info(" Warning last node still not optimized fully")
        self.gsm.xyz_writer(
            f"grown_string_{self.gsm.ID:03}.xyz",
            self.gsm.geometries,
            self.gsm.energies,
            self.gsm.gradrmss,
            self.gsm.dEs,
        )
        logging.info(" SSM growth phase over")
        self.gsm.done_growing = True

    def post_process(self):
        gsm_plot(self.gsm.energies, x=range(len(self.gsm.energies)), title=0)
        ICs = [self.gsm.nodes[0].primitive_internal_coordinates]

    def post_processing(self, analyze_ICs=False, have_TS=True):
        gsm_plot(
            fx=self.gsm.energies, x=range(len(self.gsm.energies)), title=self.gsm.ID
        )

        ICs = [self.gsm.nodes[0].primitive_internal_coordinates]
        # TS energy
        if have_TS:
            minnodeR = np.argmin(self.gsm.energies[: self.gsm.TSnode])
            TSenergy = self.gsm.energies[self.gsm.TSnode] - self.gsm.energies[minnodeR]
            logging.info(f" TS energy: {TSenergy:5.4f}")
            logging.info(
                f" absolute energy TS node {self.gsm.nodes[self.gsm.TSnode].energy:5.4f}"
            )
            minnodeP = self.gsm.TSnode + np.argmin(self.gsm.energies[self.gsm.TSnode :])
            logging.info(
                " min reactant node: %i, min product node %i, TS node is %i"
                % (minnodeR, minnodeP, self.gsm.TSnode)
            )
            # write TS node
            self.gsm.xyz_writer(
                f"TS_{self.gsm.ID:03d}.xyz",
                [self.gsm.geometries[self.gsm.TSnode]],
                [self.gsm.energies[self.gsm.TSnode]],
                [self.gsm.gradrmss[self.gsm.TSnode]],
                [self.gsm.dEs[self.gsm.TSnode]],
            )

            ICs.extend(
                (
                    self.gsm.nodes[minnodeR].primitive_internal_values,
                    self.gsm.nodes[self.gsm.TSnode].primitive_internal_values,
                    self.gsm.nodes[minnodeP].primitive_internal_values,
                )
            )
            with Path.open(Path(f"IC_data_{self.gsm.ID:04d}.txt"), "w") as f:
                f.write(
                    f"Internals \t minnodeR: {minnodeR} \t TSnode: {self.gsm.TSnode} \t minnodeP: {minnodeP}\n"
                )
                for x in zip(*ICs):
                    f.write("{}\t{}\t{}\t{}\n".format(*x))

        else:
            minnodeR = 0
            minnodeP = self.gsm.nR
            logging.info(
                f" absolute energy end node {self.gsm.nodes[self.gsm.nR].energy:5.4f}"
            )
            logging.info(
                f" difference energy end node {self.gsm.nodes[self.gsm.nR].difference_energy:5.4f}"
            )
            ICs.extend(
                (
                    self.gsm.nodes[minnodeR].primitive_internal_values,
                    self.gsm.nodes[minnodeP].primitive_internal_values,
                )
            )
            with Path.open(Path(f"IC_data_{self.gsm.ID}.txt"), "w") as f:
                f.write(f"Internals \t Beginning: {minnodeR} \t End: {self.gsm.TSnode}")
                for x in zip(*ICs):
                    f.write("{}\t{}\t{}\n".format(*x))

        # Delta E
        deltaE = self.gsm.energies[minnodeP] - self.gsm.energies[minnodeR]
        logging.info(f" Delta E is {deltaE:5.4f}")

    def clean_scratch(self):
        if self.cleanup_scratch:
            cmd = f"rm scratch/growth_iters_{self.gsm.ID:03d}_*.xyz"
            os.system(cmd)
            cmd = f"rm scratch/opt_iters_{self.gsm.ID:03d}_*.xyz"
            os.system(cmd)

    def run(self):
        self.atoms, self.xyz, self.geom = self.atom2geom()
        self.build_lot()
        self.build_pes()
        self.build_topology()
        self.build_primitives()
        self.build_delocalized_coords()
        self.build_molecule()
        self.create_optimizer()
        self.optimize_reactant()
        self.run_gsm()
        self.post_processing()
        self.clean_scratch()
