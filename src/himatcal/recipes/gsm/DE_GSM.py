from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

from ase.io import read, write
from pydantic import field_validator
from pyGSM.bin.gsm import cleanup_scratch, post_processing  # type: ignore
from pyGSM.coordinate_systems.delocalized_coordinates import (  # type: ignore
    DelocalizedInternalCoordinates,
)
from pyGSM.coordinate_systems.primitive_internals import (  # type: ignore
    PrimitiveInternalCoordinates,
)
from pyGSM.coordinate_systems.topology import Topology  # type: ignore
from pyGSM.growing_string_methods import DE_GSM  # type: ignore
from pyGSM.level_of_theories.ase import ASELoT  # type: ignore
from pyGSM.molecule import Molecule  # type: ignore
from pyGSM.optimizers.eigenvector_follow import eigenvector_follow  # type: ignore
from pyGSM.optimizers.lbfgs import lbfgs  # type: ignore
from pyGSM.potential_energy_surfaces import PES  # type: ignore
from pyGSM.utilities.manage_xyz import XYZ_WRITERS  # type: ignore

from himatcal.recipes.gsm.core import atoms2geom, gsm2atoms

if TYPE_CHECKING:
    from ase import Atoms
    from ase.calculators.calculator import Calculator


class ASE_DE_GSM:
    def __init__(self,
                    reactant: Atoms | str | Any = None,
                    product: Atoms | str | Any = None,
                    calculator: Calculator | Any = None,
                    multiplicity: int = 1,
                    fixed_reactant: bool | Any = False,
                    fixed_product: bool | Any = False,
                    coordinate_type: Literal["TRIC", "DLC", "HDLC"] | Any = "TRIC",
                    optimizer_method: Literal["eigenvector_follow", "lbfgs"] | Any = "eigenvector_follow",
                    num_of_nodes: int | Any = 11,
                    line_search: Literal["NoLineSearch", "backtrack"] | Any = "NoLineSearch",
                    conv_Ediff: float | Any = 100.0,
                    conv_gmax: float | Any = 100.0,
                    DMAX: float | Any = 0.1,
                    ID: int = 0,
                    r_type: Literal[0, 1, 2] | Any = 1,
                    max_gsm_iterations: int = 50,
                    max_opt_steps: int = 3):
            self.reactant = reactant
            self.product = product
            self.calculator = calculator
            self.multiplicity = multiplicity
            self.fixed_reactant = fixed_reactant
            self.fixed_product = fixed_product
            self.coordinate_type = coordinate_type
            self.optimizer_method = optimizer_method
            self.num_of_nodes = num_of_nodes
            self.line_search = line_search
            self.conv_Ediff = conv_Ediff
            self.conv_gmax = conv_gmax
            self.DMAX = DMAX
            self.ID = ID
            self.r_type = r_type
            self.max_gsm_iterations = max_gsm_iterations
            self.max_opt_steps = max_opt_steps

    # if reactant or product is str ,read atoms from file
    @field_validator("reactant", "product")
    def check_atoms(cls, values):
        reactant, product = values
        if isinstance(reactant, str):
            values["reactant"] = read(reactant)
        if isinstance(product, str):
            values["product"] = read(product)
        return values

    def union_bonds(self):
        for bond in self.product_topo.edges():
            if (
                bond in self.reactant_topo.edges()
                or (bond[1], bond[0]) in self.reactant_topo.edges()
            ):
                continue
            logging.info(f" Adding bond {bond} to reactant topology")
            if bond[0] > bond[1]:
                self.reactant_topo.add_edge(bond[0], bond[1])
            else:
                self.reactant_topo.add_edge(bond[1], bond[0])

    def run(self):
        # * 0. convert atoms to geom
        self.reactant_gemos = atoms2geom(self.reactant)
        self.product_gemos = atoms2geom(self.product)

        # * 1. build the LoT
        logging.info("Building the LoT")
        self.lot = ASELoT.from_options(self.calculator, geom=self.reactant_gemos[2])

        # * 2. build the PES
        logging.info("Building the PES")
        self.pes = PES.from_options(
            lot=self.lot, ad_idx=0, multiplicity=self.multiplicity
        )

        # * 3. build the topology
        logging.info("Building the topologies")
        self.reactant_topo = Topology.build_topology(
            xyz=self.reactant_gemos[1], atoms=self.reactant_gemos[0]
        )
        self.product_topo = Topology.build_topology(
            xyz=self.product_gemos[1], atoms=self.product_gemos[0]
        )
        self.union_bonds()

        # * 4. build the primitice internal coordinates
        logging.info("Building Primitive Internal Coordinates")
        self.reactant_prim = PrimitiveInternalCoordinates.from_options(
            xyz=self.reactant_gemos[1],
            atoms=self.reactant_gemos[0],
            topology=self.reactant_topo,
            connect=self.coordinate_type == "DLC",
            addtr=self.coordinate_type == "TRIC",
            addcart=self.coordinate_type == "HDLC",
        )
        self.product_prim = PrimitiveInternalCoordinates.from_options(
            xyz=self.product_gemos[1],
            atoms=self.product_gemos[0],
            topology=self.product_topo,
            connect=self.coordinate_type == "DLC",
            addtr=self.coordinate_type == "TRIC",
            addcart=self.coordinate_type == "HDLC",
        )
        # * 4.1. add product coords to reactant coords
        self.reactant_prim.add_union_primitives(self.product_prim)

        # * 5. build the delocalized internal coordinates
        logging.info("Building Delocalized Internal Coordinates")
        self.deloc_coords_reactant = DelocalizedInternalCoordinates.from_options(
            xyz=self.reactant_gemos[1],
            atoms=self.reactant_gemos[0],
            connect=self.coordinate_type == "DLC",
            addtr=self.coordinate_type == "TRIC",
            addcart=self.coordinate_type == "HDLC",
            primitives=self.reactant_prim,
        )

        # * 6. build the molecule
        logging.info("Building Molecule")
        self.molecule_reactant = Molecule.from_options(
            geom=self.reactant_gemos[2],
            PES=self.pes,
            coord_obj=self.deloc_coords_reactant,
            Form_Hessian=self.optimizer_method == "eigenvector_follow",
        )
        self.molecule_product = Molecule.copy_from_options(
            self.molecule_reactant,
            xyz=self.product_gemos[1],
            new_node_id=self.num_of_nodes - 1,
            copy_wavefunction=False,
        )

        # * 7. create the optimizer
        logging.info("Creating optimizer")
        opt_options = {
            "print_level": 1,
            "Linesearch": self.line_search,
            "update_hess_in_bg": False,
            "conv_Ediff": self.conv_Ediff,
            "conv_gmax": self.conv_gmax,
            "DMAX": self.DMAX,
            "opt_climb": self.r_type in [1, 2],
        }
        if self.optimizer_method == "eigenvector_follow":
            self.optimizer = eigenvector_follow.from_options(**opt_options)
        elif self.optimizer_method == "lbfgs":
            self.optimizer = lbfgs.from_options(**opt_options)
        else:
            raise NotImplementedError

        # * 7.1 optimize reactant and product if needed
        if not self.fixed_reactant:
            path = str(Path.cwd() / "scratch" / f"{self.ID:03}" / "0")
            self.optimizer.optimize(
                molecule=self.molecule_reactant,
                refE=self.molecule_reactant.energy,
                opt_steps=100,
                path=path,
            )
        if not self.fixed_product:
            path = str(
                Path.cwd() / "scratch" / f"{self.ID:03}" / str(self.num_of_nodes - 1)
            )
            self.optimizer.optimize(
                molecule=self.molecule_product,
                refE=self.molecule_product.energy,
                opt_steps=100,
                path=path,
            )

        # * 8. build the GSM
        logging.info("Building the GSM object")
        self.gsm = DE_GSM.from_options(
            reactant=self.molecule_reactant,
            product=self.molecule_product,
            nnodes=self.num_of_nodes,
            CONV_TOL=0.0005,
            CONV_gmax=self.conv_gmax,
            CONV_Ediff=self.conv_Ediff,
            ADD_NODE_TOL=0.1,
            growth_direction=0,
            optimizer=self.optimizer,
            ID=self.ID,
            print_level=1,
            mp_cores=1,
            interp_method="DLC",
            xyz_writer=XYZ_WRITERS["multixyz"]
        )

        # * 9. run the GSM
        logging.info("Main GSM Calculation")
        self.gsm.go_gsm(
            max_iters=self.max_gsm_iterations,
            opt_steps=self.max_opt_steps,
            rtype=self.r_type,
        )

        # * 10. write the results into an extended xyz file
        string_ase, ts_ase = gsm2atoms(self.gsm)
        write(f"opt_converged_{self.gsm.ID:03d}_ase.xyz", string_ase)
        write(f"TSnode_{self.gsm.ID}.xyz", ts_ase)

        # * 11. post processing
        logging.info("Post processing")
        post_processing(self.gsm, have_TS=True)

        # * 12. cleanup
        logging.info("Cleaning up")
        cleanup_scratch(self.gsm.ID)
