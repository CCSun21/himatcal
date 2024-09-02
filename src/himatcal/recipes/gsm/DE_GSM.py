from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Annotated, Literal, Optional

from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.io import read, write
from pydantic import BaseModel, Field, SkipValidation, field_validator
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

from himatcal.recipes.gsm.core import atoms2geom, gsm2atoms


class ASE_DE_GSM(BaseModel):
    """
    Represents a DE_GSM (Double-ended Geometric Surface Mapping) model for simulating molecular transformations. This class encapsulates the necessary attributes for defining the reactant and product molecules, along with their associated computational settings.

    The DE_GSM class inherits from BaseModel and includes fields for specifying the reactant and product as either Atoms objects or file paths. It also allows for the inclusion of an ASE calculator and options to fix the geometries of the reactant and product during simulations.

    Attributes:
        reactant (Atoms | str): Reactant Atoms object or path to a file.
        product (Atoms | str): Product Atoms object or path to a file.
        calculator (Calculator | None): ASE calculator for performing calculations.
        fixed_reactant (bool): Indicates whether to fix the reactant geometry.
        fixed_product (bool): Indicates whether to fix the product geometry.

    Examples:
        To create a DE_GSM model, instantiate the class with the required parameters:
        >>> model = DE_GSM(reactant='path/to/reactant.xyz', product='path/to/product.xyz')
    """


    reactant: SkipValidation[Atoms | str] = Field(
        None, description="Reactant Atoms object or path to a file"
    )
    product: Atoms | str = Field(
        None, description="Product Atoms object or path to a file"
    )
    calculator: Calculator | None = Field(None, description="ASE calculator")
    multiplicity: int = Field(1, ge=1, description="Multiplicity")

    # * GSM options
    fixed_reactant: bool | None = Field(False, description="Fix the reactant geometry")
    fixed_product: bool | None = Field(False, description="Fix the product geometry")
    coordinate_type: Literal["TRIC", "DLC", "HDLC"] | None = Field(
        "TRIC", description="Coordinate type"
    )
    optimizer_method: Literal["eigenvector_follow", "lbfgs"] | None = Field(
        "eigenvector_follow", description="Optimizer method"
    )
    num_of_nodes: int | None = Field(11, ge=2, description="Number of nodes")
    line_search: Literal["NoLineSearch", "backtrack"] | None = Field(
        "NoLineSearch", description="Line search method"
    )
    conv_Ediff: float | None = Field(100.0, description="Energy difference convergence")
    conv_gmax: float | None = Field(100.0, description="Max grad rms threshold")
    DMAX: float | None = Field(0.1, description="Step size cap")
    ID: int = Field(0, description="ID")
    r_type: Literal[0, 1, 2] = Field(
        1, description="0 for no climb, 1 for only climb, 2 for climb and find TS"
    )
    max_gsm_iterations: int = Field(100, description="Max GSM iterations")
    max_opt_steps: int = Field(3, description="Max optimization steps")

    class Config:
        arbitrary_types_allowed = True
        validate_assignment = False
        validate_all = False

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
        )

        # * 9. run the GSM
        logging.info("Main GSM Calculation")
        self.gsm.go_gsm(
            max_gsm_iterations=self.max_gsm_iterations,
            max_opt_steps=self.max_opt_steps,
            rtype=self.r_type,
        )

        # * 10. write the results into an extended xyz file
        string_ase, ts_ase = gsm2atoms(self.gsm)
        write(f"opt_converged_{self.gsm.ID:03d}_ase.xyz", string_ase)
        write(f"TSnode_{self.gsm.ID}.xyz", string_ase)

        # * 11. post processing
        logging.info("Post processing")
        post_processing(self.gsm, have_TS=True)

        # * 12. cleanup
        logging.info("Cleaning up")
        cleanup_scratch(self.gsm.ID)

ASE_DE_GSM.model_rebuild()
