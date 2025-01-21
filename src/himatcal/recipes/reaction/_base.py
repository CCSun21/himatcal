"""reaction calculations"""

from __future__ import annotations

import logging

from ase import Atoms
from pydantic import BaseModel
from pymatgen.io.ase import MSONAtoms

Atoms.as_dict = MSONAtoms.as_dict  # type: ignore[attr-defined]
Atoms.from_dict = MSONAtoms.from_dict  # type: ignore[attr-defined]


class MolGraph(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    atoms: Atoms | None = None  # ASE Atoms object
    smiles: str | None = None
    energy: float | None = None
    state: str = "init"  # init, opt, ts, final
    label: str = "mol"

    def __init__(
        self,
        smiles: str | None = None,
        atoms: Atoms | None = None,
        energy: float | None = None,
        state: str = "init",
        label: str = "mol",
    ):
        super().__init__(
            smiles=smiles, atoms=atoms, energy=energy, state=state, label=label
        )
        if self.smiles and self.atoms:
            logging.info("Both smiles and atoms provided, using atoms.")
        elif self.smiles:
            self.smiles2atoms()
        elif self.atoms:
            self.atoms2smiles()

    def smiles2atoms(self):
        if self.smiles and self.atoms is None:
            from himatcal.utils.rdkit.core import smiles2atoms

            self.atoms = smiles2atoms(self.smiles)
        return self.atoms

    def atoms2smiles(self):
        if self.atoms and self.smiles is None:
            from himatcal.utils.rdkit.core import atoms2smiles

            self.smiles = atoms2smiles(self.atoms)
        return self.smiles

    def show_atoms(self):
        from nglview import show_ase

        return show_ase(self.atoms)

    def show_mol(self):
        from rdkit import Chem
        from rdkit.Chem import Draw

        mol = Chem.MolFromSmiles(self.smiles)
        return Draw.MolToImage(mol)

    def to_json(self):
        from monty.json import jsanitize

        return jsanitize(self.model_dump(), enum_values=True, recursive_msonable=True)

    def to_file(self, filename):
        from monty.serialization import dumpfn

        dumpfn(self.to_json(), filename)

    @classmethod
    def from_json(cls, data):
        from monty.json import MontyDecoder

        return cls(**MontyDecoder().process_decoded(data))

    @classmethod
    def from_file(cls, filename, smiles=None):
        from monty.serialization import loadfn

        data = loadfn(filename)
        if not smiles:
            return cls.from_json(data)
        try:
            moldata = next(item for item in data if item["smiles"] == smiles)
            return cls.from_json(moldata)
        except StopIteration as e:
            raise ValueError(
                f"No molecule with smiles '{smiles}' found in the file."
            ) from e


class Reaction(BaseModel):
    reactant: MolGraph | list[MolGraph]
    product: MolGraph | list[MolGraph]
    ts: MolGraph | None = None

    @property
    def reactant_energy(self):
        if not isinstance(self.reactant, list):
            return self.reactant.energy
        energies = []
        for mol in self.reactant:
            if mol.energy is None:
                logging.warning(f"Reactant molecule {mol.label} has no energy.")
            else:
                energies.append(mol.energy)
        return sum(energies)

    @property
    def product_energy(self):
        if not isinstance(self.product, list):
            return self.product.energy
        energies = []
        for mol in self.product:
            if mol.energy is None:
                logging.warning(f"Product molecule {mol.label} has no energy.")
            else:
                energies.append(mol.energy)
        return sum(energies)

    @property
    def barrier(self):
        if self.ts is None:
            return None
        else:
            return (
                self.ts.energy - self.reactant_energy
                if self.ts.energy is not None and self.reactant_energy is not None
                else None
            )

    @property
    def enthalpy(self):
        if self.product_energy is not None and self.reactant_energy is not None:
            return self.product_energy - self.reactant_energy
        else:
            return None

    def reverse(self):
        return Reaction(reactant=self.product, product=self.reactant, ts=self.ts)

    @property
    def reactant_smiles(self):
        if not isinstance(self.reactant, list):
            return self.reactant.smiles
        smiles = []
        for mol in self.reactant:
            if mol.smiles is None:
                logging.warning(f"Reactant molecule {mol.label} has no smiles.")
            else:
                smiles.append(mol.smiles)
        return ".".join(smiles)

    @property
    def product_smiles(self):
        if not isinstance(self.product, list):
            return self.product.smiles
        smiles = []
        for mol in self.product:
            if mol.smiles is None:
                logging.warning(f"Product molecule {mol.label} has no smiles.")
            else:
                smiles.append(mol.smiles)
        return ".".join(smiles)

    @property
    def reaction_results(self):
        return ReactionResults(
            rsmi=self.reactant_smiles or "",
            psmi=self.product_smiles or "",
            ea=self.barrier,
            dh=self.enthalpy,
        )


class ReactionResults(BaseModel):
    rsmi: str
    psmi: str
    ea: float | None = None
    dh: float | None = None

    def __str__(self):
        return f"{self.rsmi}>>{self.psmi}"

    def reverse(self):
        return ReactionResults(
            rsmi=self.psmi, psmi=self.rsmi, ea=self.ea - self.dh, dh=-self.dh
        )

    @property
    def rxn_smi(self):
        return str(self)

    @property
    def reactants(self):
        # spilt the reactants by the . character
        return self.rsmi.split(".")

    @property
    def products(self):
        # spilt the products by the . character
        return self.psmi.split(".")
