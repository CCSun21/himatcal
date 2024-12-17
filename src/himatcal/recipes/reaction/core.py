from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from himatcal.recipes.orca.core import neb_job

if TYPE_CHECKING:
    from himatcal.recipes.reaction._base import MolGraph


def orca_neb(
    reactant: MolGraph,
    product: MolGraph,
    chg: int,
    mult: int,
    ts: MolGraph | None = None,
    n_images: int = 5,
    preopt: bool = False,
    guessTS: bool = False,
):
    Product_dir = Path("Product")
    Product_dir.mkdir(exist_ok=True)
    if product and product.atoms is not None:
        product.atoms.write(Product_dir / "product.xyz")
    if ts and ts.atoms is not None:
        ts.atoms.write(Product_dir / "guessts.xyz")

    return neb_job(
        atoms=reactant.atoms,
        charge=chg,
        spin_multiplicity=mult,
        xc="b97-3c",
        basis="def2-SVP",
        n_images=n_images,
        preopt=preopt,
        guessTS=guessTS,
        copy_files=f"{Product_dir.resolve()}",
    )


def plot_rxn_smiles(reaction_smi, useSmiles=True):
    from rdkit.Chem import Draw, rdChemReactions

    reaction = rdChemReactions.ReactionFromSmarts(reaction_smi, useSmiles=useSmiles)
    logging.info(f"input reaction: {reaction_smi}")
    logging.info(
        f"Rdkit formmated reaction: {rdChemReactions.ReactionToSmiles(reaction)}"
    )
    return Draw.ReactionToImage(reaction)
