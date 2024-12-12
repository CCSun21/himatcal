from __future__ import annotations

import logging


def plot_rxn_smiles(reaction_smi, useSmiles=True):
    from rdkit.Chem import Draw, rdChemReactions

    reaction = rdChemReactions.ReactionFromSmarts(reaction_smi, useSmiles=useSmiles)
    logging.info(f"input reaction: {reaction_smi}")
    logging.info(
        f"Rdkit formmated reaction: {rdChemReactions.ReactionToSmiles(reaction)}"
    )
    return Draw.ReactionToImage(reaction)
