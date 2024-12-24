from __future__ import annotations

import json
import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from himatcal.recipes.reaction import MolGraph


def update_molgraph(molgraph: MolGraph, filename: str = "molgraph.json"):
    """
    update the molgraph in the JSON file
    """
    with open(filename) as json_file:
        content = json_file.read()
        molgraph_list = [] if content == "" else json.loads(content)

    # update the molgraph with same smiles
    mol_json = molgraph.to_json()
    for i, mol in enumerate(molgraph_list):
        if mol["smiles"] == molgraph.smiles:
            molgraph_list[i] = mol_json
            break
    else:
        molgraph_list.append(mol_json)

    # Save the list of dictionaries to a JSON file
    with open(filename, "w") as json_file:
        json.dump(molgraph_list, json_file, indent=4)


def get_charge_and_spin(smiles):
    from rdkit import Chem

    # 从SMILES字符串创建分子对象
    mol = Chem.MolFromSmiles(smiles)

    # 计算分子的电荷
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())

    # 计算分子的自旋多重度
    num_radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    spin_multiplicity = num_radicals + 1

    return charge, spin_multiplicity


def molgraph_relax(molgraph: MolGraph, charge: int | None, mult: int | None, method: str="orca"):

    if charge is None or mult is None:
        charge, mult = get_charge_and_spin(molgraph.smiles)
    if method == "orca":
        from quacc.recipes.orca.core import relax_job
        result = relax_job(
            atoms=molgraph.atoms,
            charge=charge,
            spin_multiplicity=mult,
            xc="b97-3c",
            basis="def2-tzvp",
        )
    elif method == "aimnet2":
        from himatcal.calculator.aimnet import AIMNet2ASE
        from himatcal.recipes.quacc.core import relax_job
        calc = AIMNet2ASE(
            "aimnet2_b973c",
            charge=charge,
            mult=mult,
        )
        result = relax_job(
            atoms=molgraph.atoms,
            calc=calc,
        )
    logging.info(f"Relaxation of {molgraph.smiles} is done.")

    molgraph.atoms = result["atoms"]
    molgraph.energy = result["results"]["energy"]
    molgraph.state = "opt"
    molgraph.label = "b97-3c"

    return molgraph


def molgraph_spe(molgraph: MolGraph, charge: int | None, mult: int | None):
    from himatcal.recipes.gaussian.flow import calc_free_energy

    if charge is None or mult is None:
        charge, mult = get_charge_and_spin(molgraph.smiles)

    freeE = calc_free_energy(
        atoms=molgraph.atoms,
        charge=charge,
        mult=mult,
        label="molgraph",
        relax=False,
    )
    freeE.run()
    logging.info(f"SPE calculation of {molgraph.smiles} is done.")

    molgraph.energy = freeE.extract_free_energy()
    molgraph.state = "final"
    molgraph.label = "b97-3c//b3lyp/6-311+g(d)/gd3bj/acetone"

    return molgraph
