"""
uvi git+https://gitlab.com/ase/ase.git
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import pandas as pd
import yaml
from ase.io import read, write
from monty.os import cd
from pysisyphus.plot import plot_irc
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops

from himatcal.recipes.reaction import MolGraph
from himatcal.recipes.reaction._base import Reaction
from himatcal.recipes.reaction.newtonnet import geodesic_ts_hess_irc_newtonnet
from himatcal.recipes.reaction.utils import (
    molgraph_relax,
    molgraph_spe,
    update_molgraph,
)


def mol_with_atom_and_bond_indices(smiles, output_file: str | None = None):
    from rdkit import Chem
    from rdkit.Chem.Draw import IPythonConsole

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    # Add atom indices
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    # Configure drawing options
    IPythonConsole.drawOptions.addBondIndices = True
    IPythonConsole.molSize = (350, 300)

    # Generate and save image
    if output_file:
        with open(output_file, "wb") as f:
            f.write(mol._repr_png_())

    return mol


# make xyz from mol
def mol_to_xyz(mol):
    # Add hydrogens to the molecule
    mol = Chem.AddHs(mol)
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol)
    # Optimize the 3D coordinates
    AllChem.MMFFOptimizeMolecule(mol)
    # Get the atomic coordinates
    conf = mol.GetConformer()
    xyz = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    # Get the atomic symbols
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    return symbols, xyz


# write xyz file
def write_xyz(filename, symbols, xyz):
    with open(filename, "w") as f:
        f.write(f"{len(symbols)}\n")
        f.write("\n")
        for i, (symbol, coord) in enumerate(zip(symbols, xyz)):
            x, y, z = coord
            f.write(f"{symbol} {x} {y} {z}\n")


def get_fragments(mol, bond_idx):
    # 获取指定键的两个原子
    bond = mol.GetBondWithIdx(bond_idx)
    atom1 = bond.GetBeginAtomIdx()
    atom2 = bond.GetEndAtomIdx()

    # 断裂指定键
    fragments = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=False)

    # 获取断裂后的两个碎片
    frags = Chem.GetMolFrags(fragments, asMols=False, sanitizeFrags=False)

    # 获取每个碎片上原子的序号列表
    frag1_atoms = list(frags[0])
    if len(frags) == 1:
        frag2_atoms = list(frags[0])
    else:
        frag2_atoms = list(frags[1])

    return frag1_atoms, frag2_atoms, atom1, atom2


def update_rxn_results(filename, reaction_results):
    reaction_data = {
        "Reactant_SMILES": [reaction_results.rsmi],
        "Product_SMILES": [reaction_results.psmi],
        "Activation_Energy": [reaction_results.ea],
        "Enthalpy_Change": [reaction_results.dh],
    }

    # Convert the dictionary to a DataFrame
    reaction_df = pd.read_csv(filename)
    reaction_df_add = pd.DataFrame(reaction_data)
    reaction_df = pd.concat([reaction_df, reaction_df_add], ignore_index=True)
    reaction_df.drop_duplicates().to_csv(filename, index=False)


def prepare_molgraph(smiles, charge, mult, bond_break, filename, break_scale):
    mol = mol_with_atom_and_bond_indices(smiles)
    symbols, xyz = mol_to_xyz(mol)
    write_xyz(filename, symbols, xyz)
    logging.info(symbols)
    logging.info(xyz)

    frag1_atoms, frag2_atoms, atom1, atom2 = get_fragments(mol, bond_break[0])
    bond_break_config = [f"{atom1 + 1} {atom2 + 1} B"]
    with open("autoCG.com", "w") as f:
        f.write(f"{charge} {mult}\n")
        for i, (symbol, coord) in enumerate(zip(symbols, xyz)):
            x, y, z = coord
            f.write(f"{symbol} {x} {y} {z}\n")
        f.write("\n")
        for config in bond_break_config:
            f.write(f"{config}")

    cache_dir = Path("autoCG_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    frag_dir = Path("autoCG_fragments")
    frag_dir.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "/home/suncc/micromamba/envs/marimo/bin/python",
            "/home/suncc/Code/mmscd2/TAT/try-local/autoCG/autoCG/generate.py",
            "autoCG.com",
            "-sd",
            f"{frag_dir.resolve()}",
            "-wd",
            f"{cache_dir.resolve()}",
            "-bs",
            str(break_scale),
        ],
        check=True,
    )

    reactant = MolGraph(atoms=read(frag_dir / "1" / "R.xyz"))
    product = MolGraph(atoms=read(frag_dir / "1" / "P.xyz"))
    return reactant, product


def optimize_molgraph(molgraph, charge, mult):
    return molgraph_relax(molgraph, charge=charge, mult=mult, method="aimnet2")


def run_de_gsm(reactant, product, charge, mult):
    GS_wd = Path("DE-GS")
    GS_wd.mkdir(exist_ok=True)
    reactant.atoms.write(GS_wd / "initial.xyz")
    product.atoms.write(GS_wd / "final.xyz")

    gs_yaml = f"""geom:
 type: dlc
 fn: [initial.xyz, final.xyz]
calc:
 type: xtb
 pal: 16
 charge: {charge}
 mult: {mult}
 gbsa: acetone
preopt:
cos:
 type: gs
 max_nodes: 18
 climb: True
opt:
 type: string
 align: False
 stop_in_when_full: 5
 max_cycles: 20
tsopt:
 type: rsirfo
 do_hess: True
 thresh: gau
 hessian_recalc: 3
irc:
 type: eulerpc
 rms_grad_thresh: 0.0005
endopt:
barriers:
 solv_calc:
  type: xtb
  charge: {charge}
  mult: {mult}
  gbsa: acetone
  pal: 16
 do_standard_state_corr: True
"""

    with Path.open(GS_wd / "GS.yaml", "w") as f:
        f.write(gs_yaml)

    subprocess.run("pysis GS.yaml | tee pysis.log", cwd=GS_wd, shell=True, check=False)

    with cd(GS_wd):
        irc_img = plot_irc()

    try:
        ts = MolGraph(atoms=read(GS_wd / "ts_final_geometry.xyz"))
    except FileNotFoundError:
        ts = None
        logging.warning("No ts_final_geometry.xyz found by DE-GS method")

    return ts


def run_newtonnet(reactant, product):
    job1, job2, job3, job4 = geodesic_ts_hess_irc_newtonnet(reactant, product)
    ts = MolGraph(atoms=job2["atoms"])
    reactant_irc = MolGraph(atoms=job3["trajectory"][-1])
    write("ts_irc_forward.xyz", job3["trajectory"])
    product_irc = MolGraph(atoms=job4["trajectory"][-1])
    write("ts_irc_reverse.xyz", job4["trajectory"])
    return reactant_irc, product_irc, ts


def save_molgraph(molgraph, filename):
    update_molgraph(molgraph, filename=filename)


def calculate_free_energy(molgraph, charge, mult, filename):
    molgraph_free_e = molgraph_spe(molgraph, charge, mult)
    update_molgraph(molgraph_free_e, filename=filename)
    return molgraph_free_e


def load_config(config_file):
    with open(config_file) as file:
        config = yaml.safe_load(file)
    return config


def main():
    # * Load config file
    config_file = "config.yaml"
    config = load_config(config_file)
    reactant_smi = config["reactant_smi"]
    product_smi = config["product_smi"]
    chg = config["chg"]
    mult = config["mult"]
    bond_break = config["bond_break"]
    filename = config["filename"]
    break_scale = config["break_scale"]
    b97_3c_file = config["b97_3c_file"]
    molgraph_file = config["molgraph_file"]
    reaction_results_file = config["reaction_results_file"]

    # * Prepare molgraph using autoCG
    reactant, product = prepare_molgraph(
        reactant_smi, chg, mult, bond_break, filename, break_scale
    )

    # * write molgraph to file
    logging.info("Writing reactant and product to file")
    reactant.atoms.write("reactant.xyz")
    product.atoms.write("product.xyz")

    # * Optimize reactant and product
    # reactant = optimize_molgraph(reactant, chg, mult)
    # product = optimize_molgraph(product, chg, mult)

    # reactant.atoms.write("reactant_opt.xyz")
    # product.atoms.write("product_opt.xyz")
    # logging.info("Optimized reactant and product written to file")

    # * Run DE-GSM calculation to find TS
    # ts = run_de_gsm(reactant, product, chg, mult)
    reactant, product, ts = run_newtonnet(reactant.atoms, product.atoms)
    # if ts:
    #     ts.atoms.write("ts.xyz")
    #     logging.info("TS found, written to ts.xyz")
    #     result = bare_job(
    #         atoms=ts.atoms,
    #         charge=chg,
    #         spin_multiplicity=mult,
    #         xc="b97-3c",
    #         basis="def2-SVP",
    #         orcasimpleinput=["OptTS NumFreq"],
    #         orcablocks=["%geom Calc_Hess True NumHess True Recalc_Hess 5 END"],
    #     )
    #     ts = MolGraph(atoms=result["atoms"])
    #     ts.atoms.write("ts_opt.xyz")
    #     logging.info("TS relaxed, written to ts_opt.xyz")
    #     save_molgraph(ts, b97_3c_file)
    # else:
    #     logging.info("No TS found by DE-GS method")
    save_molgraph(ts, b97_3c_file)
    save_molgraph(reactant, b97_3c_file)
    save_molgraph(product, b97_3c_file)
    reactant.atoms.write("reactant_final.xyz")
    product.atoms.write("product_final.xyz")
    ts.atoms.write("ts_final.xyz")
    logging.info("TS relaxed, written to ts_opt.xyz")

    reactant_free_e = calculate_free_energy(reactant, chg, mult, molgraph_file)
    product_free_e = calculate_free_energy(product, chg, mult, molgraph_file)
    ts_free_e = calculate_free_energy(ts, chg, mult, molgraph_file)

    reaction = Reaction(
        reactant=reactant_free_e,
        product=product_free_e,
        ts=ts_free_e,
    )
    logging.info(str(reaction.barrier))
    logging.info(reaction.reaction_results.json)
    update_rxn_results(reaction_results_file, reaction.reaction_results)


if __name__ == "__main__":
    main()
