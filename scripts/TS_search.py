from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import fire
import yaml
from ase.io import read, write
from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem import AllChem

from himatcal.recipes.quacc._base import clear_quacc_cache
from himatcal.recipes.reaction import MolGraph, Reaction
from himatcal.recipes.reaction.newtonnet import geodesic_ts_hess_irc_newtonnet
from himatcal.recipes.reaction.utils import molgraph_relax, molgraph_spe
from himatcal.utils.rdkit.core import mol_with_atom_and_bond_indices


def mol_to_xyz(mol):
    """Convert RDKit molecule to xyz coordinates"""
    # Add hydrogens and generate 3D coordinates
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, randomSeed=42)
    if result == -1:
        raise ValueError("Could not embed molecule")

    # Optimize the structure
    AllChem.MMFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    xyz = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    return symbols, xyz


def write_xyz(filename, symbols, coords):
    """Write xyz coordinates to file"""
    with Path(filename).open("w") as f:
        f.write(f"{len(symbols)}\n\n")
        for symbol, coord in zip(symbols, coords):
            f.write(f"{symbol} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")


def get_fragments(mol, bond_idx):
    """Get fragments from broken bond

    Args:
        mol: RDKit molecule object
        bond_idx: Index of bond to break

    Returns:
        tuple: (atom indices of fragment 1, atom indices of fragment 2,
               begin atom index, end atom index)
    """
    # 获取指定键的两个原子
    bond = mol.GetBondWithIdx(bond_idx)
    atom1 = bond.GetBeginAtomIdx()
    atom2 = bond.GetEndAtomIdx()

    # 断裂指定键
    fragments = Chem.rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=False)

    # 获取断裂后的碎片原子索引
    frags = Chem.GetMolFrags(fragments, asMols=False, sanitizeFrags=False)

    if len(frags) != 1:
        return frags[0], frags[1], atom1, atom2
    logging.warning(f"Bond {bond_idx} breakage did not result in fragments")
    # 如果没有分成两个片段,返回同一组原子
    return frags[0], frags[0], atom1, atom2


def prepare_molgraph(smiles, charge, mult, bond_break, filename, break_scale, work_dir):
    mol = mol_with_atom_and_bond_indices(smiles, work_dir / "molecule.png")
    symbols, xyz = mol_to_xyz(mol)
    write_xyz(work_dir / filename, symbols, xyz)
    logging.info(symbols)
    logging.info(xyz)

    frag1_atoms, frag2_atoms, atom1, atom2 = get_fragments(mol, bond_break[0])
    bond_break_config = [f"{atom1 + 1} {atom2 + 1} B"]

    with (work_dir / "autoCG.com").open("w") as f:
        f.write(f"{charge} {mult}\n")
        for symbol, coord in zip(symbols, xyz):
            x, y, z = coord
            f.write(f"{symbol} {x} {y} {z}\n")
        f.write("\n")
        for config in bond_break_config:
            f.write(f"{config}")

    # 获取输入文件的绝对路径
    input_file = work_dir.resolve() / "autoCG.com"
    cache_dir = work_dir.resolve() / "autoCG_cache"
    frag_dir = work_dir.resolve() / "autoCG_fragments"

    cache_dir.mkdir(parents=True, exist_ok=True)
    frag_dir.mkdir(parents=True, exist_ok=True)

    # 在指定目录中执行命令, 使用绝对路径
    subprocess.run(
        [
            "/home/suncc/micromamba/envs/marimo/bin/python",
            "/home/suncc/Code/mmscd2/TAT/try-local/autoCG/autoCG/generate.py",
            str(input_file),  # 使用绝对路径
            "-sd",
            str(frag_dir),
            "-wd",
            str(cache_dir),
            "-bs",
            str(break_scale),
        ],
        check=True,
    )

    atoms_r = read(frag_dir / "1" / "R.xyz")
    atoms_p = read(frag_dir / "1" / "P.xyz")
    reactant = MolGraph(atoms=atoms_r[0] if isinstance(atoms_r, list) else atoms_r)
    product = MolGraph(atoms=atoms_p[0] if isinstance(atoms_p, list) else atoms_p)
    return reactant, product


def run_newtonnet(reactant, product):
    job1, job2, job3, job4 = geodesic_ts_hess_irc_newtonnet(reactant, product)
    ts = MolGraph(atoms=job2["atoms"])
    reactant_irc = MolGraph(atoms=job3["trajectory"][-1])
    write("ts_irc_forward.xyz", job3["trajectory"])
    product_irc = MolGraph(atoms=job4["trajectory"][-1])
    write("ts_irc_reverse.xyz", job4["trajectory"])
    return reactant_irc, product_irc, ts


def clear_cache(path):
    logging.info("Cleaning up QUACC cache...")
    clear_quacc_cache(path)
    logging.info("Cleaning AutoCG cache...")
    directory_patterns = ["autoCG_cache", "autoCG_fragments"]
    for pattern in directory_patterns:
        for dir_path in path.glob(pattern):
            if dir_path.is_dir():
                for file in dir_path.iterdir():
                    if file.is_file():
                        file.unlink()
                dir_path.rmdir()
            else:
                dir_path.unlink()
    logging.info("Cache cleanup complete")


class ReactionConfig(BaseModel):
    """Configuration model for reaction search parameters"""

    reactant_smi: str = Field(..., description="Reactant SMILES string")
    chg: int = Field(0, description="Molecular charge")
    mult: int = Field(1, description="Spin multiplicity")
    cleanup_cache: bool | None = Field(
        False, description="Whether to clean up QUACC cache after calculation"
    )
    spe_method: str = Field(
        "aimnet2", description="Method to use for single point energy calculations"
    )


def main(config_path="config.yaml"):
    # * 1. Load and validate config
    with Path(config_path).open() as f:
        config_data = yaml.safe_load(f)
    config = ReactionConfig(**config_data)

    # * 2. Create main output directory
    output_dir = Path("reactions_output")
    output_dir.mkdir(exist_ok=True)

    # * 3. read the reactant SMILES and generate molecule image
    mol = mol_with_atom_and_bond_indices(
        config.reactant_smi, str(output_dir / "molecule.png")
    )

    # * 4. Create a separate directory for each bond and perform calculations

    for bond_idx in range(mol.GetNumBonds()):
        # * 4.1 Create a working directory for the bond
        work_dir = output_dir / f"bond_{bond_idx}"
        work_dir.mkdir(exist_ok=True)

        try:
            # * 4.2 Prepare the reactant and product
            reactant, product = prepare_molgraph(
                config.reactant_smi,
                config.chg,
                config.mult,
                [bond_idx],
                f"molecule_{bond_idx}.xyz",
                break_scale=5.0,
                work_dir=work_dir,
            )

            # * 4.3 Relax the reactant and product
            logging.info("Relaxing reactant and product")
            reactant = molgraph_relax(
                reactant, config.chg, config.mult, method="aimnet2"
            )
            product = molgraph_relax(product, config.chg, config.mult, method="aimnet2")

            # * 4.4 Save the results in the corresponding directory
            reactant.atoms.write(work_dir / f"reactant_{bond_idx}.xyz")
            product.atoms.write(work_dir / f"product_{bond_idx}.xyz")
            logging.info(f"Successfully processed bond {bond_idx} in {work_dir}")

            # * 4.5 Run NEB calculation using the NewtonNet MLP
            reactant, product, ts = run_newtonnet(reactant.atoms, product.atoms)
            reactant.atoms.write(work_dir / f"reactant_irc_{bond_idx}.xyz")
            product.atoms.write(work_dir / f"product_irc_{bond_idx}.xyz")
            ts.atoms.write(work_dir / f"ts_irc_{bond_idx}.xyz")

            # * 4.6 Calculate free energy of reactant, product, and TS
            reactant = molgraph_spe(reactant, config.chg, config.mult)
            product = molgraph_spe(product, config.chg, config.mult)
            ts = molgraph_spe(ts, config.chg, config.mult)

            # * 4.7 save the results
            reaction = Reaction(reactant=reactant, product=product, ts=ts)
            with (work_dir / "reaction.json").open("w") as f:
                f.write(reaction.reaction_results.model_dump_json())
            logging.info(f"results: {reaction.reaction_results.json}")

        except Exception as e:
            logging.error(f"Error processing bond {bond_idx} in {work_dir}: {e!s}")
            continue

        # * 5. Clean up  cache
        if config.cleanup_cache is True:
            clear_cache(output_dir)


if __name__ == "__main__":
    fire.Fire(main)
