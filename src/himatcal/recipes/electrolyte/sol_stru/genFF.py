"""generate the GAFF force field with/without the RESP2 charge, better use the chk file form desired level of theory"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

from monty.os import cd

from himatcal import SETTINGS
from himatcal.recipes.electrolyte.sol_stru._base import (
    formchk,
    gen_resp2_chg,
)
from himatcal.recipes.gaussian.core import relax_job, static_job

if TYPE_CHECKING:
    from ase.atoms import Atoms


def sobtop_genFF(fchk_file: str, chg_file: str | None = None):
    """
    Generate a force field using the Sobtop method from a formatted check file.

    This function processes a formatted check file by converting it to various molecular formats and generating a GAFF force field with RESP2 charges. It utilizes external tools such as Multiwfn and Open Babel to perform the necessary conversions and calculations.

    Args:
        fchk_file (str): The path to the formatted check file to be processed.
        chg_file (str | None): The path to the charge file; if None, a charge file will be generated from the formatted check file.

    Returns:
        None

    Raises:
        subprocess.CalledProcessError: If any of the subprocess commands fail.

    Examples:
        sobtop_genFF("path/to/file.fchk")
    """

    sobtop_path = SETTINGS.SOBTOP_PATH
    sobtop_parent_path = Path(sobtop_path).parent
    multiwfn_path = SETTINGS.MULTIWFN_PATH
    obabelpath = SETTINGS.OBABEL_PATH
    pdbfile_path = Path(fchk_file).with_suffix(".pdb")
    mol2file_path = Path(fchk_file).with_suffix(".mol2")
    grofile_path = Path(fchk_file).with_suffix(".gro")
    topofile_path = Path(fchk_file).with_suffix(".top")
    itpfile_path = Path(fchk_file).with_suffix(".itp")
    with cd(str(sobtop_parent_path)):
        # * convert the fchk file to pdb file using Multiwfn
        subprocess.run(
            [multiwfn_path],
            input=f"{fchk_file}\n100\n2\n1\n{pdbfile_path}\n0\nq\n",
            text=True,
            check=True,
        )

        # * convert the pdb file to mol2 file using obabel
        subprocess.run(
            [
                obabelpath,
                "-ipdb",
                str(pdbfile_path),
                "-omol2",
                "-O",
                str(mol2file_path),
            ],
            check=True,
        )

        # * generate the GAFF force field with the RESP2 charge using sobtop
        chg_file = str(gen_resp2_chg(fchk_file)) if chg_file is None else str(chg_file)
        sobtop_command = f"{mol2file_path}\n7\n10\n{chg_file}\n0\n2\n{grofile_path}\n1\n2\n7\n{fchk_file}\n{topofile_path}\n{itpfile_path}\n0\n"
        subprocess.run([sobtop_path], input=sobtop_command, text=True, check=True)


def genChk(atoms: Atoms, chg: int = 0, mult: int = 1, label="mol"):
    """generate the check file with the b3lyp/def2tzvp level of theory"""
    # * relax the structure at the b3lyp/6-31G* level
    result = relax_job(
        atoms,
        charge=chg,
        spin_multiplicity=mult,
        xc="b3lyp",
        basis="6-31G*",
        label=f"{label}-relax",
    )
    relaxed_atoms = result["atoms"]
    # * run the single point calculation at the b3lyp/def2tzvp level default
    static_job_result = static_job(
        relaxed_atoms,
        charge=chg,
        spin_multiplicity=mult,
        xc="b3lyp",
        basis="def2tzvp",
        freq=True,
        label=f"{label}-static",
    )

    # * get the check file
    chkfile_path = (
        Path(static_job_result["dir_name"]) / static_job_result["parameters"]["chk"]
    )
    gzfile_path = chkfile_path.with_suffix(".chk.gz")

    if not chkfile_path.exists() and gzfile_path.exists():
        subprocess.run(["gunzip", str(gzfile_path)], check=True)

    if not chkfile_path.exists():
        raise FileNotFoundError(f"Neither {chkfile_path} nor {gzfile_path} exists.")

    return chkfile_path


def genFF(atoms: Atoms, chg: int = 0, mult: int = 1, label="mol"):
    """generate the GAFF force field with the RESP2 charge"""
    chkfile_path = genChk(atoms, chg, mult, label)
    fchkfile_path = formchk(chkfile_path)
    sobtop_genFF(str(fchkfile_path))
    logging.info(f"GAFF force field generated for {atoms.get_chemical_formula()}")
    # copy the gro,pdb,top and itp files to the current directory
    for suffix in ["gro", "pdb", "top", "itp"]:
        source = chkfile_path.with_suffix(f".{suffix}")
        target = Path(f"{label}.{suffix}")
        source.replace(target)
    return {
        "atoms": atoms,
        "chkfile_path": str(chkfile_path),
        "grofile_path": str(chkfile_path.with_suffix(".gro")),
        "topfile_path": str(chkfile_path.with_suffix(".top")),
        "itpfile_path": str(chkfile_path.with_suffix(".itp")),
    }


def genFF_noRESP2(atoms: Atoms, label="mol"):
    obabelpath = SETTINGS.OBABEL_PATH
    sobtop_path = SETTINGS.SOBTOP_PATH
    pdbfile_path = Path(f"{label}.pdb")
    topofile_path = Path(f"{label}.top")
    itpfile_path = Path(f"{label}.itp")
    mol2file_path = Path(f"{label}.mol2")
    grofile_path = Path(f"{label}.gro")

    atoms.write("tmp.pdb")
    # * convert the pdb file to mol2 file using obabel
    subprocess.run(
        [
            obabelpath,
            "-ipdb",
            str(pdbfile_path),
            "-omol2",
            "-O",
            str(mol2file_path),
        ],
        check=True,
    )
    # * generate the GAFF force field with the RESP2 charge using sobtop
    sobtop_command = f"{mol2file_path}\n2\n{grofile_path}\n1\n2\n4{topofile_path}\n{itpfile_path}\n0\n"
    subprocess.run([sobtop_path], input=sobtop_command, text=True, check=True)
