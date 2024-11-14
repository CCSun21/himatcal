"""core recipes for crest calculations"""

from __future__ import annotations

import logging
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

from ase.io import read, write
from monty.os import cd

from himatcal import SETTINGS

if TYPE_CHECKING:
    from typing import Literal

    from ase import Atoms

logger = logging.getLogger(__name__)


def relax(
    atoms: Atoms,
    chg: int = 0,
    mult: int = 1,
    gfn_level: Literal["gfn1", "gfn2", "gfnff", "gfn2//gfnff"] = "gfn2",
    alpb: str | None = None,
    threads: int = 16,
):
    """
    Relax a molecular system using the CREST optimization program.

    Args:
        atoms: The molecular system to relax.
        chg: The charge of the system.
        mult: The multiplicity of the system.
        gfn_level: The level of the GFN method to use (default is "gfn2").
        alpb: The solvent model to use (default is "acetone").
        threads: The number of threads to use for optimization (default is 4).

    Returns:
        The relaxed molecular system.

    Raises:
        FileNotFoundError: If the output file "crestopt.xyz" is not found after optimization.
    """
    scratch_folder = Path(f"crest_opt_{datetime.now(timezone.utc).strftime('%Y-%m-%d-%H-%M-%S-%f')}")
    Path.mkdir(scratch_folder)
    with cd(scratch_folder):
        atoms_name = "input.xyz"
        write(atoms_name, atoms)
        uhf = mult - 1
        protonate_cmd = f"{SETTINGS.CREST_EXE_PATH_V3} {atoms_name} --opt --{gfn_level} -chrg {chg} -uhf {uhf} --T {threads}"
        if alpb:
            protonate_cmd += f" -alpb {alpb}"
        with Path.open(Path("crest_opt.sh"), "w") as f:
            f.write(f"#!/bin/bash\n{protonate_cmd}")
        log_file_path = Path("crest_opt.log")
        with log_file_path.open("w") as log_file:
            subprocess.run(
                ["bash", "crest_opt.sh"],
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=True,
            )
        try:
            return read("crestopt.xyz")
        except FileNotFoundError:
            logger.error(
                "The relaxation did not complete successfully, please check the log file."
            )
            return None

def iMTD_GC(
    atoms: Atoms,
    chg: int = 0,
    mult: int = 1,
    gfn_level: Literal["gfn1", "gfn2", "gfnff", "gfn2//gfnff"] = "gfn2",
    alpb: str | None = None,
    topo_change: bool = False,
    threads: int = 16,
):
    scratch_folder = Path(f"crest_opt_{datetime.now(timezone.utc).strftime('%Y-%m-%d-%H-%M-%S-%f')}")
    Path.mkdir(scratch_folder)
    with cd(scratch_folder):
        atoms_name = "input.xyz"
        write(atoms_name, atoms)
        uhf = mult - 1
        if SETTINGS.CREST_EXE_PATH_V3 is None:
            logger.error("CREST_EXE_PATH_V3 is not set in the settings.")
            return None
        protonate_cmd = f"{SETTINGS.CREST_EXE_PATH_V3} {atoms_name} --{gfn_level} -chrg {chg} -uhf {uhf} --T {threads}"
        if alpb:
            protonate_cmd += f" -alpb {alpb}"
        if topo_change:
            protonate_cmd += " --noreftopo"
        with Path.open(Path("crest_opt.sh"), "w") as f:
            f.write(f"#!/bin/bash\n{protonate_cmd}")
        log_file_path = Path("crest_opt.log")
        with log_file_path.open("w") as log_file:
            subprocess.run(
                ["bash", "crest_opt.sh"],
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=True,
            )
        try:
            return read("crest_best.xyz")
        except FileNotFoundError:
            logger.error(
                "The computation did not finish successfully, please check the log file."
            )
            return None

def protonate(
    atoms: Atoms,
    ion: str = "li+",
    chg: int = 0,
    mult: int = 1,
    gfn_level: Literal["gfn1", "gfn2", "gfnff", "gfn2//gfnff"] = "gfn2",
    alpb: str | None = None,
    threads: int = 16,
):
    """
    Protonate a structure using CREST, default is to protonate with Li
    """
    scratch_folder = Path(f"crest_opt_{datetime.now(timezone.utc).strftime('%Y-%m-%d-%H-%M-%S-%f')}")
    Path.mkdir(scratch_folder)
    with cd(scratch_folder):
        atoms_name = "input.xyz"
        write(atoms_name, atoms)
        uhf = mult - 1
        protonate_cmd = f"{SETTINGS.CREST_EXE_PATH_V3} {atoms_name} --protonate --swel {ion} --{gfn_level} -chrg {chg} -uhf {uhf} --T {threads}"
        if alpb:
            protonate_cmd += f" -alpb {alpb}"
        with Path.open(Path("crest_opt.sh"), "w") as f:
            f.write(f"#!/bin/bash\n{protonate_cmd}")
        log_file_path = Path("crest_opt.log")
        with log_file_path.open("w") as log_file:
            subprocess.run(
                ["bash", "crest_opt.sh"],
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=True,
            )
        try:
            return read("protonated.xyz", index=0)
        except FileNotFoundError:
            logger.error(
                "The protonation did not complete successfully, please check the log file."
            )
            return None
