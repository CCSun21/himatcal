from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path

from monty.os import cd

from himatcal import SETTINGS

formchkpath = str(getattr(SETTINGS, "FORMCHK_PATH", ""))
if not formchkpath:
    logging.info("FORMCHK_PATH not found in SETTINGS")

multiwfnpath = str(getattr(SETTINGS, "MULTIWFN_PATH", ""))
if not multiwfnpath:
    logging.info("MULTIWFN_PATH not found in SETTINGS")


def formchk(chk_file: str):
    """
    Convert a check file to a formatted check file using an external command.

    This function takes a check file as input, executes a subprocess to convert it to a formatted check file, and returns the path of the newly created file. It ensures that the conversion process is executed correctly by checking for errors during the subprocess call.

    Args:
        chk_file (str): The path to the check file to be converted.

    Returns:
        Path: The path to the newly created formatted check file with a .fchk suffix.

    Raises:
        subprocess.CalledProcessError: If the subprocess command fails.

    Examples:
        fchk_file = formchk("/path/to/file.chk")
    """

    chk_file_path = Path(chk_file)
    subprocess.run(f"{formchkpath} {chk_file}", shell=True, check=True)
    return chk_file_path.with_suffix(".fchk")


def charges_to_json(charges, label):
    charge_data = []
    for charge in charges:
        import re

        if match := re.match(r"(\d+)\((\w+)\s*\)", charge):
            atom_index, atom_type = match.groups()
        else:
            raise ValueError(f"Unexpected charge format: {charge}")
        atom_charge = float(charge.split()[2])
        charge_data.append(
            {
                "atom_index": atom_index,
                "atom_type": atom_type,
                "atom_charge": atom_charge,
            }
        )
    return json.dumps({"label": label, "charges": charge_data}, indent=4)


def extract_resp2(fchk_file):
    label = Path(fchk_file).stem
    result = subprocess.Popen(
        multiwfnpath,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if result.stdin:
        result.stdin.write(f"{fchk_file}\n")
        result.stdin.write("7\n18\n1\ny\n0\n0\nq\n")  # * Multiwfn RESP2 calculation
    else:
        raise RuntimeError("Failed to write to stdin of the subprocess")
    out, err = result.communicate()
    charges = []
    capture = False
    for line in out.splitlines():
        if "Center       Charge" in line:
            capture = True
            continue
        if capture:
            if line.strip() == "" or "Sum of charges" in line:
                break
            charges.append(line.strip())
    return charges_to_json(charges, label)


def gen_resp2_chg(chk_file: str):
    """
    Generate a charge file from a given check file by performing a RESP calculation.

    This function converts a check file to a formatted check file, runs a RESP calculation to extract charge data, and writes the results to a JSON file. It returns the path to the generated charge file with a .chg suffix.

    Args:
        chk_file (str): The path to the check file to be processed, both .chk and .fchk files are supported.

    Returns:
        Path: The path to the generated charge file with a .chg suffix.

    Examples:
        charge_file = gen_resp2_chg("path/to/check_file.chk")
    """
    # * Convert the check files to formatted check files.
    fchk_file = chk_file if chk_file.endswith(".fchk") else str(formchk(chk_file))
    # * Run the RESP calculation using Multiwfn.
    with cd(Path(fchk_file).parent):
        charges_json = extract_resp2(fchk_file=fchk_file)
        # write the chages to a json file
        json_file_path = Path(chk_file).with_suffix(".json")
        with json_file_path.open("w") as f:
            f.write(charges_json)
        return Path(chk_file).with_suffix(".chg")


def merge_chgs(chg_file_1, chg_file_2):
    """
    Merge two charge files by averaging the charges of corresponding atoms.

    This function reads two charge files, parses their contents to extract atom data and charges, and computes the average charge for each atom. The results are then written to a new charge file, which is returned with a modified filename.

    Args:
        chg_file_1 (str): The path to the first charge file to be merged.
        chg_file_2 (str): The path to the second charge file to be merged.

    Returns:
        Path: The path to the newly created merged charge file with a .merged.chg suffix.

    Raises:
        ValueError: If there is a mismatch in atom data between the two charge files.

    Examples:
        merged_file = merge_chgs("path/to/first.chg", "path/to/second.chg")
    """

    def parse_chg_file(chg_file):
        with Path(chg_file).open() as f:
            lines = f.readlines()
        parsed_data = []
        for line in lines:
            parts = line.split()
            atom_data = parts[:-1]
            charge = float(parts[-1])
            parsed_data.append((atom_data, charge))
        return parsed_data

    def average_charges(data1, data2):
        averaged_data = []
        for (atom_data1, charge1), (atom_data2, charge2) in zip(data1, data2):
            if atom_data1 != atom_data2:
                raise ValueError("Atom data mismatch between charge files")
            averaged_charge = (charge1 + charge2) / 2
            averaged_data.append((atom_data1, averaged_charge))
        return averaged_data

    def write_chg_file(data, output_file):
        with Path(output_file).open("w") as f:
            for atom_data, charge in data:
                f.write(" ".join(atom_data) + f" {charge:.10f}\n")

    data1 = parse_chg_file(chg_file_1)
    data2 = parse_chg_file(chg_file_2)
    averaged_data = average_charges(data1, data2)
    output_file = Path(chg_file_1).with_suffix(".merged.chg")
    write_chg_file(averaged_data, output_file)
    return output_file


def write_sobtopini():
    ## * write a sobtop.ini file in the current directory
    content = f"""  nthreads= 4  // Number of threads used for parallel calculation
        iskipgendih= 0  // 1: Skip generating dihedral terms, 0: Do not skip
        ioutatminfo= 0  // 1: Output atomic coordinates and connectivities when to atminfo.txt in current folder when loading input file, 0: Do not output to file but shown on screen
        ichggeom= 1  // 1: When loading .chg file, replace current geometry with that in .chg file, 0: Do not replace
        k_method= 2  // Default method of determining k. 1: Seminario, 2: mSeminario, 3: m2Seminario, 4: DRIH
        bondcrit= 1.15  // When pdb/pqr is used as input, two atoms are considered as bonded if their distance is smaller than sum of their covalent radii multiplied by this factor. Priority is lower than the criteria defined in bondcrit.dat
        Multiwfn_cmd= "{SETTINGS.MULTIWFN_PATH}"// Path of executable file of Multiwfn
        OpenBabel_cmd= "{SETTINGS.OBABEL_PATH}" // Path of executable file of OpenBabel"""
    with Path.cwd().joinpath("sobtop.ini").open("w") as f:
        f.write(content)
