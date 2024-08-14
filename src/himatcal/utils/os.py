"""os utilities for himatcal"""

from __future__ import annotations

import os
import re
from pathlib import Path


def labeled_dir(main_workdir: Path, label: str):
    """
    Create a new folder in the main working directory with the provided label.

    Args:
        main_workdir (Path): The main working directory.
        label (str): The label of the folder.

    Returns:
        Path: The path of the new folder.
    """
    folder_names = [p.name for p in main_workdir.iterdir() if p.is_dir()]
    numbers = (
        int(re.search(r"\d+", name).group())
        for name in folder_names
        if re.search(r"\d+", name)
    )
    new_number = max(numbers, default=0) + 1
    # Create new folder
    folder_name = f"{new_number:02d}.{label}"
    folder_path = main_workdir / folder_name
    Path.mkdir(folder_path, parents=True, exist_ok=True)
    return folder_path


def get_chg_mult(molname):
    """
    Get the label, charge, and multiplicity from the name of a molecule.

    Args:
        molname (str): The name of the molecule in the format {label}-c{charge}s{multiplicity}.

    Returns:
        tuple: A tuple containing the label (str), charge (int), and multiplicity (int) of the molecule.
               If the name does not match the expected format, returns (None, None, None).
    """
    import re

    pattern = r"(.*?)-c(n?\d)s(\d+)"
    if not (match := re.match(pattern, molname)):
        return None, None, None
    label, chg, mult = match.groups()
    chg = f"-{chg[1:]}" if chg.startswith("n") else chg
    return label, int(chg), int(mult)


def write_chg_mult_label(label, chg, mult):
    """Write the label, chg and mult to a string, format: {label}-c{charge}s{mult}"""
    if chg < 0:
        chg = f"n{abs(chg)}"
    return f"{label}-c{chg}s{mult}"


def extract_fchk(label, dzip=False):
    """
    Extracts the formatted checkpoint file (.fchk) from a Gaussian checkpoint file (.chk).

    Args:
        label (str): The label to use for the extracted .fchk file.
        dzip (bool, optional): Whether to decompress the Gaussian checkpoint file if it is gzipped.

    Returns:
        None
    """
    if dzip:
        os.system("gzip -d Gaussian.chk.gz")
    chk_file = "Gaussian.chk"
    if not os.path.exists(chk_file):
        print(f"{chk_file} not found")
        return
    os.system(f"formchk {chk_file}")
    os.system(f"mv Gaussian.fchk {label}.fchk")
    print(f"fchk file extracted for {label}")


def get_homo_lumo(logfile):
    """
    Extracts HOMO, LUMO, and related energies and gaps from a computational chemistry log file.

    Args:
        logfile (str): Path to the computational chemistry log file.

    Returns:
        dict: A dictionary containing HOMO and LUMO orbitals, energies, gaps, and the minimum HOMO-LUMO gap.
    """
    import cclib
    from quacc.schemas.cclib import _get_homos_lumos

    data = cclib.io.ccread(logfile)
    HOMO = data.homos + 1
    LUMO = data.homos + 2
    homo_energies, lumo_energies, gaps = _get_homos_lumos(data.moenergies, data.homos)
    min_gap = min(gaps)
    return {
        "homo_orbital": HOMO,
        "lumo_orbital": LUMO,
        "homo_energies": homo_energies,
        "lumo_energies": lumo_energies,
        "homo_lumo_gaps": gaps,
        "min_homo_lumo_gap": min_gap,
    }


def cclib_result(log_path):
    """Extracts and reads computational chemistry log files.

    This function checks for compressed log files in the specified directory, decompresses the first found file, and reads the contents using the cclib library. It returns the parsed data from the log file.

    Args:
        log_path (Path): The directory path where log files are located.

    Returns:
        object: The parsed data from the log file.

    Raises:
        FileNotFoundError: If no log files are found and the directory is empty.

    Examples:
        result = cclib_result(Path("/path/to/logs"))
    """

    import contextlib
    import gzip

    import cclib

    with contextlib.suppress(FileNotFoundError):
        if gzip_log := list(log_path.glob("*.log.gz")):
            unzip_file = gzip.decompress(Path.open(gzip_log[0], "rb").read())
            logfile = gzip_log[0].with_suffix("")
            with Path.open(logfile, "w") as f:
                f.write(unzip_file.decode())
        log_files = list(log_path.glob("*.log"))
        return cclib.io.ccread(log_path / log_files[0])
