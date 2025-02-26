"""os utilities for himatcal"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path

from pydantic import BaseModel, field_validator

from himatcal.utils.locache import persist


class MoleculeInfo(BaseModel):
    """Represents information about a molecule, including its label, charge, and multiplicity.

    This class provides methods to parse charge values and to create an instance from a filename that follows a specific naming convention.
    """

    label: str
    charge: int
    mult: int

    @field_validator("charge", mode="before")
    def parse_charge(cls, value):
        """Parses the charge value from a string format.

        Args:
            value (str): The charge value as a string, which may start with 'n' to indicate a negative charge.

        Returns:
            int: The parsed charge as an integer.
        """
        value = str(value)
        return -int(value[1:]) if value.startswith("n") else int(value)

    @classmethod
    def from_filename(cls, filename: str):
        """Creates a MoleculeInfo instance from a filename.

        The filename must match the pattern '<label>_c<charge><mult>', where <charge> can be prefixed with 'n' for negative charges.

        Args:
            filename (str): The filename from which to extract molecule information.

        Returns:
            MoleculeInfo: An instance of MoleculeInfo with the extracted label, charge, and multiplicity.

        Raises:
            ValueError: If the filename does not match the expected format.
        """
        pattern = r"(.*?)-c(n?\d)s(\d+)"
        if not (match := re.match(pattern, filename)):
            raise ValueError(f"Filename {filename} does not match the expected format")
        label, chg, mult = match.groups()
        return cls(label=label, charge=chg, mult=int(mult))

    def write_filename(self):
        return (
            f"{self.label}-cn{-self.charge}s{self.mult}"
            if self.charge < 0
            else f"{self.label}-c{self.charge}s{self.mult}"
        )


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


def get_chg_mult(molname: str):
    """
    Get the label, charge, and multiplicity from the name of a molecule.
    Deprecated: Method moved to MoleculeInfo(BaseModel)
    """
    mol_info = MoleculeInfo.from_filename(molname)
    return mol_info.label, mol_info.charge, mol_info.mult


def write_chg_mult_label(label, chg, mult):
    """Write the label, chg and mult to a string, format: {label}-c{charge}s{mult}"""
    return MoleculeInfo(label=label, charge=chg, mult=mult).write_filename()


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
    chk_file = Path("Gaussian.chk")
    if not Path.exists(chk_file):
        logging.info(f"{chk_file} not found")
        return
    os.system(f"formchk {chk_file}")
    os.system(f"mv Gaussian.fchk {label}.fchk")
    logging.info(f"fchk file extracted for {label}")


@persist(max_age=1)
def parse_logfile(logfile: str):
    """
    Parses a Gaussian log file to extract computational chemistry data.
    This function utilizes the cclib library to read and parse the contents of the specified logfile, returning structured data for further analysis.

    Args:
        logfile (str): The path to the logfile to be parsed.

    Returns:
        object: A structured data object containing the parsed information from the logfile.

    Raises:
        FileNotFoundError: If the specified logfile does not exist.
        ValueError: If the logfile cannot be parsed due to invalid format.

    Examples:
        >>> data = parse_logfile("path/to/logfile.log")
    """

    import cclib

    return cclib.io.ccopen(logfile).parse()  # type: ignore


def parse_homo_lumo(logfile: str):
    """
    Extracts the highest occupied molecular orbital (HOMO) and lowest unoccupied molecular orbital (LUMO) energies from a specified logfile.
    This function computes the HOMO and LUMO energies along with the energy gaps between them, providing insights into the electronic structure of the molecule.

    Args:
        logfile (str): The path to the logfile containing molecular orbital energies.

    Returns:
        tuple: A tuple containing three elements:
            - List of HOMO energies.
            - List of LUMO energies or None if not available.
            - List of HOMO-LUMO energy gaps or None if not available.

    Raises:
        ValueError: If the logfile does not contain valid data.

    Examples:
        >>> homo_energies, lumo_energies, gaps = parse_homo_lumo("path/to/logfile.log")
    """
    data = parse_logfile(logfile)
    moenergies = data.moenergies  # type: ignore
    homo_indices = data.homos  # type: ignore
    homo_energies = [moenergies[i][h] for i, h in enumerate(homo_indices)]

    for i, h in enumerate(homo_indices):
        if len(moenergies[i]) < h + 2:
            return homo_energies, None, None
    lumo_energies = [moenergies[i][h + 1] for i, h in enumerate(homo_indices)]
    homo_lumo_gaps = [
        lumo_energies[i] - homo_energies[i] for i in range(len(homo_energies))
    ]
    return homo_energies, lumo_energies, homo_lumo_gaps


def get_homo_lumo(logfile: str):
    """
    Extracts HOMO, LUMO, and related energies and gaps from a computational chemistry log file.

    Args:
        logfile (str): Path to the computational chemistry log file.

    Returns:
        dict: A dictionary containing HOMO and LUMO orbitals, energies, gaps, and the minimum HOMO-LUMO gap.
    """

    data = parse_logfile(logfile)
    HOMO = data.homos + 1
    LUMO = data.homos + 2
    homo_energies, lumo_energies, gaps = parse_homo_lumo(logfile)
    min_gap = min(gaps)
    return {
        "homo_orbital": HOMO,
        "lumo_orbital": LUMO,
        "homo_energies": homo_energies,
        "lumo_energies": lumo_energies,
        "homo_lumo_gaps": gaps,
        "min_homo_lumo_gap": min_gap,
    }


def cclib_result(log_path: Path):
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
        return parse_logfile(log_path / log_files[0])
