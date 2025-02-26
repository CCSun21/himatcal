from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import psutil
from pydantic import BaseModel
from quacc import job

from himatcal.recipes.gaussian._base import run_and_summarize

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from quacc.types import Filenames, RunSchema, SourceDirectory

logger = logging.getLogger(__name__)


@job
def relax_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    label: str = "relax",
    xc: str = "wb97xd",
    basis: str = "def2tzvp",
    freq: bool = False,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Run a Gaussian relaxation calculation.

    Args:
        atoms: Atoms - The atoms for the calculation.
        charge: int - The charge of the system.
        spin_multiplicity: int - The spin multiplicity of the system.
        label: str - (Optional) The label for the calculation (default is "relax").
        xc: str - (Optional) The exchange-correlation functional to use (default is "wb97xd").
        basis: str - (Optional) The basis set to use (default is "def2tzvp").
        freq: bool - (Optional) Whether to perform frequency calculations (default is False).
        copy_files: Union[SourceDirectory, dict[SourceDirectory, Filenames], None] - (Optional) Files to copy after the calculation.
        **calc_kwargs - Additional keyword arguments for the calculation.

    Returns:
        RunSchema - The summarized result of the Gaussian relaxation calculation.
    """

    calc_defaults = {
        "mem": "64GB",
        "chk": "Gaussian.chk",
        "nprocshared": psutil.cpu_count(logical=False),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": spin_multiplicity,
        "opt": "",
        "pop": "CM5",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "ioplist": ["2/9=2000"],  # ASE issue #660
    }
    if freq:
        calc_defaults["freq"] = ""

    return run_and_summarize(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        label=label,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "Gaussian Relax"},
        copy_files=copy_files,
    )


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    label: str = "static",
    xc: str = "wb97xd",
    basis: str = "def2tzvp",
    freq: bool = True,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    xc
        Exchange-correlation functional
    basis
        Basis set
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Gaussian calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.gaussian.Gaussian][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results
    """
    calc_defaults = {
        "mem": "64GB",
        "chk": "Gaussian.chk",
        "nprocshared": psutil.cpu_count(logical=True),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": spin_multiplicity,
        # "force": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "pop": "CM5",
        "gfinput": "",
        "ioplist": ["6/7=3", "2/9=2000"],  # see ASE issue #660
    }
    if freq:
        calc_defaults["freq"] = ""

    return run_and_summarize(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        label=label,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "Gaussian Static"},
        copy_files=copy_files,
    )


class ParsedGSLogData(BaseModel):
    atoms: list[str]
    coordinates: list[list[float]]
    homo_energy: float
    lumo_energy: float
    homo_lumo_gap: float
    molecular_charge: int
    molecular_multiplicity: int
    molecular_mass: float
    molecular_scfenergy: float
    atoms_charges: dict[str, list[float]]

    @classmethod
    def from_dict(cls, data: dict) -> ParsedGSLogData:
        return cls(**data)

    def to_xyz_string(self) -> str:
        """
        Convert the molecular data to XYZ file format string.

        Returns:
            str: XYZ formatted string representation of the molecule
        """
        # First line: number of atoms
        num_atoms = len(self.atoms)

        # Second line: comment with molecular information
        comment = f"HOMO={self.homo_energy:.6f} LUMO={self.lumo_energy:.6f} Gap={self.homo_lumo_gap:.6f} Charge={self.molecular_charge} Mult={self.molecular_multiplicity}"

        # Build the XYZ content
        xyz_lines = [str(num_atoms), comment]

        # Add atom coordinates
        for i in range(num_atoms):
            atom = self.atoms[i]
            x, y, z = self.coordinates[i]
            xyz_lines.append(f"{atom:<2} {x:12.6f} {y:12.6f} {z:12.6f}")

        # 使用统一的换行符连接所有行
        return "\n".join(xyz_lines)

    def to_ase_atoms(self):
        """
        Convert the parsed data to ASE Atoms object.

        Returns:
            Atoms: ASE Atoms object representation of the molecular data
        """
        from ase import Atoms

        return Atoms(symbols=self.atoms, positions=self.coordinates)

    @classmethod
    def from_logfile(cls, logfile: str) -> ParsedGSLogData | None:
        """
        从高斯日志文件解析并创建ParsedGSLogData实例

        参数:
            logfile: 高斯日志文件路径

        返回:
            ParsedGSLogData实例或None（解析失败时）
        """
        import logging

        import cclib

        # 设置cclib的日志级别为WARNING，隐藏INFO消息
        cclib_logger = logging.getLogger("cclib")
        cclib_logger.setLevel(logging.WARNING)

        from himatcal.atoms.core import elements

        try:
            data = cclib.io.ccread(logfile)
            if data is None:
                raise ValueError("No data found in the log file.")

            logger.debug(f"Available attributes: {', '.join(dir(data))}")
            required_attrs = {
                "atomnos": "atomic numbers",
                "atomcoords": "atomic coordinates",
                "atomcharges": "atomic charges",
                "atommasses": "atomic masses",
                "homos": "HOMO indices",
                "moenergies": "molecular orbital energies",
                "charge": "molecular charge",
                "mult": "multiplicity",
                "scfenergies": "SCF energies",
            }
            for attr, desc in required_attrs.items():
                if not hasattr(data, attr):
                    raise ValueError(f"Missing {desc} in the Gaussian output file")
                logger.debug(f"{desc} ({attr}): {getattr(data, attr)}")
            atoms = [elements[atomic_num] for atomic_num in data.atomnos]
            logger.debug(f"Atoms: {atoms}")

            final_coords = data.atomcoords[-1].tolist()

            logger.debug(f"Final coordinates: {final_coords}")

            # parse homo lumo energy
            try:
                homo_energy = float(data.moenergies[-1][data.homos[-1]])
                lumo_energy = float(data.moenergies[-1][data.homos[-1] + 1])

                logger.debug(f"HOMO energy: {homo_energy}, LUMO energy: {lumo_energy}")
            except IndexError as e:
                logger.error("Index error when accessing moenergies or homos: %s", e)
                raise ValueError("Error calculating HOMO/LUMO energies.")

            parsed_data = ParsedGSLogData(
                atoms=atoms,
                coordinates=final_coords,
                homo_energy=homo_energy,
                lumo_energy=lumo_energy,
                homo_lumo_gap=lumo_energy - homo_energy,
                molecular_charge=data.charge,
                molecular_multiplicity=data.mult,
                molecular_mass=sum(data.atommasses),
                molecular_scfenergy=float(data.scfenergies[-1]),
                atoms_charges={
                    charge_type: charges.tolist()
                    for charge_type, charges in data.atomcharges.items()
                },
            )

            # 打印解析出的数据
            logger.debug(f"Parsed data: {parsed_data}")

            return parsed_data
        except Exception as e:
            logger.error(f"An error occurred while parsing the log file: {e}")
            return None  # 返回None而不是空字典，保持类型一致性
