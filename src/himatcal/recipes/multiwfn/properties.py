from __future__ import annotations

import re
import subprocess
from pathlib import Path


def get_polar(logfile: Path | str):
    """Get polarizability from Multiwfn output."""
    logfile = Path(logfile)
    if not logfile.exists():
        raise FileNotFoundError(f"Logfile {logfile} does not exist.")
    if not logfile.is_file():
        raise ValueError(f"Logfile {logfile} is not a file.")
    logfile_path = str(logfile.absolute())

    # 启动 Multiwfn 进程
    process = subprocess.Popen(
        ["Multiwfn"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    input_commands = f"{logfile_path}\n24\n1\n1\n0\n0\nq\n"
    out, err = process.communicate(input=input_commands)
    dipole_pattern = r"Dipole\s+moment.*?\n.*?X,Y,Z=\s*([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)\s+Norm=\s*([-+]?\d*\.\d+)"

    dipole_match = re.search(dipole_pattern, out, re.DOTALL)

    if dipole_match:
        return {
            "x": float(dipole_match.group(1)),
            "y": float(dipole_match.group(2)),
            "z": float(dipole_match.group(3)),
            "norm": float(dipole_match.group(4)),
        }

    return None
