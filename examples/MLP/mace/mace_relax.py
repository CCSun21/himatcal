"""Example of using the relax_job function to relax a structure using the mace-mp-0 recipe."""

from __future__ import annotations

from pathlib import Path

from ase.io import read, write
from quacc.recipes.mlp.core import relax_job

atoms_path = Path.cwd() / "LFP.vasp"
atoms = read(atoms_path)
result = relax_job(atoms, "mace-mp-0", relax_cell=True)
write(f"{atoms_path.stem}_relaxed.vasp", result["atoms"])
write(f"{atoms_path.stem}_relaxed.xyz", result["atoms"])
