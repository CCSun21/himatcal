"""Example of using the relax_job function to relax a structure using the mace-mp-0 recipe."""

from __future__ import annotations

from pathlib import Path

from ase.io import read

from himatcal.calculator.aimnet import AIMNet2ASE

atoms_path = Path.cwd() / "input.xyz"
charge = 0
mult = 1

calc = AIMNet2ASE("aimnet2", charge=charge, mult=mult)

atoms = read(atoms_path)
atoms.calc = calc
print(atoms.get_potential_energy())
