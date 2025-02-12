from __future__ import annotations

from pathlib import Path

from himatcal.recipes.electrolyte.sol_stru.build_box import ElectrolyteBuilder
from himatcal.recipes.electrolyte.sol_stru.md import (
    gmx_solvation_md,
)

yaml_file = "/home/suncc/Code/pub/himatcal/examples/EleCal/MD/test/input/input.yaml"

gmx_solvation_md(yaml_file, submit_job=False)
