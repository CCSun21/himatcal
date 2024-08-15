"""crest recipes"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal
from pathlib import Path

import toml


def write_general_input(
    input_type: str = "struc.xyz",
    runtype: Literal[
        "none",
        "ancopt",
        "optimize",
        "ancopt_ensemble",
        "optimize_ensemble",
        "md",
        "mtd",
        "dynamics",
        "metadynamics",
        "mecp",
        "mecp_search",
    ] = "none",
    threads: int = 4,
    preopt: bool = False,
    topo: bool = "False",
):
    with Path.open(Path("input.toml"),"w") as f:
        toml.dump(
            {
                "input": input_type,
                "runtype": runtype,
                "threads": threads,
                "preopt": preopt,
                "topo": topo,
                "calculation":{
                    "type" : "mecp",
                    "hess_update": "bfgs",
                    "level":{
                        "method":"xtb",
                    }
                }
            },
            f,
        )
