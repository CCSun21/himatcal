from __future__ import annotations

import logging

from ase import Atoms
from ase.build import molecule

from himatcal.recipes.crest.core import relax

logger = logging.getLogger(__name__)

def test_relax(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")
    result = relax(atoms, 0, 1, "gfn2")
    logger.info(f"Relaxed atoms: {result}")
    assert type(result) is Atoms
