"""functions using for visualization"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms


def view_atoms(atoms: Atoms, fmt="xyz"):
    """
    View atoms using a 3D molecular viewer.

    Args:
        atoms: The molecular structure to visualize.
        fmt: The format in which to save the molecular structure (default is "xyz").

    Returns:
        A 3D molecular viewer displaying the atoms.

    Raises:
        FileNotFoundError: If the temporary file "tmp_atoms" is not found.
    """

    import py3Dmol
    from ase.io import write

    Path(".cache").mkdir(exist_ok=True)
    tmp_path = Path(".cache/tmp_atoms")
    write(tmp_path, atoms, format=fmt)
    atoms_data = Path.open(tmp_path).read()
    view = py3Dmol.view(width=800, height=400)
    view.addModel(atoms_data, format)
    view.setStyle({"stick": {}})
    view.zoomTo()
    return view


def show_xyz_mol(xyz_file: Path):
    """
    Visualize a stk molecule using py3Dmol.
    """
    import py3Dmol

    mol = Path.open(xyz_file).read()
    p = py3Dmol.view(
        data=mol,
        style={"stick": {"colorscheme": "Jmol"}},
        width=400,
        height=400,
    )
    p.setBackgroundColor("white")
    p.zoomTo()
    p.show()


def xyz_to_mol(xyz_file: Path, write_mol=True):
    """
    Convert a xyz file to a mol file and block.
    """
    from openbabel import pybel as pb

    # ! openbabel is a conda package, try other packages if openbabel is not available.
    mol = next(pb.readfile("xyz", xyz_file))
    if write_mol:
        mol.write("mol", f"{xyz_file}.mol", overwrite=True)
        return Path(f"{xyz_file}.mol").open().read()
    return None


def init_style():
    """
    Use the science style for matplotlib plots.

    This function sets the style and font family for matplotlib plots to a predefined science style.

    Args:
        None

    Returns:
        None
    """

    import matplotlib.pyplot as plt
    import pkg_resources

    plt.style.use(
        pkg_resources.resource_filename("himatcal", "tools/science-1.mplstyle")
    )
    plt.rcParams["font.family"] = "Calibri, Microsoft YaHei"


def mpl_view_atoms(atoms):
    """
    view atoms using matplotlib at top and side view
    """
    import matplotlib.pyplot as plt
    from ase.visualize.plot import plot_atoms

    fig, axs = plt.subplots(1, 2, figsize=(5, 5))
    plot_atoms(atoms, axs[0])
    plot_atoms(atoms, axs[1], rotation=("-90x"))
