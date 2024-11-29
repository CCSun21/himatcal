"""functions for visualization
# TODO:

"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms


from pydantic import BaseModel, ValidationError, field_validator
from typing_extensions import Literal


class InitStyleConfig(BaseModel):
    vis: Literal["mpl", "bokeh"]

    @field_validator("vis")
    def check_vis(cls, v):
        if v not in ["mpl", "bokeh"]:
            raise ValueError("vis must be either 'mpl' or 'bokeh'")
        return v


def init_style(vis: Literal["mpl", "bokeh"] = "mpl"):
    """
    Use the science style for matplotlib plots.

    This function sets the style and font family for matplotlib plots to a predefined science style.

    Args:
        vis (Literal["mpl", "bokeh"]): The visualization library to use. Options are 'mpl' or 'bokeh'.

    Returns:
        None
    """
    try:
        config = InitStyleConfig(vis=vis)
    except ValidationError as e:
        import logging

        logging.error(e)
        return

    import pkg_resources

    if config.vis == "mpl":
        import matplotlib.pyplot as plt

        plt.style.use(
            pkg_resources.resource_filename("himatcal", "tools/science-1.mplstyle")
        )
        plt.rcParams["font.family"] = "Calibri, Microsoft YaHei"
    elif config.vis == "bokeh":
        import json

        from bokeh.io import curdoc
        from bokeh.themes import Theme

        with Path.open(
            Path(pkg_resources.resource_filename("himatcal", "tools/bokeh_theme.json"))
        ) as f:
            theme_dict = json.load(f)
            theme = Theme(json=theme_dict)

        curdoc().theme = theme
        # TODO: Test the bokeh theme


def show_xyz_mol(xyz_file: Path):
    """
    Visualize a stk molecule using py3Dmol.
    """
    import py3Dmol

    with xyz_file.open("r") as file:
        mol = file.read()

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
    from openbabel import pybel as pb  # type: ignore

    # TODO: Find another way to convert xyz to mol
    # ! openbabel is a conda package, try other packages if openbabel is not available.
    mol = next(pb.readfile("xyz", xyz_file))
    if write_mol:
        mol.write("mol", f"{xyz_file}.mol", overwrite=True)
        return Path(f"{xyz_file}.mol").open().read()
    return None


def mpl_view_atoms(atoms):
    """
    view atoms using matplotlib at top and side view
    """
    import matplotlib.pyplot as plt
    from ase.visualize.plot import plot_atoms

    fig, axs = plt.subplots(1, 2, figsize=(5, 5))
    plot_atoms(atoms, axs[0])
    plot_atoms(atoms, axs[1], rotation=("-90x"))


def view_atoms(atoms: Atoms, engine: str = "weaswidget"):
    if engine == "py3Dmol":
        return view_atoms_py3dmol(atoms)
    elif engine == "weaswidget":
        from weas_widget import WeasWidget

        return WeasWidget(from_ase=atoms)
    else:
        raise ValueError("Invalid engine")


def view_atoms_py3dmol(atoms):
    import py3Dmol
    from ase.io import write

    Path(".cache").mkdir(exist_ok=True)
    tmp_path = Path(".cache/tmp_atoms")
    write(tmp_path, atoms, format="xyz")
    atoms_data = Path.open(tmp_path).read()
    view = py3Dmol.view(width=800, height=400)
    view.addModel(atoms_data, format)
    view.setStyle({"stick": {}})
    view.zoomTo()
    return view
