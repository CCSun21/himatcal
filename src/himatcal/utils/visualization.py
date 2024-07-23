""" functions using for visualization """
from __future__ import annotations

from pathlib import Path


def view_atoms(atoms, fmt="xyz"):
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
    write("tmp_atoms", atoms, format=fmt)
    atoms_data = Path.open("tmp_atoms").read()
    view = py3Dmol.view(width=800, height=400)
    view.addModel(atoms_data, format)
    view.setStyle({"stick":{}})
    view.zoomTo()
    return view

def show_xyz_mol(xyz_file):
    """
    Visualize a stk molecule using py3Dmol.
    """
    import py3Dmol
    mol = open(xyz_file).read()
    p = py3Dmol.view(
        data = mol,
        style = {"stick": {"colorscheme": "Jmol"}},
        width=400,
        height=400,
    )
    p.setBackgroundColor("white")
    p.zoomTo()
    p.show()

def xyz_to_mol(xyz_file, write_mol=True):
    """
    Convert a xyz file to a mol file and block.
    """
    from openbabel import pybel as pb
    # ! openbabel is a conda package, try other packages if openbabel is not available.
    mol = next(pb.readfile("xyz", xyz_file))
    if write_mol:
        mol.write("mol", f"{xyz_file}.mol", overwrite=True)
        return open(f"{xyz_file}.mol").read()

# * Gasiteger charge visualization
def plot_gasteiger_charges(mol):
    """
    Plot Gasteiger charges on a molecule.
    """
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import SimilarityMaps
    AllChem.ComputeGasteigerCharges(mol)
    contribs = [float(mol.GetAtomWithIdx(i).GetProp("_GasteigerCharge")) for i in range(mol.GetNumAtoms())]
    SimilarityMaps.GetSimilarityMapFromWeights(
        mol, contribs, colorMap="jet", contourLines=10
    )

def init_style():
    """
    use the science style for matplotlib plots.
    """
    try:
        import matplotlib.pyplot as plt
        import pkg_resources
        plt.style.use(pkg_resources.resource_filename("himatcal", "tools/science-1.mplstyle"))
        plt.rcParams["font.family"] ="Calibri, Microsoft YaHei"
    except ImportError:
        print("matplotlib is not installed.")
