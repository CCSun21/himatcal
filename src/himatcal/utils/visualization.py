# view the slab using py3Dmol
import py3Dmol
from ase.io import read, write

def view_atoms(atoms, format='xyz'):
    """
    View the atoms in jupyter notebook using py3Dmol.
    """
    write('tmp_atoms', atoms, format=format)
    atoms_data = open('tmp_atoms', 'r').read()
    view = py3Dmol.view(width=800, height=400)
    view.addModel(atoms_data, format)
    view.setStyle({'stick':{}})
    view.zoomTo()
    return view

def show_xyz_mol(xyz_file):
    """
    Visualize a stk molecule using py3Dmol.
    """
    mol = open(xyz_file).read()
    p = py3Dmol.view(
        data = mol,
        style = {'stick': {'colorscheme': 'Jmol'}},
        width=400,
        height=400,
    )
    p.setBackgroundColor('white')
    p.zoomTo()
    p.show()

def xyz_to_mol(xyz_file, write_mol=True):
    """
    Convert a xyz file to a mol file and block.
    """
    from openbabel import pybel as pb
    mol = next(pb.readfile('xyz', xyz_file))
    if write_mol:
        mol.write('mol', f'{xyz_file}.mol', overwrite=True)
        return open(f'{xyz_file}.mol').read()

# * Gasiteger charge visualization
def plot_gasteiger_charges(mol):
    """
    Plot Gasteiger charges on a molecule.
    """
    from rdkit.Chem.Draw import SimilarityMaps
    from rdkit.Chem import AllChem
    AllChem.ComputeGasteigerCharges(mol)
    contribs = [float(mol.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(mol.GetNumAtoms())]
    SimilarityMaps.GetSimilarityMapFromWeights(
        mol, contribs, colorMap='jet', contourLines=10
    )