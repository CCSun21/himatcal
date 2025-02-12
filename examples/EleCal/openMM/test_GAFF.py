import openmm.app
from openff.toolkit import Molecule
from openff.units import Quantity, unit
from openmmforcefields.generators import GAFFTemplateGenerator
from rich.pretty import pprint

molecule = Molecule.from_smiles("O=S(=O)(N)c1c(Cl)cc2c(c1)S(=O)(=O)NCN2")
molecule.generate_conformers(n_conformers=1)

topology = molecule.to_topology()
topology.box_vectors = Quantity([4, 4, 4], unit.nanometer)

gaff = GAFFTemplateGenerator(molecules=molecule)

# If using this alongside another force field, we'd load that in here. But
# in this case, we're only using GAFF for this one molecule, so a "blank"
# force field is a sufficient starting point
forcefield_gaff = openmm.app.ForceField()
forcefield_gaff.registerTemplateGenerator(gaff.generator)

system_gaff = forcefield_gaff.createSystem(
    topology=topology.to_openmm(),
    nonbondedMethod=openmm.app.PME,
)
