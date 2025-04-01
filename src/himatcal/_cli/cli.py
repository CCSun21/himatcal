from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated

import typer
from ase.io import read
from monty.os import cd
from rich import print

app = typer.Typer(
    name="himatcal",
    help="A CLI for HiMatCal tool",
    no_args_is_help=True,
    rich_markup_mode="rich",
)


def machine_name() -> list:
    """
    Returns a list of machine names.

    Returns:
        list: A list of machine names.
    """
    return ["local", "sghpc1", "yeesuan", "XMU"]


@app.command("hello", help="Command to say hello :red_heart-emoji:")
def hello(name: Annotated[str | None, typer.Argument(help="test command")] = None):
    """
    test command
    """
    if name is not None:
        print(f"Hello, {name}!")
        return f"Hello, {name}!"
    else:
        print("Hello, world!")
        return "Hello, world!"


@app.command(
    "sub_gs",
    help="Submit a [bold green]Gaussian[/bold green] job",
    no_args_is_help=True,
)
def sub_gs(
    input_file: Annotated[str, typer.Argument(help="Gaussian input file.")],
    label: Annotated[str, typer.Option(help="label for the calculation")] = "Gaussian",
    machine: Annotated[
        str,
        typer.Option(
            help="Machine to run the job, select 'local' and 'sghpc1'.",
            autocompletion=machine_name,
        ),
    ] = "sghpc1",
):
    """
    Submit a Gaussian job
    """
    from himatcal.recipes.gaussian.flow import sub_gs

    dispatch_id = sub_gs(input_file, label, machine)
    print(f"Dispatched job with id: {dispatch_id}")


@app.command("gs", help="Submit a Gaussian job", no_args_is_help=True)
def gs_jobs(
    inputfile: Annotated[str, typer.Argument(help="file as input")],
    command: Annotated[
        str, typer.Argument(help="Command to submit the job")
    ] = "sub_gs",
):
    """
    Submit a Gaussian job
    """
    if command == "sub_gs":
        sub_gs(input_file=inputfile, label="Gaussian", machine="sghpc1")
    elif command == "relax":
        from ase.io import read

        from himatcal.recipes.gaussian.core import relax_job

        atoms = read(inputfile)
        relax_job(
            atoms=atoms,
            charge=0,
            spin_multiplicity=1,
            label="relax",
            xc="B3LYP",
            basis="6-31+g(d) em=GD3BJ",
            freq=False,
        )


@app.command(
    "ts_opt", help="Submit a transition state optimization job", no_args_is_help=True
)
def ts_opt(
    atoms_path: Annotated[str, typer.Argument(help="Atoms to optimize")],
    charge: Annotated[int, typer.Argument(help="Charge of the molecule")],
    mult: Annotated[int, typer.Argument(help="Spin multiplicity of the molecule")],
    label: Annotated[str, typer.Option(help="label for the calculation")] = "TS",
    machine: Annotated[
        str,
        typer.Option(
            help="Machine to run the job, select 'local' and 'sghpc1'.",
            autocompletion=machine_name,
        ),
    ] = "sghpc1",
):
    """
    Submit a transition state optimization job
    """
    from himatcal.recipes.gaussian.flow import ts_opt

    dispatch_id = ts_opt(atoms_path, charge, mult, label, machine)
    print(f"Dispatched job with id: {dispatch_id}")


@app.command("irc", help="Submit a Gaussian IRC calculation", no_args_is_help=True)
def irc(
    atoms_path: Annotated[str, typer.Argument(help="Atoms to optimize")],
    chg: Annotated[int, typer.Argument(help="Charge of the molecule")],
    mult: Annotated[int, typer.Argument(help="Spin multiplicity of the molecule")],
    label: Annotated[str, typer.Option(help="label for the calculation")] = "IRC",
    machine: Annotated[
        str,
        typer.Option(
            help="Machine to run the job, select 'local' and 'sghpc1'.",
            autocompletion=machine_name,
        ),
    ] = "sghpc1",
):
    """
    Submit a Gaussian IRC calculation
    """
    from himatcal.recipes.gaussian.flow import irc_flow

    dispatch_id = irc_flow(atoms_path, chg, mult, label, machine)
    print(f"Dispatched job with id: {dispatch_id}")


@app.command(
    "ts_flow", help="Submit a transition state optimization job", no_args_is_help=True
)
def ts_flow(
    atoms_path: Annotated[str, typer.Argument(help="Atoms to optimize")],
    chg: Annotated[int, typer.Argument(help="Charge of the molecule")],
    mult: Annotated[int, typer.Argument(help="Spin multiplicity of the molecule")],
    label: Annotated[str, typer.Option(help="label for the calculation")] = "TS-flow",
    machine: Annotated[
        str,
        typer.Option(
            help="Machine to run the job, select 'local' and 'sghpc1'.",
            autocompletion=machine_name,
        ),
    ] = "sghpc1",
):
    """
    Submit a transition state optimization job
    """
    from himatcal.recipes.gaussian.flow import ts_flow

    dispatch_id = ts_flow(atoms_path, chg, mult, label)
    print(f"Dispatched job with id: {dispatch_id}")


@app.command(
    "relax", help="Submit a Gaussian relaxation calculation", no_args_is_help=True
)
def relax(
    atoms_path: Annotated[str, typer.Argument(help="The path of atoms file to relax")],
    chg: Annotated[int, typer.Argument(help="Charge of the molecule")],
    mult: Annotated[int, typer.Argument(help="Spin multiplicity of the molecule")],
    label: Annotated[str, typer.Option(help="label for the calculation")] = "relax",
    machine: Annotated[
        str,
        typer.Option(
            help="Machine to run the job, select 'local' and 'sghpc1'.",
            autocompletion=machine_name,
        ),
    ] = "sghpc1",
):
    """
    Submit a Gaussian relaxation calculation
    """
    from himatcal.recipes.gaussian.flow import relax

    dispatch_id = relax(atoms_path, chg, mult, label, machine)
    print(f"Dispatched job with id: {dispatch_id}")


@app.command(
    "cas", help="Get molecular structure from CAS number", no_args_is_help=True
)
def cas(
    cas: Annotated[str, typer.Argument(help="CAS number")],
    write: Annotated[bool, typer.Option(help="Write the structure to file")] = True,
):
    """
    Get molecular structure from CAS number
    """
    from himatcal.recipes.mol.core import get_molecular_structure

    get_molecular_structure(molecular_cas=cas, write_mol=write)


@app.command(
    "cas2xyz",
    help="Get molecular structure from CAS number and convert to xyz file",
    no_args_is_help=True,
)
def cas2xyz(
    cas: Annotated[str, typer.Argument(help="CAS number")],
    relax: Annotated[bool, typer.Option(help="relax the molecule using crest")] = True,
):
    """
    Get molecular structure from CAS number
    """
    from himatcal.recipes.mol.core import cas2xyz

    cas2xyz(cas, relax_atoms=relax)


@app.command("GSM", no_args_is_help=True)
def GSM(
    file: Annotated[str, typer.Argument(help="The file path of the molecule")],
    dc: Annotated[
        str,
        typer.Argument(
            help='driving coordinates with format ["BREAK"/"ADD", atom1, atom2] or ["ANGLE", atom1, atom2, atom3], or ["TORSION", atom1, atom2, atom3, atom4], or ["OOP", atom1, atom2, atom3, atom4]'
        ),
    ],
    calc: Annotated[
        str | None,
        typer.Argument(
            help="The program you use for calc (xtb, gaussian and orca is supported)"
        ),
    ] = "xtb",
    chg: int = 0,
    mult: int = 1,
):
    CWD = Path.cwd()
    cache_path = (
        CWD / f"crest_opt_{datetime.now(timezone.utc).strftime('%Y-%m-%d-%H-%M-%S-%f')}"
    )
    Path.mkdir(cache_path, exist_ok=True)
    atoms = [read(Path(file))]
    if calc == "xtb":
        from himatcal.calculator.xtb import XTB
        from himatcal.recipes.gsm.SE_GSM import ASE_SE_GSM

        gsm = ASE_SE_GSM(
            atoms=atoms,
            driving_coords=eval(dc),
            calculator=XTB(
                method="gfn2-xTB", charge=chg, uhf=mult - 1, gbsa={"solvent": "acetone"}
            ),
        )
    elif calc == "gaussian":
        from ase.calculators.gaussian import Gaussian

        calc = Gaussian(
            charge=-1,
            mult=1,
            label="IMI",
            method="B3LYP",
            basis="6-31G(d)",
            scf="xqc",
            force="",  # * remember to return force
            nosymm="",
            mem="64GB",
            nprocshared=16,
        )
        gsm = ASE_SE_GSM(
            atom=atoms,
            driving_coords=eval(dc),
            calculator=calc,
        )
    elif calc == "orca":
        from ase.calculators.orca import ORCA, OrcaProfile

        profile = OrcaProfile(command="/home/suncc/orca_6_0_0/orca")

        calc = ORCA(
            profile=profile,
            charge=-1,
            mult=1,
            orcasimpleinput="B3LYP g-d3 def2-TZVP EnGrad",  # using EnGrad for force calculation
            orcablocks="%pal nprocs 16 end \n%maxcore 1000",
        )

        gsm = ASE_SE_GSM(
            atom=atoms,
            driving_coords=eval(dc),
            calculator=calc,
        )

    with cd(cache_path):
        gsm.run()


@app.command(
    "draw",
    help="Draw molecule structure with atom and bond indices",
    no_args_is_help=True,
)
def draw_molecule(
    smiles: Annotated[str, typer.Argument(help="SMILES string of the molecule")],
    output: Annotated[
        str, typer.Option("--output", "-o", help="Output file name (PNG format)")
    ] = "molecule.png",
):
    """Draw molecule structure with atom and bond indices"""
    from himatcal.utils.rdkit.core import mol_with_atom_and_bond_indices

    try:
        mol_with_atom_and_bond_indices(smiles, output_file=output)
        print(f"Molecule structure saved to {output}")
    except ValueError as e:
        print(f"Error: {e}")
        raise typer.Exit(1)


@app.command(
    "smi2conf",
    help="Convert SMILES string to XYZ file",
    no_args_is_help=True,
)
def smi2conf(
    smiles: Annotated[str, typer.Argument(help="SMILES string of the molecule")],
    output: Annotated[
        str, typer.Option("--output", "-o", help="Output file name (XYZ format)")
    ] = "molecule.xyz",
):
    from himatcal.utils.rdkit.core import smiles2atoms

    Path().mkdir(exist_ok=True)
    file_path = Path(f"{output}")
    try:
        atoms = smiles2atoms(smiles)
        atoms.write(file_path)
        return True
    except Exception as e:
        print(f"Error converting SMILES to atoms: {e}")
        return False


@app.command(
    "relax_mol",
    help="Relax a molecule using CREST",
    no_args_is_help=True,
)
def relax_mol(
    file,
    chg: int = 0,
    mult: int = 1,
    output: Annotated[
        str, typer.Option("--output", "-o", help="Output file name (XYZ format)")
    ] = "molecule_relaxed.xyz",
):
    """
    Relax a molecule using CREST
    """
    from himatcal.recipes.crest.core import relax

    try:
        atoms_path = Path(f"{file}")
        atoms = relax(read(atoms_path.absolute()), chg, mult)
        atoms.write(f"{output}")
        return True
    except Exception as e:
        print(f"Error relaxing molecule: {e}")
        return False


if __name__ == "__main__":
    app()
