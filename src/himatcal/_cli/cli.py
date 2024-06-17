from typing import Optional

import typer
from rich import print
from typing_extensions import Annotated

app = typer.Typer(
    name="himatcal",
    help="A CLI for Himatcal project",
    no_args_is_help=True,
    rich_markup_mode="rich",
)


def machine_name():
    """
    Returns a list of machine names.

    Returns:
        list: A list of machine names.
    """
    return ["local", "sghpc1", "yeesuan", "XMU"]


@app.command("hello", help="Command to say hello :red_heart-emoji:")
def hello(name: Annotated[Optional[str], typer.Argument(help="test command")] = None):
    """
    test command
    """
    if name is not None:
        print(f"Hello, {name}!")
        return f"Hello, {name}!"
    else:
        print("Hello, world!")


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


@app.command("cancel", help="Cancel a covalent job", no_args_is_help=True)
def cancel_job(job_id: Annotated[str, typer.Argument(help="Job id to cancel")]):
    """
    Cancel a covalent job
    """
    import covalent as ct

    ct.cancel(job_id)  # TODO: seems like this is not working
    print(f"Job with id {job_id} cancelled successfully")


@app.command("test_exec")

@app.command("status", help="Get the status of a covalent job", no_args_is_help=True)
def job_status(job_id: Annotated[str, typer.Argument(help="Job id to get status")]):
    """
    Get the status of a covalent job
    """
    import covalent as ct

    result = ct.get_result(job_id)
    status = result.status
    print(f"Job with id {job_id} is {status}")
    return status


@app.command("redispatch", help="Redispatch a covalent job", no_args_is_help=True)
def redispatch_job(
    job_id: Annotated[str, typer.Argument(help="Job id to redispatch")],
    reuse: Annotated[bool, typer.Option(help="Reuse the previous job")] = False,
):
    """
    Redispatch a covalent job
    """
    import covalent as ct

    redispatch_func = ct.redispatch(
        job_id,
        reuse_previous_results=reuse,
    )
    redispatch_id = redispatch_func()
    print(f"Job redispathed to {redispatch_id} (reuse_previous_results={reuse})")


@app.command(
    "extract_result", help="Extract the result of a covalent job", no_args_is_help=True
)
def extrcat_result(
    job_id: Annotated[str, typer.Argument(help="Job id to extract result")],
):
    """
    Extract the result of a covalent job
    """
    from himatcal.utils.ct import extract_result

    results = extract_result(job_id)
    return results


@app.command(
    "cas", help="Get molecular structure from CAS number", no_args_is_help=True
)
def cas(
    cas: Annotated[str, typer.Argument(help="CAS number")],
    write: Annotated[str, typer.Option(help="Write the structure to file")] = "xyz",
):
    """
    Get molecular structure from CAS number
    """
    from himatcal.recipes.mol.core import get_molecular_structure

    atoms = get_molecular_structure(molecular_cas=cas, write_to_file=write)
    print(atoms)


if __name__ == "__main__":
    app()
