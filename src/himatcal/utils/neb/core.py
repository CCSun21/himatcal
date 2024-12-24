from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ase import Atoms
from ase.filters import FrechetCellFilter
from ase.io import Trajectory
from ase.mep import NEB
from ase.optimize import MDMin
from quacc.atoms.core import copy_atoms
from quacc.runners.prep import calc_cleanup, calc_setup
from quacc.utils.dicts import recursive_dict_merge

from himatcal.utils.neb import Geodesic, redistribute

if TYPE_CHECKING:
    from typing import Any

    from ase.mep.neb import NEBOptimizer
    from ase.optimize.optimize import Optimizer
    from quacc.types import (
        Filenames,
        SourceDirectory,
    )


def run_neb(
    images: list[Atoms],
    relax_cell: bool = False,
    fmax: float = 0.01,
    max_steps: int | None = 1000,
    optimizer: NEBOptimizer | Optimizer = MDMin,
    optimizer_kwargs: dict[str, Any] | None = None,
    run_kwargs: dict[str, Any] | None = None,
    neb_kwargs: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> Optimizer:
    """
    Run NEB optimization.

    Parameters
    ----------
    images
        List of images representing the initial path.
    relax_cell
        Whether to relax the unit cell shape and volume.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        Optimizer class to use. All Optimizers except BFGSLineSearch
    optimizer_kwargs
        Dictionary of kwargs for the optimizer.
    run_kwargs
        Dictionary of kwargs for the run() method of the optimizer.
    neb_kwargs
        Dictionary of kwargs for the NEB.
    copy_files
        Files to copy before running the calculation.

    Returns
    -------
    Optimizer
        The ASE Optimizer object.
    """
    if optimizer.__name__ == "BFGSLineSearch":
        raise ValueError("BFGSLineSearch is not allowed as optimizer with NEB.")
    # Copy atoms so we don't modify it in-place
    images = copy_atoms(images)
    # settings = get_settings()

    neb = NEB(images, **neb_kwargs)

    dir_lists = []
    # Perform staging operations
    # this calc_setup function is not suited for multiple Atoms objects
    for image in images:
        tmpdir_i, job_results_dir_i = calc_setup(image, copy_files=copy_files)
        dir_lists.append([tmpdir_i, job_results_dir_i])

    # Set defaults
    optimizer_kwargs = recursive_dict_merge(
        {
            "logfile": dir_lists[0][0] / "opt.log",
            "restart": dir_lists[0][0] / "opt.json",
        },
        optimizer_kwargs,
    )
    run_kwargs = run_kwargs or {}
    traj_filename = "opt.traj"
    # Check if trajectory kwarg is specified
    if "trajectory" in optimizer_kwargs:
        msg = "Quacc does not support setting the `trajectory` kwarg."
        raise ValueError(msg)

    # Define the Trajectory object
    traj_file = dir_lists[0][0] / traj_filename
    traj = Trajectory(traj_file, "w", atoms=neb)

    # Set volume relaxation constraints, if relevant
    for i in range(len(images)):
        if relax_cell and images[i].pbc.any():
            images[i] = FrechetCellFilter(images[i])

    # Run optimization
    dyn = optimizer(neb, **optimizer_kwargs)
    dyn.attach(traj.write)

    dyn.run(fmax, max_steps)
    traj.close()

    traj.filename = traj_file
    dyn.trajectory = traj

    # Perform cleanup operations skipping the first images's directory
    # because that is where the trajectory is stored. It will get deleted
    # eventually.
    for i, image in enumerate(images[1:], start=1):
        calc_cleanup(image, dir_lists[i][0], dir_lists[i][1])

    return dyn


def _geodesic_interpolate_wrapper(
    reactant: Atoms,
    product: Atoms,
    n_images: int = 20,
    perform_sweep: bool | None = None,
    convergence_tolerance: float = 2e-3,
    max_iterations: int = 15,
    max_micro_iterations: int = 20,
    morse_scaling: float = 1.7,
    geometry_friction: float = 1e-2,
    distance_cutoff: float = 3.0,
) -> list[Atoms]:
    """
    Interpolates between two geometries and optimizes the path with the geodesic method.

    Parameters
    ----------
    reactant
        The ASE Atoms object representing the initial geometry.
    product
        The ASE Atoms object representing the final geometry.
    n_images
        Number of images for interpolation. Default is 20.
    perform_sweep
        Whether to sweep across the path optimizing one image at a time.
        Default is to perform sweeping updates if there are more than 35 atoms.
    convergence_tolerance
        Convergence tolerance. Default is 2e-3.
    max_iterations
        Maximum number of minimization iterations. Default is 15.
    max_micro_iterations
        Maximum number of micro iterations for the sweeping algorithm. Default is 20.
    morse_scaling
        Exponential parameter for the Morse potential. Default is 1.7.
    geometry_friction
        Size of friction term used to prevent very large changes in geometry. Default is 1e-2.
    distance_cutoff
        Cut-off value for the distance between a pair of atoms to be included in the coordinate system. Default is 3.0.

    Returns
    -------
    list[Atoms]
        A list of ASE Atoms objects representing the smoothed path between the reactant and product geometries.
    """
    reactant = copy_atoms(reactant)
    product = copy_atoms(product)

    # Read the initial geometries.
    chemical_symbols = reactant.get_chemical_symbols()

    # First redistribute number of images. Perform interpolation if too few and subsampling if too many images are given
    raw_interpolated_positions = redistribute(
        chemical_symbols,
        [reactant.positions, product.positions],
        n_images,
        tol=convergence_tolerance * 5,
    )

    # Perform smoothing by minimizing distance in Cartesian coordinates with redundant internal metric
    # to find the appropriate geodesic curve on the hyperspace.
    geodesic_smoother = Geodesic(
        chemical_symbols,
        raw_interpolated_positions,
        morse_scaling,
        threshold=distance_cutoff,
        friction=geometry_friction,
    )
    if perform_sweep is None:
        perform_sweep = len(chemical_symbols) > 35
    if perform_sweep:
        geodesic_smoother.sweep(
            tol=convergence_tolerance,
            max_iter=max_iterations,
            micro_iter=max_micro_iterations,
        )
    else:
        geodesic_smoother.smooth(tol=convergence_tolerance, max_iter=max_iterations)
    return [
        Atoms(symbols=chemical_symbols, positions=geom)
        for geom in geodesic_smoother.path
    ]
