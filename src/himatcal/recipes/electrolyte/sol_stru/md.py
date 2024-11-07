from __future__ import annotations

import subprocess
from pathlib import Path

import yaml

mdp_min = """
title                      = minimization
; run control
integrator                 = steep      ; steepest descents minimization
emtol                      = 100        ; stop when max force < 100.0 kJ/mol/nm
emstep                     = 0.01       ; energy step size
nsteps                     = 50000     ; max number of steps
; output control
nstxout                    = 1000       ; steps per position save
nstlog                     = 1000       ; steps per log entry
nstenergy                  = 1000       ; steps per energy file entry
; cutoffs
cutoff-scheme              = Verlet     ; Buffered neighbor searching
nstlist                    = 10         ; neighbor list update frequency
ns_type                    = grid       ; neighbor list method (simple, grid)
pbc                        = xyz        ; periodic boundary conditions
coulombtype                = PME        ; method for electrostatics
rcoulomb                   = 1.2        ; Short-range electrostatic cutoff
rvdw                       = 1.2        ; van der Waals cutoff
"""

mdp_nvt1 = """title                     = 0.01 ns NVT
; Start time and timestep in ps
integrator                = md
nsteps                    = 10000   ; 10 ps
dt                        = 0.001   ; set timestep to 1fs
nstxout	                  = 0
nstvout	                  = 0
nstenergy                 = 5000
nstlog                    = 5000
continuation              = no      ; first dynamics run
constraint_algorithm      = Lincs   ; holonomic constraints
constraints               = none    ; no constraints
lincs_iter                = 1
lincs_order               = 4
; cutoffs
cutoff-scheme             = Verlet
ns_type                   = grid
nstlist                   = 10
rcoulomb                  = 1
rvdw                      = 1
coulombtype               = PME
pme_order                 = 4
fourierspacing            = 0.16
; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling
Tcoupl                    = V-rescale
nsttcouple                = -1
nh-chain-length           = 1
; Groups to couple separately
tc-grps                   = System
; Time constant (ps) and reference temperature (K)
tau_t                     = 1
ref_t                     = 298.15
; Pressure coupling
Pcoupl                    = no
Pcoupltype                = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau_p                     = 0.5
compressibility           = 5e-5
ref_p                     = 1.
refcoord_scaling          = com
pbc                       = xyz
Dispcorr                  = Enerpres
gen_vel                   = no
"""

mdp_npt1 = """title                     = 0.1ns npt berendsen
; run control
integrator                = md
dt                        = 0.001            ; 1 fs
nsteps                    = 10000        ; steps, 0.1ns
; Ooutput control
nstxout                   = 0
nstvout                   = 0
nstfout                   = 0
nstlog                    = 5000
nstenergy                 = 1000
nstxout-compressed        = 1000
;Vel
gen_vel                   = no              ; generate velocities
;gen_temp                  = 300
;gen_seed                  = -1
; cutoffs
cutoff-scheme             = Verlet
pbc                       = xyz         ; periodic boundary conditions
coulombtype               = PME        ; method for electrostatics
rcoulomb                  = 1.0        ; Short-range electrostatic cutoff
rvdw                      = 1.0        ; van der Waals cutoff
DispCorr                  = EnerPres    ; long-distance contributions to E, P
; temperature control
tcoupl                    = v-rescale ; velocity rescaling thermostat
tc_grps                   = system    ; coupled to entire system
tau_t                     = 0.2        ; time constant
ref_t                     = 298.15        ; temperature (K)
; pressure control
pcoupl                    = C-rescale ; barostat type: Parrinello-Rahman, berendsen,
pcoupltype                = isotropic
tau_p                     = 4            ; time constant
ref_p                     = 1.0            ; pressure (bar)
compressibility           = 5e-5        ; pressure bath compressibility (of water, bar^-1)
; bond parameters
continuation              = yes            ; restart from NVT
; constraints               = hbonds            ; h-bond lengths constrained
; constraint_algorithm      = LINCS    ; constraint method
"""

mdp_npt2 = """
title                     = 1ns npt berendsen 1.8T annealing 0.25ns 500K 0.25ns 298K
; run control
integrator                = md
dt                        = 0.001          ; 1 fs
nsteps                    = 1000000        ; steps, 1 ns
; Ooutput control
nstxout                   = 0
nstvout                   = 0
nstfout                   = 0
nstlog                    = 5000
nstenergy                 = 1000
nstxout-compressed        = 1000
;Vel
gen_vel                   = no             ; generate velocities
;gen_temp                  = 300
;gen_seed                  = -1
; cutoffs
cutoff-scheme             = Verlet
pbc                       = xyz            ; periodic boundary conditions
coulombtype               = PME            ; method for electrostatics
rcoulomb                  = 1.0            ; Short-range electrostatic cutoff
rvdw                      = 1.0            ; van der Waals cutoff
DispCorr                  = EnerPres       ; long-distance contributions to E, P
; temperature control
tcoupl                    = v-rescale      ; velocity rescaling thermostat
tc_grps                   = system         ; coupled to entire system
tau_t                     = 0.2            ; time constant
ref_t                     = 298.15            ; temperature (K)
; pressure control
pcoupl                    = Parrinello-Rahman ; barostat type: Parrinello-Rahman, berendsen,
pcoupltype                = isotropic
tau_p                     = 4               ; time constant
ref_p                     = 1.0             ; pressure (bar)
compressibility           = 5e-5            ; pressure bath compressibility (of water, bar^-1)
; bond parameters
continuation              = yes             ; restart from NVT
; constraints               = hbonds        ; h-bond lengths constrained
; constraint_algorithm      = LINCS         ; constraint method

; annealing
annealing                 = single          ; single or double
annealing_npoints         = 4               ; number of points
annealing_time            = 0 250 500 750       ; time points
annealing_temp            = 298.15 500 500 298.15  ; temperature points
"""

mdp_npt3 = """
title                     = 1ns npt Parrinello-Rahman
; run control
integrator                = md
dt                        = 0.001            ; 1 fs
nsteps                    = 1000000        ; steps, 1ns
; Ooutput control
nstxout                   = 0
nstvout                   = 0
nstfout                   = 0
nstlog                    = 5000
nstenergy                 = 1000
nstxout-compressed        = 1000
;Vel
gen_vel                   = no
;gen_temp                 = 300
;gen_seed                 = -1
; cutoffs
cutoff-scheme             = Verlet
pbc                       = xyz         ; periodic boundary conditions
coulombtype               = PME        ; method for electrostatics
rcoulomb                  = 1.0        ; Short-range electrostatic cutoff
rvdw                      = 1.0        ; van der Waals cutoff
DispCorr                  = EnerPres    ; long-distance contributions to E, P
; temperature control
tcoupl                    = v-rescale ; velocity rescaling thermostat
tc_grps                   = system    ; coupled to entire system
tau_t                     = 0.2        ; time constant
ref_t                     = 298.15        ; temperature (K)
; pressure control
pcoupl                    = Parrinello-Rahman ; barostat type: Parrinello-Rahman, berendsen,
pcoupltype                = isotropic
tau_p                     = 4            ; time constant
ref_p                     = 1.0            ; pressure (bar)
compressibility           = 5e-5        ; pressure bath compressibility (of water, bar^-1)
; bond parameters
continuation              = yes            ; restart from NVT
;constraints              = hbonds            ; h-bond lengths constrained
;constraint_algorithm     = LINCS    ; constraint method

"""

slurm_title = """#!/bin/bash
#SBATCH -J gmx   ##作业名
#SBATCH -p normal   ##队列
#SBATCH -n 4  ##每节点进程数

source /usr/local/gromacs/bin/GMXRC.bash
"""

def slurm_min(pre_step: str = "init.pdb", current_step: str = "min") -> str:
    return f"""# gromacs minimization
## generate the tpr file for gromacs
gmx editconf -f {pre_step} -o conf.gro
gmx grompp -f {current_step}.mdp -o {current_step}.tpr
## run gromacs
mpirun -np 2 gmx mdrun -nt 4 -ntomp 2 -deffnm min

"""


def slurm_md(pre_step: str = "min", current_step: str = "nvt1") -> str:
    return f"""# gromacs {current_step}
## generate the tpr file for gromacs
gmx grompp -f {current_step}.mdp -c {pre_step}.gro -o {current_step}.tpr
## run gromacs
mpirun -np 2 gmx mdrun -nt 4 -ntomp 2 -v -deffnm {current_step}

"""


def gmx_min(pdb_file: str, submit_job=False, **calc_kwards) -> None:
    """
    Run Gromacs minimization
    """
    # Define file paths
    mdp_path = Path("min.mdp")
    slurm_path = Path("min.slurm")
    pdb_path = Path("init.pdb")

    # Write the mdp file
    mdp_path.write_text(mdp_min)

    # Write the slurm file
    slurm_path.write_text(slurm_title + slurm_min())

    # Write the pdb file
    pdb_path.write_text(pdb_file)

    # * check the topo file
    if not Path("topol.top").exists():
        raise FileNotFoundError(
            "topol.top file not found, please run gmx_write_topo() first"
        )

    # Run the slurm job
    if submit_job:
        subprocess.run(["sbatch", str(slurm_path)], check=True, text=True)


def gmx_write_topo(compound_dict, itp_folder: str = None) -> None:
    with Path("topol.top").open("w") as f:
        f.write("[ defaults ]\n")
        f.write("  ;nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n")
        f.write("    1        3          yes        0.5      0.5\n")
        f.write(f'#include "{itp_folder}/head.itp"\n')
        for compound in compound_dict:
            f.write(f'#include "{itp_folder}/{compound}.itp"\n')
        f.write("\n\n[ system ]")
        f.write("\n; Name\n")
        for compound in compound_dict:
            f.write(f"{compound}-")
        f.write("\n[ molecules ]\n")
        f.write("; Compound        #mols\n")
        for compound in compound_dict:
            f.write(f'{compound}           {compound_dict[compound]["mol"]}\n')


def build_box(data, compound_dict, mol_path):
    from himatcal.recipes.electrolyte.sol_stru.build_box import ElectrolyteBuilder

    density = data["density"]
    box = data["box"]
    capital = True

    builder = ElectrolyteBuilder(mol_path=mol_path)
    builder.build_box(
        box_electrolyte_composition=compound_dict,
        density=density,
        box=box,
        save_path="init.pdb",
        capital=capital,
    )

def gmx_solvation_md(yaml_file, submit_job=False):
    with Path.open(Path(yaml_file)) as f:
        data = yaml.safe_load(f)
        compound_dict = data["compound_dict"]
        mdp_files = data["mdp_files"]
        build_pdb = data["build_pdb"]
        mol_path = data["pdb_folder"]
        itp_folder = data["itp_folder"]

    # * write the pdb file
    if build_pdb:
        build_box(data, compound_dict, mol_path)
    else:
        pdb_file = data["pdb_file"]
        # * copy  pdb_file to init.pdb
        Path(pdb_file).rename("init.pdb")

    # * write the topol.top file
    gmx_write_topo(compound_dict, itp_folder)

    # * write the mdp files from the list
    for mdp_file in mdp_files:
        mdp_path = Path(f"{mdp_file}.mdp")
        if mdp_file == "min":
            mdp_path.write_text(mdp_min)
        elif mdp_file == "nvt1":
            mdp_path.write_text(mdp_nvt1)
        elif mdp_file == "npt1":
            mdp_path.write_text(mdp_npt1)
        elif mdp_file == "npt2":
            mdp_path.write_text(mdp_npt2)
        elif mdp_file == "npt3":
            mdp_path.write_text(mdp_npt3)
        else:
            raise ValueError(f"mdp file {mdp_file} not found")

    # * write the slurm file for each step in the list
    slurm_path = Path("md.slurm")
    slurm_content = slurm_title
    if "min" in mdp_files:
        slurm_content = slurm_content + slurm_min()
    for i in range(len(mdp_files) - 1):
        slurm_content = slurm_content + slurm_md(mdp_files[i], mdp_files[i + 1])
    slurm_path.write_text(slurm_content)

    # Run the slurm job
    if submit_job:
        result = subprocess.run(
            ["sbatch", str(slurm_path)], check=True, text=True, capture_output=True
        )
        # print(result.stdout)
        # print(result.stderr)  # TODO: capture the slurm job id and return it
