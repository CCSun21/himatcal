from __future__ import annotations

import json
import os
import subprocess

from ase.atoms import Atoms
from ase.io import write
from calculator.gaussian import Gaussiancal
from options_software.software import get_software_json_path

from himatcal.recipes.electrolyte.sol_stru._base import gen_resp2_chg

mdp_min = """
;define                  = -DFLEXIBLE
; run control
integrator      = steep                 ; steepest descents minimization
emtol           = 100.0                 ; stop when max force < 100.0 kJ/mol/nm
emstep          = 0.01                  ; energy step size
nsteps          = 50000                 ; max number of steps


; output control
nstxout     = 500                       ; steps per position save
nstlog      = 500                       ; steps per log entry
nstenergy   = 50                        ; steps per energy file entry

; cutoffs
cutoff-scheme   = Verlet
nstlist                 = 10            ; neighbor list update frequency
ns_type                 = grid          ; neighbor list method (simple, grid)
pbc                     = xyz           ; periodic boundary conditions

coulombtype             = PME           ; method for electrostatics
rcoulomb                = 1.2           ; Short-range electrostatic cutoff
rvdw                    = 1.2           ; van der Waals cutoff
"""
mdp_npt = """
; run control
integrator      = md
dt              = 0.001                 ; 2 fs
nsteps          = 100000                ; steps, =100 ps

; Ooutput control
nstxout     = 500                       ; steps per position save
nstlog      = 500                       ; steps per log entry
nstenergy   = 50                        ; steps per energy file entry

; cutoffs
cutoff-scheme   = Verlet
nstlist                 = 10            ; neighbor list update frequency
ns_type                 = grid          ; neighbor list method (simple, grid)
pbc                     = xyz           ; periodic boundary conditions

coulombtype             = PME           ; method for electrostatics
rcoulomb                = 1.2           ; Short-range electrostatic cutoff
rvdw                    = 1.2           ; van der Waals cutoff
DispCorr                = EnerPres      ; long-distance contributions to E, P

; temperature control
tcoupl          = v-rescale ; velocity rescaling thermostat
tc_grps         = system        ; coupled to entire system
tau_t           = 1.0           ; time constant
ref_t           = 300           ; temperature (K)

; pressure control
pcoupl          = berendsen ; barostat type
tau_p           = 10.0                  ; time constant
ref_p           = 1.0                   ; pressure (bar)
compressibility = 4.5e-5                ; pressure bath compressibility (of water, bar^-1)

; bond parameters
continuation    = yes                   ; restart from NVT
constraints     = h-bonds                       ; h-bond lengths constrained
constraint_algorithm = LINCS    ; constraint method
"""

mdp_nvt = """
; run control
integrator      = md
dt              = 0.002                 ; 2 fs
nsteps          = 50000         ; steps, = 100 ps

; Ooutput control
nstxout     = 500                       ; steps per position save
nstlog      = 500                       ; steps per log entry
nstenergy   = 50                        ; steps per energy file entry

; cutoffs
cutoff-scheme   = Verlet
nstlist                 = 10            ; neighbor list update frequency
ns_type                 = grid          ; neighbor list method (simple, grid)
pbc                     = xyz           ; periodic boundary conditions

coulombtype             = PME           ; method for electrostatics
rcoulomb                = 1.2           ; Short-range electrostatic cutoff
rvdw                    = 1.2           ; van der Waals cutoff
DispCorr                = EnerPres      ; long-distance contributions to E, P

; temperature control
tcoupl          = v-rescale ; velocity rescaling thermostat
tc_grps         = system        ; coupled to entire system
tau_t           = 1.0           ; time constant
ref_t           = 300           ; temperature (K)

; pressure control
pcoupl          = no                    ; no pressure control in NVT

; velocity generation
gen_vel         = yes                   ; generate velocities from Maxwell distribution
gen_temp        = 300                   ; temperature for Maxwell distribution

; bond parameters
continuation    = no                    ; first dynamics run
constraints     = h-bonds                       ; h-bond lengths constrained
constraint_algorithm = LINCS    ; constraint method
"""
mdp_pro = """
; Run control
integrator               = md
tinit                    = 0
dt                       = 0.002
nsteps                   = 2500000    ; 10 ns

; Output control
nstxout                  = 5000
nstlog                   = 5000
nstenergy                = 500

; cutoffs
cutoff-scheme            = verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz

coulombtype              = PME
rcoulomb                 = 1.2
vdwtype                  = cutoff
rvdw                     = 1.2
DispCorr                = EnerPres

; Temperature coupling
tcoupl                   = v-rescale
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = 300

; Pressure coupling
Pcoupl                   = berendsen
tau_p                    = 10.
compressibility          = 4.5e-05
ref_p                    = 1.0

; bond parameters
continuation    = yes                   ; restart from NPT
constraints     = h-bonds                       ; h-bond lengths constrained
constraint_algorithm = LINCS    ; constraint method
"""


class DielectricConstantCal:
    def __init__(
        self,
        label: str,
        molecular_data: Atoms,
        charge: int,
        mult: int,
        mem="64GB",
        np=64,
    ):
        self.label = label
        self.molecular_data = molecular_data
        self.charge = charge
        self.mult = mult
        self.mem = mem
        self.np = np

    def get_dielect_constant(self):
        opt_method = "b3lyp"
        opt_basis = "6-31G*"
        gas_basis = "def2tzvp"
        solvent = "Acetone"
        work_path = os.getcwd()

        # *  opt at the b3lyp/6-31G* level
        opt_atoms = Gaussiancal(
            method=opt_method, basis=opt_basis, mem=self.mem, nprocshared=self.np
        ).run(
            self.molecular_data,
            self.charge,
            self.mult,
            cal_type="em=GD3BJ opt",
            label=f"{self.label}-opt",
        )
        # gas_sp
        Gaussiancal(
            method=opt_method, basis=gas_basis, mem=self.mem, nprocshared=self.np
        ).run(
            opt_atoms,
            self.charge,
            self.mult,
            cal_type="em=GD3BJ",
            label=f"{self.label}-gas",
            chk=f"./{self.label}-gas.chk",
        )
        # sol_sp
        sol_atoms = Gaussiancal(
            method=opt_method, basis=gas_basis, mem=self.mem, nprocshared=self.np
        ).run(
            opt_atoms,
            self.charge,
            self.mult,
            cal_type="em=GD3BJ freq",
            label=f"{self.label}-sol",
            scrf=f"SMD, solvent={solvent}",
            chk=f"./{self.label}-sol.chk",
        )
        with open(get_software_json_path(), "r") as f:
            sobtoppath = json.load(f)["sobtop_path"]
        self.resp2_chg(f"./{self.label}-gas.chk", f"./{self.label}-sol.chk")
        with open(get_software_json_path(), "r") as f:
            obabelpath = json.load(f)["obabel_path"]
        write(f"./{self.label}-tmp.xyz", sol_atoms)
        subprocess.run(
            f"{obabelpath} -ixyz ./{self.label}-tmp.xyz -omol2 -O ./{self.label}.mol2",
            shell=True,
        )
        subprocess.run(
            f"cp ./{self.label}-gas.fchk ./{self.label}-sol.fchk ./{self.label}.chg ./{self.label}.mol2 {sobtoppath}",
            shell=True,
        )
        os.chdir(sobtoppath)
        sobtopexe = "./sobtop"
        process = subprocess.Popen(
            sobtopexe,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        process.stdin.write(f"./{self.label}.mol2\n")
        process.stdin.write(
            f"7\n10\n{self.label}.chg\n0\n2\n\n1\n2\n7\n{self.label}-sol.fchk\n\n\n0\n"
        )
        out, err = process.communicate()
        subprocess.run(
            f"cp ./{self.label}.gro ./{self.label}.top ./{self.label}.itp {work_path}",
            shell=True,
        )


        os.chdir(work_path)
        with open(get_software_json_path(), "r") as f:
            gromacpath = json.load(f)["gromacs_path"]
        subprocess.run(
            f"{gromacpath} insert-molecules -ci {self.label}.gro -nmol 1200 -box 5 5 5 -o {self.label}-box.gro",
            shell=True,
        )
        box_num = int(
            int(open(f"./{self.label}-box.gro", "r").readlines()[1].rstrip().split()[0])
            / len(self.molecular_data)
        )
        top_file = open(f"./{self.label}.top", "r").readlines()
        top_file_result = ""
        for i in top_file:
            if self.label in i.rstrip().split() and len(i.rstrip().split()) > 1:
                top_file_result = top_file_result + f"{self.label}\t{box_num}\n"
            else:
                top_file_result = top_file_result + i
        open(f"./{self.label}.top", "w").write(top_file_result)
        open("./min.mdp", "w").write(mdp_min)
        open("./nvt.mdp", "w").write(mdp_nvt)
        open("./npt.mdp", "w").write(mdp_npt)
        open("./pro.mdp", "w").write(mdp_pro)
        subprocess.run(
            f"mkdir min; cp {self.label}.top {self.label}.itp min.mdp {self.label}-box.gro min ; cd min ; {gromacpath} grompp -f min.mdp  -c {self.label}-box.gro  -p {self.label}.top -o min.tpr; {gromacpath} mdrun -nt {self.np} -deffnm min; cd ..",
            shell=True,
        )
        subprocess.run(
            f"mkdir nvt; cp {self.label}.top {self.label}.itp nvt.mdp nvt ; cd nvt ; cp ../min/min.gro . ; {gromacpath} grompp -f nvt.mdp  -c min.gro  -p {self.label}.top -o nvt.tpr; {gromacpath} mdrun -nt {self.np} -deffnm nvt; cd ..",
            shell=True,
        )
        subprocess.run(
            f"mkdir npt; cp {self.label}.top {self.label}.itp npt.mdp npt ;cd npt ; cp ../nvt/nvt.gro . ; {gromacpath} grompp -f npt.mdp  -c nvt.gro  -p {self.label}.top -o npt.tpr -maxwarn 1; {gromacpath} mdrun -nt {self.np} -deffnm npt; cd ..",
            shell=True,
        )
        subprocess.run(
            f"mkdir pro; cp {self.label}.top {self.label}.itp pro.mdp pro ;cd pro ; cp ../npt/npt.gro . ; {gromacpath} grompp -f pro.mdp  -c npt.gro  -p {self.label}.top -o pro.tpr -maxwarn 1; {gromacpath} mdrun -nt {self.np} -deffnm pro; cd ..",
            shell=True,
        )
        subprocess.run(f"cp pro/pro.trr pro/pro.tpr .", shell=True)
        process = subprocess.Popen(
            f"{gromacpath} current -f pro.trr -s pro.tpr",
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            shell=True,
        )
        process.stdin.write(f"0\n")
        out, err = process.communicate()
        open(f"./{self.label}-constant.out", "w").write(out + "\n" + err)
        result = open(f"./{self.label}-constant.out", "r").readlines()
        for i in range(0, len(result)):
            if "Absolute" in result[i].rstrip().split():
                constant = result[i + 1].rstrip().split("=")[1]
                break
        return float(constant)

    def resp2_chg(self, gas_chk, sol_chk):
        with open(get_software_json_path(), "r") as f:
            formchkpath = json.load(f)["formchk_path"]
        subprocess.run(f"{formchkpath} {self.label}-gas.chk", shell=True)
        subprocess.run(f"{formchkpath} {self.label}-sol.chk", shell=True)
        with open(get_software_json_path(), "r") as f:
            multiwfnpath = f"{json.load(f)['multiwfn_path']}"
        process = subprocess.Popen(
            multiwfnpath,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        process.stdin.write(f"{self.label}-gas.fchk\n")
        process.stdin.write("7\n18\n1\ny\n0\n0\nq\n")
        out, err = process.communicate()
        process = subprocess.Popen(
            multiwfnpath,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        process.stdin.write(f"{self.label}-sol.fchk\n")
        process.stdin.write("7\n18\n1\ny\n0\n0\nq\n")
        out, err = process.communicate()
        gas_chg = open(f"{self.label}-gas.chg", "r").readlines()
        sol_chg = open(f"{self.label}-sol.chg", "r").readlines()
        result = ""
        for gas_chg_line, sol_chg_line in zip(gas_chg, sol_chg):
            gas_chg = float(gas_chg_line.rstrip().split()[4])
            sol_chg = float(sol_chg_line.rstrip().split()[4])
            element = gas_chg_line.rstrip().split()[0]
            x_corr = gas_chg_line.rstrip().split()[1]
            y_corr = gas_chg_line.rstrip().split()[2]
            z_corr = gas_chg_line.rstrip().split()[3]
            chg = (gas_chg + sol_chg) / 2
            result = result + f"{element}\t{x_corr}\t{y_corr}\t{z_corr}\t{chg}\n"
        open(f"./{self.label}.chg", "w").write(result)
