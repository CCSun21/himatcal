#!/bin/bash
#SBATCH -J gmx   ##作业名
#SBATCH -p normal   ##队列
#SBATCH -n 4  ##每节点进程数

source /usr/local/gromacs/bin/GMXRC.bash
# gromacs minimization
## generate the tpr file for gromacs
gmx editconf -f init.pdb -o conf.gro
gmx grompp -f min.mdp -o min.tpr
## run gromacs
mpirun -np 2 gmx mdrun -nt 4 -ntomp 2 -deffnm min

# gromacs nvt1
## generate the tpr file for gromacs
gmx grompp -f nvt1.mdp -c min.gro -o nvt1.tpr
## run gromacs
mpirun -np 2 gmx mdrun -nt 4 -ntomp 2 -v -deffnm nvt1

# gromacs npt1
## generate the tpr file for gromacs
gmx grompp -f npt1.mdp -c nvt1.gro -o npt1.tpr
## run gromacs
mpirun -np 2 gmx mdrun -nt 4 -ntomp 2 -v -deffnm npt1

# gromacs npt2
## generate the tpr file for gromacs
gmx grompp -f npt2.mdp -c npt1.gro -o npt2.tpr
## run gromacs
mpirun -np 2 gmx mdrun -nt 4 -ntomp 2 -v -deffnm npt2

# gromacs npt3
## generate the tpr file for gromacs
gmx grompp -f npt3.mdp -c npt2.gro -o npt3.tpr
## run gromacs
mpirun -np 2 gmx mdrun -nt 4 -ntomp 2 -v -deffnm npt3

