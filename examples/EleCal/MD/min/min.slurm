#!/bin/bash
#SBATCH -J gmx   ##作业名
#SBATCH -p normal   ##队列
#SBATCH -n 4   ##每节点进程数

source /usr/local/gromacs/bin/GMXRC.bash

# generate the tpr file for gromacs
gmx editconf -f init.pdb -o conf.gro
gmx grompp -f min.mdp -o min.tpr

# run gromacs
mpirun -np 1 gmx mdrun -nt 4 -ntomp 4 -deffnm min
