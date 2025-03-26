#!/bin/bash
#SBATCH -J gmx   ##作业名
#SBATCH -p normal   ##队列
#SBATCH -n 16  ##每节点进程数
#SBATCH --cpus-per-task=2  ##每任务CPU数

source /usr/local/gromacs/bin/GMXRC.bash
export OMP_NUM_THREADS=2
# gromacs minimization
## generate the tpr file for gromacs
gmx editconf -f init.pdb -o conf.gro
gmx grompp -f min.mdp -o min.tpr
## run gromacs
gmx mdrun -nt 16 -ntomp 2 -deffnm min -s min.tpr -rdd 1.0 -pin on

# gromacs nvt1
## generate the tpr file for gromacs
gmx grompp -f nvt1.mdp -c min.gro -o nvt1.tpr
## run gromacs
gmx mdrun -nt 16 -ntomp 2 -deffnm nvt1 -s nvt1.tpr -rdd 1.0 -pin on

# gromacs npt1
## generate the tpr file for gromacs
gmx grompp -f npt1.mdp -c nvt1.gro -o npt1.tpr
## run gromacs
gmx mdrun -nt 16 -ntomp 2 -deffnm npt1 -s npt1.tpr -rdd 1.0 -pin on

# gromacs npt2
## generate the tpr file for gromacs
gmx grompp -f npt2.mdp -c npt1.gro -o npt2.tpr
## run gromacs
gmx mdrun -nt 16 -ntomp 2 -deffnm npt2 -s npt2.tpr -rdd 1.0 -pin on

# gromacs npt3
## generate the tpr file for gromacs
gmx grompp -f npt3.mdp -c npt2.gro -o npt3.tpr
## run gromacs
gmx mdrun -nt 16 -ntomp 2 -deffnm npt3 -s npt3.tpr -rdd 1.0 -pin on

