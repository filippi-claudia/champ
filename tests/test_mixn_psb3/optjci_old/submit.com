#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -N 15
#SBATCH -n 360
#SBATCH --ntasks-per-node 24
#SBATCH -J cn5-cas
#SBATCH -o vmc.%j.out
#SBATCH -e vmc.%j.err
# SBATCH --constraint=haswell
#SBATCH -p short
# SBATCH -p broadwell
##################################################################################
# Load the modules to create the proper environment
# Please note that namd needs the intel mpi implementation which is loaded via the module
# surfsara. We do  a "module purge" first, because slurm imports the list of active modules
# from the interactive environment
##################################################################################
module purge
module load surfsara
#module load g09
module list

##################################################################################
# Define testcase
inputdir=`pwd`

srun /home/cuzzocre/programs/champ-dev-fix-lind/bin/vmc.mov1 < vmc.inp  > vmc.out
