#!/bin/bash -x
#SBATCH --account=prcoe10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --output=mpi-out.%j
#SBATCH --error=mpi-err.%j
#SBATCH --time=01:00:00
#SBATCH --partition=batch

srun ../../../software/champ/bin/vmc.mov1 < vmc.inp > vmc.out_core48
