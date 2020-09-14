#!/bin/bash
#SBATCH -p short                   # partition (queue)
#SBATCH -n 5                     # number of cores
#SBATCH -t 00:05:00                 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out        # STDOUT
#SBATCH -e slurm.%N.%j.err        # STDERR

module load pre2019
module load intel/2018b mpi/impi/18.0.4

SRC=$HOME/champ/bin
JOBNAME=$(echo $1 | sed s/\.inp//g)
srun  $SRC/vmc.mov1  < $JOBNAME.inp > $JOBNAME.out
