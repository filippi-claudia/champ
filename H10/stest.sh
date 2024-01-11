#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 64
###SBATCH -c 64
#SBATCH -J C60stest 
#SBATCH -o C60stest.%j.out
#SBATCH -e C60stest.%j.err
#SBATCH --exclusive
#SBATCH -p ccp22
##################################################################################
# Load the modules to create the proper environment
# Please note that namd needs the intel mpi implementation which is loaded via the module
# surfsara. We do  a "module purge" first, because slurm imports the list of active modules
# from the interactive environment
##################################################################################
module purge
module load icc
module load mpi
module load mkl
module load hdf5
module load hdf5/1.8.12

##################################################################################

cd $PWD

##export OMP_NUM_THREADS=64

#srun --cpus-per-task 1 --cpu-bind=cores --ntasks 64 ../../../bin/dmc.mov1 -i vmc_optjas1.inp -o vmc_optjas1.out -e error_optjas1

mpirun -np 64 -ppn 64 ../bin/vmc.mov1 -i vmc_test.inp -o vmc_test.out -e error_test
