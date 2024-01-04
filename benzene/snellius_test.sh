#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --gpus-per-node=1
#SBATCH -J champ-gpu-test_2
#SBATCH -o champ-gpu-test_2.%j.out
#SBATCH -e champ-gpu-test_2.%j.err
#SBATCH --exclusive
#SBATCH --partition gpu

##################################################################################
# Load the modules to create the proper environment
# Please note that namd needs the intel mpi implementation which is loaded via the module
# surfsara. We do  a "module purge" first, because slurm imports the
# list of active modules
# from the interactive environment
##################################################################################


module purge
module restore gpu_test
module load ScaLAPACK/2.2.0-gompi-2022a-fb
module load CMake/3.23.1-GCCcore-11.3.0




#export OMPI_MCA_pml=ucx


cd $PWD

srun --cpus-per-task 1 --cpu-bind=cores --ntasks 1 ../bin/vmc.mov1 -i vmc_test.inp -o vmc_test.out -e error_test

#srun ../bin/vmc.mov1 -i vmc_test.inp -o vmc_test.out -e error_test 





