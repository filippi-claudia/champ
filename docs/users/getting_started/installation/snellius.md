---
layout: default
title: Snellius
nav_order: 3
parent: Installation
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - installation
    - snellius
---

# Installation on **Snellius** (snellius.surf.nl) Supercomputer

Here are a couple of recipes for commonly used computing facilities, which can be easily adapted.

To compile the code, first load the required modules:

```bash
module purge
module load 2022
module load intel/2022a
module load CMake/3.23.1-GCCcore-11.3.0
module load HDF5/1.12.2-iimpi-2022a
```

then set-up the build:

```bash
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifort
```

and finally build:
```bash
cmake --build build -j8 --clean-first
```

To run the code, you need to submit a job to the queue system:
```bash
sbatch job.sh
```

where `job.sh` is a SLURM job script. Here are some sample scripts:

## Sample VMC job script

```bash
#!/bin/bash
#SBATCH -t 0-12:00:00            # time in (day-hours:min:sec)
#SBATCH -N 1                     # number of nodes
#SBATCH -n 128                   # number of cores
#SBATCH --ntasks-per-node 128    # tasks per node
#SBATCH -J vmc-128-water         # name of the job
#SBATCH -o vmc.%j.out            # std output file name for slurm
#SBATCH -e vmc.%j.err            # std error file name for slurm
#SBATCH --exclusive              # specific requirements about node
#SBATCH --partition thin         # partition (queue)
#
module purge
module load 2022
module load intel/2022a
module load HDF5/1.12.2-iimpi-2022a
#
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
cd $PWD
srun champ/bin/vmc.mov1 -i vmc.inp -o vmc.out -e error
```

## Sample DMC job script

```bash
#!/bin/bash
#SBATCH -t 5-00:00:00            # time in (day-hours:min:sec)
#SBATCH -N 2                     # number of nodes
#SBATCH -n 256                   # number of cores
#SBATCH --ntasks-per-node 128    # tasks per node
#SBATCH -J dmc-256-water         # name of the job
#SBATCH -o dmc.%j.out            # std output file name for slurm
#SBATCH -e dmc.%j.err            # std error file name for slurm
#SBATCH --exclusive              # specific requirements about node
#SBATCH --partition thin         # partition (queue)
#
module purge
module load 2022
module load intel/2022a
module load HDF5/1.12.2-iimpi-2022a
#
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
cd $PWD
srun champ/bin/dmc.mov1 -i dmc.inp -o dmc.out -e error
```