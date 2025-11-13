---
layout: default
title: LUMI
nav_order: 4
parent: Installation
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - installation
    - LUMI
---

# Installation on **LUMI** (lumi.csc.fi) Supercomputer

Here are a couple of recipes for commonly used computing facilities, which can be easily adapted.

To compile the code, first load the required modules:

```bash
module swap PrgEnv-cray PrgEnv-gnu
module load LUMI/23.03
module load cray-hdf5-parallel
module load buildtools/23.03
```

then set-up the build:

```bash
export TREXIO_DIR=/users/shindera
export QMCKL_DIR=/users/shindera

cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_C_COMPILER=cc -DENABLE_TREXIO=ON -DENABLE_QMCKL=ON
```

and finally build:
```bash
cmake --build build -j --clean-first
```

To run the code, you need to submit a job to the queue system:
```bash
sbatch job.sh
```

where `job.sh` is a SLURM job script. Here are some sample scripts:

## Sample VMC job script

```bash
#!/bin/bash -l
#SBATCH --job-name=ex01         # Job name
#SBATCH --output=ex01.o%j       # Name of stdout output file
#SBATCH --error=ex01.e%j        # Name of stderr error file
#SBATCH --partition=standard    # Partition (queue) name
#SBATCH --nodes=1               # Total number of nodes
#SBATCH --ntasks=128            # Total number of mpi tasks
#SBATCH --mem=0                 # Allocate all the memory on the node
#SBATCH --time=0:10:00          # Run time (d-hh:mm:ss)
#SBATCH --mail-type=all         # Send email at begin and end of job
#SBATCH --account=              # Project for billing
#SBATCH --reservation=          # Reservation name

# Any other commands must follow the #SBATCH directives
export PMI_NO_PREINITIALIZE=y

# Launch MPI code
srun champ/bin/vmc.mov1 -i vmc_h2o_hf_trexio.inp -o vmc_h2o_hf_trexio.out  -e error
```

## Sample DMC job script

```bash
#!/bin/bash -l
#SBATCH --job-name=ex05         # Job name
#SBATCH --output=ex05.o%j       # Name of stdout output file
#SBATCH --error=ex05.e%j        # Name of stderr error file
#SBATCH --partition=standard    # Partition (queue) name
#SBATCH --nodes=1               # Total number of nodes
#SBATCH --ntasks=128            # Total number of mpi tasks
#SBATCH --mem=0                 # Allocate all the memory on the node
#SBATCH --time=0:05:00          # Run time (d-hh:mm:ss)
#SBATCH --mail-type=all         # Send email at begin and end of job
#SBATCH --account=              # Project for billing
#SBATCH --reservation=          # Reservation for training

# Any other commands must follow the #SBATCH directives
export PMI_NO_PREINITIALIZE=y

# Launch MPI code
srun champ/bin/vmc.mov1 -i vmc_h2o_hf_jas2body_trexio.inp -o vmc_h2o_hf_jas2body_trexio.out -e error

cat mc_configs_new* >> mc_configs
rm mc_configs_new*

srun champ/bin/dmc.mov1 -i dmc_h2o_hf_jas2body_trexio.inp -o dmc_h2o_hf_jas2body_trexio.out -e error
```