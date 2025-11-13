---
layout: default
title: Fugaku
nav_order: 5
parent: Installation
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - installation
    - Fugaku
---

# Installation on **Fugaku** (fugaku.r-ccs.riken.jp) Riken Supercomputer

Here are a couple of recipes for commonly used computing facilities, which can be easily adapted.

To compile the code, first load the required modules:

```bash
. /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load cmake@3.24.3%fj@4.8.1/p5qsrqc
spack load fujitsu-mpi@head%fj@4.8.1
spack load hdf5@1.12.2%fj@4.8.1/tpglq6h
spack load fujitsu-ssl2@head%fj@4.8.1/nndozbk
```

then set-up the build:

```bash
export MPIFC='mpifrt'    # Fujitsu Fortran compiler
export MPICC='mpifcc'    # Fujitsu C compiler

cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=${MPIFC} -DCMAKE_C_COMPILER=${MPICC}
```

and finally build:
```bash
cmake --build build -j --clean-first
```

To run the code, you need to submit a job to the queue system:
```bash
pjsub job.sh
```

where `job.sh` is a SLURM job script. Here are some sample scripts:

## Sample VMC job script

```bash
#!/bin/sh
#PJM -L node=1
#PJM -L rscgrp=small-s5
#PJM -L elapse=00:20:00
#PJM -g hp230349
#PJM -x PJM_LLIO_GFSCACHE=/vol0004

. /vol0004/apps/oss/spack/share/spack/setup-env.sh

module load fujitsu-mpi/head-fj-4.8.1-gncxc6a
module load fujitsu-ssl2/head-fj-4.8.1-r3hdjbl
module load hdf5/1.12.2-fj-4.8.1-tpglq6h

# Launch MPI code
mpiexec champ-main/bin/vmc.mov1 -i vmc.inp
```
