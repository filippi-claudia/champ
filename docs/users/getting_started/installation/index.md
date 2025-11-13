---
layout: default
title: CCPGate
nav_order: 2
parent: Installation
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - installation
    - ccpgate
---

# Installation on **CCPGate** (ccpgate.tnw.utwente.nl) Cluster

To build with mpiifort, load the required modules of the Intel Compiler and MPI:

```bash
module load compiler
module load compiler-rt
module load mkl
module load mpi
module load trexio/2.3.0-intel     # Optional
module load python3                # Optional
```
Setup the build:
```
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifort -DENABLE_TREXIO=yes
```

To run the code with Intel Compilers and MPI:
```bash
mpirun -np 24  champ/bin/vmc.mov1 -i input.inp -o output.out -e error
```

To build with gfortran:

Setup the build:

```bash
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=/usr/bin/mpif90
```

which will use LAPACK & BLAS from the Ubuntu repository. (Cmake should find them already if none of the Intel MKL variables are set.) Combining gfortran with the Intel MKL is possible but requires special care to work with the compiler flag `-mcmodel=large`.

To run the code on 60 processors:
```bash
mpirun -np 60 -machinefile "machinefile" -i input.inp -o output.out -e error
```

