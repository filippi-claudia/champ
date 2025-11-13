---
layout: default
title: Desktop
nav_order: 1
parent: Installation
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - installation
    - desktop
---

# Installation on Ubuntu Desktop

Ubuntu 22.04 onwards:

## Using Intel oneAPI compilers

Setup the build:
```
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifort -DENABLE_TREXIO=yes
```

To run the code with Intel Compilers and MPI:
```bash
mpirun -np 24  champ/bin/vmc.mov1 -i input.inp -o output.out -e error
```

## Using GNU compilers

Install the required packages:
```bash
sudo apt install gfortran openmpi-bin libopenmpi-dev gawk libblacs-mpi-dev liblapack-dev
```
Set-up the build:
```bash
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpifort
```
Build:
```bash
cmake --build build -- -j2
```
To run in parallel:
```bash
mpirun -n 2 path_to_CHAMP/bin/vmc.mov1 -i vmc.inp -o vmc.out -e error
```
