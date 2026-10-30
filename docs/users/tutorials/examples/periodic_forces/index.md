---
title: Periodic Systems and Forces
tags:
    - tutorial
    - VMC
    - periodic
    - forces
---

# Periodic Systems and Forces

This tutorial demonstrates a VMC calculation for a periodic system (likely LiF based on the title `Test LiF`, though the folder says `Hlattice`) computing atomic forces.

## Input File Setup

**Input File**: `vmc.inp`

```perl
%module general
    title           'Test LiF'
    basis            BASISGRID  # Grid basis
    mode            'vmc_one_mpi'
    iforce_analy     1          # Enable force analysis
    iperiodic        1          # Enable periodic boundary conditions
    np_coul          2          # Ewald summation parameter
    cutg             3.84       # G-vector cutoff
    n_images         2          # Periodic images
%endmodule

load lattice         box.txt    # Define unit cell
load molecule        geo        # Atomic positions

# ... standard loads ...
```

## Key Keywords

*   `iperiodic 1`: Activates periodic handling in the Hamiltonian.
*   `load lattice`: Reads the simulation cell vectors.
*   `iforce_analy 1`: Calculation of Hellmann-Feynman (and Pulay) forces.

Resources: [Periodic Force Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-jas1-Hlattice-force_analy)
