---
title: Water TREXIO Optimization
tags:
    - tutorial
    - VMC
    - optimization
    - TREXIO
    - H2O
---

# Water Optimization with TREXIO

This tutorial demonstrates a full wavefunction optimization (Jastrow + Orbitals) for Water ($H_2O$) using the TREXIO file format input. This represents the modern, streamlined workflow for CHAMP calculations.

## System Configuration

*   **Molecule**: Water ($H_2O$)
*   **Source Format**: TREXIO (HDF5)
*   **Method**: VMC Optimization (`sr_n`)
*   **Optimization Targets**: Jastrow factor + Orbital coefficients

## Input File Setup

**Input File**: `vmc_h2o_dft_optall.inp`

```perl
%module general
    title           'H2O DFT calculation'
    pool            './pool/'
    mode            vmc_one_mpi
%endmodule

# Single command to load all wavefunction data
load trexio          H2O_DFT.hdf5

# Load Jastrow parameters
load jastrow         jastrow.dft_optimal_2body
load jastrow_der     jastrow.der

%module electrons
    nup           4
    nelec         8
%endmodule

%module optwf
    ioptwf        1       # Enable optimization
    ioptci        0       # Keep CI coefficients fixed (single determinant)
    ioptjas       1       # Optimize Jastrow
    ioptorb       1       # Optimize Orbitals

    method        'sr_n'  # Stochastic Reconfiguration
    ncore         0
    nextorb       100     # Number of external orbitals to use for rotation
    nblk_max      4000
    nopt_iter     1       # 1 iteration for testing (increase for production)
    sr_tau        0.05
    sr_eps        0.01
    sr_adiag      0.01
%endmodule

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      1000
%endmodule
```

## Key Features

1.  **Compact Input**: Geometry, Basis, Determinants, and Orbitals are all loaded via `load trexio`.
2.  **Orbital Optimization**: `ioptorb 1` allows relaxation of the DFT orbitals within the VMC variational principle.
3.  **DFT Starting Point**: The input `H2O_DFT.hdf5` contains orbitals generated from a DFT calculation.

Resources for this test can be found in the [Water TREXIO CI Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-TREXIO-H2O-DFT-optall) folder.
