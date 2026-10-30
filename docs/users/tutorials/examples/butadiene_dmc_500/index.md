---
title: Butadiene DMC
tags:
    - tutorial
    - DMC
    - Butadiene
    - pseudopotential
---

# Butadiene Diffusion Monte Carlo

This tutorial covers a Diffusion Monte Carlo (DMC) calculation on Butadiene ($C_4H_6$) using BFD pseudopotentials and a large determinant expansion (500 determinants).

## System Configuration

*   **Molecule**: Butadiene ($C_4H_6$)
*   **Pseudopotential**: BFD (Burkatzki-Filippi-Dolg)
*   **Basis Set**: BFD-T (Triple-zeta quality)
*   **Wavefunction**: 500 Determinants (from CIPSI)
*   **Method**: DMC

## Input File Setup

**Input File**: `dmc.inp`

```perl
%module general
    title        'butadiene'
    pool         'pool/'
    pseudopot    BFD
    basis        BFD-T
    mode         'dmc_one_mpi1'
%endmodule

# Load molecular geometry
load molecule        $pool/champ_v3_butadiene.xyz

# Load basis set info
load basis_num_info  $pool/champ_v3_BFD-T_basis_pointers.bfinfo

# Load wavefunction components
load determinants    TZ_1M_500.det
load orbitals        champ_v3_trexio_order_ci1010_pVTZ_1_orbitals.lcao
load jastrow         jastrow_good_b3lyp.0
load jastrow_der     jastrow.der
load symmetry        ci1010_pVTZ_1_symmetry.sym

%module electrons
    nup           11
    nelec         22    # 22 valence electrons (C:4*4 + H:6*1)
%endmodule

%module blocking_dmc
    dmc_nstep     60    # Steps per block
    dmc_nblk      10    # Number of blocks
    dmc_nblkeq    1     # Equilibration blocks
    dmc_nconf     20    # Target population per walker file
%endmodule

%module dmc
    tau           0.05      # Component time step
    etrial        -26.3d0   # Trial energy guess (Hartree)
    icasula       -1        # Standard DMC algorithm
%endmodule
```

## Running the Calculation

This calculation uses `dmc_one_mpi1` mode, suitable for MPI parallelization. It requires pre-optimized Jastrow factors (`jastrow_good_b3lyp.0`) and a determinant file (`TZ_1M_500.det`).

Resources for this test can be found in the [Butadiene CI Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/DMC-C4H6-ci1010_pVTZ-500-dets) folder.
