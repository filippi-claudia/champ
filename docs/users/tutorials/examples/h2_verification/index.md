---
title: H2 Verification Tests
tags:
    - tutorial
    - VMC
    - optimization
    - H2
---

# H2 Verification Tests

This example demonstrates Variable Monte Carlo (VMC) verification tests on the Hydrogen molecule ($H_2$) using Slater-Type Orbitals (STOs). It covers basic energy evaluation and wavefunction optimization using different methods.

## System Configuration

*   **Molecule**: Hydrogen ($H_2$)
*   **Basis Set**: STO-CVB1 (Slater-Type Orbital)
*   **Method**: VMC
*   **Mode**: `vmc_one_mpi` (Standard VMC)

## 1. VMC Energy Evaluation

The first test evaluates the total energy of the $H_2$ molecule using a standard VMC run with Stochastic Reconfiguration parameters.

**Input File**: `revised_vmc.inp`

```perl
%module general
    title        'H2'
    pool         './pool/'
    basis        sto_cvb1
    mode         'vmc_one_mpi'
%endmodule

load molecule        $pool/champ_v3_h2.xyz
load basis_num_info  $pool/champ_v3_sto_cvb1_basis_pointers.bfinfo

load determinants    rhf.det
load orbitals        champ_v3_trexio_order_rhf.sto_cvb1.lcao
load jastrow         jastrow.0
load jastrow_der     jastrow.der

%module electrons
    nup           1
    nelec         2
%endmodule

%module optwf
    ioptwf        0       # No optimization, just evaluation
    method        'sr_n'
    nblk_max      200
    sr_tau        0.05
    sr_eps        0.001
    sr_adiag      0.01
%endmodule

%module blocking_vmc
    vmc_nstep     50
    vmc_nblk      200
    vmc_nblkeq    1
%endmodule
```

**Expected Energy**: ~ -1.0046 Ha

## 2. Linear Method Optimization

This test uses the **Linear Method** (`lin_d`) to optimize the Jastrow factor (`ioptjas 1`).

**Input File**: `revised_vmc_lin.inp`

```perl
%module optwf
    ioptwf        1       # Optimization enabled
    ioptjas       1       # Optimize Jastrow
    method        'lin_d' # Linear method with diagonal approximation
    lin_eps       0.1
    lin_adiag     0.01
    lin_jdav      1
    lin_nvec      2       # Number of vectors for subspace
%endmodule
```

## 3. Optimization Check (Correlated Sampling)

This scenario performs an optimization run checking wavefunction updates.

**Input File**: `revised_vmc_corsamp.inp`

```perl
%module optwf
    ioptwf        1
    ioptci        1        # Optimize CI coefficients
    method        'linear' # General Linear method
    multiple_adiag 1
%endmodule
```

## Input Files

The required input files (geometry, basis, orbitals) are located in the `pool/` directory and can be found in the [H2 CI Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-H2) folder.
