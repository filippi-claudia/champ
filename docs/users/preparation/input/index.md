---
title: CHAMP Input File Examples
tags:
    - input files
    - templates
    - VMC
    - DMC
---

# CHAMP Input File Examples

This section provides templates for common CHAMP input configurations. These are based on the continuous integration tests (`tests/CI_test`) and represent validated, working input structures.

## 1. VMC Energy

A basic VMC run to calculate the energy of a system without any optimization.

**File**: `vmc_energy.inp`

```perl
%module general
    title        'H2 Energy'
    pool         './pool/'
    basis        sto_cvb1
    mode         'vmc_one_mpi'
    nloc         0                # All-electron calculation
%endmodule

# Molecular System
load molecule        $pool/h2.xyz
load basis_num_info  $pool/basis_pointers.bfinfo

# Wavefunction components
load determinants    rhf.det
load orbitals        rhf.sto_cvb1.lcao
load jastrow         jastrow.0
load jastrow_der     jastrow.der

%module electrons
    nup           1
    nelec         2
%endmodule

%module blocking_vmc
    vmc_nstep     50    # Steps per block
    vmc_nblk      200   # Number of blocks
    vmc_nblkeq    1     # Equilibration blocks
%endmodule
```

## 2. VMC Wavefunction Optimization

Optimization of Jastrow factor, orbitals, and CI coefficients using the Stochastic Reconfiguration (`sr_n`) method.

**File**: `vmc_optimization.inp`

```perl
%module general
    title        'Butadiene Optimization'
    pool         'pool/'
    pseudopot    BFD
    basis        BFD-T
    mode         vmc_one_mpi
%endmodule

load molecule        $pool/butadiene.xyz
load basis_num_info  $pool/BFD-T_basis_pointers.bfinfo

load determinants    ci1010_pVTZ_determinants.det
load orbitals        ci1010_pVTZ_orbitals.lcao
load jastrow         jastrow_good_b3lyp.0
load jastrow_der     jastrow.der
load symmetry        ci1010_pVTZ_symmetry.sym

%module electrons
    nup           11
    nelec         22
%endmodule

%module optwf
    ioptwf        1       # Optimization enabled
    ioptci        1       # Optimize CI coefficients
    ioptjas       1       # Optimize Jastrow factor
    ioptorb       1       # Optimize Orbitals
    method        'sr_n'  # Stochastic Reconfiguration
    ncore         0       # Core orbitals
    nextorb       280     # Virtual orbitals included in rotation
    nblk_max      100
    nopt_iter     1       # Number of optimization steps
    sr_tau        0.025
    sr_eps        0.001
    sr_adiag      0.01
%endmodule

%module blocking_vmc
    vmc_nstep     20      # Steps per block
    vmc_nblk      100     # Number of blocks
%endmodule
```

## 3. Diffusion Monte Carlo

A standard fixed-node DMC calculation. Requires a good trial energy (`etrial`).

**File**: `dmc.inp`

```perl
%module general
    title        'Butadiene DMC'
    pool         'pool/'
    pseudopot    BFD
    basis        BFD-T
    mode         'dmc_one_mpi1'   # DMC mode
%endmodule

load molecule        $pool/butadiene.xyz
load basis_num_info  $pool/BFD-T_basis_pointers.bfinfo

load determinants    TZ_1M_500.det
load orbitals        ci1010_pVTZ_orbitals.lcao
load jastrow         jastrow_good_b3lyp.0
load jastrow_der     jastrow.der
load symmetry        ci1010_pVTZ_symmetry.sym

%module electrons
    nup           11
    nelec         22
%endmodule

%module blocking_dmc
    dmc_nstep     60
    dmc_nblk      10
    dmc_nblkeq    1
    dmc_nconf     20      # Target walkers per thread
%endmodule

%module dmc
    tau           0.05    # Time step $\tau$ (in a.u.)
    etrial      -26.3d0   # Trial Energy
    icasula      -1       # Flag for T-move algorithm
%endmodule
```

## 4. Input with TREXIO

Using the HDF5-based TREXIO format to load wavefunction data significantly simplifies input files.

**File**: `vmc_trexio.inp`

```perl
%module general
    title           'H2O TREXIO'
    mode            vmc_one_mpi
%endmodule

# Main data source
load trexio          H2O_DFT.hdf5

# Additional components can still be loaded explicitly
load jastrow         jastrow.dft_optimal_2body
load jastrow_der     jastrow.der

%module optwf
    ioptwf        1
    ioptci        0
    ioptjas       1
    ioptorb       1
    method        'sr_n'
    nextorb       100
    nblk_max      4000
    sr_tau        0.05
%endmodule

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      1000
%endmodule
```

## 5. State-Specific Optimization

Optimizing for a specific excited state (or multiple states) typically involves defining state weights and norms.

**File**: `vmc_state_specific.inp`

```perl
%module general
    title           'HNO State-Specific'
    pool            './pool/'
    pseudopot       BFD
    basis           BFD
    mode            'vmc_one_mpi'
    
    # State-Specific Parameters
    nstates         2
    weights         [ 1.0d0,  1.0d0 ]
    weights_guiding [ 1.0d0,  1.0d0 ]
    sr_lambda       [ 1.0d0 ]
    anorm           [ 1.0d0, 0.32 ]
%endmodule

%module mstates
    iguiding        2
    iefficiency     1
%endmodule

load molecule        $pool/hno.xyz
load basis_num_info  $pool/BFD-Dan_basis_pointers.bfinfo

load determinants    start.det
load orbitals        start.orb
load jastrow         start.jas
load jastrow_der     jastrow.der

%module electrons
    nup           6
    nelec         12
%endmodule

%module optwf
    ioptwf        1 
    ioptci        1
    ioptjas       1
    ioptorb       1
    method        'sr_n'
    nblk_max      250
%endmodule

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      250
%endmodule
```
