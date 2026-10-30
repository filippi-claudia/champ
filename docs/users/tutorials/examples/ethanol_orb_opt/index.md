---
title: Ethanol Orbital Optimization
tags:
    - tutorial
    - VMC
    - optimization
    - QMCKL
    - TREXIO
---

# Ethanol Orbital Optimization

This tutorial focuses on optimizing orbital coefficients for Ethanol ($C_2H_6O$) using VMC. It is based on a test setup designed to verify **QMCkl** integration with TREXIO.

## Input File Setup

**Input File**: `vmc_ethanol_opt.inp`

```perl
%module general
    title        'Ethanol vmc calculation jastrow and orbitals'
    mode         'vmc_one_mpi'
%endmodule

load trexio          ethanol.hdf5

load determinants    ethanol_rhf.det 
load jastrow         jastrow_optimized.jas
load jastrow_der     jastrow.der

%module optwf
    ioptwf        1 
    ioptci        0 
    ioptjas       0       # Jastrow fixed (previously optimized)
    ioptorb       1       # Optimize Orbitals
    method        'sr_n'
    ncore         0
    nextorb       250 
    isample_cmat  0       # Sampling controls
%endmodule
```

## Description

Orbital optimization allows the nodal surface of the wavefunction to change, potentially lowering the fixed-node energy in DMC. Here, we optimize the orbitals starting from an RHF guess.

Resources: [Ethanol OptOrb Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-TREXIO-QMCKL-C2H6O-optorb)
