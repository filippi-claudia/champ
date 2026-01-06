---
title: Butadiene CSF Optimization
tags:
    - tutorial
    - VMC
    - optimization
    - CSFs
---

# Butadiene CSF Optimization

This tutorial demonstrates the simultaneous optimization of Jastrow, Orbitals, and Configuration State Functions (CSFs) for Butadiene.

## System Configuration

*   **Molecule**: Butadiene
*   **Expansion**: 20 Determinants + 12 CSFs
*   **Basis**: BFD-Q (Quadruple-zeta)

## Input File Setup

**Input File**: `vmc_optall_ci44.inp`

```perl
%module general
    title        'butadiene ci44 vmc calculation with 20 dets and 12 csf'
    pseudopot    BFD
    basis        BFD-Q
    mode         'vmc_one_mpi'
%endmodule

load determinants    champ_v2_gamess_ci44_pVQZ_cart.det
load orbitals        champ_v3_trexio_order_champ_v2_gamess_ci44_pVQZ_cart.lcao
load symmetry        champ_v2_gamess_ci44_pVQZ_cart.sym
load jastrow         jastrow_good_b3lyp.0

%module optwf
    ioptwf        1
    ioptci        1       # Optimize CSF coefficients
    ioptjas       1       # Optimize Jastrow
    ioptorb       1       # Optimize Orbitals
    method        'sr_n'
    ncore         0
    nextorb       280
%endmodule
```

## Note

CSFs are linear combinations of determinants that respect the spin and spatial symmetry of the state. Optimizing CSFs instead of individual determinants preserves these symmetries.

Resources: [Butadiene CSF Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-C4H6-ci44_pVQZ_20_dets_12_csf_optall)
