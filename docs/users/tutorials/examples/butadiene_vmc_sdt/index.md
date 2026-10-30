---
title: Butadiene VMC SDT
tags:
    - tutorial
    - VMC
    - optimization
    - excitations
---

# Butadiene VMC with SDT Excitations

This tutorial involves a VMC optimization calculation on Butadiene using a wavefunction that includes Singles, Doubles, and Triples (SDT) excitations (1522 determinants).

## Input File Setup

**Input File**: `vmc_optimization_1522.inp`

```perl
%module general
    title        'butadiene'
    mode         vmc_one_mpi
%endmodule

load molecule        $pool/champ_v3_butadiene.xyz
load basis_num_info  $pool/champ_v3_BFD-T_basis_pointers.bfinfo

load determinants    champ_v2_gamess_ci1010_SDT_pVTZ_determinants.det
load orbitals        champ_v3_trexio_order_champ_v2_gamess_ci1010_SDT_pVTZ_orbitals.lcao
load jastrow         jastrow_good_b3lyp.0
load jastrow_der     jastrow.der
load symmetry        champ_v2_gamess_ci1010_SDT_pVTZ_symmetry.sym

%module optwf
    ioptwf        1
    ioptci        1       # Optimize CI coefs
    ioptjas       1       # Optimize Jastrow
    ioptorb       1       # Optimize Orbitals
    method        'sr_n'
    nopt_iter     1
%endmodule
```

## Significance

This test verifies the stability and performance of the optimization algorithm (`sr_n`) when handling higher-order excitations (Triples) in the determinant expansion.

Resources: [Butadiene SDT Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-C4H6-ci1010_SDT_pVTZ-1522-dets)
