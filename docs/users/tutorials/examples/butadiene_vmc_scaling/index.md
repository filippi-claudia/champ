---
title: Butadiene VMC Scaling
tags:
    - tutorial
    - VMC
    - scaling
    - large systems
---

# Butadiene VMC Determinant Scaling

This set of tutorials/tests demonstrates VMC calculations on Butadiene with varying sizes of determinant expansions: 500, 5000, and 15000 determinants. This is useful for benchmarking and verifying code performance with increasing wavefunction complexity.

## Configurations

All tests use BFD pseudopotentials and the BFD-T basis set.

### 500 Determinants
A standard test case.
[Folder: VMC-C4H6-ci1010_pVTZ-00500-dets](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-C4H6-ci1010_pVTZ-00500-dets)

### 5000 Determinants
Intermediate scaling.
[Folder: VMC-C4H6-ci1010_pVTZ-05000-dets](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-C4H6-ci1010_pVTZ-05000-dets)

### 15000 Determinants
Large scale test.
[Folder: VMC-C4H6-ci1010_pVTZ-15000-dets](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-C4H6-ci1010_pVTZ-15000-dets)

## Typical Input (`vmc.inp`)

Example for 15,000 determinants:

```perl
%module general
    title        'butadiene'
    pool         'pool/'
    pseudopot    BFD
    basis        BFD-T
    mode         vmc_one_mpi
%endmodule

load molecule        $pool/champ_v3_butadiene.xyz
load basis_num_info  $pool/champ_v3_BFD-T_basis_pointers.bfinfo

load determinants    TZ_1M_15k.det
load orbitals        champ_v3_trexio_order_ci1010_pVTZ_1_orbitals.lcao
load jastrow         jastrow_good_b3lyp.0
load jastrow_der     jastrow.der
load symmetry        ci1010_pVTZ_1_symmetry.sym

%module electrons
    nup           11
    nelec         22
%endmodule

%module optwf
    ioptwf        1
    ioptci        1
    ioptjas       1
    ioptorb       1
    method        'sr_n'
    ncore         0
    nextorb       280
    no_active     0
    nblk_max      100
    nopt_iter     1
    sr_tau        0.025
    sr_eps        0.001
    sr_adiag      0.01
%endmodule

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      100
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```
