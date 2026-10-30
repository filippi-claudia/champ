---
title: Water DMC with TREXIO
tags:
    - tutorial
    - DMC
    - TREXIO
    - H2O
---

# Water DMC with TREXIO

This tutorial demonstrates a Diffusion Monte Carlo (DMC) calculation for the Water molecule using the TREXIO input format and a 2-body Jastrow factor. This complements the VMC optimization tutorial.

## System Configuration

*   **Molecule**: Water ($H_2O$)
*   **Format**: TREXIO (HDF5)
*   **Method**: DMC (`dmc_one_mpi1`)
*   **Jastrow**: 2-body

## Input File Setup

**Input File**: `dmc_h2o_dft_jas2body_tau0.05.inp`

```perl
%module general
    title           'H2O DFT calculation'
    pool            './pool/'
    mode            'dmc_one_mpi1'
%endmodule

load trexio          H2O_DFT.hdf5
load orbitals        champ_v3_trexio_order_orbitals.lcao
load basis_num_info  $pool/champ_v3_champ_v2_trexio_H2O_DFT_with_g_basis_pointers.bfinfo
load jastrow         jastrow.final_optall

%module electrons
    nup           4
    nelec         8
%endmodule

%module blocking_dmc
    dmc_nstep     60
    dmc_nblk      40
    dmc_nblkeq    1
    dmc_nconf     100
%endmodule

%module dmc
    tau           0.05d0
    etrial      -17.26d0
    icasula      -1
%endmodule
```

## Description

This test runs a standard DMC simulation. The `load trexio` command is used along with explicit loading of orbitals and basis info (demonstrating how to override or supplement TREXIO data if needed, though often TREXIO is self-contained).

Resources: [DMC TREXIO Water Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/DMC-TREXIO-water-DFT-jas2body_tau0.05)
