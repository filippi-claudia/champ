---
title: Water 3-Body Jastrow
tags:
    - tutorial
    - VMC
    - jastrow
    - TREXIO
---

# Water with 3-Body Jastrow

This tutorial demonstrates the use of a **3-Body Jastrow factor** (Electron-Electron-Nucleus correlation) in a VMC calculation for Water.

## Input File Setup

**Input File**: `vmc.inp`

```perl
%module general
    title           'H2O VMC energy-only calculation with 3-body Jastrow'
    mode            vmc_one_mpi
    pseudopot       ccECP
%endmodule

load trexio         water_c2v_LDA.hdf5
load jastrow        jastrow_3body.jas    # Explicit 3-body Jastrow file
load determinants   single.det

%module blocking_vmc
    vmc_nstep     60
    vmc_nblk      1000
%endmodule
```

## Jastrow Files

The `jastrow_3body.jas` file typically contains sections for all three correlation types:
1.  e-n (1-body)
2.  e-e (2-body)
3.  e-e-n (3-body)

Including 3-body terms significantly improves the quality of the wavefunction and reduces the variance.

Resources: [Water 3-Body Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-TREXIO-H2O-3body-Jastrow)
