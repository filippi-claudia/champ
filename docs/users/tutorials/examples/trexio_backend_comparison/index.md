---
title: TREXIO Backend Comparison
tags:
    - tutorial
    - TREXIO
    - technical
---

# TREXIO Backend Comparison

CHAMP supports both the **HDF5** (binary, high performance) and **Text** (human-readable) backends of the TREXIO library. This example demonstrates that calculations are invariant to the choice of backend.

## Test Workflow

The calculation performs a short VMC run on Butadiene using identical wavefunctions stored in HDF5 and Text formats.

### 1. HDF5 Run

**Input File**: `vmc_opt_ci1010_pVTZ_1522_hdf5.inp`

```perl
%module general
    title           'butadiene'
    pool            'pool/'
    backend         hdf5
    mode            vmc_one_mpi
    pseudopot	    BFD
%endmodule

load trexio          gamess_butadiene_ci1010_pVTZ.hdf5
load jastrow         jastrow_good_b3lyp.0
load jastrow_der     jastrow.der

%module electrons
   nup     11
   nelec   22
%endmodule 

%module blocking_vmc
    vmc_nstep     10
    vmc_nblk      10
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

### 2. Text Run

**Input File**: `vmc_opt_ci1010_pVTZ_1522_text.inp`

The input is identical except for the `backend` keyword and the file extension.

```perl
%module general
    title           'butadiene'
    backend         text     # text backend
    # ...
%endmodule

load trexio         gamess_butadiene_ci1010_pVTZ.trexio  # Directory name (text format)
# ...
```

## Validation

The total energy and its error bars should check out to be statistically consistent between the two runs.

Resources: [Backend Comparison Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-TREXIO-text-hdf5-backend-comparison)
