---
tags:
    - setup
    - input files
---

# Preparing the input files

CHAMP needs the following input files to describe a system

1. Molecular coordinates
1. ECP / Pseudopotentials
1. Basis Set (Radial Grid files)
1. Basis pointers
1. MO coefficients
1. Determinants and/or CSF files
1. Molecular orbital symmetries (Optional)
1. Molecular orbital eigenvalues (Optional)
1. Jastrow parameters file
1. Jastrow derivative parameters file (Optional)


CHAMP input file itself has a modular structure. For example,

```
1. general
2. electrons
3. blocking_vmc
4. blocking_dmc
5. optwf
6. ...
```

A typical input file looks like:

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

load determinants    TZ_1M_5k.det
load orbitals        champ_v3_ci1010_pVTZ_1_orbitals.lcao
load symmetry        champ_v3_ci1010_pVTZ_1_symmetry.sym
load jastrow         jastrow_good_b3lyp.0
load jastrow_der     jastrow.der



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
    multiple_adiag 0
    ncore         0
    nextorb       280
    no_active     0
    nblk_max      100
    nopt_iter     1
    sr_tau        0.025
    sr_eps        0.001
    sr_adiag      0.01
    isample_cmat  0
    energy_tol    0.0
%endmodule

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      100
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```
