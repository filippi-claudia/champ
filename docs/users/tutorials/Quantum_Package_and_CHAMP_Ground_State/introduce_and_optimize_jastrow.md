---
title: QMC runs (b) Introduce and optimize Jastrow
tags:
    - tutorial
    - CIPSI
    - ground state
    - jastrow
    - QP
---

# Introduce and optimize Jastrow

The Jastrow factor accounts for dynamic electron correlation. We will now introduce a parameterized Jastrow factor and optimize it using the Variation Monte Carlo (VMC) method in CHAMP.

## Jastrow Factor Definition

The Jastrow factor contains three types of terms:

- **Electron-Nucleus ($f_{en}$)**: Describes correlation between electrons and nuclei (1-body).
- **Electron-Electron ($f_{ee}$)**: Describes correlation between electrons (2-body).
- **Electron-Electron-Nucleus ($f_{een}$)**: Describes 3-body correlations.

Since we use pseudopotentials for H and O, the electron-nucleus cusp conditions are less critical (no singularity), but we optimize the polynomial part. For electron-electron terms, we enforce the cusp condition ($b_1$ parameter).

## Starting Jastrow Factor

Create `jastrow.start` with minimal parameters:

```python
jastrow_parameter   1
  5  5  0           norda,nordb,nordc
   0.60000000         scalek
   0.00000000   0.00000000 0. 0. 0. 0. (a(iparmj),iparmj=1,nparma) ! e-n O
   0.00000000   0.00000000 0. 0. 0. 0. (a(iparmj),iparmj=1,nparma) ! e-n H
   0.50000000   1. 0. 0. 0. 0. (b(iparmj),iparmj=1,nparmb) ! e-e
 (c(iparmj),iparmj=1,nparmc) ! e-e-n O
 (c(iparmj),iparmj=1,nparmc) ! e-e-n H
end
```

## Optimization Setup

We need to tell CHAMP which parameters to optimize using a derivative mapping file, `jastrow.der`:

```python
jasderiv
4 4 5 0 0 0 0 nparma,nparmb,nparmc,nparmf
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma) ! e-n O
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma) ! e-n H
2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb) ! e-e
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c(iparmj),iparmj=1,nparmc)
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c(iparmj),iparmj=1,nparmc)
end
```

Update your CHAMP input file layout to include optimization blocks:

```python
%module general
    title           'Jastrow Optimization'
    mode            'vmc_one_mpi1'
%endmodule

load trexio          h2o_hf.trexio
load jastrow         jastrow.start
load jastrow_der     jastrow.der

%module electrons
    nup           4
    nelec         8
%endmodule

%module optwf
    ioptwf        1       # Enable wavefunction optimization
    ioptci        0       # Do not optimize CI coefs yet
    ioptjas       1       # Optimize Jastrow parameters
    ioptorb       0       # Do not optimize orbitals

    method        'sr_n'  # Stochastic Reconfiguration
    nopt_iter     20      # Number of optimization steps
    nblk_max      4000    # Max blocks per step
    ncore         0
    nextorb       100

    sr_tau        0.05
    sr_eps        0.001
    sr_adiag      0.01
%endmodule

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      100     # Fewer blocks sufficient for start
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

Run CHAMP. The code will generate updated Jastrow files (e.g., `jastrow_optimal.1.iter20`). Use the last optimized file for subsequent calculations.