---
title: QMC wave function optimizations
tags:
    - tutorial
    - CIPSI
    - excited state
    - QP
---

# QMC Wavefunction Optimizations

We will optimize the Jastrow factor and then the Coupled Cluster (CI) coefficients for both the Ground State (GS) and Excited State (ES).

## 1. Ground State Optimization

Create a directory for the Ground State calculation and link the TREXIO file:

```bash
mkdir GS_Calc
cd GS_Calc
ln -s ../COH2_GS.trexio .
```

### Initial Jastrow

Create a `jastrow.start` file suitable for $COH_2$ (C, O, H atoms):

```python
jastrow_parameter   1
  5  5  0           norda,nordb,nordc
   0.60000000   0.00000000     scalek,a21
   0.00000000   0.00000000  0. 0. 0. 0. (a(iparmj),iparmj=1,nparma) ! e-n C
   0.00000000   0.00000000  0. 0. 0. 0. (a(iparmj),iparmj=1,nparma) ! e-n O
   0.00000000   0.00000000  0. 0. 0. 0. (a(iparmj),iparmj=1,nparma) ! e-n H
   0.50000000   1.00000000  0. 0. 0. 0. (b(iparmj),iparmj=1,nparmb)
 (c(iparmj),iparmj=1,nparmc) ! e-e-n C
 (c(iparmj),iparmj=1,nparmc) ! e-e-n O
 (c(iparmj),iparmj=1,nparmc) ! e-e-n H
end
```

Prepare a `jastrow.der` file matching this structure (see previous tutorials).

### Optimization Input

Create `vmc_opt_gs.inp`:

```python
%module general
    title           'GS Optimization'
    mode            'vmc_one_mpi1'
%endmodule

load trexio         COH2_GS.trexio
load jastrow        jastrow.start
load jastrow_der    jastrow.der

%module electrons
    nup           6
    nelec         12
%endmodule

%module optwf
    ioptwf        1
    ioptci        0      # Optimize Jastrow first
    ioptjas       1
    ioptorb       0
    method        'sr_n'
    nopt_iter     10
    nblk_max      5000
    ncore         0
    nextorb       600    # Ensure large enough for virtuals
    sr_tau        0.05
    sr_eps        0.01
    sr_adiag      0.01
%endmodule

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      100
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

Run CHAMP. Save the resulting optimized Jastrow (e.g., `cp jastrow_optimal.1.iter10 jastrow_gs.opt`).

### Optimize CI Coefficients

Create a new input or modify the existing one to optimize CI coefficients as well:

```python
    ioptci        1      # Enable CI optimization
    load jastrow  jastrow_gs.opt
```

Run CHAMP again to get the fully optimized GS wavefunction.

## 2. Excited State Optimization

Create a directory for the Excited State calculation:

```bash
cd ..
mkdir ES_Calc
cd ES_Calc
ln -s ../COH2_ES.trexio .
```

### reuse GS Jastrow

We can use the optimized Jastrow from the Ground State as a good starting point for the Excited State.

```bash
cp ../GS_Calc/jastrow_gs.opt jastrow.start
cp ../GS_Calc/jastrow.der .
```

### ES Optimization Input

Create `vmc_opt_es.inp`. Use the same settings as GS, but load the ES TREXIO file.

```python
load trexio         COH2_ES.trexio
load jastrow        jastrow.start
load jastrow_der    jastrow.der
```

Run the optimization (Jastrow first, then Jastrow+CI).

## 3. Calculate Excitation Energy

Run DMC calculations for both optimized GS and ES wavefunctions (see Ground State tutorial for DMC setup).

The excitation energy is:
$$ \Delta E = E_{\text{DMC}}^{\text{ES}} - E_{\text{DMC}}^{\text{GS}} $$

Compute the error bars using standard error propagation:
$$ \delta (\Delta E) = \sqrt{\delta E_{\text{GS}}^2 + \delta E_{\text{ES}}^2} $$
