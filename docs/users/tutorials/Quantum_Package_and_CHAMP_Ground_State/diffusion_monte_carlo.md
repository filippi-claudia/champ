---
title: QMC runs (c) Diffusion Monte Carlo
tags:
    - tutorial
    - CIPSI
    - ground state
    - DMC
    - QP
---

# Diffusion Monte Carlo

After optimizing the Jastrow factor, we can perform a Diffusion Monte Carlo (DMC) simulation. DMC provides a more accurate energy by projecting out the ground state component of the trial wavefunction.

## 1. Generate Walkers

First, perform a short VMC run using the optimized Jastrow factor to generate a population of walkers distributed according to $|\Psi_T|^2$.

Create an input file `vmc_generate_walkers.inp`:

```python
%module general
    title           'Generate Walkers'
    mode            'vmc_one_mpi1'
%endmodule

load trexio         h2o_hf.trexio
load jastrow        jastrow_optimal.1.iter20  # Use your best Jastrow

%module electrons
    nup           4
    nelec         8
%endmodule

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      200
    vmc_nblkeq    1
    vmc_nconf_new 100   # Generate 100 walkers per core
%endmodule
```

Run this. It will produce `mc_configs_newX` files. Concatenate them into a single `mc_configs` file:

```bash
cat mc_configs_new* >> mc_configs
rm mc_configs_new*
```

## 2. Run DMC

Create the DMC input file `dmc_h2o.inp`. Note the change in `mode` and the inclusion of the `dmc` and `blocking_dmc` modules.

```python
%module general
    title           'DMC Run'
    mode            'dmc_one_mpi1'  # Switch to DMC mode
%endmodule

load trexio         h2o_hf.trexio
load jastrow        jastrow_optimal.1.iter20
# Walker file 'mc_configs' is read automatically if present

%module electrons
    nup           4
    nelec         8
%endmodule

%module blocking_dmc
    dmc_nstep     60       # Steps per block
    dmc_nblk      40       # Number of blocks
    dmc_nblkeq    1        # Equilibration blocks
    dmc_nconf     100      # Target population per core
%endmodule

%module dmc
    tau           0.05     # Time step (try 0.05, 0.02, 0.01)
    etrial        -17.240  # Estimate from VMC
    icasula       -1       # Algorithm choice (standard)
%endmodule
```

Run the DMC calculation:

```bash
champ -i dmc_h2o.inp
```

## 3. Analyze Results

Check the output file for the average DMC energy:

```bash
grep 'total E' dmc_h2o.out
```

Or look for the summary at the end of the blocks.

!!! warning "Time Step Bias"
    DMC has a time-step error. Run calculations with different `tau` values (e.g., 0.05, 0.02, 0.01) and extrapolate to $\tau \to 0$ for high precision results.

## Orbital Optimization (Optional)

You can also optimize the orbital coefficients (`ioptorb 1`) along with the Jastrow factor to further improve the node surface and variational energy.
