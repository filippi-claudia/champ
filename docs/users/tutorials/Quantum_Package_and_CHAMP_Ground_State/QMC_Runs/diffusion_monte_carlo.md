---
layout: default
title: QMC runs (c) Diffusion Monte Carlo
nav_order: 7
grand_parent: Tutorials
parent: '01. QP and CHAMP : Ground State Calculation'
mathjax: true
authors:
    - Ravindra Shinde
    - Claudia Filippi
    - Anthony Scemama
tags:
    - CHAMP
    - tutorial
    - CIPSI
    - ground state
    - DMC
    - QP
---

# Diffusion Monte Carlo

Let us start to run a DMC simulation with the HF orbitals and the
optimal Jastrow factor you have just generated.

Create a new directory and copy the wave function TREXIO info and the
optimal Jastrow factor (for simplicity, pick the last one).

First, generate an input file as before where you read the wave function
files (careful to load the new Jastrow factor) and perform a short VMC
calculation to generate the walkers for DMC.

To shorten the VMC run, you can choose a small `vmc_nblk` in the main
input file and modify `vmc_nconf_new` to be the number of walkers per
core you wish. Here, we use the same values as for the starting
iterations of the Jastrow factor optimization:

```python
%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      200
    vmc_nblkeq    1
    vmc_nconf_new 100
%endmodule
```

This will generate 100 walkers per core (`vmc_nconf_new`) by writing the
coordinates of a walker every $$20 \times 200 / 100$$ steps. Since the
correlation time is less than 2 step in VMC, your walkers will be
decorrelated.

A bunch of `mc_configs_newX` files will appear in your directory, each
containing 100 walkers.

```python
cat mc_configs_new* >> mc_configs
rm mc_configs_new*
```

`mc_configs` contains now all walkers.

Generate a DMC input

```python
%module blocking_dmc
    dmc_nstep     60
    dmc_nblk      40
    dmc_nblkeq    1
    dmc_nconf     100
%endmodule

%module dmc
    tau           0.05
    etrial      -17.240
    icasula      -1
%endmodule
```

You also need to change the `mode` keyword in the input file:

```python
mode            'dmc_one_mpi1'
```

within the general module.

Some debug files are being created, that you can just erase.

```bash
rm problem*
rm mc_configs_new*
```

To look at the energy, you can do

```bash
grep '( 100) =' dmc*out
```

In the last column, you have the correlation time.

{: .warning}
Make sure that you have chosen `dmc_nstep` about two times larger.


Also perform another calculation with a smaller time step.

{: .warning}
Make sure that you increase `dmc_nstep` by as much as you have decreased
$$\tau$$.


{: .warning}
You do not need to regenerate the file `mc_configs` containing the
walkers.

Repeat the optimization and DMC calculation with the DFT orbitals.
Compare the VMC and DMC energies.


## Optimal one-determinant Jastrow-Slater wave function


Finally, starting from the DFT orbitals and the optimal two-body Jastrow
optimize the full wave function (Jastrow and orbitals).

To this aim, set `ioptorb` to `1` in the `optwf` module.

```python
ioptorb     1
```
---

## More examples to play with

Multiple files with geometries, basis sets and pseudopotentials can be
downloaded here:
[Examples](https://github.com/TREX-CoE/school-slovakia-2022/tree/master/docs/examples)
