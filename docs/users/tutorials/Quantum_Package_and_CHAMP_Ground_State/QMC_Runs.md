---
title: QMC runs (a) Check setup
tags:
    - tutorial
    - CIPSI
    - QP
---

# QMC runs: Check setup and Initial VMC

## Check QP Energies

First, verify the energies of the wavefunctions generated in Quantum Package:

```bash
qp set_file h2o_hf
qp run print_energy

qp set_file h2o_dft
qp run print_energy
```

You should observe energies close to:

- **HF**: -16.9503842 Ha
- **DFT**: -16.9465884 Ha

These serve as reference values for your VMC calculations.

## Setup CHAMP Calculation

We will use the TREXIO files exported from Quantum Package directly.

1.  Create a directory for the run (e.g., `H2O_HF`) and enter it.
2.  Link or copy the `h2o_hf.trexio` directory/file here.
3.  Create a `pool` directory if you wish to organize files, though with TREXIO it's less critical.

### Input File

Create a CHAMP input file `vmc_h2o_hf.inp`:

```python
%module general
    title           'H2O HF calculation'
    mode            'vmc_one_mpi1'
%endmodule

# Load all wavefunction data from TREXIO
load trexio        h2o_hf.trexio

# Load Jastrow factor (initially empty/simple)
load jastrow       jastrow.start

%module electrons
    nup           4
    nelec         8
%endmodule

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      20000
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

### Initial Jastrow Factor

Create a simple starting Jastrow file `jastrow.start`. Since we use TREXIO, we need to ensure the Jastrow format matches the atoms in the TREXIO file.

```python
jastrow_parameter   1
  0  0  0           norda,nordb,nordc
   0.60000000   0.00000000     scalek,a21
   0.00000000   0.00000000   (a(iparmj),iparmj=1,nparma)
   0.00000000   0.00000000   (a(iparmj),iparmj=1,nparma)
   0.00000000   1.00000000   (b(iparmj),iparmj=1,nparmb)
 (c(iparmj),iparmj=1,nparmc)
 (c(iparmj),iparmj=1,nparmc)
end
```

This defines a trivial Jastrow factor ($e^J = 1$) effectively running VMC with the bare HF/DFT wavefunction.

### Run Calculation

Submit the job using your system's scheduler or run locally:

```bash
/path/bin/vmc.mov1 -i vmc_h2o_hf.inp -o vmc_h2o_hf.out -e error
```

Compare the resulting total energy with the QP reference. They should match closely (within statistical error), confirming the interface works correctly.

Repeat the process for the DFT wavefunction (`H2O_DFT`).
