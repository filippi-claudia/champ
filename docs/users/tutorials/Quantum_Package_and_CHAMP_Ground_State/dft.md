---
layout: default
title: DFT Calculation
nav_order: 4
grand_parent: Tutorials
parent: '01. QP and CHAMP : Ground State Calculation'
authors:
    - Ravindra Shinde
    - Claudia Filippi
    - Anthony Scemama
tags:
    - CHAMP
    - tutorial
    - CIPSI
    - ground state
    - QP
---

# DFT calculation

Create the EZFIO directory with the geometry, basis and pseudopotential
parameters:

```bash
qp create_ezfio --pseudo=PSEUDO --basis=BASIS h2o.xyz --output=h2o_dft
```

Specify that you want to use the PBE functional.

```bash
qp set dft_keywords exchange_functional pbe
qp set dft_keywords correlation_functional pbe
```

The default DFT grid is very fine. We can specify we want a coarser grid
to accelerate the calculations:

```bash
qp set becke_numerical_grid grid_type_sgn 1
```

Run the Kohn-Sham calculation

```bash
qp run ks_scf | tee h2o_dft.out
```

Export the wave function into TREXIO format

```bash
qp set trexio trexio_file h2o_dft.trexio
qp run export_trexio
```

