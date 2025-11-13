---
layout: default
title: Hartree-Fock calculation
nav_order: 3
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

# Hartree-Fock calculation

Create the EZFIO directory with the geometry, basis and pseudopotential
parameters:

```bash
qp create_ezfio --pseudo=PSEUDO --basis=BASIS h2o.xyz --output=h2o_hf
```

Run the Hartree-Fock calculation

```bash
qp run scf | tee h2o_hf.out
```

Export the wave function into TREXIO format

```bash
qp set trexio trexio_file h2o_hf.trexio
qp run export_trexio
```
