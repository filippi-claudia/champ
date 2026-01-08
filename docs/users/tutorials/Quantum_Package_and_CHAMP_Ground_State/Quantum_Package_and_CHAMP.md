---
title: 'Ground State Calculation'
tags:
    - tutorial
    - CIPSI
    - ground state
    - QP
---

# Ground State Calculation

This tutorial guides you through the process of calculating the ground state of a water molecule using Quantum Package (QP) for wavefunction generation and CHAMP for Quantum Monte Carlo calculations.

We will:

1.  Use Quantum Package to generate single-determinant wavefunctions (Hartree-Fock and DFT-PBE).
2.  Export these wavefunctions to the [TREXIO](https://github.com/trex-coe/trexio) format.
3.  Perform VMC and DMC calculations in CHAMP using the TREXIO file directly.
