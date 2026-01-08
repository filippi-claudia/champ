---
title: CHAMP and Quantum Package
tags:
    - tutorial
    - overview
    - QP
---

# CHAMP and Quantum Package

This set of tutorials demonstrates the interoperability between **Quantum Package (QP)** and **CHAMP** using the **TREXIO** file format.

## Overview

Quantum Package is used as the **wavefunction generator**, providing high-quality trial wavefunctions (Hartree-Fock, DFT, or CIPSI expansions). CHAMP is then used to perform **Quantum Monte Carlo (QMC)** calculations (VMC and DMC) using these wavefunctions.

## Tutorials in this Series

1.  **[Ground State Calculation](Quantum_Package_and_CHAMP_Ground_State/Quantum_Package_and_CHAMP.md)**
    *   Learn the basic workflow: Generating a single-determinant wavefunction (HF/DFT) for a water molecule in QP, exporting it to TREXIO, and running VMC/DMC in CHAMP.
    *   Covers basis set setup, geometry optimization (in QP), and Jastrow optimization (in CHAMP).

2.  **[Excited States Calculation](Quantum_Package_and_CHAMP_Excited_States/Quantum_Package_and_CHAMP_Excited_States.md)**
    *   Advanced workflow for excited states: Using CIPSI in QP to generate multi-determinant wavefunctions for ground and excited states ($COH_2$).
    *   Covers state-specific export to TREXIO and optimizing wavefunctions for specific electronic states in CHAMP.
