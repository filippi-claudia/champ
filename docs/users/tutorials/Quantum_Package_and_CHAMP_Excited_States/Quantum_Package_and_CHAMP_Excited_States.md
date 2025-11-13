---
title: '02. QP and CHAMP : Excited State Calculation'
mathjax: true
tags:
    - tutorial
    - CIPSI
    - excited state
    - QP
---

# Quantum Package and CHAMP : Excited States Calculation

We will import a Hartree-Fock wavefunction for the formaldehyde ($$COH_2$$)
molecule from a TREXIO file into Quantum Package (QP), and run a
two-state CIPSI calculation with these orbitals. The wavefunctions for
the 2 states will be stored in the TREXIO file, and we will run wave
function optimization in CHAMP, followed by a Diffusion Monte Carlo
calculation.

!!! note
    The theoretical best estimate of the excitation energy (complete basis
    set extrapolation from coupled cluster calculations) is 3.97 eV.

