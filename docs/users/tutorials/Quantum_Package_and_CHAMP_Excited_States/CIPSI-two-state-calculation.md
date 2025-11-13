---
layout: default
title: CIPSI 2-state calculations
nav_order: 2
grand_parent: Tutorials
parent: '02. QP and CHAMP : Excited State Calculation'
authors:
    - Ravindra Shinde
    - Claudia Filippi
    - Anthony Scemama
tags:
    - CHAMP
    - tutorial
    - CIPSI
    - excited state
    - QP
---

# CIPSI 2-state calculations

You can import the single-determinant wave function from the provided

[TREXIO file
COH2.trexio](https://github.com/TREX-CoE/school-slovakia-2022/raw/master/docs/TrexioFiles/COH2.trexio){: .btn .btn-blue }

as:

```bash
qp_import_trexio.py COH2.trexio -o COH2
qp set_file COH2
```

Specify that you want to run a two-state calculation:

```bash
qp set determinants n_states 2
```

Tell QP to stop when the number of determinants is larger than 2000

```bash
qp set determinants n_det_max 2000
```

and run the CIPSI in the Full-CI space:

```bash
qp run fci | tee COH2.fci.out
```

{: .note}
>The extrapolated excitation energy is around 8 eV and we expect 4 ev, so
we did not catch the correct state. This is because the orbitals in the
TREXIO file are symmetry adapted, so it is impossible to make a
determinant from another symmetry enter in the determinant space.
>
>To obtain a solution from another symmetry, we need to put at least one
determinant of each symmetry.


{: .new-title}
> Tip
>
>The simplest practical solution is to first perform a CIS, and then
continue with a CIPSI in the FCI space.


By default, at every iteration QP tries to double the size of the wave
function. In QMC, we will use a small number of determinants, so we
should tell QP to add less determinants at each iteration to have a
finer selection.

```bash
qp set determinants selection_factor 0.5
```

```bash
qp run cis | tee COH2.cis.out
qp set determinants read_wf true
qp run fci | tee COH2.fci.out
```

`read_wf = true` specifies that the wave function stored in the EZFIO
database should be used as a starting point for the the CI calculation.

Now, we have obtained a more reasonable excitation energy, around 4 eV.
We are now ready to export the data for CHAMP.