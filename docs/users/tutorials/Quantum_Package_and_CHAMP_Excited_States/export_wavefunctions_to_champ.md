---
layout: default
title: Export wave functions to CHAMP
nav_order: 3
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

# Export wave functions to CHAMP

The excited states are of different symmetries, so we will generate two
different setups in CHAMP, one for each state. To do that, we will save
two different files, one for each state, and containing only the
non-zero determinants.

First, copy the `COH2` directory into `COH2_GS` and `COH2_ES`, one
directory for each state:

```bash
cp -r COH2 COH2_GS
cp -r COH2 COH2_ES
```

Then, we will use `qp_edit` to extract one state in each EZFIO
directory:

```bash
qp set_file COH2_GS
qp edit --state=1

qp set_file COH2_ES
qp edit --state=2
```

The states have been extracted, but the EZFIO databases still contain
the determinants with almost zero coefficients. We can remove them by
running

```bash
qp set_file COH2_GS
qp run truncate_wf
```

This last program is interactive and asks for the minimum weight of the
kept configurations. Answer `1.d-10` to this question.

Similarly, remove the negligible determinants from the excited state:

```bash
qp set_file COH2_ES
qp run truncate_wf
```

We can now export the wave functions in two different TREXIO files. To
do that, for each state we copy the initial TREXIO file and add the
determinants information:

```bash
cp COH2.trexio COH2_GS.trexio
qp set_file COH2_GS
qp set trexio trexio_file  COH2_GS.trexio
qp run export_trexio
```

```bash
cp COH2.trexio COH2_ES.trexio
qp set_file COH2_ES
qp set trexio trexio_file  COH2_ES.trexio
qp run export_trexio
```

Now, we are ready to run the QMC calculations for each state.
