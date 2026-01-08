---
title: CIPSI 2-state calculations
tags:
    - tutorial
    - CIPSI
    - excited state
    - QP
---

# CIPSI 2-state calculations

## Import Wavefunction

Import the single-determinant wavefunction from the TREXIO file:

[Download COH2.trexio](https://github.com/TREX-CoE/school-slovakia-2022/raw/master/docs/TrexioFiles/COH2.trexio){: .btn .btn-blue }

```bash
qp_import_trexio.py COH2.trexio -o COH2
qp set_file COH2
```

## Setup Two-State Calculation

Configure QP to seek two states:

```bash
qp set determinants n_states 2
qp set determinants n_det_max 2000
```

## Running the Calculation

First, perform a CIS calculation to ensure we capture states of different symmetries (Ground State and Excited State):

```bash
qp run cis | tee COH2.cis.out
```

Now, refine the wavefunction using CIPSI (FCI in truncated space). We use a selection factor to carefully grow the determinant space:

```bash
qp set determinants selection_factor 0.5
qp set determinants read_wf true
qp run fci | tee COH2.fci.out
```

Check `COH2.fci.out`. You should see two states with an excitation energy around 4 eV.