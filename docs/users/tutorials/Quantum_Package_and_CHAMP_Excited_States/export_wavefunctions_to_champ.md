---
title: Export wave functions to CHAMP
tags:
    - tutorial
    - CIPSI
    - excited state
    - QP
---

# Export Wavefunctions to CHAMP

Since the ground and excited states have different symmetries, we handle them separately in CHAMP. We need to create two separate TREXIO files, one for each state.

## Prepare State Directories

Copy the QP data directory `COH2` for each state:

```bash
cp -r COH2 COH2_GS
cp -r COH2 COH2_ES
```

## Extract States

Use `qp_edit` to isolate state 1 (Ground State) and state 2 (Excited State) in their respective directories:

```bash
qp set_file COH2_GS
qp edit --state=1

qp set_file COH2_ES
qp edit --state=2
```

## Truncate Negligible Determinants

Remove determinants with very small coefficients to keep the QMC calculation efficient:

```bash
qp set_file COH2_GS
qp run truncate_wf  # Answer 1.d-10

qp set_file COH2_ES
qp run truncate_wf  # Answer 1.d-10
```

## Export to TREXIO

Now create specific TREXIO files for each state. We start by copying the original file (to keep basis/geometry) and then update the determinants.

**Ground State:**
```bash
cp COH2.trexio COH2_GS.trexio
qp set_file COH2_GS
qp set trexio trexio_file COH2_GS.trexio
qp run export_trexio
```

**Excited State:**
```bash
cp COH2.trexio COH2_ES.trexio
qp set_file COH2_ES
qp set trexio trexio_file COH2_ES.trexio
qp run export_trexio
```

We now have `COH2_GS.trexio` and `COH2_ES.trexio` ready for CHAMP.
