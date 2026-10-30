---
tags:
    - tutorial
    - CIPSI
    - ground state
    - QP
---

# Hartree-Fock calculation

Create the EZFIO directory using the `h2o.xyz` geometry, `BASIS`, and `PSEUDO` files created in previous steps:

```bash
qp create_ezfio --pseudo=PSEUDO --basis=BASIS h2o.xyz --output=h2o_hf
```

Run the Hartree-Fock calculation:

```bash
qp run scf | tee h2o_hf.out
```

Export the results to a TREXIO file (`h2o_hf.trexio`). This single file will contain all necessary information (geometry, basis set, orbitals, determinants, pseudopotentials) for the QMC calculation.

```bash
qp set trexio trexio_file h2o_hf.trexio
qp run export_trexio
```
