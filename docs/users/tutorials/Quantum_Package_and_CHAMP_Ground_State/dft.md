---
tags:
    - tutorial
    - CIPSI
    - ground state
    - QP
---

# DFT calculation

Create a new EZFIO directory for the DFT calculation:

```bash
qp create_ezfio --pseudo=PSEUDO --basis=BASIS h2o.xyz --output=h2o_dft
```

Set the calculation to use the PBE exchange-correlation functional:

```bash
qp set dft_keywords exchange_functional pbe
qp set dft_keywords correlation_functional pbe
```

Set a coarser grid to accelerate the calculation (default is very fine):

```bash
qp set becke_numerical_grid grid_type_sgn 1
```

Run the Kohn-Sham DFT calculation:

```bash
qp run ks_scf | tee h2o_dft.out
```

Export the results to a TREXIO file (`h2o_dft.trexio`):

```bash
qp set trexio trexio_file h2o_dft.trexio
qp run export_trexio
```
