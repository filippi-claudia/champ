# VMC → DMC Workflow

(Reference: `tests/CI_test/DMC-TREXIO-water-DFT-jas2body_tau0.05`)

This is the standard production workflow: optimize the wavefunction in VMC, then run DMC to get the final energy.

## Steps

1.  **VMC Optimization**:
    -   Optimize the wavefunction as described in [VMC Optimization](vmc_optimization.md).
    -   Save the optimized parameters (usually done automatically to `jastrow_optimal` and `orbitals_optimal`).

2.  **VMC Walker Generation**:
    -   Run a short VMC calculation with the *optimized* wavefunction.
    -   Set `vmc_nconf_new` to the desired number of walkers (e.g., 100 per core).
    -   This generates `mc_configs` files.

3.  **DMC Run**:
    -   Create a new input file (or modify the existing one).
    -   Load the optimized Jastrow and orbitals.
    -   Set `mode = 'dmc_one_mpi1'` (or similar DMC mode).
    -   Configure the `dmc` module (time step `tau`, `etrial`).

## Example Input (DMC Step)

```python
load jastrow      jastrow_optimal
load orbitals     orbitals_optimal

%module blocking_dmc
    dmc_nstep     100
    dmc_nblk      200
    dmc_nconf     100
%endmodule

%module dmc
    tau           0.01
    etrial        -76.40  # From VMC
%endmodule
```
