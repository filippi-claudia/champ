# VMC Optimization

(Reference: `tests/CI_test/VMC-TREXIO-H2O-DFT-optall`)

This workflow describes how to optimize all parameters of the wavefunction (Jastrow, Orbitals, CI coefficients) using VMC.

## Steps

1.  **Input Preparation**:
    -   Load the TREXIO file containing the starting orbitals (e.g., from DFT).
    -   Define the Jastrow factor (or load a starting one).

2.  **Input Configuration**:
    -   Set `ioptwf = 1` to enable optimization.
    -   Set `ioptjas = 1`, `ioptorb = 1`, `ioptci = 1` to optimize all sets of parameters.
    -   Choose an optimization method (e.g., `method = 'sr_n'` for Stochastic Reconfiguration or `method = 'linear'` for the Linear Method).

3.  **Execution**:
    -   Run CHAMP.
    -   Monitor the energy and gradient norm in the output.

## Example Input

```python
%module optwf
    ioptwf        1
    ioptjas       1
    ioptorb       1
    ioptci        0  # Set to 1 if using multi-determinant WF
    method        'sr_n'
    nopt_iter     20
    sr_tau        0.05
%endmodule
```

## Tips

-   Start by optimizing Jastrow only (`ioptjas=1`, `ioptorb=0`).
-   Then optimize Jastrow + Orbitals together.
-   Use `vmc_nblk_max` to increase the precision as the optimization proceeds.
