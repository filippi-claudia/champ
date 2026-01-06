# Excitation Energy

(Reference: `tests/CI_test/VMC-TREXIO-CH2O-excited-8724-dets-BFD-aug-cc-pVDZ`)

Calculating excitation energies involves comparing the energy of the excited state to the ground state.

## Strategies

1.  **Separate Energy Calculations**:
    -   Run VMC/DMC for the ground state.
    -   Run VMC/DMC for the excited state.
    -   Calculate the energy difference.

2.  **State-Specific Optimization**:
    -   Optimize the ground state wavefunction.
    -   Separately optimize the excited state wavefunction.
    -   Calculate the energy difference.
    -   *Challenge*: Ensuring the excited state doesn't collapse to the ground state (variational collapse).

3.  **State-Average Optimization**:
    -   Optimize parameters to minimize the average energy of multiple states.
    -   Useful when states are degenerate or close in energy.


## Example (state-specific optimization)

```python
%module general
    ...
    nstates       2                         # Number of states in the calculation
    weights         [ 1.0d0,  1.0d0 ]       # State weights
    weights_guiding [ 1.0d0,  1.0d0 ]       # Guiding wavefunction weights
    sr_lambda       [ 1.0d0 ]               # SR lambda
    anorm           [ 1.0d0, 0.3215d0 ]     # State amplitudes
    ...
%endmodule

%module mstates
    iguiding        2
    iefficiency     1
%endmodule

```
