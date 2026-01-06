# Energy Minimization

Energy minimization is the most common optimization strategy in modern QMC. It aims to find the parameters that minimize the expectation value of the Hamiltonian.

## Stochastic Reconfiguration (SR)

Stochastic Reconfiguration (SR) is a robust energy minimization method similar to the natural gradient method. It is selected with `method = 'sr_n'`.

### Key Parameters

- `sr_tau`: The time step (step size) for the parameter update. Typical values are 0.01 - 0.1.
- `sr_eps`: A small regularization parameter to stabilize the inversion of the overlap matrix.
- `sr_adiag`: Diagonal shift for the overlap matrix (similar to `sr_eps`).

### Usage Example

```python
%module optwf
    ioptwf      1
    method      'sr_n'
    nopt_iter   20
    sr_tau      0.05
%endmodule
```

## Comparison with Variance Minimization

While variance minimization is robust, there are strong motivations for optimizing the energy directly:

1.  **Primary Goal**: One typically seeks the lowest energy in a VMC or DMC calculation, rather than the lowest variance.
2.  **Parameters**: Variance minimization is highly effective for **Jastrow coefficients**. However, for **determinantal coefficients** (coefficients of determinants, orbital expansions, csf coefficients), it can take many iterations and get stuck in local minima. Most authors use variance minimization primarily for Jastrow parameters.
3.  **Observables**: For a given trial wave function form, energy-minimized wave functions on average yield more accurate values of other expectation values.
4.  **Forces**: The Hellmann-Feynman theorem can be better exploited with energy-minimized wave functions to compute forces on nuclei.

Despite these points, variance minimization remains a cornerstone technique, particularly for the stable optimization of Jastrow factors.
