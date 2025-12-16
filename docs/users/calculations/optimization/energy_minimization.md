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

