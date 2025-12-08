# Diffusion Monte Carlo (DMC)

Diffusion Monte Carlo (DMC) is a projector method that extracts the ground state energy from a trial wavefunction by evolving it in imaginary time.

## Overview

DMC solves the Schrödinger equation in imaginary time $\tau = it$:

$$ -\frac{\partial \Psi}{\partial \tau} = (\hat{H} - E_T) \Psi $$

As $\tau \to \infty$, the wavefunction projects onto the ground state $\Phi_0$, provided the trial wavefunction $\Psi_T$ has a non-zero overlap with $\Phi_0$.

## Key Features

- **Fixed-Node Approximation**: To maintain the fermionic nature of the wavefunction (antisymmetry), the nodal surface is fixed to that of the trial wavefunction $\Psi_T$. This provides a variational upper bound to the ground state energy.
- **Importance Sampling**: The random walk is guided by the trial wavefunction to improve efficiency.
- **Time Step Error**: The simulation uses a finite time step $\tau$. Results should be extrapolated to $\tau \to 0$.

## Running DMC

DMC calculations are controlled by the `blocking_dmc` and `dmc` modules.

```python
%module blocking_dmc
    dmc_nstep     100
    dmc_nblk      200
%endmodule

%module dmc
    tau           0.01
    etrial        -15.8
%endmodule
```

See [Input Keywords](input_keywords.md) for details.
