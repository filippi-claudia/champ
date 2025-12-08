# Variational Monte Carlo (VMC)

Variational Monte Carlo (VMC) is a quantum Monte Carlo method that evaluates the expectation value of the Hamiltonian (energy) and other properties for a given trial wavefunction.

## Overview

In VMC, the energy is calculated as:

$$ E_V = \frac{\int \Psi_T^*(\mathbf{R}) \hat{H} \Psi_T(\mathbf{R}) d\mathbf{R}}{\int |\Psi_T(\mathbf{R})|^2 d\mathbf{R}} = \int P(\mathbf{R}) E_L(\mathbf{R}) d\mathbf{R} $$

where:
- $\Psi_T(\mathbf{R})$ is the trial wavefunction.
- $P(\mathbf{R}) = \frac{|\Psi_T(\mathbf{R})|^2}{\int |\Psi_T(\mathbf{R})|^2 d\mathbf{R}}$ is the probability distribution sampled by the Metropolis algorithm.
- $E_L(\mathbf{R}) = \frac{\hat{H} \Psi_T(\mathbf{R})}{\Psi_T(\mathbf{R})}$ is the local energy.

## Key Features

- **Zero-Variance Principle**: If $\Psi_T$ is the exact eigenstate, $E_L(\mathbf{R})$ is constant (the eigenenergy), and the statistical error vanishes.
- **Optimization**: VMC is used to optimize the parameters of the trial wavefunction (Jastrow, orbitals) by minimizing the energy or variance.
- **Walker Generation**: VMC generates configurations (walkers) distributed according to $|\Psi_T|^2$, which are used as starting points for DMC calculations.

## Running VMC

VMC calculations are controlled by the `blocking_vmc` module in the input file.

```python
%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      200
    vmc_nblkeq    1
    vmc_nconf_new 100
%endmodule
```

See [Input Keywords](input_keywords.md) for details on configuration.
