# Wavefunction Optimization

Wavefunction optimization is a critical step in QMC calculations. It involves adjusting the parameters of the trial wavefunction (Jastrow factors, orbital coefficients, CI coefficients) to minimize the energy or variance.

## Overview

CHAMP supports several optimization algorithms, controlled by the `optwf` module.

## Key Concepts

- **Parameters**: The variables being optimized (e.g., Jastrow coefficients, orbital rotation parameters).
- **Cost Function**: The quantity being minimized (usually energy, variance, or a combination).
- **Algorithm**: The numerical method used to update parameters (e.g., Stochastic Reconfiguration, Linear Method).

## Input Keywords

The `optwf` module controls the optimization process.

```python
%module optwf
    ioptwf        1       # Enable optimization
    ioptjas       1       # Optimize Jastrow parameters
    ioptorb       0       # Optimize orbitals (0=no, 1=yes)
    ioptci        0       # Optimize CI coefficients
    
    method        'sr_n'  # Optimization method (sr_n, lin_d, linear)
    nopt_iter     20      # Number of optimization iterations
    nblk_max      4000    # Max blocks per iteration
    
    sr_tau        0.05    # Step size for SR
    sr_eps        0.001   # Regularization for SR
%endmodule
```

## Available Methods

- **[Energy Minimization](energy_minimization.md)**: Minimizing the energy using Stochastic Reconfiguration (SR) or other gradient-based methods.
- **[Linear Method](linear_method.md)**: A robust method that can minimize energy, variance, or a mix.
- **[Variance Minimization](variance_minimization.md)**: Minimizing the variance of the local energy (often used for initial optimization).
