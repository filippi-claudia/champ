# Calculations

This section describes the different types of calculations that can be performed using CHAMP.

## Available Calculation Types

CHAMP supports the following main calculation types:

- **[Variational Monte Carlo (VMC)](vmc/index.md)**:  
  Perform VMC calculations to evaluate the energy and other properties of a trial wavefunction.

- **[Wavefunction Optimization](optimization/index.md)**:  
  Optimize the parameters of the trial wavefunction (Jastrow factors, orbitals, coefficients) to minimize energy or variance.

- **[Diffusion Monte Carlo (DMC)](dmc/index.md)**:  
  Perform DMC calculations to project out the ground state component of the trial wavefunction, recovering correlation energy.

## General Workflow

A typical workflow involves:

1.  **Preparation**: Generating the initial wavefunction (orbitals, basis sets) using external tools or internal converters.
2.  **VMC**: Evaluating the quality of the initial wavefunction.
3.  **Optimization**: Optimizing the Jastrow factor and potentially other parameters.
4.  **DMC**: Performing a production DMC run for high-accuracy results.

Please refer to the specific subsections for detailed instructions on each calculation type.
