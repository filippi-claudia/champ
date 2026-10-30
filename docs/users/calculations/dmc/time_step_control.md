# Time Step Control

The time step $\tau$ is a crucial parameter in DMC.

## Time Step Error

The DMC algorithm involves a short-time approximation for the Green's function, which introduces a bias proportional to $\tau$ (or $\tau^2$ depending on the algorithm).

To obtain accurate results, one should:
1.  Perform calculations at multiple time steps (e.g., $\tau = 0.01, 0.005, 0.002$).
2.  Extrapolate the energy to $\tau \to 0$.

## Choosing $\tau$

- **Large $\tau$**: Faster equilibration, smaller statistical error for same CPU time, but larger bias.
- **Small $\tau$**: Smaller bias, but requires more steps to reduce statistical error.

A good starting point is often $\tau \approx 0.01$ - $0.05$ a.u. for all-electron calculations of light atoms, or larger for pseudopotential calculations.

## Acceptance Ratio

The acceptance ratio of the Metropolis step in DMC is a good indicator. It should be high (> 99%) for the time step error to be small. If it drops significantly, $\tau$ is likely too large.
