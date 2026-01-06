# VMC Theory

## The Metropolis Algorithm

VMC uses the Metropolis-Hastings algorithm to sample the probability distribution $P(\mathbf{R}) \propto |\Psi_T(\mathbf{R})|^2$.

1.  **Initialization**: Start with an initial configuration of electrons $\mathbf{R}$.
2.  **Proposal**: Propose a move to a new configuration $\mathbf{R}'$ (e.g., by moving one electron).
3.  **Acceptance/Rejection**: Accept the move with probability:
    $$ A(\mathbf{R} \to \mathbf{R}') = \min \left( 1, \frac{|\Psi_T(\mathbf{R}')|^2 T(\mathbf{R}' \to \mathbf{R})}{|\Psi_T(\mathbf{R})|^2 T(\mathbf{R} \to \mathbf{R}')} \right) $$
    where $T(\mathbf{R} \to \mathbf{R}')$ is the transition probability (often symmetric).
4.  **Accumulation**: Accumulate averages of local observables (e.g., local energy) over the sampled configurations.

## Local Energy

The local energy $E_L(\mathbf{R})$ is a central quantity in VMC.

$$ E_L(\mathbf{R}) = -\frac{1}{2} \sum_i \frac{\nabla_i^2 \Psi_T}{\Psi_T} + V(\mathbf{R}) $$

- The kinetic energy term involves the Laplacian of the wavefunction.
- The potential energy term $V(\mathbf{R})$ includes electron-electron, electron-nucleus, and nucleus-nucleus interactions.

## Statistical Error

The error in the VMC energy estimate decreases as $1/\sqrt{N}$, where $N$ is the number of independent samples. Because samples in a Markov chain are correlated, the effective number of samples is $N_{eff} = N / \tau_{corr}$, where $\tau_{corr}$ is the autocorrelation time.

Blocking analysis is used to estimate the true statistical error by grouping steps into blocks to remove correlation.
