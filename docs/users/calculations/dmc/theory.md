# DMC Theory

## Imaginary Time Evolution

The imaginary time Schrödinger equation is a diffusion equation with a source/sink term.

$$ -\frac{\partial \Psi}{\partial \tau} = -\frac{1}{2} \nabla^2 \Psi + (V - E_T) \Psi $$

By simulating a population of walkers that diffuse (kinetic energy) and branch (potential energy), DMC samples the ground state distribution.

## Importance Sampling

To improve efficiency, we sample the distribution $f(\mathbf{R}, \tau) = \Psi_T(\mathbf{R}) \Psi(\mathbf{R}, \tau)$. The evolution equation for $f$ includes a drift term guided by the trial wavefunction.

## Fixed-Node Approximation

The fixed-node approximation enforces the nodes (zeros) of the wavefunction to be the same as those of the trial wavefunction $\Psi_T$. This prevents the "fermion sign problem" where the wavefunction would otherwise collapse to the bosonic ground state.

The fixed-node energy is variational: $E_{FN} \ge E_0$. The error depends on the quality of the nodal surface of $\Psi_T$.

## Non-Local Pseudopotentials

When using non-local pseudopotentials, the locality approximation or T-moves (Casula moves) are used to handle the non-local integration. In CHAMP, `icasula` controls this behavior.
