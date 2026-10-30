# DMC Equation of State

Calculating the Equation of State (EOS) involves computing the energy at different volumes.

## Workflow

1.  **Generate Structures**: Create supercells with different lattice constants (volumes).
2.  **VMC & Optimization**: For each volume, optimize the wavefunction (or at least the Jastrow).
3.  **DMC**: Run DMC for each volume.
4.  **Fit**: Fit the Energy vs. Volume data to an EOS (e.g., Murnaghan or Birch-Murnaghan) to extract the equilibrium volume and bulk modulus.

