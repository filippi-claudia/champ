# Workflows

This section provides step-by-step guides for common calculation workflows in CHAMP.

## Categories

- **[Input Preparation](input_preparation/index.md)**:  
  Guides on preparing input files, including geometry conversion, basis set generation, and Jastrow tuning.

- **[Molecular Systems](molecular/index.md)**:  
  Workflows specific to isolated molecules, such as VMC optimization, DMC production runs, and excitation energy calculations.

- **[Solid-State Systems](solid_state/index.md)**:  
  Workflows for periodic systems, including supercell construction and equation of state calculations.

- **[Benchmark & Testing](benchmarks/index.md)**:  
  Standard benchmark calculations to verify installation and performance (e.g., H2 dimer, LiH).

## General Advice

Before running large-scale production calculations, it is highly recommended to:
1.  Run a small test case (like the [H2 dimer benchmark](benchmarks/h2_dimer.md)) to ensure the code is working correctly.
2.  Start with a coarse VMC optimization to get reasonable parameters before refining them.
3.  Always check the [Output Analysis](../calculations/vmc/output_analysis.md) to verify convergence.
