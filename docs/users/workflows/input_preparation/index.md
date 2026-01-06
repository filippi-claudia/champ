# Input Preparation

Proper input preparation is essential for a successful QMC calculation. This section covers the steps to generate the necessary input files for CHAMP.

## Key Steps

1.  **[Geometry Conversion](geometry_conversion.md)**:  
    Convert molecular geometries from standard formats (XYZ, PDB) or quantum chemistry codes (TREXIO) into CHAMP format.

2.  **[Basis & MO Generation](basis_mo_generation.md)**:  
    Generate basis set and molecular orbital files. This often involves running a Hartree-Fock or DFT calculation using an external code (e.g., Quantum Package, PySCF, GAMESS) and converting the output.

3.  **[Jastrow Tuning](jastrow_tuning.md)**:  
    Create and tune the initial Jastrow factor to ensure numerical stability before starting the main optimization.

## Tools

CHAMP provides several scripts and tools in the `tools/` directory to assist with these tasks.
