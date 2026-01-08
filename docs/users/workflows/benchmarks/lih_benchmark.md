# LiH Benchmark

Lithium Hydride (LiH) is a standard system for testing all-electron QMC calculations.

## System

-   Molecule: LiH
-   Bond length: 3.015 a.u.
-   All-electron (no ECPs).

## Workflow

1.  **Cusp Correction**: Since this is an all-electron calculation, ensure that the electron-nucleus cusp conditions are satisfied.
2.  **VMC & DMC**: Run both VMC and DMC.
3.  **Accuracy**: Compare the DMC energy to the exact non-relativistic energy (-8.0705 a.u.).
