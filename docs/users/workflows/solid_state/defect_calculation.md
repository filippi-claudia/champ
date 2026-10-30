# Defect Calculation

Calculating the formation energy of a point defect (e.g., vacancy, interstitial).

## Workflow

1.  **Perfect Crystal**: Calculate the energy of the perfect supercell ($E_{perf}$).
2.  **Defect Supercell**: Create a supercell with the defect (remove/add atom).
3.  **Relaxation**: Relax the geometry of the defect supercell using DFT.
4.  **QMC**: Calculate the energy of the relaxed defect supercell ($E_{def}$).
5.  **Formation Energy**:
    $$ E_f = E_{def} - E_{perf} \pm \mu $$
    where $\mu$ is the chemical potential of the removed/added species.

## Finite-Size Corrections

Defect calculations in QMC are subject to significant finite-size errors (image interactions). Corrections (e.g., Makov-Payne) are often necessary.
