---
title: Tutorials
tags:
    - tutorial
    - learning path
---

# Tutorials

Welcome to the CHAMP tutorials section. Here you will find step-by-step guides ranging from basic workflows to advanced optimization techniques.

## [CHAMP and Quantum Package](champ_and_qp.md)

This series is the core learning path for new users. It demonstrates the seamless integration between **Quantum Package (QP)** (for wavefunction generation) and **CHAMP** (for QMC), utilizing the **TREXIO** file format.

*   **[Ground State Calculation](Quantum_Package_and_CHAMP_Ground_State/Quantum_Package_and_CHAMP.md)**:
    *   Learn to generate HF/DFT wavefunctions in QP.
    *   Export to TREXIO.
    *   Perform VMC and DMC calculations (energy and optimization) on a water molecule.

*   **[Excited States Calculation](Quantum_Package_and_CHAMP_Excited_States/Quantum_Package_and_CHAMP_Excited_States.md)**:
    *   Generate multi-determinant wavefunctions (CIPSI) for ground and excited states ($COH_2$).
    *   Optimize wavefunctions specifically for excited states.

## [Validation Examples](examples/index.md)

This section provides a comprehensive collection of diverse examples derived from the CHAMP **Continuous Integration (CI) test suite**. These are excellent for understanding specific features or specific input configurations.

*   **Basic Verification**: [H2](examples/h2_verification/index.md), [Butadiene DMC](examples/butadiene_dmc_500/index.md), [Water DMC](examples/water_trexio_dmc/index.md).
*   **Advanced Wavefunctions**: [SDT Excitations](examples/butadiene_vmc_sdt/index.md), [Determinant Scaling](examples/butadiene_vmc_scaling/index.md), [3-Body Jastrow](examples/water_3body_jastrow/index.md).
*   **Optimization Techniques**:
    *   State-Specific & Multi-State: [HNO](examples/hno_state_specific/index.md), [PSB](examples/psb_multistate/index.md).
    *   Orbital & CSF Optimization: [Ethanol](examples/ethanol_orb_opt/index.md), [Butadiene CSF](examples/butadiene_csf_opt/index.md).
    *   Geometry: [Force Analysis & Optimization](examples/butadiene_geo_opt/index.md).
*   **Advanced Features**: [Periodic Forces](examples/periodic_forces/index.md), [TREXIO Backends](examples/trexio_backend_comparison/index.md).
