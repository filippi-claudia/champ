---
title: HNO State-Specific Optimization
tags:
    - tutorial
    - VMC
    - optimization
    - CSFs
    - HNO
---

# Nitroxyl (HNO) State-Specific Optimization

This advanced tutorial demonstrates the **State-Specific** optimization capabilities of CHAMP. It optimizes a Jastrow factor, Orbital coefficients, and Configuration State Function (CSF) coefficients simultaneously for the Nitroxyl ($HNO$) molecule.

## System Configuration

*   **Molecule**: Nitroxyl ($HNO$)
*   **Initial Orbitals**: CAS(12,9) from GAMESS
*   **Expansion**: 2-state CIPSI (322 Determinants mapped to 143 CSFs)
*   **Method**: VMC State-Specific `sr_n` optimization

## The Challenge

Standard optimization often targets the ground state. State-specific optimization allows for a unique Jastrow and set of orbitals to be optimized for a specific electronic state, which is crucial for balanced excited state descriptions.

## Input File Setup

**Input File**: `vmc.inp`

```perl
%module general
    title           'hno-cas129-cipsi322'
    pool            './pool/'
    pseudopot       BFD
    basis           BFD
    
    # State-Specific Optimization Parameters
    mode            'vmc_one_mpi'
    nstates         2                       # Number of states in calculation
    weights         [ 1.0d0,  1.0d0 ]       # Weights for state averaging
    weights_guiding [ 1.0d0,  1.0d0 ]       # Weights for guiding function
    sr_lambda       [ 1.0d0 ]               # State-specific shift (?)
    anorm           [ 1.0d0, 0.32154638 ]   # Norms for states
%endmodule

%module mstates
    iguiding        2       # State index used for guiding (or mode?)
    iefficiency     1       # Efficiency flag
%endmodule

load molecule        $pool/champ_v2_cas129_geom.xyz
load basis_num_info  $pool/champ_v3_cas129_basis_pointers.bfinfo

# Load CAS(12,9) derived data
load determinants    start.det
load orbitals        start.orb
load jastrow         start.jas
load jastrow_der     jastrow.der
load symmetry        champ_v2_cas129_symmetry.sym

%module electrons
    nup           6
    nelec         12
%endmodule

%module optwf
    ioptwf        1 
    ioptci        1       # Optimize CSF coefficients
    ioptjas       1       # Optimize Jastrow
    ioptorb       1       # Optimize Orbitals
    method        'sr_n'
    ncore         0
    nextorb       600
    nblk_max      250
    no_active     0
    nopt_iter     5
    sr_tau        0.025
    sr_eps        0.0005
    sr_adiag      0.01
%endmodule

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      250
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

## Important Parameters

For state-specific optimization, several specific parameters are crucial:

*   **`nstates 2`**: Explicitly tells CHAMP to handle two electronic states.
*   **`weights [1.0, 1.0]`**: Defines the weighting of states in the optimization target function (usually the state-averaged energy or variance).
*   **`weights_guiding [1.0, 1.0]`**: Defines the weights used to construct the guiding wavefunction $\Psi_G$ for importance sampling. Typically, $\Psi_G^2 = \sum_i w_i \Psi_i^2$.
*   **`anorm`**: Normalization factors for the states, ensuring balanced sampling.
*   **`sr_lambda`**: Parameters related to the Stochastic Reconfiguration shift or stabilization in the multi-state context.

## Workflow Context

The input wavefunctions were generated via a multi-step process:

1.  **GAMESS**: CAS(12,9) calculation.
2.  **Quantum Package**: 2-state CIPSI using MCSCF orbitals.
3.  **CHAMP**: Extensive pre-optimization using state-specific `sr_n` method.

Resources and description for this test can be found in the [HNO CI Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-HNO-cipsi_322_dets_143_csfs_optall) folder.
