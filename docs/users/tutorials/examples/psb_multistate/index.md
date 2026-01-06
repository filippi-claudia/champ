---
title: PSB Multi-State Optimization
tags:
    - tutorial
    - VMC
    - optimization
    - multi-state
    - CSFs
---

# PSB Multi-State Optimization

This tutorial demonstrates **Multi-State Optimization** using the `mix_n` method for PSB (Penta-2,4-dieniminium) molecules.

## System Configuration

*   **Molecule**: PSB2 / PSB3
*   **Method**: `mix_n` (mixtures of normal distributions? or mixed state optimization)
*   **Target**: Optimizing multiple states simultaneously (Ground + Excited).

## Input File Setup

**Input File**: `revised_psb2_mix_n.inp`

```perl
%module general
    title           'psb-20'
    weights         [ 1.0d0,  1.0d0 ]  # Weights for state averaging
    weights_guiding [ 1.0d0,  1.0d0 ]
    mode            'vmc_one_mpi'
%endmodule

%module mstates
    iguiding        2      # Number of guiding states / states to optimize
    iefficiency     1
%endmodule

load molecules       ...
load determinants    champ_v2_gamess_PSB2_casci44.det
# ...

%module optwf
    ioptwf        1
    ioptci        1
    ioptjas       1
    ioptorb       1
    method        'mix_n'  # Multi-state optimization method
    ncore         0
    nextorb       86
%endmodule
```

## Significance

Multi-state optimization is critical for correctly describing regions of potential energy surfaces where states utilize, such as conical intersections. The `weights` keyword controls the balance between states in the penalty function.

Resources: [PSB2 Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-PSB2-36_dets-20_csfs_ci44_BFD-Da)
