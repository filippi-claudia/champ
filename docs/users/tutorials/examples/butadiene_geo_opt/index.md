---
title: Butadiene Geometry Optimization
tags:
    - tutorial
    - VMC
    - optimization
    - geometry
---

# Butadiene Geometry Optimization

This tutorial demonstrates a VMC calculation with a large determinant expansion (64,000 determinants) that includes geometry force analysis settings, preparing for geometry optimization.

## System Configuration

*   **Molecule**: Butadiene
*   **Expansion**: ~64,000 Determinants (RAS1022 SDT)
*   **Basis**: BFD-Q

## Input File Setup

**Input File**: `vmc_opt_ras1022_pVQZ_65000.inp`

```perl
%module general
    title        'butadiene vmc optimization 65000 dets'
    mode         vmc_one_mpi
%endmodule

load determinants    ras1022_SDT_pVQZ_determinants_state1.det
# ... load orbitals, jastrow ...

%module optwf
    ioptwf        1
    ioptci        1
    ioptjas       1
    ioptorb       1
    method        'sr_n'
    nextorb       500
%endmodule

%module optgeo
    iforce_analy  1      # Enable force analysis
    iuse_zmat     0      # Use Cartesian coordinates
    alfgeo        0.25d0 # Step size factor for geometry update
%endmodule
```

## Description

The `optgeo` module controls geometry optimization parameters. `iforce_analy 1` computes atomic forces. `alfgeo` sets the scaling for the geometry update step if optimization were fully enabled (though often this is done via external drivers or specific modes).

Resources: [Butadiene Geo Opt Test](https://github.com/filippi-claudia/champ/tree/master/tests/CI_test/VMC-C4H6-ras1022_Q_optWF+geo-64000dets)
