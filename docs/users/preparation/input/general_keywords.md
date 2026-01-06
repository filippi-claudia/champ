---
title: General Keywords
tags:
    - input
    - general
    - load
    - keywords
---

# General Keywords and Load Commands

This section describes the common keywords used in the `%module general` block and the various `load` commands used to define the system and wavefunction components.

## The `general` module

The `general` module sets up the fundamental parameters for the calculation.

```perl
%module general
    title        'butadiene'
    pool         'pool/'
    pseudopot    BFD
    basis        BFD-T
    mode         'dmc_one_mpi1'
%endmodule
```


### Basic Setup

<div class="grid cards single-col" markdown>

-   __Title__

    ---

    Defines a descriptive title for the run, which is printed to the output.

    ```perl
    title 'Butadiene Optimization'
    ```

</div>

<div class="grid cards single-col" markdown>

-   __Pool__

    ---

    Sets the default directory prefix for file paths accessed via the `$pool` variable.

    ```perl
    pool './pool/'
    ```

</div>

<div class="grid cards single-col" markdown>

-   __Mode__

    ---

    Specifies the type of calculation to run.

    ```perl
    mode 'vmc_one_mpi'
    ```
    
    [:octicons-arrow-right-24: View Calculations](../../calculations/index.md)

</div>

### System Configuration

<div class="grid cards single-col" markdown>

-   __Pseudopotential__

    ---

    Specifies the type of pseudopotential being used.

    ```perl
    pseudopot 'BFD'
    ```

    [:octicons-arrow-right-24: Effective Core Potentials](../ecp.md)

</div>

<div class="grid cards single-col" markdown>

-   __Basis__

    ---

    Specifies the basis set type.

    ```perl
    basis 'BFD'
    ```

    [:octicons-arrow-right-24: Basis Sets](../basis.md)

</div>

<div class="grid cards single-col" markdown>

-   __nloc__

    ---

    Number of localized electrons (pseudo-electrons). Set to 0 for all-electron or standard ECP.

    ```perl
    nloc 0
    ```

</div>

### State-Specific / Multi-State Options

<div class="grid cards single-col" markdown>

-   __nstates__

    ---

    Number of electronic states involved in the calculation.

    ```perl
    nstates 2
    ```

</div>

<div class="grid cards single-col" markdown>

-   __weights__

    ---

    Array of weights for each state in the optimization target function.

    ```perl
    weights [ 0.5, 0.5 ]
    ```

</div>

<div class="grid cards single-col" markdown>

-   __weights_guiding__

    ---

    Array of weights defining the guiding wavefunction $\Psi_G$.

    ```perl
    weights_guiding [ 1.0, 0.0 ]
    ```

</div>

<div class="grid cards single-col" markdown>

-   __sr_lambda__

    ---

    Array of state-specific shift (?) parameters for importance sampling.

    ```perl
    sr_lambda [ 1.0d0 ]
    ```

</div>

<div class="grid cards single-col" markdown>

-   __anorm__

    ---

    Array of normalization constants for each state.

    ```perl
    anorm [ 1.0, 0.32 ]
    ```

</div>

## File Loading Commands

CHAMP uses `load` commands to import data from external files. These commands are placed outside of `%module` blocks.

### System Geometry & Basis

<div class="grid cards single-col" markdown>

-   __Molecule__

    ---

    Reads atomic coordinates (typically in CHAMP-XYZ format).

    ```perl
    load molecule butadiene.xyz
    ```

    [:octicons-arrow-right-24: Molecular Geometry](../molecule.md)

</div>

<div class="grid cards single-col" markdown>

-   __Basis Info__

    ---

    Reads basis set definitions and pointer information.

    ```perl
    load basis_num_info BFD-T.bfinfo
    ```

    [:octicons-arrow-right-24: Basis Set Pointers](../basis_pointers.md)

</div>

<div class="grid cards single-col" markdown>

-   __Lattice__

    ---

    Reads lattice vectors for periodic systems.

    ```perl
    load lattice box.txt
    ```

</div>

### Wavefunction Components

<div class="grid cards single-col" markdown>

-   __Determinants__

    ---

    Reads the determinant expansion (CI) and CSF definitions.

    ```perl
    load determinants ci.det
    ```

    [:octicons-arrow-right-24: Determinants](../determinants.md)

</div>

<div class="grid cards single-col" markdown>

-   __Orbitals__

    ---

    Reads molecular orbital coefficients.

    ```perl
    load orbitals orbitals.lcao
    ```

    [:octicons-arrow-right-24: Molecular Orbitals](../orbitals.md)

</div>

<div class="grid cards single-col" markdown>

-   __Jastrow__

    ---

    Reads Jastrow factor parameters.

    ```perl
    load jastrow jastrow.Start
    ```

    [:octicons-arrow-right-24: Jastrow Factors](../jastrow.md)

</div>

<div class="grid cards single-col" markdown>

-   __Jastrow Deriv.__

    ---

    Reads Jastrow derivative mappings.

    ```perl
    load jastrow_der jastrow.der
    ```

</div>

<div class="grid cards single-col" markdown>

-   __Symmetry__

    ---

    Reads symmetry information for the system.

    ```perl
    load symmetry molecule.sym
    ```

</div>

### TREXIO and QMCkl

<div class="grid cards single-col" markdown>

-   __TREXIO__

    ---

    Loads geometry, basis, orbitals, and determinants from a single [TREXIO](../trexio.md) file.
    
    ```perl
    load trexio molecule.hdf5
    ```

    [:octicons-arrow-right-24: TREXIO Guide](../trexio.md)

</div>
