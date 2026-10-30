---
title: Effective Core Potentials
tags:
    - ECP
    - pseudopotential
    - core electrons
---

# Effective Core Potentials (ECPs)

Effective Core Potentials (ECPs), also known as pseudopotentials, replace core electrons with a smooth potential, significantly reducing computational cost while maintaining chemical accuracy. This is especially beneficial for heavy elements where relativistic effects and large numbers of core electrons make all-electron calculations prohibitively expensive.

## Overview

ECPs approximate the effect of core electrons on valence electrons through an analytical potential, allowing you to:

- **Reduce computational cost** - Treat only valence electrons explicitly
- **Include relativistic effects** - Incorporate scalar relativistic corrections implicitly
- **Improve efficiency** - Fewer electrons means faster QMC calculations
- **Maintain accuracy** - Well-designed ECPs preserve chemical properties

**Common ECP families**:

- **BFD** (Burkatzki-Filippi-Dolg) - Optimized for QMC, available in CHAMP
- **ccECP** (correlation-consistent ECP) - Energy-consistent with cc-pVnZ basis sets

## File Format

ECP files in CHAMP follow a specific format. CHAMP includes BFD ECPs in the `champ/pool/BFD/ECP_champ/` directory.

### ECP File Structure

Each element requires a separate ECP file with the naming convention: `{FAMILY}.gauss_ecp.dat.{Element}`

**Example**: `BFD.gauss_ecp.dat.Si` (BFD ECP for Silicon)

```perl
BFD Si pseudo
3
3
4.00000000 1 1.80721061
7.22884246 3 9.99633089
-13.06725590 2 2.50043232
1
21.20531613 2 2.26686403
1
15.43693603 2 2.11659661
```

### Format Specification

The file format follows the order: **local component first**, then **non-local projectors** for l = 0, 1, 2, ... up to lmax-1.

**Line 1**: Arbitrary label/comment (free format, written to log file)
```
BFD Si pseudo
```

**Line 2**: Number of projectors + 1 (total number of components = lpot)
```
3
```
This means there are 3 total components: local potential (stored at index lpot=3), plus non-local projectors for l=0 (s-channel) and l=1 (p-channel). The local potential corresponds to l=2 (d-channel).

**Reading order in file:**
1. First component → Local potential (stored at index lpot = 3, highest l)
2. Second component → l = 0 (s-channel, non-local, stored at index 1)
3. Third component → l = 1 (p-channel, non-local, stored at index 2)

**For each component:**

First line of component: Number of Gaussian terms
```
3
```

Subsequent lines: One term per line with three values
```perl
coefficient  power  exponent
4.00000000   1      1.80721061
7.22884246   3      9.99633089
-13.06725590 2      2.50043232
```

**Understanding the power notation:**

!!! important "Power Convention"
    The power value `n` in the file represents $r^{n-2}$ in the actual formula. So:
    
    - Power 1 in file → $r^{-1}$ in formula
    - Power 2 in file → $r^{0} = 1$ in formula  
    - Power 3 in file → $r^{1} = r$ in formula

Each term represents: $C \times r^{n-2} \times e^{-\alpha r^2}$

**Example component blocks:**

```perl
# Local potential (3 terms)
3
4.00000000   1      1.80721061      # Term: 4.0 × r^(-1) × exp(-1.807 r²)
7.22884246   3      9.99633089      # Term: 7.229 × r × exp(-9.996 r²)
-13.06725590 2      2.50043232      # Term: -13.067 × exp(-2.500 r²)

# s-channel non-local projector (1 term)
1
21.20531613  2      2.26686403      # Term: 21.205 × exp(-2.267 r²)

# p-channel non-local projector (1 term)
1
15.43693603  2      2.11659661      # Term: 15.437 × exp(-2.117 r²)
```

### Mathematical Form

The ECP potential is expressed as:

$$V_{\text{ECP}}(r) = V_{\text{local}}(r) + \sum_{l=0}^{l_{\max}-1} \sum_{m=-l}^{l} |Y_{lm}\rangle \Delta V_l(r) \langle Y_{lm}|$$

Where:

- $V_{\text{local}}(r)$ is the local potential (highest l channel, here l=2 for d)
- $\Delta V_l(r) = V_l(r) - V_{\text{local}}(r)$ are the non-local projectors (l = 0, 1, ...)
- $Y_{lm}$ are spherical harmonics

Each potential component is a sum of Gaussians with the power convention:

$$V(r) = \sum_i C_i \, r^{n_i-2} \, e^{-\alpha_i r^2}$$

where $n_i$ is the power value read from the file.

## Using ECPs in CHAMP

### Specifying ECPs in Input File

ECPs are specified in the `%module general` section with the `pseudopot` keyword:

```perl
%module general
    title           'VMC with BFD ECP'
    pool            './pool/'
    mode            'vmc_one_mpi'
    pseudopot       BFD
    basis           cc-pVTZ
%endmodule

load molecule        $pool/molecule.xyz
load basis_num_info  $pool/BFD_basis_pointers.bfinfo
# ... other load commands ...
```

!!! info "Key Points"
    `pseudopot BFD` tells CHAMP to look for `BFD.gauss_ecp.dat.{Element}` files

!!! Warning "Placement"
    Files must exist in the `pool/` directory for each element type
    
!!! danger "Naming convention"    
    Naming must match exactly: `{FAMILY}.gauss_ecp.dat.{Element}`

### File Organization

Recommended directory structure:

```
project/
├── vmc.inp
├── pool/
│   ├── molecule.xyz
│   ├── BFD.gauss_ecp.dat.Si
│   ├── BFD.gauss_ecp.dat.O
│   ├── BFD.gauss_ecp.dat.H
│   └── BFD_basis_pointers.bfinfo
```

### Geometry Requirements with ECPs

When using ECPs, you **must specify valence electron counts** in the geometry file:

**Correct format** (`molecule.xyz`):
```perl
3
SiO2 unit with BFD ECP
  Si   0.00000000   0.00000000   0.00000000    4.0
  O    2.83459300   0.00000000   0.00000000    6.0
  O   -2.83459300   0.00000000   0.00000000    6.0
```

The fourth column (4.0, 6.0) specifies the number of valence electrons:
- Si: 4 valence (BFD treats [Ne] core as frozen)
- O: 6 valence (BFD treats [He] core as frozen)

See [Molecular Geometry](molecule.md) for details on explicit valence specification.

## Available ECP Libraries in CHAMP

### BFD (Burkatzki-Filippi-Dolg)

**Location**: `champ/pool/BFD/ECP_champ/`

**Available elements**:

- First row: H, He
- Second row: Li-Ne
- Third row: Na-Ar
- Fourth row: K-Kr (includes 3d transition metals)
- Fifth row: Rb-Xe (includes 4d transition metals)
- Sixth row: Cs-Rn (excludes lanthanides, 5d transition metals)

### From TREXIO Files

TREXIO files can contain ECP information, but be cautious:

!!! danger "GAMESS ECP Precision Loss"
    
    GAMESS output files truncate ECP coefficients to a limited precision. TREXIO files generated from GAMESS will have insufficient digits, causing errors in QMC calculations. **Always supply ECP files separately when using TREXIO from GAMESS.**

**Extracting ECPs from TREXIO**:
```bash
python3 trex2champ.py \
            --trex benzene.hdf5 \       # Specify the trexio file
			--geom \                    # This flag is required
			--ecp                       # Get the ecp files
```

## ECP and Basis Set Compatibility

ECPs must be used with compatible basis sets:

### Basis Set Requirements

- **BFD ECPs** → BFD basis sets or BFD-augmented correlation-consistent (BFD-aug-cc-pVnZ)
- **ccECP** → cc-pVnZ or aug-cc-pVnZ basis sets

**Mismatched combinations lead to**:
- Unbalanced core-valence treatment
- Poor energy convergence
- Incorrect chemical properties

## Example

### Calculating Total Electrons

With ECPs, total electrons = sum of valence electrons:

**Example: H₂O with BFD**

- H: 1 valence × 2 atoms = 2 electrons
- O: 6 valence × 1 atom = 6 electrons
- **Total: 8 electrons**

**All-electron H₂O would have**:

- H: 1 electron × 2 = 2
- O: 8 electrons × 1 = 8
- **Total: 10 electrons**

### Specifying in CHAMP Input

**An ECP calculation**

```perl
%module general
    title       'H2O with BFD'
    pool        './pool/'
    pseudopot   BFD
    basis       ccVTZ
%endmodule

%module electrons
    nup           4     # Alpha electrons
    nelec         8     # Total electrons (nup + ndown)
%endmodule

# Other modules
```

**An all-electron calculation**

```perl
%module general
    title       'H2O with BFD'
    pool        './pool/'
    basis       ccVTZ
    nloc          0     # Note this flag
%endmodule

%module electrons
    nup           5     # Alpha electrons
    nelec        10     # Total electrons (nup + ndown)
%endmodule

# Other modules
```


For closed-shell: `nup = nelec/2`

For open-shell: Specify alpha spin count explicitly

## Related Topics

- [Molecular Geometry](molecule.md) - Specifying valence electrons with ECPs
- [Basis Sets](basis.md) - Compatible basis sets for ECPs
- [TREXIO Files](using_trexio_file.md) - ECP handling in TREXIO format
- [Solid State Systems](../workflows/solid_state/index.md) - ECPs in periodic calculations

## Additional Resources

- [BFD ECP Database](http://www.burkatzki.com/pseudos/index.html)
- [ccECP Database](https://pseudopotentiallibrary.org/)
- [Pseudopotential Library](https://www.pseudopotentiallibrary.org/)
- [EMSL Basis Set Exchange](https://www.basissetexchange.org/) - Includes ECPs

## Getting Help

- Verify ECP file format matches specification exactly
- Check precision of coefficients (should have 8+ significant digits)
- Ensure valence counts in geometry match ECP definitions
- Consult [CHAMP GitHub](https://github.com/filippi-claudia/champ) for issues

