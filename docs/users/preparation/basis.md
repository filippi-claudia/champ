---
title: Basis Sets on Radial Grids
tags:
    - basis
    - orbitals
    - wavefunction
---

# Basis Sets on Radial Grids

Basis set files in CHAMP define atomic orbitals on numerical radial grids. These files store the radial parts of the basis functions evaluated on a discrete grid, which CHAMP uses to construct molecular orbitals for the trial wavefunction. The radial grid representation provides high accuracy for the orbital representation while maintaining computational efficiency.

## Overview

CHAMP represents atomic orbitals by evaluating their radial parts on a logarithmic grid and combining them with angular functions (spherical harmonics). This approach:

- Provides accurate representation of orbital shapes, including near-nucleus behavior
- Enables efficient calculation of orbital values and derivatives
- Supports various basis set types (Gaussian, Slater, numerical)
- Works seamlessly with pseudopotentials

### Key Concepts

- **Radial shells**: Groups of basis functions with the same angular momentum quantum number
- **Grid points**: Discrete radial distances where basis functions are evaluated
- **Logarithmic spacing**: Grid points distributed to capture both short and long-range behavior

## File Format

Basis set files have a fixed format and are typically stored in the `pool/` directory. Files generated from the `trex2champ` converter can be used directly.

### Naming Convention

Basis files follow the naming pattern:

```
<basis_name>.basis.<element_symbol>
```

**Examples**:
- `BFD-T.basis.C` - BFD-T basis for carbon
- `ccpVTZ.basis.Si` - cc-pVTZ basis for silicon
- `VDZ.basis.O` - VDZ basis for oxygen

### File Structure

```python
9 3 2000 1.003000 20.000000 0
 0.000000000000e+00  5.469976184517e-01  2.376319920758e+00  5.557936498748e-01  3.412818210005e+00  2.206803021951e-01  8.610719484857e-01  3.738901952004e-01  3.289926074834e+00  1.106692909826e+00
 1.508957441883e-04  5.469976454488e-01  2.376319870895e+00  5.557936481942e-01  3.412817957941e+00  2.206803015581e-01  8.610719410992e-01  3.738901923954e-01  3.289925989316e+00  1.106692890335e+00
 3.018821756935e-04  5.469976724431e-01  2.376319821040e+00  5.557936465134e-01  3.412817705877e+00  2.206803009210e-01  8.610719337127e-01  3.738901895905e-01  3.289925903799e+00  1.106692870844e+00
 ... (additional grid points)
```

### Format Specification

**Line 1**: Header information
```python
9 3 2000 1.003000 20.000000 0
```

- `9` = Number of radial shells (basis functions with different angular momenta)
- `3` = Number of basis function types (typically 3: s, p, d or similar)
- `2000` = Number of radial grid points
- `1.003000` = First grid point radius (in bohr)
- `20.000000` = Maximum radius of grid (in bohr)
- `0` = Flag for grid type (0 = logarithmic)

**Lines 2+**: Radial function values

Each subsequent line contains values for all radial shells at a specific grid point:
```python
r_i  phi_1(r_i)  phi_2(r_i)  phi_3(r_i)  ...  phi_N(r_i)
```

where:
- `r_i` = Radial distance for this grid point (in bohr)
- `phi_j(r_i)` = Value of the j-th radial shell at distance r_i
- Number of values per line = 1 + number of shells

The grid points are typically logarithmically spaced to provide high resolution near the nucleus while covering large distances efficiently.

## Examples

### Example 1: Carbon with BFD-T Basis

**Basis file** (`BFD-T.basis.C`):
```python
9 3 2000 1.003000 20.000000 0
 0.000000000000e+00  5.469976184517e-01  2.376319920758e+00  5.557936498748e-01  3.412818210005e+00  2.206803021951e-01  8.610719484857e-01  3.738901952004e-01  3.289926074834e+00  1.106692909826e+00
 1.508957441883e-04  5.469976454488e-01  2.376319870895e+00  5.557936481942e-01  3.412817957941e+00  2.206803015581e-01  8.610719410992e-01  3.738901923954e-01  3.289925989316e+00  1.106692890335e+00
 3.018821756935e-04  5.469976724431e-01  2.376319821040e+00  5.557936465134e-01  3.412817705877e+00  2.206803009210e-01  8.610719337127e-01  3.738901895905e-01  3.289925903799e+00  1.106692870844e+00
 ... (1997 more lines)
```

**Explanation**:

- 9 radial shells for carbon
- 2000 grid points from ~0 to 20 bohr
- Logarithmic spacing captures both core and valence regions
- Compatible with BFD pseudopotential

### Example 2: Multi-element System (Water)

For a water molecule (H₂O), you need two basis files:

**Directory structure**:
```
pool/
  ├── ccpVTZ.basis.O
  └── ccpVTZ.basis.H
```

**Oxygen basis** (`ccpVTZ.basis.O`):
```python
14 4 2000 1.003000 20.000000 0
 0.000000000000e+00  1.234567890123e+00  ...  (14 shell values)
 1.508957441883e-04  1.234567890123e+00  ...  (14 shell values)
 ... (remaining grid points)
```

**Hydrogen basis** (`ccpVTZ.basis.H`):
```python
5 2 1500 1.003000 15.000000 0
 0.000000000000e+00  9.876543210987e-01  ...  (5 shell values)
 2.011929889229e-04  9.876543210987e-01  ...  (5 shell values)
 ... (remaining grid points)
```

**Explanation**:

- Oxygen has 14 shells (larger basis), hydrogen has 5 shells
- Different grid sizes appropriate for each element
- Both use logarithmic grids
- CHAMP loads both files automatically based on atom types

### Example 3: Pseudopotential Compatibility

When using pseudopotentials, ensure basis sets match:

**BFD pseudopotential system**:
```
pool/
  ├── BFD-T.basis.C      # Valence basis for carbon
  ├── BFD-T.basis.H      # All-electron basis for hydrogen
  └── BFD.gauss_ecp.dat.C  # Carbon pseudopotential
```

The basis must be consistent with the pseudopotential (same valence space).

## Loading Basis Sets in CHAMP Input

Specify the basis set name in the `%module general` section:

```perl
%module general
    title  'VMC calculation for molecule'
    pool   './pool/'
    basis  'BFD-T'        # Basis set prefix
%endmodule

load trexio   $pool/molecule.hdf5
load jastrow  $pool/jastrow.jas

%module electrons
    nup    5
    nelec  10
%endmodule
```

CHAMP will automatically load:
- `BFD-T.basis.C` for carbon atoms
- `BFD-T.basis.H` for hydrogen atoms
- etc., based on atom types in geometry

### Alternative: Explicit Basis Loading

You can also load basis files explicitly:

```perl
load basis_num_info  $pool/BFD-T.basis.C   # For carbon
load basis_num_info  $pool/BFD-T.basis.H   # For hydrogen
```

This is useful when:
- Using non-standard naming conventions
- Debugging basis set issues
- Mixing different basis sets for different atoms

## Generating Basis Files

### From TREXIO Files

The `trex2champ` converter automatically generates basis files:

```bash
trex2champ molecule.hdf5
```

This creates:
- Basis files for each element: `<basis_name>.basis.<element>`
- Properly formatted for CHAMP
- Consistent with molecular orbitals in TREXIO file

### From Quantum Chemistry Calculations

Most quantum chemistry codes export to TREXIO format:

**GAMESS/TurboRVB workflow**:
1. Run quantum chemistry calculation
2. Convert to TREXIO format
3. Use `trex2champ` to generate CHAMP input files

**ORCA/QMC=CHEM workflow**:
1. Run ORCA calculation with appropriate keywords
2. Use QMC=CHEM tools to generate TREXIO
3. Convert to CHAMP format

### Manual Construction

For custom basis sets, construct files following the format:

1. Define grid: Choose number of points and maximum radius
2. Generate logarithmic grid: $r_i = r_0 \times e^{(i-1) \times \text{spacing}}$
3. Evaluate basis functions at each grid point
4. Write header and grid data to file

## Grid Specifications

### Logarithmic Grid

The default grid uses logarithmic spacing:

$$r_i = r_{\text{min}} \times e^{\alpha \cdot i}$$

where:
- $r_{\text{min}}$ ≈ 10⁻³ bohr (first grid point)
- $\alpha$ determines spacing
- Chosen to reach $r_{\text{max}}$ in N points

**Advantages**:
- High resolution near nucleus (important for cusps)
- Good coverage of valence region
- Extends to long-range tail
- Efficient for most systems

### Grid Parameters

Typical values:

| Element Type | Grid Points | Max Radius | Min Radius |
|-------------|-------------|------------|------------|
| Light (H, He) | 1500-2000 | 15-20 bohr | 0.001 bohr |
| First row (C, N, O) | 2000-2500 | 20-25 bohr | 0.001 bohr |
| Second row (Si, P, S) | 2500-3000 | 25-30 bohr | 0.001 bohr |
| Transition metals | 3000-4000 | 30-40 bohr | 0.0005 bohr |

### Grid Quality Check

Verify basis quality:

1. **Completeness**: Grid extends beyond significant orbital amplitude
2. **Resolution**: Sufficient points to resolve oscillations
3. **Tail behavior**: Orbitals decay to ~10⁻⁸ at $r_{\text{max}}$

## Common Basis Sets

### Gaussian Basis Sets

**Pople-style**: 6-31G, 6-311G
- Moderate size
- Good for initial calculations
- May need augmentation for anions/excited states

**Correlation-consistent**: cc-pVXZ (X = D, T, Q, 5)
- Systematically improvable
- cc-pVDZ: Minimal for QMC
- cc-pVTZ: Standard choice
- cc-pVQZ: High accuracy

**Augmented**: aug-cc-pVXZ
- Added diffuse functions
- Essential for anions, weak interactions
- Larger memory/time requirements

### Pseudopotential Basis Sets

**BFD (Burkatzki-Filippi-Dolg)**:
- Optimized for QMC with BFD pseudopotentials
- Small cores for first/second row elements
- Excellent accuracy/efficiency balance
- Recommended for most production calculations

**ccECP**:
- Correlation-consistent effective core potentials
- Designed for coupled cluster and QMC
- Systematic basis set series
- Good for heavy elements

## Best Practices

### Basis Set Selection

1. **Match pseudopotential**: Use basis sets designed for your ECP
2. **Balance accuracy/cost**: Larger basis → better accuracy but slower
3. **System size consideration**: Use smaller basis for large systems
4. **Property-dependent**: Some properties need specific basis features

### File Management

- Keep all basis files in `pool/` directory
- Use consistent naming conventions
- Document basis set sources and versions
- Verify basis files match geometry atom types

### Common Issues

**Missing basis file**:
```
ERROR: Cannot find basis file: BFD-T.basis.Si
```
- Solution: Ensure file exists in `pool/` directory with exact naming

**Grid mismatch**:
```
WARNING: Orbital tail not converged at r_max
```
- Solution: Increase maximum radius or number of grid points

**Angular momentum issues**:
```
ERROR: Insufficient angular momentum in basis
```
- Solution: Use higher angular momentum basis or check molecular orbital requirements

## Related Topics

- [TREXIO Files](using_trexio_file.md) - Source of basis set data
- [Pseudopotentials](pseudopotential.md) - Consistent pseudopotential/basis pairs
- [Molecular Orbitals](orbitals.md) - How basis sets build molecular orbitals
- [Wavefunction Optimization](../calculations/optimization/index.md) - Impact of basis quality

## Getting Help

- Verify basis file format matches specification exactly
- Check that grid extends sufficiently for all orbitals
- Ensure one basis file exists for each unique atom type
- Use `trex2champ` for reliable basis file generation
- Test with standard basis sets before using custom ones
- Consult [Preparation Guide](index.md) for complete workflow

!!! warning "Required Files"
    Basis set files must be present in the `pool/` directory for all unique atom types in your system. CHAMP cannot proceed without the appropriate basis files.
