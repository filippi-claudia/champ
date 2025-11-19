---
title: Molecular Geometry
icon: material/molecule-co2
tags:
    - geometry
    - molecule
    - coordinates
---

# Molecular Geometry

Molecular geometry specifies the atomic positions and nuclear charges for your quantum Monte Carlo calculation. CHAMP supports multiple formats for providing geometry information, either as separate XYZ files or embedded within the input file.

!!! warning "Important: Units are in Bohr (atomic units)"
    
    All atomic coordinates in CHAMP must be specified in **Bohr** (atomic units), not Ångströms. To convert: 1 Ångström = 1.8897259886 Bohr.

## Overview

CHAMP accepts geometry in several formats:

1. **External XYZ file with automatic valence** - Most common, lets CHAMP determine valence electrons
2. **External XYZ file with explicit valence** - Specify valence electrons per atom (useful for ECPs)
3. **Inline geometry block** - Embed coordinates directly in input file
4. **TREXIO file** - Modern format containing geometry and other wavefunction data

## Loading Geometry from External File

The standard way to specify geometry is with the `load molecule` command in your CHAMP input file:

```perl
# Direct path
load molecule  molecule.xyz

# Using pool variable
load molecule  $pool/molecule.xyz

# Absolute path
load molecule  /full/path/to/molecule.xyz
```

**Where to place this**:
```perl
%module general
    title  'My calculation'
    pool   './pool/'
%endmodule

load molecule  $pool/molecule.xyz
# ... other load commands ...

%module electrons
    nup    5
    nelec  10
%endmodule
```

## Geometry File Formats

### Format 1: Automatic Valence

CHAMP automatically determines the number of valence electrons based on the element symbol. This is the simplest format for all-electron calculations.

**File format** (`molecule.xyz`):

```perl
10
# molecular complex (Symbol, X, Y, Z in Bohr)
  Si  -0.59659972  0.06162019  0.21100680
  S   -2.60025162 -2.54807062 -2.52884266
  S    2.14594449  2.17606672 -2.44253887
  S    1.75703132 -2.78062975  2.53564756
  S   -1.40663455  3.06742023  3.14712509
  H   -3.50597461  0.49044059  0.39864337
  H    0.96753971  3.57914102  3.86259992
  H   -0.57825615 -3.70197321 -3.52433897
  H    0.37416575  3.66039924 -3.47898554
  H   -0.21164931 -3.70953211  3.82669513
```

**Format specification**:

- **Line 1**: Number of atoms (integer)
- **Line 2**: Comment line (can be any text, blank or system description)
- **Lines 3+**: `Element  X  Y  Z`
  - `Element`: Chemical symbol (H, He, Li, C, O, Si, etc.)
  - `X, Y, Z`: Cartesian coordinates in Bohr

**Automatic valence electrons**:

- H: 1, He: 2
- C: 4, N: 5, O: 6, F: 7
- Si: 4, P: 5, S: 6, Cl: 7
- (Elements upto Radon Z=86 are currently supported)

**Example: Water molecule**

```perl
3
H2O molecule in Bohr
  O   0.00000000   0.00000000   0.22143139
  H   0.00000000   1.43042809  -0.88572555
  H   0.00000000  -1.43042809  -0.88572555
```

### Format 2: Explicit Valence (For ECPs)

When using effective core potentials (ECPs), you must specify the number of valence electrons per atom. These numbers should also be consistent with the number of valence electrons in the ECP. This format also allows using different labels for the same element.

**File format** (`molecule.xyz`):

```perl
10
# molecular complex (Symbol, X, Y, Z in Bohr, Zvalence)
  Si   -0.59659972  0.06162019  0.21100680    4.0
  S    -2.60025162 -2.54807062 -2.52884266    6.0
  S     2.14594449  2.17606672 -2.44253887    6.0
  S     1.75703132 -2.78062975  2.53564756    6.0
  S    -1.40663455  3.06742023  3.14712509    6.0
  H1   -3.50597461  0.49044059  0.39864337    1.0
  H2    0.96753971  3.57914102  3.86259992    1.0
  H2   -0.57825615 -3.70197321 -3.52433897    1.0
  H2    0.37416575  3.66039924 -3.47898554    1.0
  H2   -0.21164931 -3.70953211  3.82669513    1.0
```

**Format specification**:

- **Line 1**: Number of atoms
- **Line 2**: Comment line
- **Lines 3+**: `Label  X  Y  Z  Zvalence`
  - `Label`: Atom label (can be element symbol with numbers, e.g., H1, H2, O1)
  - `X, Y, Z`: Coordinates in Bohr
  - `Zvalence`: Number of valence electrons (float)

**Use cases**:

- Required when using ECPs (e.g., BFD, ccECP pseudopotentials)
- Allows distinguishing atoms of the same element with different basis sets / Jastrow
- Useful for studying isotope effects or core excitations

**Example: Water with ECP on oxygen**

```perl
3
H2O with ECP on Oxygen (6 valence electrons, 2 core electrons removed)
  O   0.00000000   0.00000000   0.22143139    6.0
  H   0.00000000   1.43042809  -0.88572555    1.0
  H   0.00000000  -1.43042809  -0.88572555    1.0
```

### Format 3: Inline Geometry Block

You can embed the geometry directly in the input file using a block syntax:

```perl
%module general
    title  'Inline geometry example'
%endmodule

%block molecule
3
H2O inline
  O   0.00000000   0.00000000   0.22143139
  H   0.00000000   1.43042809  -0.88572555
  H   0.00000000  -1.43042809  -0.88572555
%endblock
```

Alternatively, read from file using block syntax:

```perl
%block molecule < molecule.xyz
```

**When to use inline geometry**:

- Small molecules (< 10 atoms)
- Self-contained input files for testing
- Quick calculations without managing multiple files

**When to use external files**:

- Large molecules or clusters
- When geometry is shared across multiple calculations
- Better organization and reusability

## Unit Conversion

CHAMP requires coordinates in **Bohr (atomic units)**, but most quantum chemistry codes output in Ångströms.

### Converting Ångströms to Bohr

**Conversion factor**: 1 Ångström = 1.88972612546 Bohr

**Manual conversion**:
```
X_bohr = X_angstrom × 1.88972612546
```

**Python conversion**:
```python
angstrom_to_bohr = 1.88972612546

# Ångström coordinates
coords_ang = [
    [0.0000, 0.0000, 0.1172],  # O
    [0.0000, 0.7572, -0.4689],  # H
    [0.0000, -0.7572, -0.4689]  # H
]

# Convert to Bohr
coords_bohr = [[x * angstrom_to_bohr for x in atom] for atom in coords_ang]
```

### Common Pitfalls

!!! danger "Wrong units will give incorrect results"
    
    Using Ångströms instead of Bohr will make your molecule ~1.89 times larger, drastically affecting energies and wavefunctions. Always verify units!

**Symptoms of wrong units**:

- Energies off by orders of magnitude
- Wavefunction optimization fails
- Unrealistic bond lengths in output

## Validating Geometry

### Check Coordinates

Before running CHAMP, verify your geometry:

1. **Visualize the structure**
   ```bash
   # Using Avogadro
   avogadro molecule.xyz
   
   # Using VMD
   vmd molecule.xyz
   ```

2. **Check units** - Typical bond lengths:
   - C-C single bond: ~2.9 Bohr (1.54 Å)
   - C=C double bond: ~2.5 Bohr (1.34 Å)
   - H-H in H₂: ~1.4 Bohr (0.74 Å)
   - O-H in water: ~1.8 Bohr (0.96 Å)

3. **Verify atom count**
   ```bash
   # First line should match number of atoms
   head -n 1 molecule.xyz
   ```

4. **Check formatting**
   - No missing columns
   - Consistent spacing (spaces or tabs)
   - Valid element symbols

## Examples

### Example 1: Methane (CH₄)

```perl
5
Methane in Bohr
  C   0.00000000   0.00000000   0.00000000
  H   1.18886981   1.18886981   1.18886981
  H  -1.18886981  -1.18886981   1.18886981
  H  -1.18886981   1.18886981  -1.18886981
  H   1.18886981  -1.18886981  -1.18886981
```

### Example 2: Benzene with ECP

```perl
12
Benzene with BFD ECP on C (4 valence electrons each)
  C1   0.00000000   2.63650000   0.00000000   4.0
  C2   2.28350000   1.31825000   0.00000000   4.0
  C3   2.28350000  -1.31825000   0.00000000   4.0
  C4   0.00000000  -2.63650000   0.00000000   4.0
  C5  -2.28350000  -1.31825000   0.00000000   4.0
  C6  -2.28350000   1.31825000   0.00000000   4.0
  H    0.00000000   4.68650000   0.00000000   1.0
  H    4.06000000   2.34325000   0.00000000   1.0
  H    4.06000000  -2.34325000   0.00000000   1.0
  H    0.00000000  -4.68650000   0.00000000   1.0
  H   -4.06000000  -2.34325000   0.00000000   1.0
  H   -4.06000000   2.34325000   0.00000000   1.0
```

### Example 3: Crystal Structure (Periodic)

For periodic systems, geometry defines the unit cell. See [Solid State Workflows](../workflows/solid_state/index.md) for details on supercell construction.


## Related Topics
- [TREXIO Files](using_trexio_file.md) - Modern format containing geometry and other data
- [Effective Core Potentials](ecp.md) - Pseudopotential definitions requiring explicit valence
- [Basis Sets](basis.md) - Must match atom types and ordering in geometry
- [Solid State Workflows](../workflows/solid_state/index.md) - Periodic systems and supercells

## Getting Help

- Review [Troubleshooting Guide](../troubleshooting/index.md)
- Check geometry with visualization tools
- Verify units and formatting carefully
- Open an issue on [GitHub](https://github.com/filippi-claudia/champ)
