---
title: Jastrow Derivatives
tags:
    - jastrow
    - derivatives
    - optimization
    - wavefunction
---

# Jastrow Derivatives

Jastrow derivative files specify which Jastrow parameters should be optimized during wavefunction optimization. This file controls the degrees of freedom in the optimization process by identifying which parameters in the Jastrow factor can vary during energy or variance minimization.

## Overview

When optimizing a trial wavefunction, you typically want to optimize only a subset of the Jastrow parameters while keeping others fixed. The Jastrow derivative file (`.der` file) specifies which parameters are active in the optimization by listing their indices.

### Key Concepts

- **Parameter indices**: Each Jastrow parameter has an index corresponding to its position in the parameter arrays
- **Selective optimization**: Only parameters listed in the derivative file will be optimized
- **Coordinate with Jastrow file**: The derivative file must match the structure of the corresponding Jastrow file

## File Format

Jastrow derivative parameters are provided in a text file that CHAMP reads during optimization.

### Basic Structure

```python
jasderiv
4 4 5 15 15 0 0 nparma,nparmb,nparmc,nparmf
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)
2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb)
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c(iparmj),iparmj=1,nparmc)
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c(iparmj),iparmj=1,nparmc)
end
```

### Format Specification

**Line 1**: Header
```python
jasderiv
```
Identifies this as a Jastrow derivative specification file.

**Line 2**: Parameter counts
```python
4 4 5 15 15 0 0 nparma,nparmb,nparmc,nparmf
```

- `nparma` = 4: Number of `a` parameters to optimize per atom type
- `nparmb` = 4: Number of `b` parameters to optimize
- `nparmc` = 15: Number of `c` parameters to optimize per atom type
- `nparmf` = 0: Number of `f` parameters (typically 0)
- Additional zeros are placeholders for future parameter types

**Lines 3+**: Parameter indices

**a-parameter indices (electron-nucleus)**: One line per unique atom type
```python
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)
```

- Number of `a` lines = number of unique atom types
- Order must match atom type order in Jastrow file
- Values are 1-indexed positions of parameters to optimize
- Example: `3 4 5 6` means optimize the 3rd, 4th, 5th, and 6th `a` parameters
- First two parameters (indices 1 and 2) are typically fixed at 0.0

**b-parameter indices (electron-electron)**: Single line
```python
2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb)
```

- Only one `b` line for entire system
- Values are 1-indexed positions of parameters to optimize
- Example: `2 3 4 5 6` means optimize the 2nd through 6th `b` parameters
- First parameter (index 1) is often fixed

**c-parameter indices (three-body)**: One line per unique atom type (if nordc > 0)
```python
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c(iparmj),iparmj=1,nparmc)
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c(iparmj),iparmj=1,nparmc)
```

- Number of `c` lines = number of unique atom types
- Order must match atom type order in Jastrow file
- Values are 1-indexed positions of parameters to optimize
- If no `c` parameters in Jastrow (`nordc = 0`), omit these lines

**Final line**:
```
end
```

## Examples

### Example 1: Water Molecule - Simple Optimization

**Jastrow file** (`jastrow.jas`):
```python
jastrow_parameter   1
  5  5  0           norda,nordb,nordc
   0.60000000         scalek
   0.00000000   0.00000000  -0.39803661  -0.19727356  -0.04112912   0.08477303 (a for O)
   0.00000000   0.00000000  -0.25378471   0.03693657   0.02169054  -0.00707375 (a for H)
   0.50000000   0.48551602   0.09924779   0.00590014  -0.00626481  -0.00347973 (b)
end
```

**Derivative file** (`jastrow.der`):
```python
jasderiv
4 5 0 0 nparma,nparmb,nparmc,nparmf
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)
2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb)
end
```

**Explanation**:

- Two atom types (O and H), so two `a` lines
- For each atom type: optimize parameters 3-6 (4 parameters)
- For `b` parameters: optimize parameters 2-6 (5 parameters)
- No `c` parameters since `nordc = 0` in Jastrow file
- First two `a` parameters and first `b` parameter remain fixed

### Example 2: Water with 3-body Terms

**Jastrow file** (`jastrow.jas`):
```python
jastrow_parameter   1
  5  5  5           norda,nordb,nordc
0.4 0.0  scalek a21
0. 0. -4.730090305534E-02 -0.396897598319 -0.159437822975 9.708712655798E-02 (a for O)
0. 0. -1.799730417040E-02 0.134738007823 -0.112317815798 8.510888256315E-03 (a for H)
0.500000000000 0.946548937514 0.149132774876 -0.132156435317 9.053314515581E-02 -1.728006822161E-02 (b)
1.69202811660 1.24656170850 -1.88213809691 3.91512631662 -0.743838508334 -8.64443822088 -3.45929875267 7.94536776357 -0.617594991576 -7.28154417029 -3.20829812008 10.1065053694 2.84682886326 3.11253608381 -1.82059538913 -4.12168239496 -2.73374637202 7.97866751721 2.63741371261 0.589927732147 -3.22734144476 -3.43776654324 -0.805248495319 (c for O)
-0.184402482938 -2.24997275301 -1.73421689881 0.472513797867 1.54643138671 6.06904433460 3.29714804065 1.00667381009 -3.97339434983 -1.90821969292 2.48432116561 -1.92995303583 -5.13867842975 -1.55965196034 -3.12079269533 4.42562466029 3.87653442445 -4.16458516079 0.111329001065 -3.79462067595 3.68329167488 0.667608476506 1.96643077097 (c for H)
end
```

**Derivative file** (`jastrow.der`):
```python
jasderiv
4 4 15 0 nparma,nparmb,nparmc,nparmf
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)
2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb)
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c for O)
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c for H)
end
```

**Explanation**:

- Two atom types (O and H)
- Optimize 4 `a` parameters per atom type (indices 3-6)
- Optimize 5 `b` parameters (indices 2-6)
- Optimize 15 `c` parameters per atom type (non-consecutive indices shown)
- Some `c` parameters are not optimized (e.g., indices 1, 2, 4, 6, 10, 12, 19, 22)

### Example 3: Selective Optimization Strategy

For staged optimization, you might start with fewer parameters:

**Initial stage** - optimize only 2-body terms:
```python
jasderiv
2 3 0 0 nparma,nparmb,nparmc,nparmf
  3 4 (iwjasa for O)
  3 4 (iwjasa for H)
2 3 4 (iwjasb)
end
```

**Later stage** - add 1-body terms:
```python
jasderiv
4 5 0 0 nparma,nparmb,nparmc,nparmf
  3 4 5 6 (iwjasa for O)
  3 4 5 6 (iwjasa for H)
2 3 4 5 6 (iwjasb)
end
```

## Determining Parameter Indices

The parameter indices correspond to positions in the Jastrow parameter arrays:

### a-parameters
For `norda = 5`, there are `nparmja = 6` parameters per atom type:
```
Index:  1    2    3         4         5         6
Value: 0.0  0.0  -0.39804  -0.19727  -0.04113   0.08477
```
Typically optimize indices 3-6 (sometimes 4-6 if cusp is pre-set).

### b-parameters
For `nordb = 5`, there are `nparmjb = 6` parameters:
```
Index:  1    2         3         4         5          6
Value: 0.5  0.48552   0.09925   0.00590  -0.00626   -0.00348
```
Typically optimize indices 2-6 (index 1 often fixed at 0.5).

### c-parameters
For `nordc = 5`, there are 23 parameters per atom type. The selection depends on which three-body terms are important for your system.

## Loading Jastrow Derivatives in CHAMP Input

Specify both the Jastrow file and derivative file:

```perl
%module general
    title  'VMC optimization with Jastrow'
    pool   './pool/'
%endmodule

load trexio          $pool/molecule.hdf5
load jastrow         $pool/jastrow.jas
load jastrow_der     $pool/jastrow.der

%module electrons
    nup    5
    nelec  10
%endmodule

%module optwf
    nopt_iter   10
    method      'lin_d'
%endmodule
```

## Related Topics

- [Jastrow Factors](jastrow.md) - Understanding Jastrow parameter structure
- [Wavefunction Optimization](../calculations/optimization/index.md) - Using derivative files in optimization
- [VMC Calculations](../calculations/vmc/index.md) - Testing optimized wavefunctions
- [TREXIO Files](using_trexio_file.md) - Initial wavefunction setup

## Getting Help

- Verify parameter counts match between Jastrow and derivative files
- Check that all parameter indices are valid (not exceeding array sizes)
- Start with fewer parameters and gradually increase complexity
- Monitor optimization progress to ensure parameters are improving the wavefunction
- Consult [Optimization Guide](../calculations/optimization/index.md) for detailed optimization strategies