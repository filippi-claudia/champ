---
title: Jastrow Factors
tags:
    - jastrow
    - correlation
    - wavefunction
---

# Jastrow Factors

Jastrow factors are correlation functions that multiply the Slater determinant wavefunction to explicitly account for electron-electron and electron-nucleus correlations. They dramatically improve the accuracy of quantum Monte Carlo calculations by capturing short-range correlation effects that are poorly described by single-particle wavefunctions.

## Overview

The Jastrow factor $J$ modifies the trial wavefunction:

$$\Psi_T = e^J \Psi_{SD}$$

where $\Psi_{SD}$ is the Slater determinant wavefunction from quantum chemistry calculations.

### Jastrow Components

The total Jastrow factor typically includes:

1. **Electron-nucleus terms (1-body)** - `a` parameters
      - Correlation between electrons and nuclei
      - Different parameters for each atom type
      - Captures electron-nucleus cusp

2. **Electron-electron terms (2-body)** - `b` parameters
      - Correlation between electron pairs (same and opposite spin)
      - Single set of parameters for all electron pairs

3. **Electron-electron-nucleus terms (3-body)** - `c` parameters
      - Three-body correlation effects
      - Different parameters for each atom type
      - Optional but improves accuracy


The Jastrow factor depends on the electronic ($\mathbf{r}$) and nuclear ($\mathbf{R}$) coordinates. Its defined as
$\exp(J(\mathbf{r},\mathbf{R}))$, where

$$
 J = f_{en} + f_{ee} + f_{een}
$$

Electron-nucleus and electron-electron:
$R={1-e^{-\kappa r} \over \kappa}$

$$
 f_{en} = \sum_{i=1}^{N_{\rm elec}} \sum_{\alpha=1}^{N_{\rm nuc}}
 \left( {a_1 R_{i\alpha} \over 1+a_2R_{i\alpha}} + \sum_{p=2}^{N^a_{\rm ord}} a_{p+1} R_{i\alpha}^p \right)
$$

$$
 f_{ee} = \sum_{i=2}^{N_{\rm elec}} \sum_{j=1}^{i-1} \left( {b_1 R_{ij} \over 1+b_2R_{ij}} + \sum_{p=2}^{N^b_{\rm ord}} b_{p+1} R_{ij}^p \right)
$$

Electron-electron-nucleus: $R=\exp\left(-\kappa r \right)$

$$
 f_{een} = \sum_{i=2}^{N_{\rm elec}} \sum_{j=1}^{i-1} \sum_{\alpha=1}^{N_{\rm nuc}} \sum_{p=2}^{N^c_{\rm ord}} \sum_{k=p-1}^0 \sum_{l=l_{\rm max}}^0 c_n R_{ij}^k (R_{i\alpha}^l+R_{j\alpha}^l) (R_{i\alpha}R_{j\alpha})^m
$$

where $m={p-k-l \over 2}$

-   Typically $N^a_{\rm ord}=N^b_{\rm ord}=5$. If $f_{een}$ is included,
    $N^c_{\rm ord}=5$.
-   Dependence among $\{c_n\} \rightarrow f_{een}$ does not
    contribute to cusp-conditions
-   $f_{en}$ and $f_{een}$: different $\{a_n\}$ and $\{c_n\}$ for
    different atom types


## File Format

Jastrow parameters are provided in a text file that CHAMP reads during initialization.

### Basic Structure

```perl
jastrow_parameter   1
  5  5  0           norda,nordb,nordc
   0.60000000         scalek
   0.00000000   0.00000000  -0.41907755  -0.22916790  -0.04194614   0.08371252 (a(iparmj),iparmj=1,nparma)
   0.00000000   0.00000000  -0.09956809  -0.00598089   0.00503028   0.00600649 (a(iparmj),iparmj=1,nparma)
   0.50000000   0.36987319   0.06971895   0.00745636  -0.00306208  -0.00246314 (b(iparmj),iparmj=1,nparmb)
 (c(iparmj),iparmj=1,nparmc)
 (c(iparmj),iparmj=1,nparmc)
end
```

### Format Specification

**Line 1**: Header
```perl
jastrow_parameter   1
```
The `1` indicates number of types of wavefunctions (typically always 1).

**Line 2**: Expansion orders
```perl
  5  5  0           norda,nordb,nordc
```

- `norda` = 5: $N^a_{\rm ord}$ Order of electron-nucleus expansion (1-body terms)

    if we are using pseudopotentials (no e-n cusps), we always leave
    $a_1=a_2=0$ and add
    $a_3 (r_{i\alpha}^2), \ldots, a_6 (r_{i\alpha}^5)$ equal to zero,
    which we then optimize. We do so for each atom type.

- `nordb` = 5: $N^b_{\rm ord}$ Order of electron-electron expansion (2-body terms)

    We set $b_1=0.5$ (for up-down e-e cusp condition), and add $b_3$
    ($r_{ij}^2$), $\ldots$, $b_6$ ($r_{ij}^5$) equal to zero, which we
    then optimize. $b_1$ is modified to 0.25 for up-up and down-down
    electrons.

- `nordc` = 0: $N^c_{\rm ord}$ Order of electron-electron-nucleus expansion (3-body terms, 0 = not present)

**Line 3**: Scaling parameter

```perl
   0.60000000         scalek
```
`scalek` scales the electron-electron distance in the Jastrow: commonly 0.5 or 0.6.

**Lines 4+**: Parameter sets

**a-parameters (electron-nucleus)**: One line per unique atom type

```perl
   0.00000000   0.00000000  -0.41907755  -0.22916790  -0.04194614   0.08371252
   0.00000000   0.00000000  -0.09956809  -0.00598089   0.00503028   0.00600649
```

- Number of `a` lines = number of unique atom types in geometry
- Order must match atom type order in `.xyz` file
- Number of values per line = `nparmja = 2 + max(0, norda - 1)`
  - For `norda = 5`: nparmja = 2 + 4 = 6 values
  - For `norda = 3`: nparmja = 2 + 2 = 4 values
  - For `norda = 0`: nparmja = 2 values
- Comments in parentheses are optional documentation

**b-parameters (electron-electron)**: Single line (if treating up and down spin electrons equally)
```perl
   0.50000000   0.36987319   0.06971895   0.00745636  -0.00306208  -0.00246314
```

- Only one `b` line for entire system (if treating up and down spin electrons equally)
- Number of values = `nparmjb = 2 + max(0, nordb - 1)`
  - For `nordb = 5`: nparmjb = 2 + 4 = 6 values
  - For `nordb = 3`: nparmjb = 2 + 2 = 4 values
  - For `nordb = 0`: nparmjb = 2 values
- Applies to all electron pairs (parallel and antiparallel spins)

**c-parameters (three-body)**: One line per unique atom type (if nordc > 0)
```perl
 (c(iparmj),iparmj=1,nparmc)
 (c(iparmj),iparmj=1,nparmc)
```

- Number of `c` lines = number of unique atom types
- Order must match atom type order in `.xyz` file
- Number of values per line = `nparmjc = nterms4(nordc)`
- If `nordc = 0`, c-parameter lines can be blank or omitted

**Final line**:
```
end
```

## Examples

### Example 1: Water Molecule (2 atom types: O, H)

**Geometry** (`molecule.xyz`):
```perl
3
Water molecule
  O   0.00000000   0.00000000   0.22143139
  H   0.00000000   1.43042809  -0.88572555
  H   0.00000000  -1.43042809  -0.88572555
```

**Jastrow file** (`jastrow.jas`):
```perl
jastrow_parameter   1
  5  5  0           norda,nordb,nordc
   0.60000000         scalek
   0.00000000   0.00000000  -0.39803661  -0.19727356  -0.04112912   0.08477303 (a for O)
   0.00000000   0.00000000  -0.25378471   0.03693657   0.02169054  -0.00707375 (a for H)
   0.50000000   0.48551602   0.09924779   0.00590014  -0.00626481  -0.00347973 (b)
end
```

**Explanation**:

- Two `a` lines: first for O (atom type 1), second for H (atom type 2)
- One `b` line for electron-electron correlation
- `nordc = 0`, so no `c` parameters needed

### Example 2: Water with 3-body Terms

**Geometry** (`molecule.xyz`):
```perl
3
Water molecule
  O   0.00000000   0.00000000   0.22143139
  H   0.00000000   1.43042809  -0.88572555
  H   0.00000000  -1.43042809  -0.88572555
```

**Jastrow file** (`jastrow.jas`):
```perl

jastrow_parameter   1
  5  5  5           norda,nordb,nordc
 0.40000  0.00000  scalek a21
 0.00000  0.00000 -0.04730 -0.39690 -0.15944  0.09709 (a(iparmj),iparmj=1,nparma)
 0.00000  0.00000 -0.01800  0.13474 -0.11232  0.00851 (a(iparmj),iparmj=1,nparma)
 0.50000  0.94655  0.14913 -0.13216  0.09053 -0.01728 (b(iparmj),iparmj=1,nparmb)
 1.69203  1.24656 -1.88214  3.91513 -0.74384 -8.64444 -3.45930  7.94537 -0.61759 -7.28154 -3.20830 10.10651  2.84683  3.11254 -1.82060 -4.12168 -2.73375  7.97867  2.63741  0.58993 -3.22734 -3.43777 -0.80525 (c(iparmj),iparmj=1,nparmc)
-0.18440 -2.24997 -1.73422  0.47251  1.54643  6.06904  3.29715  1.00667 -3.97339 -1.90822  2.48432 -1.92995 -5.13868 -1.55965 -3.12079  4.42562  3.87653 -4.16459  0.11133 -3.79462  3.68329  0.66761  1.96643 (c(iparmj),iparmj=1,nparmc)
end
```

**Explanation**:

- Two atom types: O and H
- `norda = 5`, `nordb = 5`, `nordc = 5` (includes fully optimized 3-body terms)
- Two `a` lines (O and H), each with 6 parameters
- One `b` line with 6 parameters
- Two `c` lines (O and H), each with 23 parameters (from `nterms4(5)`)

### Example 3: Initial Guess (Simple)

For starting VMC calculations, a minimal Jastrow:

```perl
jastrow_parameter   1
  3  3  0           norda,nordb,nordc
   0.60000000         scalek
   0.0   0.0  -0.25  -0.10 (a - adjust for each atom type)
   0.5   0.2   0.05   0.0  (b)
end
```

This provides a reasonable starting point for optimization.

## Loading Jastrow in CHAMP Input

Specify the Jastrow file using the `load jastrow` command:

```perl
%module general
    title  'VMC with Jastrow'
    pool   './pool/'
%endmodule

load trexio          $pool/molecule.hdf5
load jastrow         jastrow.jas
load jastrow_der     jastrow.der  # Optional, for optimization

%module electrons
    nup    5
    nelec  10
%endmodule
```

## Common Jastrow Forms

### Standard Padé-Jastrow

Most common form in CHAMP:

$$J = \sum_I \sum_i \frac{a_{I}(r_{iI})}{1 + a_0 r_{iI}} + \sum_{i<j} \frac{b(r_{ij})}{1 + b_0 \cdot \text{scalek} \cdot r_{ij}}$$

where:

- $a_I(r)$ are polynomials in electron-nucleus distance $r_{iI}$
- $b(r)$ is a polynomial in electron-electron distance $r_{ij}$
- The denominator ensures proper asymptotic behavior

### Cusp Conditions

Jastrow factors must satisfy:

- **Electron-nucleus cusp**: $\frac{\partial J}{\partial r_{iI}}\bigg|_{r=0} = -Z_I$ for electron $i$ at nucleus $I$
- **Electron-electron cusp**: $\frac{\partial J}{\partial r_{ij}}\bigg|_{r=0} = \pm\frac{1}{2}$ for parallel/antiparallel spins


## Related Topics

- [Jastrow Derivatives](jastrow_derivatives.md) - Required for wavefunction optimization
- [Wavefunction Optimization](../calculations/optimization/index.md) - Optimizing Jastrow parameters
- [VMC Calculations](../calculations/vmc/index.md) - Using optimized Jastrows
- [TREXIO Files](using_trexio_file.md) - Jastrows must be supplied separately

## Tips

- Start with simple Jastrow forms before adding complexity
- Check that atom type ordering matches geometry exactly
- Verify parameter counts match declared orders (norda, nordb, nordc)
- Use VMC to test Jastrow quality (variance should decrease with optimization)
- Consult [Optimization Guide](../calculations/optimization/index.md) for parameter tuning
