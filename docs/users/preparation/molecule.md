---
title: Molecular Geometry
tags:
    - geometry
    - molecule
---

Molecular coordinates can be provided in the vmc or dmc input files a text file or a separate .xyz file in one of the following ways:

## 1. Geometry in the (XYZ in Bohr units) format to be read from a separate .xyz file.

```bash
load molecule  molecule.xyz
```

```bash
load molecule  $pool/molecule.xyz
```

The following are the valid examples (molecule.xyz)

## 2. Geometry in the (XYZ in Bohr units) format with automatic Zvalence
```perl
10
# molecular complex (Symbol, X,Y,Z in Bohr)
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

!!! danger "Mind the units"

    The atomic coordinates mentioned above are still in Bohr units.


## 3. Geometry in the (XYZ in Bohr units) format with explicit Zvalence. This also allows different labels for the same element.
```perl
10
# molecular complex (Symbol, X,Y,Z in Bohr, Zvalence)
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

## 4. Geometry in the (XYZ in Bohr units) format to be read from a separate .xyz file.

```bash
%block molecule < molecule.xyz
```
