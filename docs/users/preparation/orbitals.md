---
layout: default
title: Molecular Orbitals
nav_order: 6
parent: Input files
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - orbitals
---

# Molecular Orbitals file

This file contains the molecular orbital coefficients. These are arranged as [num_ao, num_mo] array. This file is obtained automatically from the `trex2champ.py` converter. Please note that the AOs in this file follow the trexio convention of AO ordering.

For example,
Four p-type shells of AOs will be arranged alphabetically as

`X Y Z   X Y Z   X Y Z   X Y Z`

Two d-type shells of AOs will be arranged alphabetically as

`XX XY XZ YY YZ ZZ   XX XY XZ YY YZ ZZ`

and so on.

The `.lcao` or `.orb` file has the following format.

```python
lcao  226 200  1
...
...

end
```

The number 226 will be number of AOs, 200 will be number of orbitals, 1 will be number of types of orbitals.
