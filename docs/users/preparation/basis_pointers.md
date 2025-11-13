---
layout: default
title: Basis pointers
nav_order: 5
parent: Input files
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - basis info
    - basis pointers
---

# Basis pointers (formerly bfinfo) files

The new format of the basis pointers file is given below. This file should be kept in the `pool` directory.
This file is generated automatically by the `trex2champ.py` converter.

```python
# Format of the new basis information file champ_v3
# num_ao_per_center, n(s), n(p), n(d), n(f), n(g)
# Index of Slm (Range 1 to 35)
# Index of column from numerical basis file
qmc_bf_info 1
54 4 4 3 2 0
1 1 1 1 2 3 4 2 3 4 2 3 4 2 3 4 5 6 7 8 9 10 5 6 7 8 9 10 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 11 12 13 14 15 16 17 18 19 20
1 2 3 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 9 9 9 10 10 10 10 10 10 11 11 11 11 11 11 12 12 12 12 12 12 12 12 12 12 13 13 13 13 13 13 13 13 13 13
35 4 3 2 1 0
1 1 1 1 2 3 4 2 3 4 2 3 4 5 6 7 8 9 10 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
1 2 3 4 5 5 5 6 6 6 7 7 7 8 8 8 8 8 8 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10
end
```

Each unique type of atom will have a pair of lines in the basis pointers file.

The first line after the comments `qmc_bf_info 1` is a specification line to make sure that we are reading basis function information file.

The second line is for the first unique atom in the system. It contains the number of atomic orbitals for that atom, the number of s-type functions, number of p-type functions, number of d-type functions, number of f-type functions, and number of g-type functions.
`num_ao_per_center, n(s), n(p), n(d), n(f), n(g)`

The third line gives the index of Slm (or real Ylm). The numbers depend on how many radial shells are there in the basis set.

The fourth line tells which column of the radial grid file to be read for the construction of MO from the AOs.


{: .warning }
Each unique type of atom will have a corresponding set of basis pointers in the file.
