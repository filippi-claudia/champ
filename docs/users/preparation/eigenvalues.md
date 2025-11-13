---
layout: default
title: MO Eigenvalues
nav_order: 9
parent: Input files
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - eigenvalues
---

# Molecular orbital eigenvalues file [Optional]
This file is also generated using the `trex2champ.py` converter if the parent .hdf5 file contains the orbital eigenvalues.

A typical file looks like this:

```python
# File created using the trex2champ converter https://github.com/TREX-CoE/trexio_tools
# Eigenvalues correspond to the RHF orbitals
eigenvalues 64
-1.3659 -0.7150 -0.5814 -0.5081 0.1201 0.1798 0.4846 0.5148 0.5767 0.6085 0.7153 0.7820 0.8691 0.8699 0.9642 1.2029 1.4091 1.4388 1.6082 1.6342 2.0787 2.1179 2.1776 2.2739 2.4123 2.5591 2.8217 3.3480 3.3840 3.4544 3.4607 3.6199 3.6237 3.9628 3.9661 4.0439 4.0481 4.2212 4.3500 4.4225 4.4577 4.5747 4.7271 4.8382 5.0086 5.5800 5.8020 6.0317 6.3754 6.5827 6.6970 6.7474 6.9245 7.0790 7.1820 7.2121 7.3257 7.3865 7.8607 8.4146 8.4733 9.0201 16.4980 27.1462
end

```

The first line contains a keyword `eigenvalues` followed by the number of orbitals. The following line contains eigenvalues as they appear in GAMESS or similar output. The file ends with keyword `end`.
