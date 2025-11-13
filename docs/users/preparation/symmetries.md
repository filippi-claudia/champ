---
layout: default
title: MO symmetries
nav_order: 8
parent: Input files
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - symmetries
---

# Molecular orbital symmetries file [Optional; useful when doing orbital optimization]
This file is also generated using the `trex2champ.py` converter if the parent .hdf5 file contains the orbital symmetries.

A typical file looks like this:

```python
sym_labels 4 226
 1 AG 2 AU 3 BG 4 BU
1 4 4 1 1 4 1 4 1 2 3 2 3 4 1 4 1 4 4 1 1 4 4 1 4 1 3 1 4 1 2 4 1 2 4 1 3 2 4 3 2 1 4 4 3 1 1 4 4 4 2 1 3 1 4 1 1 4 1 4 3 1 4 2 2 3 1 4 1 4 1 1 4 2 3 4 1 4 2 1 3 1 4 1 4 2 4 4 1 3 4 1 3 4 2 1 2 3 4 1 2 4 1 3 4 2 3 1 1 4 4 1 2 1 3 1 4 1 4 2 3 4 1 4 2 1 4 3 1 4 2 3 2 3 4 1 2 3 1 2 4 2 3 4 1 4 3 2 1 1 3 4 4 1 4 1 2 4 1 3 1 2 4 4 4 3 1 1 3 1 1 2 2 4 4 2 1 4 3 1 1 4 3 4 2 1 1 2 4 3 4 3 2 1 3 4 1 3 1 4 4 2 1 4 1 4 1 1 4 4 4 1 1 1 4 1 4 4 1 4 1 4 1 4 1 4
end
```

The numbers in front of irreducible representations are used as correspondence to identify the symmetry type of each orbital. Here in this case there are 226 molecular orbitals with 4 irreps.
