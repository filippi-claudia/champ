---
layout: default
title: Jastrow
nav_order: 10
parent: Input files
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - jastrow
---

# Jastrow parameters file
The Jastrow parameters can be provided using this file. It has the following format [Example: water].

```python
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

The set `a`should appear for each unique atom type (in the same order as in the .xyz file).

The set `b` should appear once.

The three-body Jastrow terms `c` should appear for each unique atom type (in the same order as in the .xyz file)
