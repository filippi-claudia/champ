---
layout: default
title: Basis on the grid
nav_order: 4
parent: Input files
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - basis
---

# Basis set (Basis on the radial grid) files

Basis files have a fixed format. The files generated from the trex2champ converter can also be used as they are.
These files are generally kept in the `pool` directory of the calculation folder. You just need to specify the name of the basis file (say, `ccpVTZ`) in the general module of CHAMP input file under the keyword `basis`. This will read the file `ccpVTZ.basis.Si` for the element `Si`.

The top few lines of `BFD-T.basis.C` look like

```python
9 3 2000 1.003000 20.000000 0
 0.000000000000e+00  5.469976184517e-01  2.376319920758e+00  5.557936498748e-01  3.412818210005e+00  2.206803021951e-01  8.610719484857e-01  3.738901952004e-01  3.289926074834e+00  1.106692909826e+00
 1.508957441883e-04  5.469976454488e-01  2.376319870895e+00  5.557936481942e-01  3.412817957941e+00  2.206803015581e-01  8.610719410992e-01  3.738901923954e-01  3.289925989316e+00  1.106692890335e+00
 ...
 ```
This means there are 9 radial shells in the basis set of carbon put on a radial grid of 2000 points (upto 20 bohr).


{: .warning }
Basis set files are sourced from the pool directory. All unique atoms should have a corresponding basis file.
