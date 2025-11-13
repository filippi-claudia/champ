---
layout: default
title: Pseudopotentials
nav_order: 3
parent: Input files
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - ECP
    - pseudopotential
---

# ECP / Pseudopotential files

ECP or pseudopotential files have a fixed format. Most of the BFD ECP files can be found in the `champ/pool/BFD/ECP_champ` folder. The files generated from the trexio file can also be used (except if it is coming from GAMESS. In this case, GAMESS truncates the digits of ECP information in its output, so the trexio file will not have all the digits stored.)

File format: BFD ECP for Silicon

`BFD.gauss_ecp.dat.Si`

```perl
BFD Si pseudo
3
3
4.00000000 1 1.80721061
7.22884246 3 9.99633089
-13.06725590 2 2.50043232
1
21.20531613 2 2.26686403
1
15.43693603 2 2.11659661
```
These files are generally kept in the `pool` directory of the calculation folder. You just need to specify the the name `BFD` in the general module of CHAMP input file under the keyword `pseudopot`. There should be a file for each type of an atom.

```python
%module general
    title           'VMC Calculation for a molecule'
    pool            './pool/'
    mode            'vmc'
    seed            1138139413245321
    pseudopot       BFD
    basis           ccpVTZ
    ipr             -1
%endmodule
```

{: .warning }
Note that GAMESS output file does not contain all the digits of the ECP information. So, the trexio file generated from GAMESS output will not have all the digits stored. This will cause a problem in the calculation. So, it is recommended to supply separately the ECP files if you have a trexio as input file.

