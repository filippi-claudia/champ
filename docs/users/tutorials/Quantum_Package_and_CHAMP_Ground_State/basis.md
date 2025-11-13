---
layout: default
title: Basis Sets and Pseudopotentials
nav_order: 1
grand_parent: Tutorials
parent: '01. QP and CHAMP : Ground State Calculation'
authors:
    - Ravindra Shinde
    - Claudia Filippi
    - Anthony Scemama
tags:
    - CHAMP
    - tutorial
    - CIPSI
    - ground state
    - QP
---

# Basis Sets and Pseudopotentials

For QMC calculations, we need to use pseudopotentials optimized
specifically for QMC, and basis sets optimized to be used with these
pseudopotentials. Here, we use the
[Burkatzki-Filippi-Dolg](http://burkatzki.com/pseudos/index.2.html)
(BFD) ones except for hydrogen (the hydrogen pseudo on the website is
too soft and not sufficiently accurate).

QP can read basis sets and pseudopotentials from files in GAMESS format,
if the files exist in the current directory. Otherwise, it will try to
look into its own database of basis sets and pseudopotentials.

### BFD Pseudopotential


Store the pseudopotential parameters in a file named `PSEUDO`:

```python
H GEN 0 0
3
 1.000000000000 1 25.000000000000
25.000000000000 3 10.821821902641
-8.228005709676 2  9.368618758833

O GEN 2 1
3
6.00000000 1 9.29793903
55.78763416 3 8.86492204
-38.81978498 2 8.62925665
1
38.41914135 2 8.71924452
```

### Double-Zeta basis set


Store the basis set parameters in a file named `BASIS`:

```python
HYDROGEN
s 3
1 6.46417546   0.063649375945
2 1.13891461   0.339233210576
3 0.28003249   0.702654522063
s 1
1 0.05908405   1.00000000
p 1
1 0.51368060   1.00000000

OXYGEN
s 9
1 0.125346     0.055741
2 0.268022     0.304848
3 0.573098     0.453752
4 1.225429     0.295926
5 2.620277     0.019567
6 5.602818     -0.128627
7 11.980245     0.012024
8 25.616801     0.000407
9 54.775216     -0.000076
s 1
1 0.258551     1.000000
p 9
1 0.083598     0.044958
2 0.167017     0.150175
3 0.333673     0.255999
4 0.666627     0.281879
5 1.331816     0.242835
6 2.660761     0.161134
7 5.315785     0.082308
8 10.620108     0.039899
9 21.217318     0.004679
p 1
1 0.267865     1.000000
d 1
1 1.232753     1.000000
```