---
layout: default
title: QMC runs (b) Introduce and optimize Jastrow
nav_order: 6
grand_parent: Tutorials
parent: '01. QP and CHAMP : Ground State Calculation'
mathjax: true
authors:
    - Ravindra Shinde
    - Claudia Filippi
    - Anthony Scemama
tags:
    - CHAMP
    - tutorial
    - CIPSI
    - ground state
    - jastrow
    - QP
---

# Introduce and optimize Jastrow

The Jastrow factor depends on the electronic ($$\mathbf{r}$$) and nuclear ($$\mathbf{R}$$) coordinates. Its defined as
$$\exp(J(\mathbf{r},\mathbf{R}))$$, where

$$
 J = f_{en} + f_{ee} + f_{een}
$$

Electron-nucleus and electron-electron:
$$R={1-e^{-\kappa r} \over \kappa}$$

$$
 f_{en} = \sum_{i=1}^{N_{\rm elec}} \sum_{\alpha=1}^{N_{\rm nuc}}
 \left( {a_1 R_{i\alpha} \over 1+a_2R_{i\alpha}} + \sum_{p=2}^{N^a_{\rm ord}} a_{p+1} R_{i\alpha}^p \right)
$$

$$
 f_{ee} = \sum_{i=2}^{N_{\rm elec}} \sum_{j=1}^{i-1} \left( {b_1 R_{ij} \over 1+b_2R_{ij}} + \sum_{p=2}^{N^b_{\rm ord}} b_{p+1} R_{ij}^p \right)
$$

Electron-electron-nucleus: $$R=\exp\left(-\kappa r \right)$$

$$
 f_{een} = \sum_{i=2}^{N_{\rm elec}} \sum_{j=1}^{i-1} \sum_{\alpha=1}^{N_{\rm nuc}} \sum_{p=2}^{N^c_{\rm ord}} \sum_{k=p-1}^0 \sum_{l=l_{\rm max}}^0 c_n R_{ij}^k (R_{i\alpha}^l+R_{j\alpha}^l) (R_{i\alpha}R_{j\alpha})^m
$$

where $$m={p-k-l \over 2}$$

-   Typically $$N^a_{\rm ord}=N^b_{\rm ord}=5$$. If $$f_{een}$$ is included,
    $$N^c_{\rm ord}=5$$.
-   Dependence among $$\{c_n\} \rightarrow f_{een}$$ does not
    contribute to cusp-conditions
-   $$f_{en}$$ and $$f_{een}$$: different $$\{a_n\}$$ and $$\{c_n\}$$ for
    different atom types

## Add a simple e-e and e-n Jastrow factor

-   $$N^a_{\rm ord}=5$$

    Since we are using pseudopotentials (no e-n cusps), we always leave
    $$a_1=a_2=0$$ and add
    $$a_3 (r_{i\alpha}^2), \ldots, a_6 (r_{i\alpha}^5)$$ equal to zero,
    which we then optimize. We do so for each atom type.

-   $$N^b_{\rm ord}=5$$

    We set $$b_1=0.5$$ (for up-down e-e cusp condition), and add $$b_3$$
    ($$r_{ij}^2$$), $$\ldots$$, $$b_6$$ ($$r_{ij}^5$$) equal to zero, which we
    then optimize. $$b_1$$ is modified to 0.25 for up-up and down-down
    electrons.

    The following file is your starting Jastrow factor `jastrow.start`:

    ```python
    jastrow_parameter   1
      5  5  0           norda,nordb,nordc
       0.60000000         scalek
       0.00000000   0.00000000 0. 0. 0. 0. (a(iparmj),iparmj=1,nparma) ! e-n O
       0.00000000   0.00000000 0. 0. 0. 0. (a(iparmj),iparmj=1,nparma) ! e-n H
       0.50000000   1. 0. 0. 0. 0. (b(iparmj),iparmj=1,nparmb) ! e-e
     (c(iparmj),iparmj=1,nparmc) ! e-e-n O
     (c(iparmj),iparmj=1,nparmc) ! e-e-n H
    end
    ```

## Optimize the Jastrow factor

Create the file `jastrow.der`:

```python
jasderiv
4 4 5 0 0 0 0 nparma,nparmb,nparmc,nparmf
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma) ! e-n O
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma) ! e-n H
2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb) ! e-e
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c(iparmj),iparmj=1,nparmc)
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c(iparmj),iparmj=1,nparmc)
end
```

where you are telling CHAMP to optimize $$a_i, 3\le i \le 6$$ for e-n of O
and H (4 parameters for both O and H), and $$b_i, 2 \le i \le 6$$ (5
parameters in total).

Now, specify the name of the info of the derivatives of the Jastrow in
the input file, below the line where the `jastrow.start` file is
specified. You also need to add a block with different options for the
optimizer as follows.

```python
load jastrow         jastrow.start
load jastrow_der     jastrow.der

%module optwf
    ioptwf        1
    ioptci        0
    ioptjas       1
    ioptorb       0

    method        'sr_n'
    nopt_iter     20
    nblk_max      4000

    ncore         0
    nextorb       100

    sr_tau        0.05
    sr_eps        0.001
    sr_adiag      0.01
%endmodule
```

Optimization doesn\'t require a long QMC simulation in the first SR
steps. You can reduce the number of blocks in `blocking_vmc` to 100, and
the code will slowly increase the number of blocks to `nblk_max` in the
`optwf` module.

```python
%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      100
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

If you `grep 'total E' output`, you will see the optimization
progressing and generating new Jastrow factors in
`jastrow_optimal.1.iterX`.

If you `grep nblk output` you will see that the code automatically
increases the maximum number of blocks.