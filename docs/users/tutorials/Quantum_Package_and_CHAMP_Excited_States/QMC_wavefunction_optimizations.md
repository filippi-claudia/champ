---
layout: default
title: QMC wave function optimizations
nav_order: 4
grand_parent: Tutorials
parent: '02. QP and CHAMP : Excited State Calculation'
mathjax: true
authors:
    - Ravindra Shinde
    - Claudia Filippi
    - Anthony Scemama
tags:
    - CHAMP
    - tutorial
    - CIPSI
    - excited state
    - QP
---

# QMC wave function optimizations

In this section, we will optimize a Jastrow factor for each state, and
we will the re-optimize the CI coefficients in the presence of the
Jastrow. The setup of the CHAMP files is similar to what we have done in the [Setup]({{ site.baseurl }}{% link docs/Tutorials/Quantum_Package_and_CHAMP_GS/basis.md %}#basis-sets-and-pseudopotentials)

{: .warning}
Here, we have 12 electrons, 6 up and 6 down.


## Optimization of the ground state

Create a new directory, and copy the `COH2_GS.trexio` TREXIO file inside
it. Convert the TREXIO file into CHAMP files:

```bash
python3 /path/champ/tools/trex-tools/trex2champ.py \
                        --trex          "COH2_GS.trexio" \
                        --motype        "Canonical" \
                        --backend       "HDF5" \
                        --basis_prefix  "BFD-aug-cc-pVDZ" \
                        --lcao \
                        --geom \
                        --basis \
                        --ecp \
                        --sym \
                        --det
```

$$COH_2$$ has three different atom types, so the Jastrow factor file will
be slightly different from the file for water with one extra line for
$a$ parameters. You can start by creating a file called `jastrow.start`:

```python
jastrow_parameter   1
  5  5  0           norda,nordb,nordc
   0.60000000   0.00000000     scalek,a21
   0.00000000   0.00000000  0. 0. 0. 0. (a(iparmj),iparmj=1,nparma) ! e-n C
   0.00000000   0.00000000  0. 0. 0. 0. (a(iparmj),iparmj=1,nparma) ! e-n O
   0.00000000   0.00000000  0. 0. 0. 0. (a(iparmj),iparmj=1,nparma) ! e-n H
   0.50000000   1.00000000  0. 0. 0. 0. (b(iparmj),iparmj=1,nparmb)
 (c(iparmj),iparmj=1,nparmc) ! e-e-n C
 (c(iparmj),iparmj=1,nparmc) ! e-e-n O
 (c(iparmj),iparmj=1,nparmc) ! e-e-n H
end
```

Start by optimizing the Jastrow factor and perform a \"quick\"
optimization. The following champ input file (`vmc_quick.inp`) contains
the parameters for such a \"quick\" optimization.

```python
%module optwf
    ioptwf        1
    ioptci        0
    ioptjas       1
    ioptorb       0
    method        'sr_n'

    ncore         0
    nextorb       600
    nblk_max      5000

    nopt_iter     10
    sr_tau        0.05
    sr_eps        0.01
    sr_adiag      0.01
%endmodule


%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      20
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

Move the `jastrow_optimal.1.iter10` file to `jastrow_optimal.rough_GS`
and load this optimized Jastrow factor. You can now optimize also the CI
coefficients together with the Jastrow factor by setting:

```python
ioptci        1
```

Use some more Monte Carlo steps to perform a more strict optimization.

```python
%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      500
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

In your directory, you will now have `jastrow_optimal.1.iterX` and
`det_optimal.1.iterX` files.

Set up a DMC calculation where you use the optimal Jastrow and CI
coefficients. Adjust the `etrial` to be a bit below the VMC energy.
Recall that you will have to generate the `mc_configs` file.

{: .new-title}
> Tip
>
> You could have also optimized the orbitals but we did not do it
here to keep the calculations short. If you are setting `optorb=1`, load
also the symmetry file.
>
>```python
>load symmetry champ_v2_COH2_GS_symmetry.sym
>```


## Optimization of the excited state


Create a new directory, and copy the `COH2_ES.trexio` TREXIO file inside
it. Apply the same procedure as for the ground state.

Repeat what you have done for the ground state. Start to perform a quick
optimization of the Jastrow factor but do not start from zero\'s: start
from the Jastrow factor you have for the ground state, namely,
`jastrow_optimal.rough_GS`.

Do all step until when you obtain the DMC energy.

Compute the VMC and DMC excitation energies. Recall that if your
energies are
$$
  E_{\rm GS}+\delta E_{\rm GS}  \text{ and } E_{\rm ES}+\delta E_{\rm ES},
$$

where $$\delta E$$ is the statistical error, the error on

$$\Delta E= E_{\rm ES}-E_{\rm GS}$$
is given by

$$
  \sqrt{\delta E_{\text{GS}}^2+\delta E_{\text{ES}}^2}.
$$
