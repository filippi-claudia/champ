---
layout: default
title: QMC runs (a) Check setup
nav_order: 5
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
    - QP
---

# QMC runs : Check that the setup is OK

First, we can compute with QP the energies of the single-determinant
wave functions with the 2 different sets of MOs.

```bash
qp set_file h2o_hf
qp run print_energy

qp set_file h2o_dft
qp run print_energy
```

These commands return the energy of the wavefunction contained in the
EZFIO database. These values will be useful for checking that the QMC
setup is OK. You should obtain the energies:

| Molecular orbitals type | Total Energy |
|:-------------:|:-----------------:|
| HF            | -16.9503842       |
| DFT           | -16.9465884       |


We will now convert the TREXIO files into input files suitable for
CHAMP:

{: .important}
You need the `resultsFile` and `trexio` Python packages. They can be
installed with pip as described in []().


Create a new directory named `H2O_HF` and copy the TREXIO file
`h2o_hf.trexio` into it. Go inside this directory and run

```bash
python3 ~filippi/Tutorial-QMC-School/trex2champ.py --trex "h2o_hf.trexio" \
                       --motype  "Canonical" \
                       --backend "HDF5" \
                       --basis_prefix "BFD-cc-pVDZ" \
                       --lcao \
                       --geom \
                       --basis \
                       --ecp \
                       --det
```

Many files were created. Now, create a directory named `pool`, and move
some files into the pool:

```bash
mkdir pool
mv *.xyz *bfinfo BFD-* ECP* pool
```

You can now create an input file for CHAMP `vmc_h2o_hf.inp` :

```python
%module general
    title           'H2O HF calculation'
    pool            './pool/'
    pseudopot       ECP
    basis           BFD-cc-pVDZ
    mode            'vmc_one_mpi1'
%endmodule


load molecule        $pool/champ_v2_h2o_hf_geom.xyz
load basis_num_info  $pool/champ_v2_h2o_hf_with_g.bfinfo

load orbitals        champ_v2_h2o_hf_orbitals.lcao
load determinants    champ_v2_h2o_hf_determinants.det
load jastrow         jastrow.start

%module electrons
    nup           4
    nelec         8
%endmodule


%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      20000
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

Create the file for the Jastrow factor as follows, and save it as
`jastrow.start`:

```python
jastrow_parameter   1
  0  0  0           norda,nordb,nordc
   0.60000000   0.00000000     scalek,a21
   0.00000000   0.00000000   (a(iparmj),iparmj=1,nparma)
   0.00000000   0.00000000   (a(iparmj),iparmj=1,nparma)
   0.00000000   1.00000000   (b(iparmj),iparmj=1,nparmb)
 (c(iparmj),iparmj=1,nparmc)
 (c(iparmj),iparmj=1,nparmc)
end
```

This files implies that there is no Jastrow factor $$\exp(J)=1$$.

Create the submission script as presented in [](), and submit the job. You should obtain the
Hartree-Fock energy.

Now reproduce the same steps for the TREXIO file containing the DFT
orbitals in directory `H2O_DFT`.

The energies obtained with VMC without the Jastrow factor should be the
same as those computed by QP at the beginning of this section.

