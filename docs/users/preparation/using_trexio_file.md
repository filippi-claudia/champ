---
layout: default
title: Using TREXIO as input
nav_order: 1
parent: Input files
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - TREXIO
---

# Using a single trexio file as input

We can use trexio file (in hdf5 or text backend format) to specify all the inputs (except Jastrow and Jastrow derivatives)

A sample input file would look like:

```perl
%module general
    title           'VMC Calculation for a molecule'
    pool            './pool/'
    mode            'vmc'
    seed            1138139413245321
    ipr -1
%endmodule

load trexio          molecule.hdf5
load jastrow         jastrow.jas
load jastrow_der     jastrow.der

%module electrons
    nup           20
    nelec         40
%endmodule


%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      100000
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

## Obtaining a trexio file from GAMESS-US output

{: .warning }
Make sure that the recent version of `trexio_tools` has been installed.


```bash
pip install trexio_tools
```

This will provide `trexio` executable in the path. Use the following command to generate a trexio file.

```bash
trexio convert-from --type gamess --input gamess.out --motype "RHF"  --back_end=HDF5 sample.hdf5
```

Allowed values of MOtype are `'RHF', 'ROHF', 'MCSCF', 'NATURAL', 'GUGA' ...`


{: .important }
Use `trexio --help` for a verbose list of options.


## Converting TREXIO file into text inputs


The trexio file can be converted into several text files to be used with CHAMP. The python converter is provided in the CHAMP's repository in the `champ/tools/trex_tools` folder.

A sample script is given below:

```python
python3 /home/user/champ/tools/trex_tools/trex2champ.py \
	--trex 	"COH2_GS.trexio" \
	--backend	"HDF5" \
	--basis_prefix  "BFD-aug-cc-pVDZ" \
	--lcao \
	--ecp \
	--sym \
	--geom \
	--basis \
	--det
```

{: .important }
Use `python3 trex2champ.py --help` for a verbose list of options.