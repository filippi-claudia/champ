**Pseudopotentials**

`http://www.burkatzki.com/pseudos/index.2.html`

Do not use the H pseudopotential/basis (in principle, avoid basis sets of 1- and 2-electron atoms) from this database.

For H pseudopotential, we still need to check the double basis (extended fit with weight(H)=1 or 10). From VTZ up, the basis sets are good.

Check under `pool/BFD/BASIS_gamess/`

**Gamess**

If we do not use external basis file, we need to repeat basis set for each atom (cumbersome for large molecules).

External basis file was created by Omar and modified by Riccardo, and available as `champ/pool/BFD/BASIS_gamess/BFD_Basis.EXTBAS`

**QMC codes**

make clean, clean_all, seq, mpi, all, vmc.mov1, etc.

Check Makefile to see what is available. 

dimensions in `~/champ/include/<...>.h`

For compiler commands, check make_config directory and file settings.make (must be symbolically linked to the correct settings.make.ifort etc.)

**Setup of input files of QMC**

If you type:

`~/champ/interface/gamess2qmc`

you can see all options.

*Example:*

`~/champ/interface/gamess2qmc -g -t rhf -n 6 -r -b BFD.pVDZ rhf.out`

g -> geometry
t -> type orbitals
n -> number of orbitals
r -> tabulate radial atomic basis (radial AO) on a grid
b -> name of the basis

It generates:

BFD.pVDZ.lcao -> LCAO coefficients for MO orbitals
BFD.pVDZ.geometry  -> geometry
BFD.pVDZ.basis.C.1 -> radial AO basis for C
BFD.pVDZ.basis.H.2 -> radial AO basis for H
BFD.pVDZ.bfinfo -> Information on the basis (relation radial -> Ylm: s, p etc?)

* If you change the wave function (e.g. rhf -> b3lyp) using the same 
basis (for the same system) and run again a QMC run, you do not need 
to regenerate the atomic basis/geometry.

`~/champ/interface/gamess2qmc -t rhf -n 6 b3lyp.out`

* Generate file for optimization of orbitals

`~/programs/champ_source/interface/gamess2qmc -t rhf -s rhf_pVDZ.out`

Gets all orbitals and symmetry (-s) file.

* After MCSCF/CASSF in Gamess, you need to run a CI starting from the MCSCF orbitals
to have the correct CI coefficients.

Then, you will get the orbitals as the INITIAL of the CI run (put PRTMO=.TRUE.);
include the MCSCF orbitals in the CI input as the VEC field (from the mcscf.dat).
You may choose to use natural orbitals from the VEC file (if you truncate in 
QMC, the surviving determinants might be more representative). If you want to
optimize and start from natural orbitals, recall to include the virtual MCSCF 
orbitals in the VEC file of the CI run.

* To get the CI coefs from the CI output:

`~/programs/champ_source/interface/gamess2qmc -t initial -d 0.0 -s ci22_pVDZ.out`

`-d 0.0` -> the threshold for the coefficients is set to zero (get them all)

* If you are running more that one state (after an SA-MCSCF run), you will need
to specify which states you want as 

`~/programs/champ_source/interface/gamess2qmc -t initial -d 0.1 -s -w 1,2 ci66_pVDZ.out`

where you are getting (for example) states 1 and 2.

* To run vmc etc.

There are different excecutables to run vmc, dmc, serial and parallel.

`~/codes/champ/bin/vmc.mov1 < vmc.inp > vmc.out &`

Use the vmc.mov1, dmc.mov1. vmc.mov1.mpi etc., namely, the single-electron move versions.

* MPI champ

Recall to create a file 'filename' with the name of the inputfile.
vmc.mov1.mpi etc. will open 'filename' to read the name of the input.

* From VMC (1 walker) to DMC (nconf 100 or more)

You need to generate the walkers in a VMC run by setting for instance
(nconf_new 100). The VMC run generates 'mc_configs_new' (for parallel run,
mc_configs_newIPROC). Careful that nstep\*nblk > Tcorr\*10\*nconf_new.

Then, you create a file `mc_configs` which will be read by DMC.
On 1 proc: 

`mv mc_configs_new mc_configs`

on N procs:

`cat mc_configs_new* >> mc_configs`
`rm mc_configs_new*`

The file contains 100\*NPROC 3N coordinates.

These coordinates will be the starting population of walkers for DMC.

