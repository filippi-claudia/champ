# Problematic Common blocks

Common blocks that were not deleted in the second round. 

- [x] ***#3 common block /Bloc/***\
**PROBLEM:** \
Variables xmatd, xmatu are shared with common block /scratch/.
\
**SOLUTION:** \
The common block scratch was used only by optjas.f. Solution is simple: delete xmatd, xmatu from scratch module and introduce a "use Bloc, only: xmatd, xmatu" statement.\
FYI: \
common /Bloc/ \
b(MORB,MELEC), xmat(MELEC * * 2,2), xmatu(MELEC* * 2), xmatd(MELEC * * 2)

- [x] ***#17 common block /contrl/***\
**PROBLEM:** \
contrl shares variable nconf with common block /sr_mat_n/.\
**SOLUTION:** \
Variable nconf was changed to n_conf (only used in read_input.f for dmc calculations) before doing the python refactoring, that is: \
common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar\
To:\
common /contrl/ nstep,nblk,nblkeq,n_conf,nconf_new,isite,idump,irstar

- [ ] ***#21 common block /cuspmat4/***\
**PROBLEM:** \
It requires to move the parameters NEQSX and MTERMS to an include file (done!) and to solve icusp(NEQSX) array issue: this variable also appears as an integer and an array with dimensions icusp(MCTYPE) throughout the code.\ 
**SOLUTION:** \
TODO

- [x] ***#22 common block /da_energy_ave/***\
**PROBLEM:** \
Common block force_fin uses array with same name as da_energy_ave.\
**SOLUTION:** \
To change the name of the (common block) module from da_energy_ave to da_energy_ave_m.

- [x] ***#27  common block /deloc_dj/***\
**PROBLEM:** \
optjas.f allocates an array with the same name as deloc_dj.\
To change the name of the (common block) module from deloc_dj to deloc_dj_m.

- [x] ***#28 common block /denergy_det/***\
**PROBLEM:** \
Common block and variable have the same name, e.i:\
common /denergy_det/ denergy_det(MDET,2)\
**SOLUTION:** \
To change the name of the (common block) module from denergy_det to denergy_det_m. 

- [ ] ***#33 common block /distance/***\
**PROBLEM:** \
After substituting the common block it fails to reproduce numbers of butadiene-3wfs/TZ-1M-128 test:\
total E = -26.1345117 +- 0.0114277 1.14277 0.82670 0.82670 1.91\
total E = -26.1668044 +- 0.0105012 1.05012 0.78215 0.78215 1.80\
total E = -26.1868742 +- 0.0090494 0.99132 0.74752 0.74752 1.76\
total E = -26.1825821 +- 0.0080614 0.96737 0.73212 0.73212 1.75\
and should be:\
totalE=-26.1345117\
totalE=-26.1560896\
totalE=-26.1885011\
totalE=-26.2094842\
running in 5 processors.\
After multiple checks and several trials, number differences keep appearing. I could not find the reason why. I am blaming the fact that variables r_ee, r_en, rshift, rvec_ee, rvec_en are quite spread on the code, sometimes even redefined instead of calling the /distance/ common block. I pass to the next one...\
**SOLUTION:** \
TODO

- [x] ***#35 common block /dorb/***\
**PROBLEM:** \
 Common block /orbval/ has an array with same name as dorb common block.
**SOLUTION:** \
To change the name of the (common block) module from denergy_det to dorb_m. 

- [ ] ***#54 common block /gradhess_all/***\
**PROBLEM:** \
MPARMALL parameter defined in optorb.f should be in a include file.\
**SOLUTION:** \
TODO

- [ ] ***#67 common block /icount_ci/***\
**PROBLEM:** \
Common blocks /icount_ci/, /icount_orb/ and /icount_prop/ are initialized in optwf_matrix_corsamp.f as:

      block data optprt_count

      implicit real*8(a-h,o-z)

      common /icount_ci/ icount_ci
      common /icount_orb/ icount_orb
      common /icount_prop/ icount_prop
      data icount_ci /1/
      data icount_orb /1/
      data icount_prop /1/

      end
      
However, /icount_ci/ in optci.f is used as:\
common /icount_ci/ icount 

and /icount_prop/ in properties.f as:\
common /icount_prop/ icount

Besides, icount variable is used in other subroutines like pcm.f without calling the block, and in pcm_reduce.f, giving even dimensions to the variable:
dimension icount(0:NPROCX),idispl(0:NPROCX)\
- [ ] **SOLUTION:** \
TODO

- [ ] ***#68 common block /icount_orb/***
**PROBLEM:** \
Same as #67.\
**SOLUTION:** \
TODO

- [ ] ***#69 common block /icount_prop/***\
**PROBLEM:** \
Same as #67.\
**SOLUTION:** \
TODO

- [ ] ***#88 common block /multimat/***\
**PROBLEM:** \
MEXCIT parameter must be included in a include file.\
**SOLUTION:** \
TODO

- [ ] ***#89 common block /multimatn/***\
**PROBLEM:** \
MEXCIT parameter must be included in a include file.\
**SOLUTION:** \
TODO

- [ ] ***#90 common block /multislater/***\
**PROBLEM:** \
Unknown.\
**SOLUTION:** \
TODO

- [ ] ***#91 common block /multislatern/***\
**PROBLEM:** \
It has equivalent arguments to common block /orbval/.\
**SOLUTION:** \
TODO

- [ ] ***#96 common block /numexp/***\
**PROBLEM:** \
NCOEF must be include in a include file.\
**SOLUTION:** \
TODO

- [ ] ***#104 common block /orbital_num_spl2/***
**PROBLEM:** \
Variables bc and wk are used in 3dgrid_orbitals.f with different dimensions.\
**SOLUTION:** \
TODO

- [ ] ***#105 common block /orbval/***
**PROBLEM:** \
It shares variables with common block /multislatern/, they should be condensed in just one.\
**SOLUTION:** \
TODO

- [ ] ***#118 common block /qua/***\
**PROBLEM:** \
Unclear why but tests numbers are not reproducible.\
**SOLUTION:** \
TODO

- [ ] ***#120 common block /rnyucm/***\
**PROBLEM:** \
m(4),l(4) variables are incompatible with implicit-real statement. We need to change the variable names.\
**SOLUTION:** \
TODO

- [ ] ***#125 common block /slater/***\
**PROBLEM:** \
Unknown. \
**SOLUTION:** \
TODO

- [x] ***#135 common block /vardep/***\
**PROBLEM:** \
NEQSX parameter must be include in an input file.\
**SOLUTION:** \
NEQSX and MTERMS were included in the vmc.h file.

- [x] ***#142 common block /zmatrix_grad/***\
**PROBLEM:** \
Parameter MCENT3 defined in misc_grdnts.f should be included in an include file.\
**SOLUTION:** \
To add:\
parameter (MCENT3=3*MCENT)

In the vmc.h file.
