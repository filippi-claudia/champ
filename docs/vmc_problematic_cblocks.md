# Problematic Common blocks

Common blocks that were not deleted in the second round. 

***#3 common block /Bloc/***

**PROBLEM**: \
Variables xmatd, xmatu are shared with common block /scratch/.
\
**SOLUTION:**
CB scratch block was used only by optjas.f. Delete xmatd, xmatu from scratch module and introduce a "use Bloc, only: xmatd, xmatu" statement. 

FYI: \
common /Bloc/ \
								b(MORB,MELEC), xmat(MELEC * * 2,2), xmatu(MELEC* * 2), xmatd(MELEC * * 2) 

#3 common block /contrl/ 
**PROBLEMATIC** \
contrl shares variable nconf with common block /sr_mat_n/.

**SOLUTION:**
Variable nconf was changed to n_conf (only used in read_input.f for dmc calculations) before doing the python refactoring: \
common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
To:
common /contrl/ nstep,nblk,nblkeq,n_conf,nconf_new,isite,idump,irstar

cuspmat4,5,Todo
da_energy_ave,7,Todo
deloc_dj,6,Todo
denergy_det,3,Todo
distance,14,Todo
dorb,24,Todo
gradhess_all,4,Todo
icount_ci,2,Todo
icount_orb,1,Todo
icount_prop,2,Todo
multimat,10,Todo
multimatn,3,Todo
multislater,23,Todo
multislatern,8,Todo
numexp,2,Todo
orbital_num_spl2,1,Todo
orbval,28,Todo
qua,14,Todo
rnyucm,6,Todo
slater,17,Todo
vardep,3,Todo
zmatrix_grad,1,Todo
