mpirun -np 1  ../../../../bin/vmc.mov1 -i vmc_noopt_cas44.inp  -o vmc_noopt_cas44_single.out  -e error
mpirun -np 2  ../../../../bin/vmc.mov1 -i vmc_noopt_cas44.inp  -o vmc_noopt_cas44_double.out  -e error
