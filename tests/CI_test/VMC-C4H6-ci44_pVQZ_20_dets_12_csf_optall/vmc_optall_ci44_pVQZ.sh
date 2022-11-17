mpirun -np 1  ../../../bin/vmc.mov1 -i vmc_optall_ci44.inp  -o vmc_optall_ci44_single.out  -e error
mpirun -np 2  ../../../bin/vmc.mov1 -i vmc_optall_ci44.inp  -o vmc_optall_ci44_double.out  -e error
