mpirun -np 1  ../../../bin/vmc.mov1 -i vmc_opt_ras1022_pVQZ_65000.inp  -o vmc_opt_ras1022_pVQZ_65000.out -e error_single
mpirun -np 2  ../../../bin/vmc.mov1 -i vmc_opt_ras1022_pVQZ_65000.inp  -o vmc_opt_ras1022_pVQZ_65000.out -e error_double
