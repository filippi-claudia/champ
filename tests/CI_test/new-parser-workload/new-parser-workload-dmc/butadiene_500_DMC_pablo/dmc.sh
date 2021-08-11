mpirun -np 1 ../../../../../bin/dmc.mov1 -i dmc.inp  -o dmc_single.out  -e error_dmc_single
mpirun -np 2 ../../../../../bin/dmc.mov1 -i dmc.inp  -o dmc_double.out  -e error_dmc_double
