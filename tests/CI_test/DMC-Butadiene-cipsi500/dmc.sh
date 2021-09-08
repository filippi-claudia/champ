mpirun -np 1 ../../../../bin/dmc.mov1 -i dmc.inp  -o dmc_workload_single.out  -e error_dmc_workload  
mpirun -np 2 ../../../../bin/dmc.mov1 -i dmc.inp  -o dmc_workload_double.out  -e error_dmc_workload  
