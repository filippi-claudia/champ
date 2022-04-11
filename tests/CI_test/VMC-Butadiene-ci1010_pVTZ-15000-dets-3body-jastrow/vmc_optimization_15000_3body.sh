mpirun -np 1  ../../../bin/vmc.mov1 -i vmc_optimization_15000_3body.inp  -o vmc_optimization_15000_3body_single.out  -e error_vmc_optimization_15000_3body
mpirun -np 2  ../../../bin/vmc.mov1 -i vmc_optimization_15000_3body.inp  -o vmc_optimization_15000_3body_double.out  -e error_vmc_optimization_15000_3body
