mpirun -np 1  ../../../../bin/vmc.mov1 -i vmc_optimization_5000.inp  -o vmc_optimization_5000_single.out  -e error_vmc_optimization_5000
mpirun -np 2  ../../../../bin/vmc.mov1 -i vmc_optimization_5000.inp  -o vmc_optimization_5000_double.out  -e error_vmc_optimization_5000
