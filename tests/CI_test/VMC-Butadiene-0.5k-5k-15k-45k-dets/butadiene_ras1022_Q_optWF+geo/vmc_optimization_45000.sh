mpirun -np 1  ../../../../bin/vmc.mov1 -i vmc_optimization_45000.inp  -o vmc_optimization_45000_single.out -e error 
mpirun -np 2  ../../../../bin/vmc.mov1 -i vmc_optimization_45000.inp  -o vmc_optimization_45000_double.out -e error
