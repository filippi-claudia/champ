#mpirun -np 1  ../../../bin/vmc.mov1 -i vmc_optimization_500.inp  -o vmc_optimization_500_single.out  -e error_vmc_optimization_500
mpirun -np 2  ../../../bin/vmc.mov1 -i vmc_optimization_500.inp  -o vmc_optimization_500_double.out  -e error_vmc_optimization_500
#mpirun -np 1  ../../../bin/vmc.mov1 -i vmc_no_optimization_500.inp  -o vmc_no_optimization_500_single.out  -e error_vmc_no_optimization_500

