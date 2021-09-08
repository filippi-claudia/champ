mpirun -np 1 ../../../bin/vmc.mov1 -i revised_vmc_lin.inp -o revised_vmc_lin_single.out -e error 
mpirun -np 2 ../../../bin/vmc.mov1 -i revised_vmc_lin.inp -o revised_vmc_lin_double.out -e error 
