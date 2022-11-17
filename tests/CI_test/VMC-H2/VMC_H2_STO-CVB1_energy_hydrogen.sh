mpirun -np 1 ../../../bin/vmc.mov1 -i revised_vmc.inp -o revised_vmc_single.out -e error 
mpirun -np 2 ../../../bin/vmc.mov1 -i revised_vmc.inp -o revised_vmc_double.out -e error 
