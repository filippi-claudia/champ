mpirun -np 1 ../../../bin/vmc.mov1 -i revised_vmc_davidson_check.inp -o revised_vmc_davidson_check_single.out -e error 
mpirun -np 2 ../../../bin/vmc.mov1 -i revised_vmc_davidson_check.inp -o revised_vmc_davidson_check_double.out -e error 
