mpirun -np 1 ../../../bin/vmc.mov1 -i revised_vmc_corsamp.inp -o revised_vmc_corsamp_single.out -e error
mpirun -np 2 ../../../bin/vmc.mov1 -i revised_vmc_corsamp.inp -o revised_vmc_corsamp_double.out -e error
