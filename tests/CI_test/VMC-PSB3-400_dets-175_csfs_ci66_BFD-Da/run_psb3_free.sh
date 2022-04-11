mpirun -np 1 ../../../bin/vmc.mov1 -i revised_free.inp -o revised_free_single.out -e error 
mpirun -np 2 ../../../bin/vmc.mov1 -i revised_free.inp -o revised_free_double.out -e error 
