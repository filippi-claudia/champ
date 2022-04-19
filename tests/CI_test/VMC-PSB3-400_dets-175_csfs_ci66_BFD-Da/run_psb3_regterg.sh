mpirun -np 1 ../../../bin/vmc.mov1 -i revised_regterg.inp -o revised_regterg_single.out -e error 
mpirun -np 2 ../../../bin/vmc.mov1 -i revised_regterg.inp -o revised_regterg_double.out -e error 
