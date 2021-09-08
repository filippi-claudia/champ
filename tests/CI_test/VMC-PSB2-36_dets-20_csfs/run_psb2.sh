mpirun -np 1 ../../../bin/vmc.mov1 -i revised_psb2_mix_n.inp -o revised_psb2_mix_n_single.out -e error 
mpirun -np 2 ../../../bin/vmc.mov1 -i revised_psb2_mix_n.inp -o revised_psb2_mix_n_double.out -e error 
