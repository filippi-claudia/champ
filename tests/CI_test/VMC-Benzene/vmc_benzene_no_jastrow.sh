mpirun -np 1 ../../../bin/vmc.mov1  -i  vmc_benzene_hf_no_jastrow.inp -o vmc_benzene_hf_no_jastrow_core1_test.out  -e error
mpirun -np 8 ../../../bin/vmc.mov1  -i  vmc_benzene_hf_no_jastrow.inp -o vmc_benzene_hf_no_jastrow_core8_test.out  -e error
