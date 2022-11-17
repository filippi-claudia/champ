#Testing HDF5 backend
mpirun -np 1  ../../../bin/vmc.mov1 -i vmc_opt_ci1010_pVTZ_1522_hdf5.inp  -o vmc_opt_ci1010_pVTZ_1522_hdf5_single.out  -e error_hdf5_single

#Testing TEXT backend
mpirun -np 1  ../../../bin/vmc.mov1 -i vmc_opt_ci1010_pVTZ_1522_text.inp  -o vmc_opt_ci1010_pVTZ_1522_text_single.out  -e error_text_single
