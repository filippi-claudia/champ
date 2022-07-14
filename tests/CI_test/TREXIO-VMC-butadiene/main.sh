#Testing HDF5 backend
mpirun -np 1  /home/ravindra/softwares/champ-main/bin/vmc.mov1 -i vmc_optimization_500_hdf5.inp  -o vmc_optimization_500_hdf5_single_main.out  -e error
#mpirun -np 2  ../../../bin/vmc.mov1 -i vmc_optimization_500_hdf5.inp  -o vmc_optimization_500_hdf5_double.out  -e error_double_hdf5

#Testing TEXT backend
#mpirun -np 1  ../../../bin/vmc.mov1 -i vmc_optimization_500_text.inp  -o vmc_optimization_500_text_single.out  -e error_single_text
#mpirun -np 2  ../../../bin/vmc.mov1 -i vmc_optimization_500_text.inp  -o vmc_optimization_500_text_double.out  -e error_double_text
