echo "TREXIO backend comparison: HDF5"
input="vmc_opt_ci1010_pVTZ_1522_hdf5.inp"
output="vmc_opt_ci1010_pVTZ_1522_hdf5"

# unicore test
N=1
mpirun -np $N ../../../bin/vmc.mov1 -i ${input} -o ${output}_core_${N}.out -e error



