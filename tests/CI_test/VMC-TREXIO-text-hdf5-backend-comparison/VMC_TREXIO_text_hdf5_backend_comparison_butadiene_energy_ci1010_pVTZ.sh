echo "TREXIO backend comparison: HDF5"
input="vmc_opt_ci1010_pVTZ_1522_hdf5.inp"
output="vmc_opt_ci1010_pVTZ_1522_hdf5"

# unicore test
N=1
ReferenceEnergyhdf5=-24.1199689
ReferenceError=0.0545703
mpirun -np $N ../../../bin/vmc.mov1 -i ${input} -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=${N}           (total E = $ReferenceEnergyhdf5 +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergyhdf5     $ReferenceError


echo "TREXIO backend comparison : TEXT"
input="vmc_opt_ci1010_pVTZ_1522_text.inp"
output="vmc_opt_ci1010_pVTZ_1522_text"

# unicore test
N=1
ReferenceEnergytext=-24.1199689
ReferenceError=0.0545703
mpirun -np $N ../../../bin/vmc.mov1 -i ${input} -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=${N}           (total E = $ReferenceEnergytext +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergytext     $ReferenceError

