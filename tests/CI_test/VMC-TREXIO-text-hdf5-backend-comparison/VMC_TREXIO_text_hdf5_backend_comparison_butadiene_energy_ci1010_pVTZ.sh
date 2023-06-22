echo "TREXIO backend comparison: HDF5"
input="vmc_opt_ci1010_pVTZ_1522_hdf5.inp"
output="vmc_opt_ci1010_pVTZ_1522_hdf5"

# unicore test
N=1
ReferenceEnergyhdf5=-25.8360399
ReferenceError=0.2856974
mpirun -np $N ../../../bin/vmc.mov1 -i ${input} -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=${N}           (total E = $ReferenceEnergyhdf5 +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergyhdf5     $ReferenceError

energy_hdf5=$(grep 'total E' ${output}_core_${N}.out | awk '{print $4}')



echo "TREXIO backend comparison : TEXT"
input="vmc_opt_ci1010_pVTZ_1522_text.inp"
output="vmc_opt_ci1010_pVTZ_1522_text"

# unicore test
N=1
ReferenceEnergytext=-25.8360399
ReferenceError=0.2856974
mpirun -np $N ../../../bin/vmc.mov1 -i ${input} -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=${N}           (total E = $ReferenceEnergytext +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergytext     $ReferenceError

energy_text=$(grep 'total E' ${output}_core_${N}.out | awk '{print $4}')

echo "Comparing energy of hdf5 vs text run "
if [ "$energy_hdf5" = "$energy_text" ]; then
    echo "Total energies match in both backends "
else
    echo "Total energies do not match in both backends. Exiting!"
    exit -1
fi


