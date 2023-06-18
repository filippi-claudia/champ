echo "VMC butadiene ci44 pVQZ energy only with HDF5 restart feature"

input_dump="vmc_noopt_ci44_idump_10.inp"
input_rest="vmc_noopt_ci44_irstar_10.inp"
input_longer="vmc_noopt_ci44_longer.inp"

output_dump="vmc_noopt_ci44_idump_10"
output_rest="vmc_noopt_ci44_irstar_10"
output_longer="vmc_noopt_ci44_longer"

# HDF5 restart VMC test
ReferenceEnergyRestart=-26.2256925
ReferenceEnergyLonger=-26.2256925
ReferenceError=0.0000001


#-26.2815545 +-  0.0496094
N=1

echo "Run VMC for 10 blocks 10 steps"
mpirun -np $N ../../../bin/vmc.mov1 -i $input_dump -o ${output_dump}.out -e error

echo "Total energy after 10 blocks 10 steps"
../../../tools/compare_value.py ${output_dump}.out    "total E"  -26.2044945 0.0786157

echo "Rename the HDF5 file for restarting"
mv -v restart_vmc_*.hdf5 restart_vmc.hdf5


echo "Run VMC for another 10 blocks 10 steps"
mpirun -np $N ../../../bin/vmc.mov1 -i $input_rest -o ${output_rest}.out -e error

../../../tools/compare_value.py ${output_rest}.out    "total E"  $ReferenceEnergyRestart  $ReferenceError



echo "Run longer VMC for 20 blocks 10 steps"
mpirun -np $N ../../../bin/vmc.mov1 -i $input_longer -o ${output_longer}.out -e error

../../../tools/compare_value.py ${output_longer}.out     "total E"  $ReferenceEnergyLonger   $ReferenceError

echo "Comparing energy of restarted vs longer run "
if [ "$ReferenceEnergyRestart" == "$ReferenceEnergyLonger" ]; then
    echo "Total energies match after restart"
else
    echo "Total energies do not match after restart. Exiting!"
    exit -1
fi

echo "Removing the hdf5 files"
rm -fv *.hdf5


