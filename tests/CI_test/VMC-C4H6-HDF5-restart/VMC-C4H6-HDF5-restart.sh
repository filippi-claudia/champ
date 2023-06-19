echo "VMC butadiene ci44 pVQZ energy only with HDF5 restart feature"

input_dump="vmc_noopt_ci44_idump_10.inp"
input_rest="vmc_noopt_ci44_irstar_10.inp"
input_longer="vmc_noopt_ci44_longer.inp"

output_dump="vmc_noopt_ci44_idump_10"
output_rest="vmc_noopt_ci44_irstar_10"
output_longer="vmc_noopt_ci44_longer"

# HDF5 restart VMC test
ReferenceEnergyLonger=-26.1788661
ReferenceError=0.0487149

N=1

echo "Run VMC for 10 blocks 10 steps"
mpirun -np $N ../../../bin/vmc.mov1 -i $input_dump -o ${output_dump}.out -e error

echo "Total energy after 10 blocks 10 steps"
../../../tools/compare_value.py ${output_dump}.out    "total E"  -26.2044945 0.0786157

echo "Rename the HDF5 file for restarting"
mv -v restart_vmc_*.hdf5 restart_vmc.hdf5

echo "Run VMC for another 10 blocks 10 steps"
mpirun -np $N ../../../bin/vmc.mov1 -i $input_rest -o ${output_rest}.out -e error

after_restart=$(grep 'total E' ${output_rest}.out | awk '{print $4}')
echo $after_restart

echo "Run longer VMC for 20 blocks 10 steps"
mpirun -np $N ../../../bin/vmc.mov1 -i $input_longer -o ${output_longer}.out -e error

after_longer=$(grep 'total E' ${output_longer}.out | awk '{print $4}')
echo $after_longer

echo "Compare the energies obtained via restart calculation with binary files"
../../../tools/compare_value.py ${output_longer}.out  "total E"  $ReferenceEnergyLonger   $ReferenceError

echo "Comparing energy of restarted vs longer run "
if [ "$after_restart" = "$after_longer" ]; then
    echo "Total energies match after restart"
else
    echo "Total energies do not match after restart. Exiting!"
    exit -1
fi

echo "Removing the hdf5 files"
rm -fv *.hdf5


