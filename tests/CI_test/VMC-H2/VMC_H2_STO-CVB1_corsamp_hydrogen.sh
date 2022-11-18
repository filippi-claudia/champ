echo "H2 energy with corsamp method"

input="revised_vmc_corsamp.inp"
output="revised_vmc_corsamp"

# Unicore test
N=1
ReferenceEnergy=-0.9926131
ReferenceError=0.0120065
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o $output_core_$N.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py $output_core_$N.out     "total E"  $ReferenceEnergy     $ReferenceError

# Multicore test
N=2
ReferenceEnergy=-1.0135052
ReferenceError=0.0069582
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o $output_core_$N.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py $output_core_$N.out     "total E"  $ReferenceEnergy     $ReferenceError
