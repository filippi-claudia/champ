echo "PSB3 energy with lin free method"

input="revised_free.inp"
output="revised_free"

# Unicore test
N=1
ReferenceEnergy=-40.8932882
ReferenceError=0.0505664
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o $output_core_$N.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py $output_core_$N.out     "total E"  $ReferenceEnergy     $ReferenceError

# Multicore test
N=2
ReferenceEnergy=-41.1767726
ReferenceError=0.0289563
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o $output_core_$N.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py $output_core_$N.out     "total E"  $ReferenceEnergy     $ReferenceError
