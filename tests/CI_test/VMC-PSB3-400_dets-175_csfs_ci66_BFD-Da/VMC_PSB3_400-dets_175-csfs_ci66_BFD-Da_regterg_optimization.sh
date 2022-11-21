echo "PSB3 energy with regterg method"

input="revised_regterg.inp"
output="revised_regterg"

# Unicore test
N=1
ReferenceEnergy=-41.4894680
ReferenceError=0.0432582
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError --no_assert

# Multicore test
N=2
ReferenceEnergy=-42.1632583
ReferenceError=0.0282363
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError --no_assert
