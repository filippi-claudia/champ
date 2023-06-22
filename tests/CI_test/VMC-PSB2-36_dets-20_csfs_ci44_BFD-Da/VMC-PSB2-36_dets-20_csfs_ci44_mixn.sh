echo "PSB2 with mix_n method"

input="revised_psb2_mix_n.inp"
output="revised_psb2_mix_n"

# Unicore test
N=1
ReferenceEnergy=-29.9175911
ReferenceError=0.0385569
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError --no_assert

# Multicore test
N=2
ReferenceEnergy=-29.5422301
ReferenceError=0.036741
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError --no_assert
