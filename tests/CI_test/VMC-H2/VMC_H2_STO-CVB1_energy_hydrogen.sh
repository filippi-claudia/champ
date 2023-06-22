echo "H2 energy with sr_n method"

input="revised_vmc.inp"
output="revised_vmc_energy"

# Unicore test
N=1
ReferenceEnergy=-1.0046458
ReferenceError=0.0084583
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out         "total E"  $ReferenceEnergy     $ReferenceError

# Multicore test
N=2
ReferenceEnergy=-1.0090406
ReferenceError=0.0054465
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out         "total E"  $ReferenceEnergy     $ReferenceError
