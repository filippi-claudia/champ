echo "H2 energy with sr_n method"

input="revised_vmc.inp"
output="revised_vmc_energy"

# Unicore test
N=1
ReferenceEnergy=-1.0035409
ReferenceError=0.0067232
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out         "total E"  $ReferenceEnergy     $ReferenceError

# Multicore test
N=2
ReferenceEnergy=-1.0097594
ReferenceError=0.0049086
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out         "total E"  $ReferenceEnergy     $ReferenceError
