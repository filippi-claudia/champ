echo "H2 energy with corsamp method"

input="revised_vmc_corsamp.inp"
output="revised_vmc_corsamp"

# Unicore test
N=1
ReferenceEnergy=-0.9923642
ReferenceError=0.0102988
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out         "total E"  $ReferenceEnergy     $ReferenceError

# Multicore test
N=2
ReferenceEnergy=-1.0005539
ReferenceError=0.007121
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out         "total E"  $ReferenceEnergy     $ReferenceError
