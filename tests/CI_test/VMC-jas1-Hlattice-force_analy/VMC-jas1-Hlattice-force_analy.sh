echo "periodic H lattice Jastrow 1"

input="vmc.inp"
output="vmc.out"

# Multicore test
N=2
ReferenceEnergy=-4.2118242
ReferenceError=0.0073473
mpirun -np $N ../../../bin/vmc.mov1 -i ${input} -o ${output} -e error
echo "Comparing energy with reference Core=${N}           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}     "total E"  $ReferenceEnergy     $ReferenceError
echo "Comparing analytical forces with reference Core=${N} using reference force file "
../../../tools/compare_forces.py force_analytic force_reference
