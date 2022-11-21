echo "Butadiene DMC 500 dets "

input="dmc.inp"
output="dmc_butadiene"

# Unicore test
N=1
ReferenceEnergy=-26.2876913
ReferenceError=0.0178542
mpirun -np $N ../../../bin/dmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out "total energy ( 100) "  $ReferenceEnergy     $ReferenceError

# Multicore test
N=2
ReferenceEnergy=-26.2846795
ReferenceError=0.0077869
mpirun -np $N ../../../bin/dmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out "total energy ( 100) "  $ReferenceEnergy     $ReferenceError
