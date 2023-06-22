echo "Butadiene DMC 500 dets "

input="dmc.inp"
output="dmc_butadiene"

# Unicore test
N=1
ReferenceEnergy=-26.2980283 
ReferenceError=0.0222300
mpirun -np $N ../../../bin/dmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out "total energy ( 100) "  $ReferenceEnergy     $ReferenceError

# Multicore test
N=2
ReferenceEnergy=-26.2891738
ReferenceError=0.0131776
mpirun -np $N ../../../bin/dmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out "total energy ( 100) "  $ReferenceEnergy     $ReferenceError
