
echo "H2 energy with lin_d method"

input="revised_vmc_lin.inp"
output="revised_vmc_lin"

# Unicore test
N=1
ReferenceEnergy=-1.0555878
ReferenceError=0.0032887  
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error 
echo "Comparing energy with reference Core=$N		(total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out 	"total E"  $ReferenceEnergy     $ReferenceError

# Multicore test
N=2
ReferenceEnergy=-1.0604310 
ReferenceError=0.0021645 
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error 
echo "Comparing energy with reference Core=$N		(total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out 	"total E"  $ReferenceEnergy     $ReferenceError
