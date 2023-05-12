echo "VMC nitroxyl cipsi 322 dets full optimization"

input="vmc.inp"
output="vmc"

# unicore test
N=1
ReferenceEnergy_1=-26.4733476
ReferenceError_1=0.0131403
ReferenceEnergy_2=-26.3137309
ReferenceError_2=0.0116616

ReferenceGap=4.3433938
ReferenceGapError=0.4779437

echo " "
echo "Core=$N"

mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error

grep "total E" ${output}_core_${N}.out | tail -2 | head -1 > temp_state1_E_err
grep "total E" ${output}_core_${N}.out | tail -1 > temp_state2_E_err
E1=$(awk '{print $4}' temp_state1_E_err)
err1=$(awk '{print $6}' temp_state1_E_err)
E2=$(awk '{print $4}' temp_state2_E_err)
err2=$(awk '{print $6}' temp_state2_E_err)
D=$(echo "$E2 - $E1" | bc)
DE=$(echo "$D * 27.2114" | bc)
D=$(echo "sqrt($err1*$err1 + $err2*$err2)" | bc)
Derr=$(echo "$D * 27.2114" | bc)
echo "total E =       $DE +- $Derr" > temp_gap_E_err

echo "Comparing state 1 energy with reference        (total E = $ReferenceEnergy_1 +-  $ReferenceError_1 Ha) "
../../../tools/compare_value.py temp_state1_E_err     "total E"  $ReferenceEnergy_1     $ReferenceError_1
echo " "
echo "Comparing state 2 energy with reference        (total E = $ReferenceEnergy_2 +-  $ReferenceError_2 Ha) "
../../../tools/compare_value.py temp_state2_E_err     "total E"  $ReferenceEnergy_2     $ReferenceError_2    
echo " "
echo "Comparing excitation energy with reference        (Delta E = $ReferenceGap +-  $ReferenceGapError eV) "
echo "Your excitation energy and error                  (Delta E = $DE +-  $Derr eV)"
../../../tools/compare_value.py temp_gap_E_err     "total E"  $ReferenceGap     $ReferenceGapError    

rm -f temp_state1_E_err temp_state2_E_err temp_gap_E_err

# Multicore test
N=2
ReferenceEnergy_1=-26.4975023
ReferenceError_1=0.0081915
ReferenceEnergy_2=-26.3449000
ReferenceError_2=0.0084472

ReferenceGap=4.1525222
ReferenceGapError=0.3201230

echo " "
echo "Core=$N"

mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error

grep "total E" ${output}_core_${N}.out | tail -2 | head -1 > temp_state1_E_err
grep "total E" ${output}_core_${N}.out | tail -1 > temp_state2_E_err
E1=$(awk '{print $4}' temp_state1_E_err)
err1=$(awk '{print $6}' temp_state1_E_err)
E2=$(awk '{print $4}' temp_state2_E_err)
err2=$(awk '{print $6}' temp_state2_E_err)
D=$(echo "$E2 - $E1" | bc)
DE=$(echo "$D * 27.2114" | bc)
D=$(echo "sqrt($err1*$err1 + $err2*$err2)" | bc)
Derr=$(echo "$D * 27.2114" | bc)
echo "total E =       $DE +- $Derr" > temp_gap_E_err

echo "Comparing state 1 energy with reference        (total E = $ReferenceEnergy_1 +-  $ReferenceError_1 Ha) "
../../../tools/compare_value.py temp_state1_E_err     "total E"  $ReferenceEnergy_1     $ReferenceError_1
echo " "
echo "Comparing state 2 energy with reference        (total E = $ReferenceEnergy_2 +-  $ReferenceError_2 Ha) "
../../../tools/compare_value.py temp_state2_E_err     "total E"  $ReferenceEnergy_2     $ReferenceError_2
echo " "
echo "Comparing excitation energy with reference        (Delta E = $ReferenceGap +-  $ReferenceGapError eV) "
echo "Your excitation energy and error                  (Delta E = $DE +-  $Derr eV)"
../../../tools/compare_value.py temp_gap_E_err     "total E"  $ReferenceGap     $ReferenceGapError    
rm -f temp_state1_E_err temp_state2_E_err temp_gap_E_err

