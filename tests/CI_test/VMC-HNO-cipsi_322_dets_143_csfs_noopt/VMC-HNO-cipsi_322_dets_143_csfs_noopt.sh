echo "VMC nitroxyl cipsi 322 dets energy only"

input="vmc.inp"
output="vmc"

# unicore test
N=1
ReferenceEnergy_1=-26.4897187
ReferenceError_1=0.0029363
ReferenceEnergy_2=-26.3218125
ReferenceError_2=0.0033607

mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error

grep "total E" ${output}_core_${N}.out | head -1 > temp_state1_E_err
grep "total E" ${output}_core_${N}.out | tail -1 > temp_state2_E_err
echo " "
echo "Core=$N"
echo "Comparing state 1 energy with reference        (total E = $ReferenceEnergy_1 +-  $ReferenceError_1 ) "
../../../tools/compare_value.py temp_state1_E_err     "total E"  $ReferenceEnergy_1     $ReferenceError_1


#if [ $? -ne 0 ]; then
#    echo "Assertion Error"
#    exit -1
#fi


echo " "
echo "Comparing state 2 energy with reference        (total E = $ReferenceEnergy_2 +-  $ReferenceError_2 ) "
../../../tools/compare_value.py temp_state2_E_err     "total E"  $ReferenceEnergy_2     $ReferenceError_2

#if [ $? -ne 0 ]; then
#    echo "Assertion Error"
#    exit -1
#fi

rm -f temp_state1_E_err temp_state2_E_err

# Multicore test
N=2
ReferenceEnergy_1=-26.4897228
ReferenceError_1=0.00208
ReferenceEnergy_2=-26.3247255
ReferenceError_2=0.0022513

mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error

grep "total E" ${output}_core_${N}.out | head -1 > temp_state1_E_err
grep "total E" ${output}_core_${N}.out | tail -1 > temp_state2_E_err
echo " "
echo "Core=$N"
echo "Comparing state 1 energy with reference        (total E = $ReferenceEnergy_1 +-  $ReferenceError_1 ) "
../../../tools/compare_value.py temp_state1_E_err     "total E"  $ReferenceEnergy_1     $ReferenceError_1

#if [ $? -ne 0 ]; then
#    echo "Assertion Error"
#    exit -1
#fi

echo " "
echo "Comparing state 2 energy with reference        (total E = $ReferenceEnergy_2 +-  $ReferenceError_2 ) "
../../../tools/compare_value.py temp_state2_E_err     "total E"  $ReferenceEnergy_2     $ReferenceError_2
#if [ $? -ne 0 ]; then
#    echo "Assertion Error"
#    exit -1
#fi

rm -f temp_state1_E_err temp_state2_E_err
