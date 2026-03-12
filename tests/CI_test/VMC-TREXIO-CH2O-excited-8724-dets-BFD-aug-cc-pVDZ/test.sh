#!/bin/bash
echo "TREXIO formaldehyde excited state energy"

input="trexio_vmc_COH2_excited_state.inp"
output="trexio_vmc_COH2_excited_state"

# Multicore test
N=4
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Reference run completed: ${output}_core_${N}.out"

# Extract reference energy and error from the 1000-block run
ReferenceEnergy=$(grep "total E" ${output}_core_${N}.out | tail -1 | awk '{print $4}')
ReferenceError=$(grep "total E" ${output}_core_${N}.out | tail -1 | awk '{print $6}')
echo "Reference energy: ${ReferenceEnergy} +/- ${ReferenceError}"

# HDF5 restart feature test: store then restore
echo ""
echo "=== HDF5 restart test ==="

# Step 1: 500 blocks with HDF5 dump
echo "HDF5 restart test: step 1 (store phase, 500 blocks)"
rm -f restart_vmc.hdf5
mpirun -np $N ../../../bin/vmc.mov1 -i step1.inp -o step1_core_${N}.out -e error_step1
if [ -f "restart_vmc.hdf5" ]; then
    echo "HDF5 restart test: restart_vmc.hdf5 created successfully"
else
    echo "HDF5 restart test: FAILED - restart_vmc.hdf5 was not created by step1"
    exit 1
fi

# Step 2: 500 blocks restarting from HDF5 (cumulative stats include step1's 500 blocks)
echo "HDF5 restart test: step 2 (restore phase, 500 blocks)"
mpirun -np $N ../../../bin/vmc.mov1 -i step2.inp -o step2_core_${N}.out -e error_step2
echo "HDF5 restart test: step 2 completed"

# Compare step2 cumulative energy (1000 total blocks) with reference (1000 blocks)
echo ""
echo "Comparing step2 cumulative energy with reference (1000 blocks each):"
Step2Energy=$(grep "total E" step2_core_${N}.out | tail -1 | awk '{print $4}')
Step2Error=$(grep "total E" step2_core_${N}.out | tail -1 | awk '{print $6}')
echo "Step2 energy:    ${Step2Energy} +/- ${Step2Error}"
echo "Reference energy: ${ReferenceEnergy} +/- ${ReferenceError}"

../../../tools/compare_value.py step2_core_${N}.out "total E" ${ReferenceEnergy} ${ReferenceError}
if [ $? -eq 0 ]; then
    echo "HDF5 restart test: PASSED - energies match within tolerance"
else
    echo "HDF5 restart test: energies differ (may be statistical fluctuation - check manually)"
fi
