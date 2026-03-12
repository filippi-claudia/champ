#!/bin/bash
echo "HDF5 VMC restart test: CH2O excited state (BFD/aug-cc-pVDZ)"
# Validates the HDF5 store/restore feature for VMC:
#   step1 runs 500 blocks and dumps the full checkpoint (walker positions,
#   accumulators, RNG state) to restart_vmc.hdf5.
#   step2 restores from that file and runs 500 more blocks.
# Because the accumulated statistics are preserved, step2 reports the
# cumulative energy over all 1000 blocks, which must match the reference.

N=4
ReferenceEnergy=-22.6180281
ReferenceError=0.0042991

# Step 1: 500 blocks with HDF5 dump
echo "HDF5 restart: step 1 (store phase, 500 blocks, ${N} MPI ranks)"
rm -f restart_vmc.hdf5
mpirun -np $N ../../../bin/vmc.mov1 -i step1.inp -o step1_core_${N}.out -e error_step1
if [ ! -f "restart_vmc.hdf5" ]; then
    echo "FAILED: restart_vmc.hdf5 was not created by step1"
    exit 1
fi
echo "HDF5 file created: restart_vmc.hdf5"

# Step 2: 500 blocks restoring from HDF5 (cumulative stats cover 1000 blocks)
echo "HDF5 restart: step 2 (restore phase, 500 blocks, ${N} MPI ranks)"
mpirun -np $N ../../../bin/vmc.mov1 -i step2.inp -o step2_core_${N}.out -e error_step2

# Compare cumulative step2 energy (1000 blocks total) with reference
echo "Comparing step2 cumulative energy with reference (total E = $ReferenceEnergy +- $ReferenceError)"
../../../tools/compare_value.py step2_core_${N}.out "total E" $ReferenceEnergy $ReferenceError
