echo "VMC nitroxyl cipsi 322 dets full optimization"

input="vmc.inp"
output="vmc"

# Multicore test
N=2
ReferenceEnergy_1=-26.4970628
ReferenceError_1=0.0083161
ReferenceEnergy_2=-26.3313189
ReferenceError_2=0.0093368

ReferenceGap=4.5101235
ReferenceGapError=0.3400853

echo " "
echo "Core=$N"

mpirun -np $N ../../../bin/vmc.mov1 -i $input 
