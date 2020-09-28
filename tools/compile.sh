#! /bin/bash
rm -rf build bin
module load pre2019
module load cmake/3.7.2 lapack/mkl blas/mkl intel/2018b mpi/impi/18.0.4
/usr/bin/cp ./tests/psb3_cas66/include/* ./src/include ; \
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifort
cmake --build build  --target all -- -j4 # vmc.mov1 -- -j4

CHAMP=$HOME/champ

cd $CHAMP/tests/psb3_cas66 ; \
cp $CHAMP/tools/test_champ.sh ./ ; \

########### Test 1
echo "" ; \
echo "Running test 1: free.inp" ; \
sbatch test_champ.sh free.inp ; sleep 40 ; \
grep "eigenvalue            1  :" free.out ; \
echo "That should be:" ; \
echo " eigenvalue            1  :   -43.4420580939622" ; \
echo "" ; \

########### Test 2
echo "Running test 2: regterg.inp" ; \
sbatch test_champ.sh regterg.inp ; sleep 20 ; \
grep "LIN_D: state, vec, energy   1   1" regterg.out ; \
echo "That should be:" ; \
echo "LIN_D: state, vec, energy   1   1   -43.51855     0.10377" ; \
echo "" ; \

########### Test 3
echo "Running test 3:" ; \
sbatch test_champ.sh sr_n.inp ; sleep 20 ; \
grep "delta0" sr_n.out ; \
echo "That should be:" ; \
echo "delta0 =   0.888765908731205" ; \
echo "" ; \

########### Test 4
echo "Running test 4:" ; \
cd $CHAMP/tests/butadiene-3wfs/TZ-1M-128 ; \
cp $CHAMP/tools/test_champ.sh ./ ; \
sbatch test_champ.cmd vmc_optall.inp ; sleep 180 ; \
grep "total E =" vmc_optall.out | cut -f 1-11 -d  " " ; \
echo "That should be:" ; \
echo "totalE=-26.1345117" ; \
echo "totalE=-26.1560896" ; \
echo "totalE=-26.1885011" ; \
echo "totalE=-26.2094842" 
echo "" ; \

########### Test 5
echo "Running test 5:" ; \
cd $CHAMP/tests/butadiene_sr; \
cp $CHAMP/tools/test_champ.sh ./ ; \
cat vmc_sr.inp | sed s/"iforce_analy 0"/"iforce_analy 1"/g > tmp
rm -rf vmc_sr.inp ; mv tmp vmc_sr.inp ;
sbatch test_champ.cmd vmc_sr.inp ; sleep 180 ;
echo "Last geometry is:"
grep -10 "Norm of parm variation" vmc_sr.out | tail -10
echo "That won't concide with but should be something like:"
echo "CENT   -1.30663316596991      -0.506629269169022      -2.976475069589144E-002"
echo "CENT    1.46036019041281       0.594626857099687       0.106276790635788"
echo "CENT   -3.06356735543629       0.189040601049533       3.238462977910284E-002"
echo "CENT    3.06009819332530      -0.185735861220545       8.854917614288282E-002"
echo "CENT   -3.56774277382982        2.14566042209734      -7.624619640162496E-003"
echo "CENT    3.53353534662322       -2.10055754828399       5.140113108229803E-002"
echo "CENT   -1.01913133794436       -2.65794762655668      -4.069822941517589E-003"
echo "CENT    1.01306705446793        2.71810460820081       2.477132418401577E-002"
echo "CENT   -4.88478665443138      -0.712524793620838      -4.055343477810950E-003"
echo "CENT    4.80507470363003       0.773553855145416       5.887308044885128E-002"
echo "" ;
