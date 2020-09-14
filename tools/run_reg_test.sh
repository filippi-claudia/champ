#! /bin/bash

CHAMP=$HOME/champ 

cd $CHAMP
rm -rf build bin
module load pre2019
module load cmake/3.7.2 lapack/mkl blas/mkl intel/2018b mpi/impi/18.0.4

# compile the current branch
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifort
cmake --build build  --target vmc.mov1 -- -j4 # vmc.mov1 -- -j4
mv ./bin ./bin_cur

# compile the ref branch
git checkout $1
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifort
cmake --build build  --target vmc.mov1 -- -j4 # vmc.mov1 -- -j4
mv ./bin ./bin_ref

# restore the original 
cp -r ./bin_cur ./bin

# go to the test folder
cd $CHAMP/tests/psb3_cas66 ; \
cp $CHAMP/tools/test_branch_champ.sh ./ ; \

########### Test 1
echo "" ; \
echo "Running test 1: free.inp" ; \
sbatch test_branch_champ.sh free.inp $CHAMP/bin_cur $CHAMP/bin_ref; sleep 40 ; \
echo "current branch" ; \
grep "eigenvalue            1  :" free.cur.out ; \
echo "reference branch" ; \
grep "eigenvalue            1  :" free.ref.out ; \
echo "" ; \

########### Test 2
echo "Running test 2: regterg.inp" ; \
sbatch test_branch_champ.sh regterg.inp $CHAMP/bin_cur $CHAMP/bin_ref; sleep 20 ; \
echo "current branch" ; \
grep "LIN_D: state, vec, energy   1   1" regterg.cur.out ; \
grep "LIN_D: state, vec, energy   1   1" regterg.ref.out ; \


########### Test 3
echo "Running test 3:" ; \
sbatch test_branch_champ.sh sr_n.inp $CHAMP/bin_cur $CHAMP/bin_ref; sleep 20 ; \
echo "current branch" ; \
grep "delta0" sr_n.cur.out ; \
echo "reference branch" ; \
grep "delta0" sr_n.ref.out ; \


########### Test 4
echo "Running test 4:" ; \
cd $CHAMP/tests/butadiene-3wfs/TZ-1M-128 ; \
cp $CHAMP/tools/test_branch_champ.sh ./ ; \
sbatch test_branch_champ.sh vmc_optall.inp $CHAMP/bin_cur $CHAMP/bin_ref; sleep 180 ; \
echo "current branch" ; \
grep "total E =" vmc_optall.cur.out | cut -f 1-11 -d  " " ; \
echo "reference branch" ; \
grep "total E =" vmc_optall.ref.out | cut -f 1-11 -d  " " ; \


