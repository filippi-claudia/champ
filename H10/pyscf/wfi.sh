#!/bin/bash

#activate last trextools environment

source ~/Codes/venvtrex/bin/activate

python /home/landinez/Codes/trexio_tools/src/trexio_tools/converters/pyscf_to_trexio.py -c RHF.dump -o RHF_sph.hdf5

#since k-points are regarded the name of the textio file change a bit by the converter by default (not sure how to enforce the name of the output explivitely)

mv k0_RHF_sph.hdf5 RHF_sph.hdf5


# second step trexio spherical to cartesian 
trexio convert-to -t cartesian -o RHF_cart.hdf5 RHF_sph.hdf5

    ##to hack trex2champ
    # create a false gamess input file
touch fake.out


#third step from trexio to champ (regarding local converter inside current directory)

#python trex2champ.py --trex RHF_cart.hdf5 --gamess fake.out --motype RHF --backend hdf5 --lcao  --geometry  --basis  --pseudo   --symmetry   --determinants

python trex2champ.py --trex RHF_cart.hdf5 --gamess fake.out --motype RHF --backend hdf5 --lcao  --geometry  --basis  --pseudo --determinants




deactivate 
