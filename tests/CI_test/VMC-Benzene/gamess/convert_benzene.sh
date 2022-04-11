trexio convert-from --type gamess --input hf_cc-VDZ.out --motype "RHF" gamess_benzene.hdf5 --back_end=HDF5
trexio convert2champ gamess_benzene.hdf5  --input hf_cc-VDZ.out --motype "RHF" --back_end=HDF5
