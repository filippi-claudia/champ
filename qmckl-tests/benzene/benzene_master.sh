#python3 -m pip install .
#pip install .
trexio convert-from --type gamess --input hf_cc-VDZ.out --motype="RHF" gamess_benzene_hf_cc-VDZ.hdf5 --back_end=HDF5
trexio convert2champ gamess_benzene_hf_cc-VDZ.hdf5 --input hf_cc-VDZ.out --motype="RHF" --back_end=HDF5
