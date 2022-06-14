#!/bin/bash
python trex2champ.py 	--trex 		benzene.hdf5 \
			--gamess 	benzene_hf_cc-VDZ.out \
			--motype 	"RHF" \
			--backend	"HDF5" \
			--basis_prefix  "cc-VDZ" \
			--lcao \
			--geom \
			--basis \
			--ecp \
			--sym \
			--det
