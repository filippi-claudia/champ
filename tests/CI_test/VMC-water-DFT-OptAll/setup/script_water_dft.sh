#!/bin/bash
python3 ../../../../tools/trex_tools/trex2champ.py \
			--trex 	"../H2O_DFT.hdf5" \
			--motype 	"Canonical" \
			--backend	"HDF5" \
			--basis_prefix  "BFD-aug-cc-pVDZ" \
			--lcao \
			--geom \
			--basis \
			--ecp \
			--det

