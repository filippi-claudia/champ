#!/bin/bash
python3 trex2champ.py \
	--trex 	"COH2_GS.trexio" \
	--backend	"HDF5" \
	--basis_prefix  "BFD-aug-cc-pVDZ" \
	--lcao \
	--ecp \
	--sym \
	--geom \
	--basis \
	--det

