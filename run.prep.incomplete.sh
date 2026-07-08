#!/bin/bash

rootdir=${ROOTDIR}
sysdir=${rootdir}/zzz.testset_files/
rawdir=${rootdir}/zzz.master/

# loop through the unprocessed receptor files in the zzz.master dir and
# if the corresponding processed directory grid does not exist then run
# prep on the system
for rec in `ls ${rawdir}/????.rec.noch.pdb`; do
	sys=${rec: -17}	# assumes standard name ????.rec.noch.pdb where ???? is PDB code
	sys=${sys:0:4}
	if [ ! -s ${sysdir}/${sys}/003.gridsph/*bmp ]; then
		bash run.prep.allsteps.sh -O -s ${sys} &
	fi
done

wait
