#!/bin/bash

# prepare a given system (4-character PDB code) for covalent docking
#
# Options:
#  -s: REQ'D, 4-character PDB code for system
#  -O: OPT, if passed will overwrite any data on the system
#
# assumes:
#  - ligand file called ${system}.lig.ch.mol2 exists in zzz.master/ directory 
#     protonated and charged with Gasteiger charges (via Chimera)
#     and possessing CYS attachment in addition to ligand
#  - receptor file called ${system}.rec.noch.pdb exists in zzz.master/
#     unprotonated and uncharged
#     and modified such that the covalent residue is a GLY instead of CYS or whatever else
#  - 000.setup.sh has been run already and all modules loaded/variables set successfully

# See individual scripts for documentation on actions taken
# Output from this script should be directed to a verbose file for analysis
# in the event of errors

# sample use: bash run.prep.allsteps.sh -s 1A2B >& 1A2B.v.txt &

dockdir=${DOCKHOME}
amberdir=${AMBERHOME}
chimeradir=${CHIMERAHOME}
rootdir=${ROOTDIR}
testsetdir=${rootdir}/zzz.testset_files/
rawfiledir=${rootdir}/zzz.master/
scriptdir=${rootdir}/zzz.scripts/
system=${1}

overwrite=false;
system="UNSET"

while getopts 'Os:' flag; do
  case "${flag}" in
    O) overwrite=true ;;
    s) system="${OPTARG}";;
    *) print_usage
       exit 1 ;;
  esac
done

if [ ${system} = "UNSET" ]; then echo "Please specify a system with -s! Exiting. "; exit 1; fi

if [ -e ${testsetdir}/${system} ]; then 
	if [ ${overwrite} = true ]; then 
		rm -r ${testsetdir}/${system}; 
	else
		echo "System data already exists! Specify -O to overwrite and try again."
		exit 1
	fi
fi
mkdir ${testsetdir}/${system}

# Ligand preparation
bash ${rootdir}/001.preplig.sh ${system} >& ${system}.ligprep.v.txt && mv ${system}.ligprep.v.txt ${testsetdir}/${system}/001.ligprep/

# Protein preparation
bash ${rootdir}/002.preprec.sh ${system} >& ${system}.recprep.v.txt && mv ${system}.recprep.v.txt ${testsetdir}/${system}/002.recprep/

# Grid and sphere prep
bash ${rootdir}/003.gridsph.sh ${system} >& ${system}.gridsph.v.txt && mv ${system}.gridsph.v.txt ${testsetdir}/${system}/003.gridsph/
