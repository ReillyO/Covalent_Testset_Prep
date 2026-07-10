#!/bin/bash

# gets RMSDs for all runs of a given experiment and writes to analysis directory

# REQUIRES ARGUMENT of experiment name - eg 004.poserep

dockdir=${DOCKHOME}
amberdir=${AMBERHOME}
rootdir=${ROOTDIR}
chimeradir=${CHIMERAHOME}
testsetdir=${rootdir}/zzz.testset_files/
rawfiledir=${rootdir}/zzz.master/
scriptdir=${rootdir}/zzz.scripts/
paramdir=${rootdir}/zzz.parameters/
analysisdir=${rootdir}/zzz.analysis

prefix=${1}
outfile=${analysisdir}/${prefix}_rmsds.csv

if [ ! -e ${analysisdir} ]; then mkdir ${analysisdir}; fi

if [ -e ${outfile} ]; then mv ${outfile} ${outfile//_rmsds/_old_rmsds}; fi

for f in `ls -d ${testsetdir}/*`; do
  system=${f: -4}
  sysdir=${testsetdir}/${system}
  
  cd ${testsetdir}/${system}/${prefix}
  
  if grep -qi error ${sysdir}/${prefix}/cov.out; then 
    echo "DOCK encountered an error on ${system}! Check cov.out for error details."
    echo "${system},1000" >> ${outfile}
  else
    echo "${system}"`grep "HA_RMSDh" ${system}_out_scored.mol2` >> ${outfile}
  fi
  sed -i "s/ $//g" ${outfile}
  sed -i "s/ ########## HA_RMSDh: /,/g" ${outfile}
  sed -i "s/########## HA_RMSDh: /,/g" ${outfile}
  sed -i "s/ //" ${outfile}

done

