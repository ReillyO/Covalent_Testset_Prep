#!/bin/bash

# will check every system present in the master rec/lig directory and see
# if preparatory steps 1-3 have been run on them; if a system has not been
# suitably prepared, the step that needs to be run will be printed

step=${1}
ncomplete=0
nincomplete=0

if [ -e "complete.txt" ]; then rm "complete.txt"; fi

for f in zzz.master/*.pdb; do 
  f=${f//zzz\.master\//}
  f=${f//\.rec\.noch\.pdb/}
  if [ ! -e zzz.testset_files/${f} ]; then
    echo $f 000
    nincomplete=$((nincomplete+1))
  elif [ ! -e zzz.testset_files/${f}/001.ligprep/03.${f}.lig.noch.pdb ]; then
    echo $f 001
    nincomplete=$((nincomplete+1))
  elif [ ! -e zzz.testset_files/${f}/002.recprep/05.${f}.lig.clean.mol2 ]; then
    echo $f 002
    nincomplete=$((nincomplete+1))
  elif [ ! -e zzz.testset_files/${f}/003.gridsph/${f}.rec.bmp ]; then
    echo $f 003
    nincomplete=$((nincomplete+1))
  else
    echo $f >> "complete.txt"
    ncomplete=$((ncomplete+1))
  fi
done

echo "${ncomplete} complete"
echo "${nincomplete} incomplete"
