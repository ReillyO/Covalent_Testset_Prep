#!/bin/bash

# Argument: 4-character PDB code

# Uses code to find a stripped receptor PDB file with ligand removed and covalent
# residue replaced with GLY (via Chimera swapaa)

# Creates a surface file for the unprotonated receptor,
# Charges the receptor using Chimera (Gasteiger),
# Protonates also using Chimera
# and subsequently charges again with Antechamber (AM1BCC)


# Assumes the ligand file is in zzz.master/ 

dockdir=${DOCKHOME}
amberbin=${AMBERBIN}
rootdir=${ROOTDIR}
chimerabin=${CHIMERABIN}
testsetdir=${rootdir}/zzz.testset_files/
rawfiledir=${rootdir}/zzz.master/
scriptdir=${rootdir}/zzz.scripts/
paramdir=${rootdir}/zzz.parameters/
system=${1}
has_cofactor=false

echo "Setting up receptor prep directory..."

cd ${testsetdir}
if [ ! -e ./${system}/ ]; then
  mkdir ./${system}
fi
cd ./${system}/
if [ ! -e ./002.recprep/ ]; then
  mkdir 002.recprep
fi
cd ${testsetdir}/${system}/002.recprep


# check for cofactor and update variable accordingly
if [ -e ${rawfiledir}/${system}.cof.ch.mol2 ]; then
  has_cofactor=true
fi

# Check for and copy in ligand
if [ ! -e ${testsetdir}/${system}/001.ligprep/03.${system}.lig.am1bcc.mol2 ]; then
	echo "Ligand file 03.${system}.lig.am1bcc.mol2 could not be found! Please ensure ligand prep occurred."
	exit 1
else
	cp ${testsetdir}/${system}/001.ligprep/03.${system}.lig.am1bcc.mol2 ./00.${system}.lig.am1bcc.mol2
fi

# Check to see that if a cofactor exists, a corresponding am1bcc.mol2 also exists
if [ $has_cofactor = true ]; then
  if  [ ! -e ${testsetdir}/${system}/001.ligprep/03.${system}.cof.am1bcc.mol2 ] ; then
    echo "Cofactor expected for ${system} could not be found. Exiting."
    exit 1
  else
    cp ${testsetdir}/${system}/001.ligprep/03.${system}.cof.am1bcc.mol2 ./00.${system}.cof.am1bcc.mol2
  fi
fi

# check for and copy in receptor file
if [ ! -e ${rawfiledir}/${system}.rec.noch.pdb ]; then
  echo "Cannot find receptor file for ${system}. Please ensure it is in zzz.master/ and named SYSTEM.rec.noch.pdb where SYSTEM is the 4-character PDB code that the receptor belongs to."
  exit 1
else
  cp ${rawfiledir}/${system}.rec.noch.pdb ./00.${system}.rec.noch.pdb
fi


# copy over LEAP parameter files
cp -fr ${paramdir}/ions.frcmod ./
cp -fr ${paramdir}/ions.lib ./
cp -fr ${paramdir}/gaff_cz_mass_fix.frcmod ./
cp -fr ${paramdir}/heme.frcmod ./
cp -fr ${paramdir}/heme.prep ./
cp -fr ${paramdir}/y2p.frcmod ./
cp -fr ${paramdir}/y2p.off ./
cp -fr ${paramdir}/vdw_AMBER_parm99.defn ./vdw.defn
cp -fr ${paramdir}/chem.defn ./



# Replace \r\n newlines
perl -pe 's/\r\n/\n/g' ./00.${system}.rec.noch.pdb > temp1.pdb
# Cleanup file and only retain essential records
grep '^ATOM  \|^TER  \|^END  ' ./temp1.pdb > 01.${system}.rec.noch.pdb

rm temp1.pdb

# initial TLeap processing to normalize receptor residue numbering
# and assign formal charges
##################################################

cat << EOF > 01.rec.leap.in
set default PBradii mbondi2
source oldff/leaprc.ff14SB
source leaprc.water.tip3p
loadoff ions.lib
REC = loadpdb 01.${system}.rec.noch.pdb 
check REC
saveamberparm REC ${system}.rec.preproc.parm7 01.${system}.rec.ori.crd
charge REC
quit
EOF
##################################################

${amberbin}/tleap -s -f 01.rec.leap.in >& 01.rec.tleap.v.out

if grep FATAL 01.rec.tleap.v.out; then echo "Tleap had some errors in ${system}, exiting..."; exit 1; fi

# strip Hs off of the receptor now that it has been renumbered/processed
${amberbin}/ambpdb -p ${system}.rec.preproc.parm7 -c 01.${system}.rec.ori.crd > temp1.pdb
${chimerabin}/chimera --nogui --script "${scriptdir}/removehs.chim.py \
                                        temp1.pdb \
                                        02.${system}.rec.noch.pdb" >& 02.chim2.v.out

# may be necessary to include more steps for SS/long bonds here but leaving that alone for the moment

rm temp1.pdb


### Prepare the ligand file with antechamber / CPC Specify gaff2 at as junmei suggested for COF
echo "Creating ligand prep file with antechamber"
${amberbin}/antechamber -i 00.${system}.lig.am1bcc.mol2 -fi mol2  -o 01.${system}.lig.ante.mol2 -fo mol2 -at gaff2 -j 5 -rn LIG -dr n > 00.${system}.lig.ante.out
echo "Creating ligand prep file with parmchk2"
${amberbin}/parmchk2 -i 01.${system}.lig.ante.mol2 -f mol2 -o ${system}.lig.ante.frcmod

### Prepare the cofactor file with antechamber, if it exists 
if [ $has_cofactor = true ]; then
  echo "Creating cofactor prep file with antechamber" 
  ${amberbin}/antechamber -fi mol2 -fo prepi -i 00.${system}.cof.am1bcc.mol2 -o $system.cof.ante.prep -at gaff2 -j 5 -rn COF -dr no > 01.${system}.cof.ante.out
  echo "Creating cofactor pdb file with antechamber" 
  ${amberbin}/antechamber -fi mol2 -fo pdb -i 00.${system}.cof.am1bcc.mol2 -o 01.${system}.cof.ante.pdb -j 5 -rn COF -dr no > 01.${system}.cof.ante2.out

  echo "Creating cofactor frcmod file with parmchk2" 
  ${amberbin}/parmchk2 -i ${system}.cof.ante.prep -f prepi -o ${system}.cof.ante.frcmod
fi

# use TLeap to prepare complex using receptor, ligand, cofactor
##################################################
cat > 02.com.leap.in<<EOF
set default PBradii mbondi2
source oldff/leaprc.ff14SB
source leaprc.gaff
loadoff ions.lib
loadamberparams ions.frcmod
loadamberparams frcmod.ions234lm_126_tip3p
loadamberparams frcmod.ions1lm_126_tip3p
loadamberparams frcmod.tip3p
loadamberparams heme.frcmod
loadamberprep heme.prep
loadoff y2p.off
loadamberparams y2p.frcmod
loadamberparams gaff_cz_mass_fix.frcmod
PRO = loadpdb 02.${system}.rec.noch.pdb
loadamberparams ${system}.lig.ante.frcmod
LIG = loadmol2 01.${system}.lig.ante.mol2
EOF
##################################################

# cofactor should be added to complex if present
if [ $has_cofactor = true ];then
        echo "Generating complex = pro+lig+cof"
        echo "loadamberparams ${system}.cof.ante.frcmod" >> 02.com.leap.in
        echo "loadamberprep ${system}.cof.ante.prep" >> 02.com.leap.in
        echo "COF = loadpdb 01.${system}.cof.ante.pdb" >> 02.com.leap.in
        echo "REC = combine { PRO COF }" >> 02.com.leap.in
        echo "saveamberparm COF ${system}.cof.parm 03.${system}.cof.ori.crd"  >> 02.com.leap.in
else
        echo "Generating complex = pro+lig (no cof)"
        echo "REC = combine { PRO }" >> 02.com.leap.in
fi
echo "COM = combine { REC LIG }" >> 02.com.leap.in
echo "saveamberparm LIG ${system}.lig.parm 03.${system}.lig.ori.crd"  >> 02.com.leap.in
echo "saveamberparm PRO ${system}.pro.parm 03.${system}.pro.ori.crd"  >> 02.com.leap.in
echo "saveamberparm REC ${system}.rec.parm 03.${system}.rec.ori.crd"  >> 02.com.leap.in
echo "saveamberparm COM ${system}.com.parm 03.${system}.com.ori.crd"  >> 02.com.leap.in
echo "quit" >> 02.com.leap.in

### Use leap to generate complex
echo "------------ LEAP RUN_002 SUMMARY -------------"
echo "Purpose: Generate complex with ssbonds"
${amberbin}/tleap -s -f 02.com.leap.in >& 03.${system}.com.leap.log
${amberbin}/ambpdb -p ${system}.lig.parm -tit "lig" -c 03.${system}.lig.ori.crd > 03.${system}.lig.ori.pdb
${amberbin}/ambpdb -p ${system}.pro.parm -tit "pro" -c 03.${system}.pro.ori.crd > 03.${system}.pro.ori.pdb
${amberbin}/ambpdb -p ${system}.rec.parm -tit "rec" -c 03.${system}.rec.ori.crd > 03.${system}.rec.ori.pdb
${amberbin}/ambpdb -p ${system}.com.parm -tit "com" -c 03.${system}.com.ori.crd > 03.${system}.com.ori.pdb
echo -n "atoms in 03.${system}.lig.ori.pdb = "
grep -c ATOM 03.${system}.lig.ori.pdb
echo -n "atoms in 03.${system}.pro.ori.pdb = "
grep -c ATOM 03.${system}.pro.ori.pdb
echo -n "atoms in 03.${system}.rec.ori.pdb = "
grep -c ATOM 03.${system}.rec.ori.pdb
echo -n "atoms in 03.${system}.com.ori.pdb = "
grep -c ATOM 03.${system}.com.ori.pdb

### Run sander to minimize hydrogen positions
echo "Creating ori.mol2 files before minimization"

echo "lig"
${amberbin}/ambpdb -p ${system}.lig.parm -c 03.${system}.lig.ori.crd -mol2 > ${system}.lig.ori.0.mol2
${amberbin}/antechamber -s 2 -i ${system}.lig.ori.0.mol2 -fi mol2 -at sybyl -j 5 -o 03.${system}.lig.ori.mol2 -fo mol2 -dr n > 03.lig.ori.ante.out
echo "pro"
${amberbin}/ambpdb -p ${system}.pro.parm -c 03.${system}.pro.ori.crd -mol2 > ${system}.pro.ori.0.mol2
${amberbin}/antechamber -s 2 -i ${system}.pro.ori.0.mol2 -fi mol2 -at sybyl -j 5 -o 03.${system}.pro.ori.mol2 -fo mol2 -dr n > 03.pro.ori.ante.out
echo "rec"
${amberbin}/ambpdb -p ${system}.rec.parm -c 03.${system}.rec.ori.crd -mol2 > ${system}.rec.ori.0.mol2
${amberbin}/antechamber -s 2 -i ${system}.rec.ori.0.mol2 -fi mol2 -at sybyl -j 5 -o 03.${system}.rec.ori.mol2 -fo mol2 -dr n > 03.rec.ori.ante.out
echo "com"
${amberbin}/ambpdb -p ${system}.com.parm -c 03.${system}.com.ori.crd -mol2 > ${system}.com.ori.0.mol2
${amberbin}/antechamber -s 2 -i ${system}.com.ori.0.mol2 -fi mol2 -at sybyl -j 5 -o 03.${system}.com.ori.mol2 -fo mol2 -dr n > 03.com.ori.ante.out
if [ $has_cofactor = true ]; then 
  ${amberbin}/ambpdb -p ${system}.cof.parm -c 03.${system}.cof.ori.crd -mol2 > ${system}.cof.ori.0.mol2
  ${amberbin}/antechamber -s 2 -i ${system}.cof.ori.0.mol2 -fi mol2 -at sybyl -j 5 -o 03.${system}.cof.ori.mol2 -fo mol2 -dr n > 03.cof.ori.ante.out
  rm ${system}.cof.ori.0.mol2
fi

#rm ${system}.lig.ori.0.mol2 ${system}.pro.ori.0.mol2 ${system}.rec.ori.0.mol2 ${system}.com.ori.0.mol2


# Determine ligand and receptor atoms to freeze based on atom types and bonding
# Necessary because ligand CYS component and covalent residue GLY overlap/clash/explode violently
# See script for more documentation on which residues are frozen
# If the system explodes, check sander.in to see if mask is accurate 
mask=`python ${scriptdir}/gen_com_rest_mask.py 03.${system}.com.ori.mol2`

##################################################
cat  >03.sander.in<<EOF
01mi minimization
 &cntrl
    ibelly = 1, bellymask = "${mask}"
    imin = 1, maxcyc = 100,
    ntpr = 10, ntx=1,
    ntb = 0, cut = 10.0,
    ntr = 1, drms=0.1,
    restraintmask = "!@H=",
    restraint_wt  = 1000.0
&end
EOF
##################################################

echo "---------------------------------------------------------"
echo "Minimizing complex with sander"
${amberbin}/sander -O -i 03.sander.in -o 03.sander.out -p ${system}.com.parm -c 03.${system}.com.ori.crd -ref 03.${system}.com.ori.crd -r 04.${system}.com.min.rst
${amberbin}/ambpdb -p ${system}.com.parm -tit "04.${system}.com.min" -c 04.${system}.com.min.rst > 04.${system}.com.min.pdb

# sanity check for failure
if [ ! -e 04.${system}.com.min.rst ]; then echo "Complex minimizaton failed! Terminating."; exit; fi

## Run sander on ligand alone to see if gaff screwed up anything
## (check RMSD between minimized complex ligand and this one)
echo "---------------------------------------------------------"

##################################################
cat  >03.sander.lig.in <<EOF
01mi minimization
 &cntrl
    imin = 1, maxcyc = 1000,
    ntpr = 10, ntx=1,
    ntb = 0, cut = 10.0,
    ntr = 0, drms=0.1,
&end
EOF
##################################################

echo "Minimizing unrestrained gas-phase ligand alone with sander"
${amberbin}/sander -O -i 03.sander.lig.in -o 03.sander.lig.out -p ${system}.lig.parm -c 03.${system}.lig.ori.crd -r 04.${system}.lig.only.min.rst
${amberbin}/ambpdb -p ${system}.lig.parm -c 04.${system}.lig.only.min.rst -mol2 > ${system}.lig.only.min.0.mol2
${amberbin}/antechamber -s 2 -i ${system}.lig.only.min.0.mol2 -fi mol2 -at sybyl -j 5 -o 04.${system}.lig.only.min.mol2 -fo mol2 -dr n > 04.lig.only.ante.out

rm ${system}.lig.only.min.0.mol2

# sanity check and report RMSD
grep "SANDER BOMB" 03.sander.lig.out
grep -A1 NSTEP 03.sander.lig.out | tail -2
echo -n "Minimizing Ligand 1000 steps alone rmsd "
python ${scriptdir}/calc_rmsd_mol2.py 03.${system}.lig.ori.mol2 04.${system}.lig.only.min.mol2

### Run sander on cofactor alone (if it exists) to see if gaff screwed up anything
if [ $has_cofactor = true ];then
        echo "---------------------------------------------------------"
        echo "Minimizing unrestrained gas-phase cofactor alone with sander"
        cp 03.sander.lig.in 03.sander.cof.in
        ${amberbin}/sander -O -i 03.sander.cof.in -o 03.sander.cof.out -p ${system}.cof.parm -c 03.${system}.cof.ori.crd -r 04.${system}.cof.only.min.rst
        ${amberbin}/ambpdb -p ${system}.cof.parm -c 04.${system}.cof.only.min.rst -mol2 > ${system}.cof.only.min.0.mol2
	${amberbin}/antechamber -s 2 -i ${system}.cof.only.min.0.mol2 -fi mol2 -at sybyl -j 5 -o 04.${system}.cof.only.min.mol2 -fo mol2 -dr n > 04.cof.only.ante.out

        grep "SANDER BOMB" 03.sander.cof.out
        grep -A1 NSTEP 03.sander.cof.out | tail -2
        echo -n "Minimizing Cofactor 1000 steps alone rmsd "
        python ${scriptdir}/calc_rmsd_mol2.py 03.${system}.cof.ori.mol2 04.${system}.cof.only.min.mol2

	rm ${system}.cof.only.min.0.mol2
else
        echo "No cofactor present to minimize"
fi

### Extract some files from the minimized complex
echo "---------------------------------------------------------"
echo "Extracting receptor with cpptraj"
echo "trajin 04.${system}.com.min.rst" > 04.rec.ptraj.in
echo "strip :LIG" >> 04.rec.ptraj.in
echo "trajout 04.${system}.rec.min.rst restart"  >> 04.rec.ptraj.in
${amberbin}/cpptraj ${system}.com.parm 04.rec.ptraj.in >& 04.rec.ptraj.out
grep STRIP 04.rec.ptraj.out
echo "Writing receptor mol2"

#Yuzhang modification made for multiple ligands aka waters etc.
${amberbin}/ambpdb -p ${system}.rec.parm -c 04.${system}.rec.min.rst -mol2 > ${system}.rec.min.0.mol2
${amberbin}/antechamber -i ${system}.rec.min.0.mol2 -fi mol2 -at sybyl -j 5 -o 04.${system}.rec.min.mol2 -fo mol2 -dr n >& 04.${system}.rec.ante.out

rm ${system}.rec.min.0.mol2

# extract ligand from minimized complex
echo "Creating ligand mol2 file"
echo "trajin 04.${system}.com.min.rst" > 04.lig.ptraj.in
echo "strip !(:LIG)" >> 04.lig.ptraj.in
echo "trajout 04.${system}.lig.min.rst restart"  >> 04.lig.ptraj.in
${amberbin}/cpptraj ${system}.com.parm 04.lig.ptraj.in >& 04.lig.ptraj.out
grep STRIP 04.lig.ptraj.out
${amberbin}/ambpdb -p ${system}.lig.parm -c 04.${system}.lig.min.rst -mol2 > ${system}.lig.min.0.mol2
${amberbin}/antechamber -i ${system}.lig.min.0.mol2 -fi mol2 -at sybyl -j 5 -o 04.${system}.lig_CYS.min.mol2 -fo mol2 -dr n >& 04.lig.min.ante.out

rm ${system}.lig.min.0.mol2

${amberbin}/ambpdb -p ${system}.com.parm -c 04.${system}.com.min.rst -mol2 > ${system}.com.min.0.mol2
${amberbin}/antechamber -i ${system}.com.min.0.mol2 -fi mol2 -at sybyl -j 5 -o 04.${system}.com.min.mol2 -fo mol2 -dr n >& 04.com.min.ante.out

rm ${system}.com.min.0.mol2

# calculate RMSDs for sanity checking (should be < 0.1)
echo -n "Minimized Ligand rmsd: "
python ${scriptdir}/calc_rmsd_mol2.py 03.${system}.lig.ori.mol2 04.${system}.lig_CYS.min.mol2
echo -n "Minimized Receptor rmsd: "
python ${scriptdir}/calc_rmsd_mol2.py 03.${system}.rec.ori.mol2 04.${system}.rec.min.mol2
echo -n "Minimized Complex rmsd: "
python ${scriptdir}/calc_rmsd_mol2.py 03.${system}.com.ori.mol2 04.${system}.com.min.mol2
python ${scriptdir}/clean_mol2.py 04.${system}.rec.min.mol2 05.${system}.rec.clean.mol2
#python ${scriptdir}/clean_mol2.py 04.${system}.lig_CYS.min.mol2 05.${system}.lig_CYS.clean.mol2

${amberbin}/ambpdb -p ${system}.rec.parm -c 04.${system}.rec.min.rst > 05.${system}.rec.clean.pdb

# create PDB file of minimized ligand for use in sphere gen
${amberbin}/antechamber -fi mol2 -i 04.${system}.lig_CYS.min.mol2 -fo pdb -o 05.${system}.lig.clean.pdb -dr no > 04.lig.ante.v.out

# remove CYS from ligand and add dummy atoms
python ${scriptdir}/remove_CYS.py 04.${system}.lig_CYS.min.mol2 05.${system}.lig.clean.mol2

echo "Receptor preparation and complex minimization completed."

# remove some excess files

# antechamber verbose
declare -a to_remove=(ANTECHAMBER_AC.AC ANTECHAMBER_AC.AC0 ANTECHAMBER_BOND_TYPE.AC ANTECHAMBER_BOND_TYPE.AC0 ATOMTYPE.INF ANTECHAMBER_PREP.AC ANTECHAMBER_PREP.AC0 NEWPDB.PDB PREP.INF)
for file in "${to_remove[@]}"; do
	if [ -e ${file} ]; then
		rm ${file}
	fi
done

# param files
declare -a to_remove=(chem.defn gaff_cz_mass_fix.frcmod heme.frcmod heme.prep ions.frcmod ions.lib vdw.defn y2p.frcmod y2p.off)
for file in "${to_remove[@]}"; do
        if [ -e ${file} ]; then
                rm ${file}
        fi
done
