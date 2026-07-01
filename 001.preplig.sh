#!/bin/bash

# Currently accepts a Mol2 file of the ligand (including attached CYS residue) with hydrogens and Gasteiger charge

dockdir=${DOCKHOME}
amberdir=${AMBERHOME}
chimeradir=${CHIMERAHOME}
rootdir=${ROOTDIR}
testsetdir=${rootdir}/zzz.testset_files/
rawfiledir=${rootdir}/zzz.master/
scriptdir=${rootdir}/zzz.scripts/
system=${1}
has_cofactor=false

# check that raw file exists in zzz.master
if [ ! -e ${rawfiledir}/${system}.lig.ch.mol2 ]; then
  echo "Cannot find system; please check that corresponding ${system}.lig.ch.mol2 exists in zzz.master/ or correct the spelling"
  exit 1 
fi

# check for cofactor and update variable accordingly
if [ -e ${rawfiledir}/${system}.cof.ch.mol2 ]; then 
  has_cofactor=true
fi

# make system directory if not present
if [ ! -e ${testsetdir}/${system} ]; then mkdir ${testsetdir}/${system}; fi

# remake system ligprep directory from scratch
if [ -e ${testsetdir}/${system}/001.ligprep ]; then rm -r ${testsetdir}/${system}/001.ligprep; fi
mkdir ${testsetdir}/${system}/001.ligprep

# move the ligand file and cofactor to the ligprep directory
cp ${rawfiledir}/${system}.lig.ch.mol2 ${testsetdir}/${system}/001.ligprep/
if [ $has_cofactor = true ]; then cp ${rawfiledir}/${system}.cof.ch.mol2 ${testsetdir}/${system}/001.ligprep/; fi
cd ${testsetdir}/${system}/001.ligprep

echo "Cleaning up the ligand file..."
mv ${system}.lig.ch.mol2 00.${system}.lig.ch.mol2
perl -pe 's/\r\n/\n/g' 00.${system}.lig.ch.mol2 > temp1.mol2     # get rid of any legacy Windows newlines
python ${scriptdir}/fix_ligname_residues_mol2.py temp1.mol2 01.${system}.lig.ch.mol2 LIG  # standardize ligand name and numbering/chain

if [ $has_cofactor = true ]; then
  echo "Cleaning up cofactor file... "
  mv ${system}.cof.ch.mol2 00.${system}.cof.ch.mol2
  perl -pe 's/\r\n/\n/g' 00.${system}.cof.ch.mol2 > temp1.mol2
  python ${scriptdir}/fix_ligname_residues_mol2.py temp1.mol2 01.${system}.cof.ch.mol2 COF 
fi

rm temp1.mol2


echo "Preparing ligand with DOCK..."

##################################################
cat <<EOF >01.dock.lig.in

conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             01.${system}.lig.ch.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              no
ligand_outfile_prefix                                        02.${system}.lig
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no

EOF
##################################################

# reformat ligand Mol2 in DOCK-friendly format
dock6 -i 01.dock.lig.in -o 01.dock.lig.out

# do the same for the cofactor file if exists 
if [ $has_cofactor = true ]; then
  ##################################################
  cat <<EOF >01.dock.cof.in

conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             01.${system}.cof.ch.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              no
ligand_outfile_prefix                                        02.${system}.cof
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no

EOF
##################################################

  dock6 -i 01.dock.cof.in -o 01.dock.cof.out

fi

# check for errors in the DOCK processing
if grep -qi error 01.dock.lig.out; then echo "Error in DOCK ligand processing! Exiting..."; exit 1; fi
if [ $has_cofactor = true ] ; then 
  if grep -qi error 01.dock.cof.out; then echo "Error in DOCK cofactor processing! Exiting..."; exit 1; fi
fi

echo "Assigning AM1BCC charges to ligand with Antechamber..."
# charge the DOCK output molecule with AM1BCC empirical charges via Antechamber
${amberdir}/antechamber -fi mol2 -fo mol2 -c bcc -at sybyl -s 2 -pf y -i 02.${system}.lig_scored.mol2 -o 03.${system}.lig.am1bcc.mol2 -dr no > 02.lig.ante1.out

if grep -qi error 02.lig.ante1.out; then echo "Error in ligand Antechamber step 1! Exiting..."; exit 1; fi

# if initial SQM run doesn't converge, we run a looser parameter sim
if [ `grep "No convergence in SCF" sqm.out | wc -l` ]; then
        echo "Running second round of Antechamber due to no convergence..."
	${amberdir}/antechamber -fi mol2 -fo mol2 -c bcc -j 5 -at sybyl -s 2 -pf y -ek "itrmax=1000, qm_theory='AM1', grms_tol=0.01, tight_p_conv=0, scfconv=1.d-8" -i 02.${system}.lig_scored.mol2 -o 03.${system}.lig.am1bcc.mol2 -dr no > 02.lig.ante2.out
	if grep -qi error 02.lig.ante2.out; then echo "Error in ligand Antechamber step 2! Exiting..."; exit 1; fi
fi


# repeat the whole antechamber charging process for the cofactor if needed
if [ $has_cofactor = true ]; then
  echo "Assigning AM1BCC charges to cofactor with Antechamber..."
  ${amberdir}/antechamber -fi mol2 -fo mol2 -c bcc -at sybyl -s 2 -pf y -i 02.${system}.cof_scored.mol2 -o 03.${system}.cof.am1bcc.mol2 -dr no > 02.cof.ante1.out
  if grep -qi error 02.cof.ante1.out; then echo "Error in cofactor Antechamber step 1! Exiting..."; exit 1; fi

  if [ `grep "No convergence in SCF" sqm.out | wc -l` ]; then
    echo "Running second round of Antechamber due to no convergence..."
    ${amberdir}/antechamber -fi mol2 -fo mol2 -c bcc -j 5 -at sybyl -s 2 -pf y -ek "itrmax=1000, qm_theory='AM1', grms_tol=0.01, tight_p_conv=0, scfconv=1.d-8" -i 02.${system}.cof_scored.mol2 -o 03.${system}.cof.am1bcc.mol2 -dr no > 02.cof.ante2.out
  fi
  if grep -qi error 02.cof.ante2.out; then echo "Error in cofactor Antechamber step 2! Exiting..."; exit 1; fi
fi


# also create a PDB file of the ligand for sphere generation
${chimeradir}/chimera --nogui --script "${scriptdir}/removehs.chim.py 03.${system}.lig.am1bcc.mol2 temp1.pdb" > 03.chim2.lig.v.out
grep '^ATOM  \|^TER  \|^END  \|^HETATM' temp1.pdb > 03.${system}.lig.noch.pdb 

if grep -qi error 03.chim2.lig.v.out; then echo "Error in Chimera step 2! Exiting..."; exit 1; fi

cd ${rootdir}

echo "All ligand prep steps completed without any errors noted"

