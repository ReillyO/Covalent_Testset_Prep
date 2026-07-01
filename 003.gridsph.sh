#!/bin/bash

# requires prepared ligand and receptor files to be present in 
# the 002.recprep directory (should have run steps 001 and 002)

# will create covalent-compatible spheres and grid based on the
# files for the system provided

dockdir=${DOCKHOME}
amberdir=${AMBERHOME}
rootdir=${ROOTDIR}
chimeradir=${CHIMERAHOME}
testsetdir=${rootdir}/zzz.testset_files/
rawfiledir=${rootdir}/zzz.master/
scriptdir=${rootdir}/zzz.scripts/
paramdir=${rootdir}/zzz.parameters/
system=${1}

sph_radius=10.0
sph_count=75

# preliminary sanity checks for directory and file existence
echo "Setting up environment..."
if [ -e ${testsetdir}/${system}/003.gridsph/ ]; then
  rm -r ${testsetdir}/${system}/003.gridsph/
fi
mkdir ${testsetdir}/${system}/003.gridsph/
cd ${testsetdir}/${system}/003.gridsph/

if [ ! -e ${testsetdir}/${system}/001.ligprep ] || [ ! -e ${testsetdir}/${system}/002.recprep ]; then
  echo "please run ligand and receptor prep scripts for ${system} before creating the grid and spheres"
  exit 1
fi

# copy in pertinent files
echo "Copying pertinent files..."
# copy over requisite ligand/receptor files
if [ -e ${testsetdir}/${system}/002.recprep/05.${system}.rec.clean.mol2 ]; then
  cp ${testsetdir}/${system}/002.recprep/05.${system}.rec.clean.mol2 ./00.${system}.rec.clean.mol2
else
  echo "missing receptor file 05.${system}.rec.clean.mol2! This may indicate receptor preparation did not complete satisfactorily. exiting"
  exit 1
fi

if [ -e ${testsetdir}/${system}/002.recprep/05.${system}.rec.clean.pdb ]; then
  cp ${testsetdir}/${system}/002.recprep/05.${system}.rec.clean.pdb ./00.${system}.rec.clean.pdb
else
  echo "missing receptor file ${system}.rec.clean.pdb! This may indicate receptor preparation did not complete satisfactorily. exiting"
  exit 1
fi

if [ -e ${testsetdir}/${system}/002.recprep/05.${system}.lig.clean.mol2 ]; then
  cp ${testsetdir}/${system}/002.recprep/05.${system}.lig.clean.mol2 ./00.${system}.lig.clean.mol2
else
  echo "missing ligand file 05.${system}.lig.clean.mol2! exiting"
  exit 1
fi

if [ -e ${testsetdir}/${system}/002.recprep/05.${system}.lig.clean.pdb ]; then
  cp ${testsetdir}/${system}/002.recprep/05.${system}.lig.clean.pdb ./00.${system}.lig.clean.pdb
else
  echo "missing ligand file 05.${system}.lig.clean.pdb! exiting"
  exit 1
fi

# truncate provided receptor to atoms within 20A of ligand
${scriptdir}/keep_close_atoms.pl 00.${system}.lig.clean.mol2 00.${system}.rec.clean.pdb 15 > rec_for_dms.pdb 

# create receptor surface with no Hs using DMS
echo "Creating surface using DMS..."
/gpfs/projects/rizzo/zzz.programs/dms/bin/dms rec_for_dms.pdb -a -g ${system}.rec.dms.log -n -o ${system}.rec.dms.out > dms.v.out

if grep -iq error dms.v.out; then echo "DMS had some errors, exiting..."; exit 1; fi

# generate receptor spheres using sphgen - these will only be used in showbox,
# not for orienting the ligand
echo "Generating receptor spheres with sphgen..."
###########################################
cat << EOF > INSPH
${system}.rec.dms.out
R
X
0.0
4.0
1.4
${system}.rec.sph
EOF
###########################################

# remove old sphgen files
if [ -e temp* ]; then rm temp*; fi
if [ -e ${system}.rec.sph ]; then rm ${system}.rec.sph; fi
if [ -e OUTSPH ]; then rm OUTSPH; fi

sphgen -i INSPH -o OUTSPH > sphgen.v.out

if grep -iq error sphgen.v.out; then echo "sphgen had some errors, exiting..."; exit 1; fi

### Convert the clusters (pruned spheres, not cluster 0) to a PDB file for viewing
##################################################
cat  >showsphere.in<<EOF
${system}.rec.sph
-1
N
clustertemp
N
EOF
##################################################
${dockdir}/showsphere < showsphere.in > showsphere.v.txt

### Make a PDB file that contains all the pruned clusters for viewing
cat clustertemp* >> all.clust.pdb
rm clustertemp*

if  ! ls -l | grep -q "all.clust.pdb" ; then echo "WARNING:   Spheres were not successfully generated!"; exit 1; fi

# pare down spheres to top 75 of those within 10 Angstroms of ligand
echo "Selecting spheres near ligand..."
${scriptdir}/keep_close_spheres.pl 00.${system}.lig.clean.mol2 all.clust.pdb 10.0 75
mv temp.sph ligand_spheres_for_grid.sph

# generate orienting spheres for covalent docking from SG/CB/CA atoms
# contained in the processed ligand and receptor files
echo "Generating orienting spheres..."
python ${scriptdir}/SG_CB_from_lig.py 00.${system}.lig.clean.pdb temp1.pdb
python ${scriptdir}/pdbtosph.py temp1.pdb SG_CB.sph

rm temp1.pdb

# create showbox input file to create a box around the active
# site that grid will use
echo "Creating box around binding site..."
###########################################
cat << EOF > showbox.in
Y
8.0
ligand_spheres_for_grid.sph
1
${system}.box.pdb
EOF
###########################################
showbox < showbox.in > showbox.v.out

if grep -iq error showbox.v.out; then echo "showbox had some errors, exiting..."; exit 1; fi

# create grid files for the receptor based on the chosen box
# and receptor file mol2
echo "Generating receptor grid..."
###########################################
cat << EOF > grid.in
compute_grids                             yes
grid_spacing                              0.4
output_molecule                           no
contact_score                             no
energy_score                              yes
energy_cutoff_distance                    9999
atom_model                                a
attractive_exponent                       6
repulsive_exponent                        9
distance_dielectric                       yes
dielectric_factor                         4
allow_non_integral_charges                yes
bump_filter                               yes
bump_overlap                              0.75
receptor_file                             00.${system}.rec.clean.mol2
box_file                                  ${system}.box.pdb
vdw_definition_file                       /gpfs/projects/rizzo/zzz.programs/dock6.9_mpiv2018.0.3/parameters/vdw_AMBER_parm99.defn
chemical_definition_file                  /gpfs/projects/rizzo/zzz.programs/dock6.9_mpiv2018.0.3/parameters/chem.defn
score_grid_prefix                         ${system}.rec
EOF
##########################################

grid -i grid.in -o grid.out > grid.v.out

if grep -iq error grid.out; then echo "Grid had some errors, exiting..."; exit 1; fi
if grep -iq error grid.v.out; then echo "Grid had some errors, exiting..."; exit 1; fi

echo "Sphere and grid generation completed with no errors noted"

cd ${rootdir}
