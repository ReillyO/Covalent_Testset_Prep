#!/bin/bash

# conducts a standard virtual screening of all systems that have a 
# grid prepared using the standard covalent protocol

# sample usage: bash run.expt.poserep.sh


dockdir=${DOCKHOME}
amberdir=${AMBERHOME}
rootdir=${ROOTDIR}
chimeradir=${CHIMERAHOME}
testsetdir=${rootdir}/zzz.testset_files/
rawfiledir=${rootdir}/zzz.master/
scriptdir=${rootdir}/zzz.scripts/
paramdir=${rootdir}/zzz.parameters/
analysisdir=${rootdir}/zzz.analysis

prefix="004.poserep"

n_run=0

# runs a VS on every system in the testset directory (4 letter codes only sorry)
for f in `ls -d ${testsetdir}/*`; do
  system=${f: -4}
  sysdir=${testsetdir}/${system}
  
  echo "Copying pertinent files..."
  if [ -e ${sysdir}/${prefix} ]; then rm -rf ${sysdir}/${prefix}; fi
  mkdir ${sysdir}/${prefix}
  cd ${sysdir}/${prefix} 
  
  if [ -e ${sysdir}/002.recprep/05.${system}.lig.clean.mol2 ]; then 
    cp ${sysdir}/002.recprep/05.${system}.lig.clean.mol2 ./00.${system}.lig.clean.mol2
  else echo "ERROR: cannot find ${system} ligand file. please confirm ligand prep executed successfully."; fi

  if [ -e ${sysdir}/003.gridsph/${system}.rec.bmp ]; then
    cp ${sysdir}/003.gridsph/${system}.rec.bmp ./
    cp ${sysdir}/003.gridsph/${system}.rec.nrg ./
    cp ${sysdir}/003.gridsph/SG_CB.sph ./
  else echo "ERROR: cannot find ${system} grid files. please confirm receptor prep executed successfully."; fi

echo "Creating and running DOCK job..."
##################################################################################
  cat << EOF > cov.in
conformer_search_type                                        covalent
pruning_use_clustering                                       yes
pruning_max_orients                                          1000
pruning_clustering_cutoff                                    100
covalent_bondlength                                          1.6:0.1:2.0
covalent_bondlength2                                         -1.0
covalent_angle                                               -1.0
covalent_dihedral_step                                       10.0
pruning_conformer_score_cutoff                               100.0
use_clash_overlap                                            no
write_growth_tree                                            no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ./00.${system}.lig.clean.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes
rmsd_reference_filename                                      ./00.${system}.lig.clean.mol2
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           ./SG_CB.sph
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
grid_score_primary                                           yes
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_lig_efficiency                                          no
grid_score_grid_prefix                                       ./${system}.rec
minimize_ligand                                              yes
minimize_anchor                                              no
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
minimize_flexible_growth_ramp                                yes
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_initial_score_coverge                                5
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_grow_max_iterations                                  0
simplex_grow_tors_premin_iterations                          1000
simplex_final_min                                            no
simplex_random_seed                                          0
simplex_restraint_min                                        yes
simplex_coefficient_restraint                                10.0
atom_model                                                   all
vdw_defn_file                                                ${paramdir}/vdw_AMBER_parm99.defn
flex_defn_file                                               ${paramdir}/flex.defn
flex_drive_file                                              ${paramdir}/flex_drive.tbl
ligand_outfile_prefix                                        ${system}_out
write_mol_solvation                                          no
write_orientations                                           yes
num_final_scored_poses                                       1000
num_preclustered_conformers                                  1000
write_conformations                                          no
cluster_conformations                                        no
score_threshold                                              100.0
rank_ligands                                                 no
EOF
############################################################################
 
  cat <<EOF>submit.sh
#!/bin/bash
 
${dockdir}/dock6 -V -i cov.in -o cov.out &
EOF
  
  chmod +x submit.sh
  srun --mem-per-cpu=5000 --exclusive --ntasks-per-core=1 -N1 -n1 -W 0 submit.sh &
  
done

# wait for all docking to complete
wait

