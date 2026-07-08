# Covalent_Testset_Prep
Suite of scripts for preparing covalently binding ligand-receptor complexes for DOCKing. 


Scripts require the following programs:
 * DOCK6 compiled executable and supporting parameter files
 * Amber suite - specifically antechamber, tleap, ambpdb, parmchk2, sander, cpptraj
 * Chimera command line interface
 * DMS executable
 * Python3 available in environment

## 0: Pre-Processing

Before any automated scripts can be run, systems need to be prepared manually in a local Chimera or ChimeraX GUI. The standard test set currently uses Chimera. Guidelines for doing this are detailed in PREP.md. 

The essential files to have on hand for each system after preprocessing are:
 * `SYS.lig.ch.mol2`
 * `SYS.rec.noch.pdb`
 * `SYS.cof.ch.mol2` (if cofactor present)
These files should be copied into the `zzz.master` folder prior to automated preparation.

## 1: Automated Preparation via Scripts

The automated scripts provided in this repo can be used to:
 * Prepare the preprocessed ligand and receptor/cofactor for docking using DOCK6
 * Generate Amber-compatible coordinate and topology files
 * Run pose reproduction experiments on a local machine or cluster
 * And more to come

Before attempting any prep, `000.setup.sh` should be edited to reflect the user's home environment and directories for the specified programs. 
This script should be run using `source` whenever the environment is reloaded, for example:

```source 000.setup.sh```

### Ligand Preparation:

The `001.preplig.sh` script prepares the ligand and the cofactor by first standardizing the file with DOCK6 and subsequently running Amber's `antechamber` program to assign AM1BCC charges to the molecules. The script accepts a system name as the argument, for example:

```bash 001.preplig.sh 1ABC```

Note that charge assignment has quantum resolution and can take a significant amount of time (up to half an hour) to converge, but this is system-dependent. 

### Complex Preparation: 

The `002.preprec.sh` script first prepares the receptor, using TLeap to add hydrogens and standard charges, before creating and minimizing the hydrogens of the receptor-ligand-cofactor complex. Minimization contains some additional restraints to prevent extreme collisions of the ligand cysteine atoms with the backbone of the receptor. For example:

```bash 002.preprec.sh 1ABC```

### Grid and Sphere Generation:

DOCK6 uses a grid representation of the receptor for docking, and spheres for initially orienting the ligand to the receptor. The `003.gridsph.sh` script:
 * Truncates the receptor to residues within 20 Angstroms of the ligand
 * Creates a receptor surface using the program `dms`,
 * Generates orienting spheres using `sphgen`,
 * Retains the 75 closest spheres to the ligand
 * Generates a grid bounding box based on the spheres
 * Generates a DOCK6-compatible grid for the receptor

For example:

```bash 003.gridsph.sh 1ABC```

### Full System Preparation

Users can also use the `run.prep.allsteps.sh` script to prepare the ligand, complex, and grid in series. This also accepts the system code as the sole argument. For example:

```bash run.prep.allsteps.sh 1ABC```

### Re-Running Incomplete Preparation

The `run.prep.incomplete.sh` script performs the same actions as the `run.prep.allsteps.sh` script, for all systems which are not identified as complete (having a receptor file in zzz.master but no .bmp grid file in the system prep directory). 

Incomplete systems can be identified explicitly by the `checkincomplete.sh` script in the zzz.scripts directory. 

## 2: Benchmarking

After systems have been prepared, they can be used to benchmark the current DOCK6 implementation. 

### Pose Reproduction

The `run.expt.004.poserep.sh` script provides a template for other pose reproduction-type experiments. It expects the user to have access to an HPC-type cluster with SLURM. It creates experimental directories for each system, and then re-docks the ligand into the receptor. 

(This is still under active development and more functionality will be added in the near future)

### Data Collection

The `run.anl.getrmsd.sh` script can be used to obtain heavy-atom RMSDs for the poses resulting from a pose reproduction experiment. They are output in CSV format for ready analysis. 

(This is still under active development and more functionality will be added in the near future)
