# Manual Covalent System Prep Guidelines

1) Download system PDB using Chimera fetch or direct from RCSB PDB

### Ligand Preparation:
2) Select and isolate ligand + receptor CYS residue (ex `sel #0:216 #0:25`); save to PDB file as `SYS.lig.noch.pdb`
3) Edit `SYS.lig.noch.pdb` plaintext to:
  * remove the protein TER line,
  * change all CYS atoms to be part of the LIG (keeping the name is fine),
  * and make residue numbering/ligand name consistent
4) Add Hs and charge ligand in Chimera using Gasteiger standard, ensure formal charge matches expected vals from literature/functional groups. Sulfates need to be treated special. It is alright for the amide N to be charged as it has minimal effect on the covalently bound ligand atoms. Save the ligand as `SYS.lig.ch.mol2`

Common ligand errors tend to show up in automated antechamber script and include:
  * weird geometries leading to improper protonation
  * cov-adjacent oxygen treated as ketone instead of hydroxyl 
  * cov-adjacent Ns treated as single bonded instead of double

### Cofactor (if present):
5) Follow same protocol as ligand (minus receptor covalent residue): isolate (as `SYS.cof.noch.pdb`), protonate, and charge. Save as `SYS.cof.ch.mol2`

### Receptor:
6) Remove ligand, cofactor, HOHs, and other molecular artifacts from receptor 
7) Convert covalent residue to a GLYCINE to avoid steric effects/interactions using `swapaa GLY #0:25` if covalent residue is 25
7) Save as `SYS.rec.noch.pdb` (charging is done by the automated scripts)

Common receptor errors include:
 * missing OXT

8) After receptor, ligand, and cofactor structures have been prepared, move them to the `zzz.master` directory in the covalent script workspace for the next step of automated preparation.
