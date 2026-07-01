# Will accept a ligand-receptor complex mol2 and print out a 
# mask describing the residues to be immobilized for minimization
# as well as log them in a text file
#
# Currently this is 
# - the ligand CYS attachment including Hs 
# - the GLY of the residue where the covalent bond occurs
# - the residues directly adjacent to the GLY
# eg !:24-26&!:216@C,O,N,CA,CB,H,HN1,HN3,HN4,HB2,HB3

# assumes that the CA closest to the CA of the ligand CYS attachment
# is the residue that the ligand covalently bonds (safe assumption as
# they should overlap)

# sample: python cov_spheres_from_ligrec.py $complexfile <$outfile>
# outfile is optional, default is mask.txt

import sys, os
import numpy as np
import pdb_utils as pdb

try:
  comfile = sys.argv[1]
except:
  print("please provide a rec-lig complex file in style: python cov_spheres_from_ligrec.py $complexfile <$outfile>")
  exit()


lig_resnum = -1
in_atoms = False
in_bonds = False
CYS_atom_namelist = ['N','C','O','CA','CB','SG'] # standard
CYS_atom_idx_list = []
CYS_H_idx_list = []
lig_CA = [] 
rec_CA_list = []

# loop through file and collect:
# from atoms:
# - ligand resnum
# - ligand CYS attachment heavy atom indices
# - ligand CYS alpha carbon line
# - all CA atoms in the receptor
# from bonds:
# - indexes of Hs bonded to ligand CYS

with open(comfile, 'r') as f:
  for line in f.readlines():
    if '@<TRIPOS>ATOM' in line: in_atoms = True; continue
    elif '@<TRIPOS>BOND' in line: in_atoms = False; in_bonds = True; continue
    elif '@<TRIPOS>SUBSTRUCTURE' in line: in_bonds = False; continue
    
    if in_atoms:
      linelist = line.split()
      if linelist[7] == 'LIG':
        lig_resnum = int(linelist[6])
        if linelist[1] == 'CA':
          lig_CA = linelist
        if linelist[1] in CYS_atom_namelist:
          CYS_atom_idx_list.append(int(linelist[0]))
      else:
        if linelist[1] == 'CA':
          rec_CA_list.append(linelist)
    elif in_bonds:
      linelist = line.split()
      a1 = int(linelist[1])
      a2 = int(linelist[2])
      if a1 in CYS_atom_idx_list and a2 not in CYS_atom_idx_list:
        CYS_H_idx_list.append(a2)
      elif a2 in CYS_atom_idx_list and a1 not in CYS_atom_idx_list:
        CYS_H_idx_list.append(a1)

# go back through atoms and collect CYS-bonded H names for mask
with open(comfile, 'r') as f:
  for line in f.readlines():
    if '@<TRIPOS>ATOM' in line: in_atoms = True; continue
    elif '@<TRIPOS>BOND' in line: in_atoms = False; continue
    if in_atoms:
      linelist = line.split()
      if int(linelist[0]) in CYS_H_idx_list and linelist[5] == 'H':
        CYS_atom_namelist.append(linelist[1])


# find receptor CA with the shortest distance to ligand CA
bestdist = 100000
CA_match = rec_CA_list[0]
lig_CA_coord = np.array((float(lig_CA[2]), float(lig_CA[3]), float(lig_CA[4])))
CA_resnum = -1
for CA in rec_CA_list:
  rec_CA_coord = np.array((float(CA[2]), float(CA[3]), float(CA[4])))
  dist = np.linalg.norm(rec_CA_coord - lig_CA_coord)
  if dist < bestdist:
    CA_match = CA
    bestdist = dist
    CA_resnum = int(CA[6])


# generate atom mask describing all atoms except:
#  - ligand CYS and attached Hs
#  - covalently bound residue
mask = "".join(["!:", str(CA_resnum-1), "-", str(CA_resnum+1),
                "&!:", str(lig_CA[6]), "@", ",".join(CYS_atom_namelist)])

# print to console for direct use by script
print(mask)

