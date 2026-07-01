# Script will accept a covalent ligand file with fill CYS residue 
# covalently bound as in the receptor and write a ligand mol2 file
# with only the SG and CB atoms present

# sample usage: python remove_CYS.py $ligfile $outfile

import sys, os
import pdb_utils as pdb
from collections import defaultdict

try:
    ligfile = sys.argv[1]
    outfile = sys.argv[2]
except:
    print("Please provide a ligand file and output file as in: python remove_CYS.py $ligfile $outfile")
    exit()

# atoms always part of standard CYS residue on covalently bound ligand
reject_atoms = ["N", "C", "O", "CA", "CB"]
# list of atoms to remove
reject_atom_idxs = []
H_idxs = []
output = ""
CB_idx = -1

# state machine definition
in_header = True
in_atoms = False
in_bonds = False
in_end = False
line_num = 0

# run through file once to collect all Hs connected to the known CYS
# heavy atoms for future deletion
with open(ligfile, 'r') as f:
    for line in f.readlines():
        if "@<TRIPOS>ATOM" in line: in_header = False; in_atoms = True; continue
        elif "@<TRIPOS>BOND" in line:  in_atoms = False; in_bonds = True; continue
        elif "@<TRIPOS>SUBSTRUCTURE" in line: in_bonds = False; in_end = True; continue
       
        # collect 1) all H idxs and 2) idxs of heavy atoms to remove 
        if in_atoms:
            line_list = line.split() # num name co or ds type resnum resname charge
            a_idx = line_list[0]
            aname = line_list[1]
            atype = line_list[5]
            if aname in reject_atoms:
                reject_atom_idxs.append(a_idx) # str
                if aname == "CB": CB_idx = a_idx
            elif 'H' in atype:
                H_idxs.append(a_idx)
        # find idxs of Hs to remove based on bonding to heavy atoms 
        elif in_bonds: # after atom section
            _,a1,a2,__ = line.split() # idx atom1 atom2 type
            if a1 in reject_atom_idxs and a2 not in reject_atom_idxs and a2 in H_idxs:
                reject_atom_idxs.append(a2)
            elif a2 in reject_atom_idxs and a1 not in reject_atom_idxs and a1 in H_idxs:
                reject_atom_idxs.append(a1)

# ensure CB is not removed but attached Hs are
if CB_idx != -1:
    reject_atom_idxs.pop(reject_atom_idxs.index(CB_idx))

# track index changes as atoms are removed
atom_idx_dec = 0 
bond_idx = 1 
atoms_removed_idx = []
bond_map = defaultdict(str) # original index : new index

# state machine definition
in_header = True
in_atoms = False
in_bonds = False
in_end = False
line_num = 0

with open(ligfile, 'r') as f:
    for line in f.readlines(): 
        if "@<TRIPOS>ATOM" in line: output += line; in_header = False; in_atoms = True; continue
        elif "@<TRIPOS>BOND" in line: output += line; in_atoms = False; in_bonds = True; continue
        elif "@<TRIPOS>SUBSTRUCTURE" in line: output += line; in_bonds = False; in_end = True; continue
        
        # in header, metadata needs to be edited for # atoms and # bonds
        if in_header:
            line_list = line.split()
            if line_num == 2: # metadata line
                num_bonds = int(line_list[0]) - len(reject_atom_idxs)
                num_atoms = int(line_list[1]) - len(reject_atom_idxs)
                line_list[0] = str(num_bonds)
                line_list[1] = str(num_atoms)
                line = " ".join(line_list) + "\n"
                output += line
            else:
                output += line
        # CYS atoms will be removed and other atom indices corrected
        elif in_atoms:
            line_list = line.split()
            a_idx = line_list[0]
            aname = line_list[1]
            atype = line_list[5]
            if a_idx not in reject_atom_idxs: # keep and update idx
                # determine index based on number of already-deleted atoms
                new_idx = int(a_idx) - atom_idx_dec
                # replace index with correct one
                line = line.replace(a_idx.rjust(5), str(new_idx).rjust(5), 1)
                # add entry to bond map for future replacement
                bond_map[a_idx] = new_idx

                if aname == 'CB': # replace CB with appropriate dummy
                    line = line.replace(' CB ', ' D2 ').replace('C.3', 'Du ')
                elif aname == 'SG': # replace SG with appropriate dummy
                    line = line.replace(' SG ', ' D1 ').replace('S.3', 'Du ')

                # add to output
                output += line
            else: # increment the counter and do not write the atom
                reject_atom_idxs.pop(reject_atom_idxs.index(a_idx))
                atoms_removed_idx.append(a_idx)
                atom_idx_dec += 1
        # any bonds pertaining to CYS atoms should also be removed, and 
        # indices of remaining bonded atoms should be corrected
        elif in_bonds:
            bond_list = line.split() # idx a1 a2 type
            # only add the bond line if it doesn't participate in removed CYS atoms
            if len(bond_list) == 4 and bond_list[1] not in atoms_removed_idx and bond_list[2] not in atoms_removed_idx:
                bond_list[0] = str(bond_idx)
                bond_list[1] = str(bond_map[bond_list[1]])
                bond_list[2] = str(bond_map[bond_list[2]])
                line = "".join([bond_list[0].rjust(6),
                                bond_list[1].rjust(5), 
                                bond_list[2].rjust(5),
                                " ",
                                bond_list[3],
                                "\n"])
                output += line
                bond_idx += 1
        elif in_end:
            output += line

        # increment line counter
        line_num += 1
             
        


with open(outfile, 'w') as f:
    f.write(output)
