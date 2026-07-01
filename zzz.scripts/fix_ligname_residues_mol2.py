# accepts a ligand PDB file and:
#  - changes all residue names to LIG
#  - modifies chain of every residue to be A
# NOTE: this will incorporate any protein atoms present in the ligand
# without question, so be careful that only CB and SG are present

import sys, os

try:
  infile = sys.argv[1]
  outfile = sys.argv[2]
  ligname = sys.argv[3]
except:
  print("please provide a valid input file as argument 1 and an output filename as argument 2")
  exit(1)

outdata = ""
atom_count = 1

if "mol2" in infile:
  # Loop through each line and if line is atom:
  #  - change residue name to LIG
  #  - change residue number to 1
  #  - change chain to A
  # Also retain any lines with CONECT data
  # Otherwise discard line (cleans up header, etc) 
  with open(infile, 'r') as f:
    in_atom_section = False
    for line in f.readlines():
      if '@<TRIPOS>BOND' in line: in_atom_section = False

      if in_atom_section:
        linelist = line.split()
        linelist[6] = '1'
        linelist[7] = ligname
        line = ' '.join(linelist) + "\n"
      outdata += line
     
      if '@<TRIPOS>ATOM' in line: in_atom_section = True
      
  with open(outfile, 'w') as f:
    f.write(outdata)
else:
  print("can only process Mol2 files - please ensure the filetype is .mol2. exiting...")
  exit()
