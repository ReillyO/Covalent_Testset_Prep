# accepts a PDB file containing 1+ atoms to be converted to spheres
# and no header information
# Writes a .sph file in DOCK format where each sphere corresponds
# to one of the atoms in the original PDB (radius 0.5)
# Owen O'Reilly, August 2025

# sample: python pdbtosph.py $infile $outfile

import sys, os
import pdb_utils as pdb


try:
  infile = sys.argv[1]
  outfile = sys.argv[2]
except:
  print("issue with inputs, please follow format 'python pdbtosph.py $infile $outfile'")
  exit(1)

output = ""
count = 1

numatoms = 0
with open(infile, 'r') as f:
  for line in f.readlines():
    if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
      numatoms += 1

with open(infile, 'r') as f:
  output = "DOCK spheres generated from receptor atoms\n"
  output += "cluster     1   number of spheres in cluster "+str(numatoms).rjust(4)+"\n"
  for line in f.readlines():
    if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
      entry = pdb.aline_to_list(line)
      coord = (entry[8], entry[9], entry[10])
    
      output += str(count).rjust(5)
      output += str(coord[0]).rjust(10)
      output += str(coord[1]).rjust(10)
      output += str(coord[2]).rjust(10)
      output += "   0.500 "
      output += str(count).rjust(4)
      output += " 0  0\n"
    
      count += 1
 
with open(outfile, 'w') as f:
  f.write(output)
