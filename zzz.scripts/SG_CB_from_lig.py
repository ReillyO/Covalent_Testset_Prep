# will accept a ligand PDB file and output a 
# PDB file containing the atoms needed for covalent docking
# sphere generation

# sample: python cov_spheres_from_ligrec.py $ligfile <$outfile>
# outfile is optional, default is SG_CB.pdb 

import sys, os
import numpy as np
import pdb_utils as pdb

try:
  ligfile = sys.argv[1]
except:
  print("please provide a ligand and receptor file in style: python cov_spheres_from_ligrec.py $recfile $ligfile <$outfile>")
  exit()

try:
  outfile = sys.argv[2]
except:
  outfile = "SG_CB.pdb"

# get SG and CB atoms
SG_list = []
CB_list = []
with open(ligfile, 'r') as f:
  for line in f.readlines():
    if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
      linelist = pdb.aline_to_list(line)
      if linelist[2] == "CB":   CB_list = linelist
      elif linelist[2] == "SG": SG_list = linelist


out_txt = ""
out_txt += pdb.alist_to_line(SG_list)
out_txt += pdb.alist_to_line(CB_list)
out_txt += "TER" 

with open(outfile, 'w') as f:
  f.write(out_txt)
