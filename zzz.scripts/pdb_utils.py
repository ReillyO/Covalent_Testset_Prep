# utility functions for dealing with PDB files in python
# currently includes:
# - line_to_list()
# - list_to_line()

import os, sys

# Accepts a string representing a single atom-type row in PDB file format
# Returns a Python array of strings where each entry corresponds to a PDB parameter - NO SPACES INCLUDED
# According to documentation at https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
def aline_to_list(line):
  linelist = []
  linelist.append(line[0:6].strip())     # 0:  ATOM, HETATM, TER, CONECT, etc
  linelist.append(line[6:11].strip())    # 1:  Atom number
  linelist.append(line[12:16].strip())   # 2:  Atom name
  linelist.append(line[16].strip())      # 3:  "Alternate location indicator"
  linelist.append(line[17:20].strip())   # 4:  Residue name
  linelist.append(line[21].strip())      # 5:  Chain identifier
  linelist.append(line[22:26].strip())   # 6:  Residue number
  linelist.append(line[26].strip())      # 7:  "Code for insertions of residues"
  linelist.append(line[30:38].strip())   # 8:  X-coord (A from origin)
  linelist.append(line[38:46].strip())   # 9:  Y-coord (A from origin)
  linelist.append(line[46:54].strip())   # 10: Z-coord (A from origin)
  linelist.append(line[54:60].strip())   # 11: Occupancy (0.0 - 1.0)
  linelist.append(line[60:66].strip())   # 12: Temperature factor
  linelist.append(line[72:76].strip())   # 13: Segment identifier
  linelist.append(line[76:78].strip())   # 14: Element symbol
  linelist.append(line[78:81].strip())   # 15: Charge
  return linelist


# Accepts a list of standard PDB atom entry
# Returns a formatted PDB line string with endline char
#  - Should correct any values where extra/less spaces are necessary so that all spacing is accurate
# According to documentation at https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
#
# sidenote: rjust(1) is necessary in case the field is empty
def alist_to_line(linelist):
  line = ""
  line += linelist[0].strip().ljust(6)          # 0:  ATOM, HETATM, TER, CONECT, etc
  line += linelist[1].strip().rjust(5)          # 1:  Atom number
  line += ' '
  line += linelist[2].strip().rjust(2).ljust(4) # 2:  Atom name, special justification for 1,2-char names
  line += linelist[3].strip().rjust(1)          # 3:  "Alternate location indicator"
  line += linelist[4].strip().rjust(3)          # 4:  Residue name
  line += ' '
  line += linelist[5].strip().rjust(1)          # 5:  Chain identifier
  line += linelist[6].strip().rjust(4)          # 6:  Residue number
  line += linelist[7].strip().rjust(1)          # 7:  "Code for insertions of residues"
  line += '   '
  line += linelist[8].strip().rjust(8)          # 8:  X-coord (A from origin)
  line += linelist[9].strip().rjust(8)          # 9:  Y-coord (A from origin)
  line += linelist[10].strip().rjust(8)         # 10: Z-coord (A from origin)
  line += linelist[11].strip().rjust(6)         # 11: Occupancy (0.0 - 1.0)
  line += linelist[12].strip().rjust(6)         # 12: Temperature factor
  line += '      '
  line += linelist[13].strip().ljust(4)         # 13: Segment identifier
  line += linelist[14].strip().rjust(2)         # 14: Element symbol
  line += linelist[15].strip().rjust(2)         # 15: Charge
  line += "\n"
  return line
