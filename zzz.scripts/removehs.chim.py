import chimera
import os, sys
from chimera import runCommand

infile = sys.argv[1]
outfile = sys.argv[2]

runCommand('open '+infile)
runCommand('del #0:@H=')
runCommand('write format pdb #0 '+outfile)
