#!/bin/bash

# Sets up environment for standard covalent system generation
# RUN WITH SOURCE INSTEAD OF BASH

module load amber/22
module load chimera

# not to bin folder
export DOCKHOME="/path/to/dock/installation/"

# not to bin folder# 
export DMSHOME="/path/to/dms/installation"

# path to bin folder
export AMBERBIN="/path/to/amber/bin"

# path to bin folder
export CHIMERABIN="/path/to/chimera/bin"

export ROOTDIR=`pwd`

export VS_MPIDIR="/path/to/mpi/directory/"

if [ ! -s ${ROOTDIR}/zzz.testset_files/ ]; then
  mkdir ${ROOTDIR}/zzz.testset_files/
fi
