#!/bin/bash
#PBS -N dan10i07
#PBS -j oe
#PBS -l walltime=072:00:00
#PBS -l pmem=220mb
#PBS -l nodes=1:ppn=2
#PBS -V

# job properties
NAMD_DIR=/opt/NAMD29M/
export PATH=${NAMD_DIR}:${PATH}
export LD_LIBRARY_PATH=${NAMD_DIR}:${LD_LIBRARY_PATH}

cd $PBS_O_WORKDIR

# run job
./go.py