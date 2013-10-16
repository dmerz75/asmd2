#!/bin/bash
#PBS -N xxjobnamexx
#PBS -A TG-CTS090079
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -q batch


# job_________________________
module load python/2.7.1

NAMD_DIR=/sw/kfs/namd/2.9/centos6.2_intel12/NAMD_2.9_Linux-x86_64-ibverbs-smp-CUDA/
export PATH=${NAMD_DIR}:${PATH}
export LD_LIBRARY_PATH=${NAMD_DIR}:${LD_LIBRARY_PATH}

cd $PBS_O_WORKDIR

# run job
./go.py
