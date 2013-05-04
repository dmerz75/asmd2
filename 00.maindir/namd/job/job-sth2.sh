#!/bin/bash
#PBS -N xxjobnamexx
#PBS -q xxqueuexx
#PBS -l xxwalltimexx
#PBS -l pmem=480mb
#PBS -l xxnodesxx
#PBS -j oe
#PBS -V

# job_________________________
module load python/2.7.2

cd $PBS_O_WORKDIR

# run job
./go.py
