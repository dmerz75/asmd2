#!/bin/bash
#PBS -N xxjobnamexx
#PBS -q standby
#PBS -j oe
#PBS -l walltime=27:00
#PBS -l pmem=400mb
#PBS -l nodes=1:ppn=1
#PBS -V

# job_________________________
module load python/2.7.2

# job properties
cd $PBS_O_WORKDIR

NUM=xxnumxx

# run job
./$NUM-continue.py $NUM
