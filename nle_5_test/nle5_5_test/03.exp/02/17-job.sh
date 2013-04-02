#!/bin/bash
#PBS -N lencjob17
#PBS -j oe
#PBS -l walltime=17:00
#PBS -l pmem=310mb
#PBS -l nodes=1:ppn=1
#PBS -V

# job properties
cd $PBS_O_WORKDIR

NUM=17

# run job
./$NUM-continue.py $NUM
