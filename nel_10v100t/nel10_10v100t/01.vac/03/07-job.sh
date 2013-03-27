#!/bin/bash
#PBS -N elncjob07
#PBS -j oe
#PBS -l walltime=17:00
#PBS -l pmem=310mb
#PBS -l nodes=1:ppn=1
#PBS -V

# job properties
cd $PBS_O_WORKDIR

NUM=07

# run job
./$NUM-continue.py $NUM
