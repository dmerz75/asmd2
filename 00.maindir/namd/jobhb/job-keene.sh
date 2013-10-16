#!/bin/bash
#PBS -N xxjobnamexx
#PBS -q batch
#PBS -j oe
#PBS -l walltime=27:00
#PBS -l pmem=800mb
#PBS -l nodes=1:ppn=1
#PBS -A TG-CTS090079

# job_________________________
module load python/2.7.1

# job properties
cd $PBS_O_WORKDIR

NUM=xxnumxx

# run job
./00-hb_pkl.py $NUM
