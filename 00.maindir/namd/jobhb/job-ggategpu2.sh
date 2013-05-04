#!/bin/bash
#PBS -N xxjobnamexx
#PBS -j oe
#PBS -l walltime=27:00
#PBS -l pmem=800mb
#PBS -l nodes=1:ppn=1
#PBS -V

# job properties
cd $PBS_O_WORKDIR

NUM=xxnumxx

# run job
./00-hb_pkl.py $NUM
