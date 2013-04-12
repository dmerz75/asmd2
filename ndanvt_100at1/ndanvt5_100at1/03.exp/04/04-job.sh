#!/bin/bash
#PBS -N danvtncjob04
#PBS -j oe
#PBS -l walltime=27:00
#PBS -l pmem=310mb
#PBS -l nodes=1:ppn=1
#PBS -V

# job properties
cd $PBS_O_WORKDIR

NUM=04

# run job
./00-continue.py $NUM
