#!/bin/bash
#PBS -N dan100i01
#PBS -q tg_workq
#PBS -l walltime=072:00:00
#PBS -l pmem=220mb
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -V

# job_________________________
module load namd/2.9b2-tcp
module load python/2.7.2

cd $PBS_O_WORKDIR

# run job
./go.py
