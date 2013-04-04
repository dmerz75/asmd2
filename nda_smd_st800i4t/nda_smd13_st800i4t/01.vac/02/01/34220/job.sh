#!/bin/bash
#PBS -N dan100v01
#PBS -q tg_workq
#PBS -l walltime=03:59:00
#PBS -l pmem=220mb
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -V

# job_________________________
module load namd/2.9b2-tcp
module load python/2.7.2

cd $PBS_O_WORKDIR

# run job
./go.py
