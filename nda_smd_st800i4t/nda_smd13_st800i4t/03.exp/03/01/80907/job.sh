#!/bin/bash
#PBS -N dan10e01
#PBS -q tg_workq
#PBS -l walltime=15:00:00:00
#PBS -l pmem=220mb
#PBS -l nodes=1:ppn=3
#PBS -j oe
#PBS -V

# job_________________________
module load namd/2.9b2-tcp
module load python/2.7.2

cd $PBS_O_WORKDIR

# run job
./go.py
