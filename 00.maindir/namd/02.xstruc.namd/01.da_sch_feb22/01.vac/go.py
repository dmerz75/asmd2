#!/usr/bin/env python
import sys, os

def run_cluster():
    JOBID = os.environ['PBS_JOBID'].split('.')[0]
    os.system('/nethome/dmerz3/packages/NAMD29Lx8664m/namd2 +p3 minv.namd > run.log')  
    os.system('mv run.log run.log.%s' % JOBID)

def run_local():
    os.system('/opt/NAMD_2.9b3_Linux-x86_64-multicore/namd2 +p3 minv.namd > run.log')

run_local()
