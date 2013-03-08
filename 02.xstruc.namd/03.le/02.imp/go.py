#!/usr/bin/env python
import sys, os

JOBID = os.environ['PBS_JOBID']
JOBID1 = JOBID.split('.')[0]

os.system('/nethome/dmerz3/packages/NAMD29Lx8664m/namd2 +p4 mini.namd > run.log')  
os.system('mv run.log run.log.%s' % JOBID1)
