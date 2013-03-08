#!/usr/bin/env python
import os
import sys

JOBID = os.environ['PBS_JOBID']
JOBID1 = JOBID.split('.')[0]

os.system('/nethome/dmerz3/packages/NAMDCUDA/namd2 +idlepoll mine.namd > run.log')
#os.system('mv da_smd_tcl.out %d-tef.dat.%s' % (i,JOBID1))
#os.system('/nethome/dmerz3/packages/NAMD29Lx8664m/namd2 +p4 mine.namd > run.log')  
os.system('mv run.log run.log.%s' % (JOBID1))
