#!/usr/bin/env python
import sys, os

JOBID = os.environ['PBS_JOBID'].split('.')[0]

os.system('namd2 +p4 mini.namd > run.log')
os.system('mv run.log run.log.%s' % JOBID)
