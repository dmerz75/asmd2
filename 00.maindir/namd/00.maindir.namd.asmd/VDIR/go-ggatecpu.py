#!/usr/bin/env python
import sys,os,time
from glob import glob
import numpy as np
#import datetime
from random import randint

JOBID = os.environ['PBS_JOBID'].split('.')[0]
#now=datetime.datetime.now().strftime("%m%dt%H%M")

my_dir=os.path.abspath(os.path.dirname(__file__))
num   =my_dir.split('/')[-2]
predir1='/'.join(my_dir.split('/')[0:-1])
predir2='/'.join(my_dir.split('/')[0:-2])
cfile =os.path.join(predir2,'%s-continue.py' % num)

howmany= xxhowmanyxx
quota  = xxquotaxx
t_data = np.array([])

def run_namd(i):
    st = time.time()
    os.system('namd2 +pxxnodecountxx smd.namd > run.log')
    tt = time.time()-st
    os.system('mv smdforces.out %d-tef.dat.%s' % (i,JOBID))
    os.system('python ../hb.py %d %s' % (i,JOBID))
    return tt

for i in range(1,howmany+1):
    t1 = run_namd(i)
    t_data = np.append(t_data,t1)
    if i%howmany==0:
        np.savetxt('time.dat',t_data,fmt='%.4f')

def next_stage(dir_loc,cfile):
    os.chdir(dir_loc)
    os.system('python %s' % cfile)

count=0
for path in glob(os.path.join(predir1,'*/*-tef.dat*')):
    count+=1
if count==quota:
    next_stage(predir2,cfile)
    slp=randint(11,16)
    time.sleep(slp)
elif count!=quota: sys.exit()

#os.system('rm *.BAK')
#os.system('rm da_smd.coor')
#os.system('rm da_smd.dcd')
#os.system('rm da_smd.vel')
os.system('rm *.xsc')
os.system('rm run.log')
