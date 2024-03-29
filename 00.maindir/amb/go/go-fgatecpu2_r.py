#!/usr/bin/env python
import sys,os,time
from glob import glob
import numpy as np
#import datetime
import time
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

def run_amb(i):
    st = time.time()
    os.system('${AMBERHOME}/bin/sander -O \
                        -i smd.in \
                        -o run.log \
                        -c ../00.rst \
                        -p ../../../../xxstrucequilxx/00.prmtop \
                        -r md.rst \
                        -x md.mdcrd \
                        -v md.mdvel')
    tt = time.time()-st
    os.system("mv dist_vs_t %d-tef.dat.%s" % (i,JOBID))
    os.system("python ../hb.py %d %s" % (i,JOBID))
    return tt

for i in range(1,howmany+1):
    t1 = run_amb(i)
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
    slp=randint(12,19)
    time.sleep(slp)
    next_stage(predir2,cfile)
elif count!=quota: sys.exit()

#os.system('rm *.BAK')
#os.system('rm *.rst')
os.system('rm mden')
os.system('rm mdinfo')
#os.system('rm *.mdcrd')
os.system('rm run.log')
