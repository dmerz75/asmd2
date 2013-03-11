#!/usr/bin/env python
import sys,os,time,glob
import numpy as np
#import datetime

JOBID = os.environ['PBS_JOBID'].split('.')[0]
#now=datetime.datetime.now().strftime("%m%dt%H%M")

howmany=1
quota = 10

my_dir=os.path.abspath(os.path.dirname(__file__))
num   =my_dir.split('/')[-2]
predir1='/'.join(my_dir.split('/')[0:-1])
predir2='/'.join(my_dir.split('/')[0:-2])
cfile =os.path.join(predir2,'%s-continue.py' % num)

t_data = np.array([])

def run_namd(i):
    st = time.time()
    os.system('namd2 +p1 smd.namd > run.log')
    tt = time.time()-st
    os.system('mv tef.out %d-tef.dat.%s' % (i,JOBID))
    os.system('python hb.py %d %s' % (i,JOBID))
    return tt

for i in range(1,howmany+1):
    t1 = run_namd(i)
    t_data = np.append(t_data,t1)
    if i==1:
        os.system('mv daOut.dcd %d-daOut.dcd.%s' % (i,JOBID))
    if i%4==0 or i%howmany==0:
#        np.save('%d_time' % i,t_data)
        np.savetxt('time.dat',t_data,fmt='%.4f')

def next_stage(dir_loc,cfile):
    os.chdir(dir_loc)
    os.system('python %s' % cfile)

count=0
for path in glob.glob(os.path.join(predir1,'*/*-tef.dat*')):
    print path
    count+=1
print count

if count>=quota: next_stage(predir2,cfile)

'''
os.system('rm *.BAK')
os.system('rm da_smd.coor')
os.system('rm da_smd.dcd')
os.system('rm da_smd.vel')
os.system('rm da_smd.xsc')
os.system('rm run.log')
'''
