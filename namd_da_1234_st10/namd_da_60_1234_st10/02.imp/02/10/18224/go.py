#!/usr/bin/env python
import sys,os,re,shutil,time
from glob import glob
import numpy as np
#import datetime
from random import randint

JOBID = os.environ['PBS_JOBID'].split('.')[0]

my_dir=os.path.abspath(os.path.dirname(__file__))
num   =my_dir.split('/')[-2]
prev_num=str(int(num)-1).zfill(2)
predir1='/'.join(my_dir.split('/')[0:-1])
predir2='/'.join(my_dir.split('/')[0:-2])
cfile =os.path.join(predir2,'%s-continue.py' % prev_num)

howmany= 20
quota  = 3*20
t_data = np.array([])
slist  = []

def cp_file(f_dir,f,d_dir,d):
    shutil.copy(os.path.join(f_dir,f),os.path.join(d_dir,d))

def reg_ex(script,subout,subin,n):
    o=open(script,'r+')
    text=o.read()
    text=re.sub(subout,subin,text)
    o.close()
    o=open(script,'w+')
    o.write(text)
    o.close()

def run_namd(i):
    seed = randint(10000,99999)
    while seed in slist:
        seed = randint(10000,99999)
    slist.append(seed)
    cp_file(my_dir,'smd.namd',my_dir,'smd.namd.%s' % (seed))
    script = os.path.join(my_dir,'smd.namd.%s' % (seed))
    reg_ex(script,'xxxxx',str(seed),i)
    st = time.time()
    os.system('namd2 +p2 smd.namd.%s > run.log' % seed)
    tt = time.time()-st
    os.remove(script)
    os.rename('smdforces.out','%d-tef.dat.%s' % (i,seed))
    os.rename('daOut.coor','daOut.coor.%s' % (seed))
    os.rename('daOut.vel','daOut.vel.%s' % (seed))
    os.system('python ../hb.py %d %s' % (i,seed))
    return tt

def next_stage(dir_loc,cfile):
    os.chdir(dir_loc)
    os.system('python %s %s' % (cfile,prev_num))

def check_vel(pnum):
    if pnum=='00':
        pass
    else:
        cnt = 0
        while os.path.isfile(os.path.join(predir1,'00.vel'))==False:
            next_stage(predir2,cfile)
            cnt += 1
            print "waiting 20 seconds ..."
            print cnt
            time.sleep(20)
            if cnt >= 20:
                sys.exit()
    os.chdir(my_dir)

check_vel(prev_num)

for i in range(1,howmany+1):
    t1 = run_namd(i)
    t_data = np.append(t_data,t1)
    if i%howmany==0:
        np.savetxt('time.dat',t_data,fmt='%.4f')
