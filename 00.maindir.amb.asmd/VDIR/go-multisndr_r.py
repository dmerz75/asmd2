#!/usr/bin/env python
import sys,os,re,shutil,time
from glob import glob
import numpy as np
from random import randint

JOBID = os.environ['PBS_JOBID'].split('.')[0]

my_dir=os.path.abspath(os.path.dirname(__file__))
num   =my_dir.split('/')[-2]
prev_num=str(int(num)-1).zfill(2)
predir1='/'.join(my_dir.split('/')[0:-1])
predir2='/'.join(my_dir.split('/')[0:-2])
cfile =os.path.join(predir2,'%s-continue.py' % num)

howmany= xxhowmanyxx
quota  = xxquotaxx*xxhowmanyxx
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

def run_amb(i):
    st = time.time()
    seed = randint(10000,99999)
    while seed in slist:
        seed = randint(10000,99999)
    slist.append(seed)
    cp_file(my_dir,'smd.in',my_dir,'smd.in.%s' % (seed))
    script = os.path.join(my_dir,'smd.in.%s' % (seed))
    reg_ex(script,'xxxxx',str(seed),i)
    os.system('mpirun -np xxnodecountxx \
              -machinefile $PBS_NODEFILE ${AMBERHOME}/bin/sander -O \
              -i smd.in.%s \
              -o run.log \
              -c ../00.rst \
              -p ../../../../xxstrucequilxx/00.prmtop \
              -r md.rst.%s \
              -x md.mdcrd \
              -v md.mdvel' % (seed,seed))
    tt = time.time()-st
    os.system("mv dist_vs_t %d-tef.dat.%s" % (i,seed))
    os.system("python ../hb.py %d %s" % (i,seed))
    return tt

def next_stage(dir_loc,cfile):
    os.chdir(dir_loc)
    os.system('python %s %s' % (cfile,prev_num))

def check_vel(pnum):
    if pnum=='00':
        pass
    else:
        cnt = 0
        while os.path.isfile(os.path.join(predir1,'00.rst'))==False:
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
    t1 = run_amb(i)
    t_data = np.append(t_data,t1)
    if i%howmany==0:
        np.savetxt('time.dat',t_data,fmt='%.4f')
