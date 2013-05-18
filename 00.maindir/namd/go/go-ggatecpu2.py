#!/usr/bin/env python
import sys,os,re,shutil,time
from glob import glob
import numpy as np
from random import randint

my_dir=os.path.abspath(os.path.dirname(__file__))

num   =my_dir.split('/')[-2]
prev_num=str(int(num)-1).zfill(2)
cwd_num =int(my_dir.split('/')[-1])

predir1='/'.join(my_dir.split('/')[0:-1])
predir2='/'.join(my_dir.split('/')[0:-2])
cfile =os.path.join(predir2,'%s-continue.py' % prev_num)

howmany= xxhowmanyxx
quota  = xxquotaxx*xxhowmanyxx
t_data = np.array([])

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

def run_namd(i,st_num,c_num):
    ''' 01.txt:
             dir = the number of columns
             tpd = len(column)
        # parameters
        i     : i in howmany
        st_num: stage
        c_num : current directory
    '''
    sd_arr = np.loadtxt('../%s.txt' % st_num)
    if len(sd_arr.shape)==1:
        seed = int(sd_arr[i-1])
    else:
        seed = int(sd_arr[i-1,c_num])
    cp_file(my_dir,'smd.namd',my_dir,'smd.namd.%s' % (seed))
    script = os.path.join(my_dir,'smd.namd.%s' % (seed))
    reg_ex(script,'xxxxx',str(seed),i)
    st = time.time()
    os.system('namd2 +pxxnodecountxx smd.namd.%s > run.log' % seed)
    tt = time.time()-st
    os.remove(script)
    os.rename('smdforces.out','%d-tef.dat.%s' % (i,seed))
    os.rename('daOut.coor','daOut.coor.%s' % (seed))
    os.rename('daOut.vel','daOut.vel.%s' % (seed))
    os.system('python ../hb.py %d %s' % (i,seed))
    os.rename('daOut.dcd','daOut.dcd.%s' % (seed))
    return tt

def next_stage(dir_loc,cfile):
    os.chdir(dir_loc)
    os.system('python %s %s' % (cfile,prev_num))

def check_vel(pnum):
    if pnum=='00':
	pass
    elif os.path.isfile(os.path.join(predir1,'00.vel'))==False:
        sys.exit()
    os.chdir(my_dir)

check_vel(prev_num)

for i in range(1,howmany+1):
    t1 = run_namd(i,num,cwd_num)
    t_data = np.append(t_data,t1)
    if i%howmany==0:
        np.savetxt('time.dat',t_data,fmt='%.4f')
