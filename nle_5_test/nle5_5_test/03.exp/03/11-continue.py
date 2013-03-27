#!/usr/bin/env python
import sys,os,pickle,shutil,fnmatch,itertools
import os.path
from glob import glob
from sys import argv
import numpy as np
from random import *
#from lockfile import FileLock  # provided in 04.scripts dir ->
import datetime,time                         # export PYTHONPATH

'''
lock = FileLock(__file__)
if lock.is_locked()==True: sys.exit()
lock.acquire()
'''

#import matplotlib
#import matplotlib.pyplot as plt
#from scipy.interpolate import LSQUnivariateSpline

my_dir = os.path.abspath(os.path.dirname(__file__))
num=sys.argv[1]

class mdict(dict):
    def __setitem__(self,key,value):
        self.setdefault(key,[]).append(value)
def print_dict(dct):
    for key,val in dct.items():
        print key,val
        print ''
    return key

config = pickle.load(open('config.pkl','rb'))
key = print_dict(config)

vel  = key
dist = config[vel][0][0]
ts   = config[vel][0][1]
path_seg   = config[vel][0][2]
path_svel  = config[vel][0][3]
path_vel   = config[vel][0][4]
path_steps = config[vel][0][5]
dct        = config[vel][0][6]       # 'freq'  50*ts/1000
dt         = dct['freq']*ts/1000
path_v_aps = path_vel/ts*1000
domain     = np.cumsum(((path_steps*ts)/1000)*path_v_aps)

spos=4.0
beta=-0.6
quota=5*20

'''
# matplotlib begin
fig=plt.figure()
plt.clf()
plt.xlabel('Extension (A)')
plt.ylabel('Work (kcal/mol)')
'''

def calc_work(data,st,w_c):
    phase = int(st)-1
    if phase == 0:
        data[::,3] = np.cumsum(data[::,3]*path_v_aps[phase]*dt)
        data[::,3] = data[::,3] + w_c[phase]
        wf = data[::,3][-1]
        d = np.linspace(spos,spos+domain[phase],data.shape[0])
        #plt.plot(d,data[::,3],'c-',linewidth=0.4)
        return data,wf
    else:
        data[::,3] = np.cumsum(data[::,3]*path_v_aps[phase]*dt)
        data[::,3] = data[::,3] + w_c[phase]
        wf = data[::,3][-1]
        d = np.linspace(spos+domain[phase-1],spos+domain[phase],data.shape[0])
        #plt.plot(d,data[::,3],'c-',linewidth=0.4)
        return data,wf
def calc_pmf(data,st,w_c):
    phase = int(st)-1
    if phase == 0:
        data[::,::,3] = np.cumsum(data[::,::,3]*path_v_aps[phase]*dt,axis=1)
        data[::,::,3] = data[::,::,3] + w_c[phase]
        deltaf = np.log(np.exp(data[::,::,3]*beta).mean(axis=0))*(1/beta)
        d = np.linspace(spos,spos+domain[phase],deltaf.shape[0])
        #plt.plot(d,deltaf,'b-',linewidth=1.0)
        JA = deltaf[-1]
        return JA
    else:
        data[::,::,3] = np.cumsum(data[::,::,3]*path_v_aps[phase]*dt,axis=1)
        data[::,::,3] = data[::,::,3] + w_c[phase]
        deltaf = np.log(np.exp(data[::,::,3]*beta).mean(axis=0))*(1/beta)
        #print data.shape[0]
        d = np.linspace(spos+domain[phase-1],spos+domain[phase],deltaf.shape[0])
        #plt.plot(d,deltaf,'b-',linewidth=1.0)
        JA = deltaf[-1]
        return JA

class asmd_calcs:
    def __init__(self,wrk_pkl,w_c,d_cp):
        self.wrk = wrk_pkl
        self.w_c = w_c
        self.d_cp = d_cp
    def create_pkl(self,st):
        acc = []
        self.wrk[st]={}
        self.d_cp[st]={}
        stdir = os.path.join(my_dir,st)
        folds=[f for f in os.listdir(stdir) if os.path.isdir(os.path.join( \
                                                        stdir,f))]
        foldp=[os.path.join(stdir,f) for f in folds]
        seeds=[(p.split('/')[-2],p.split('/')[-1].split('.')[2]) for f in foldp \
               for p in glob(os.path.join(f,'*tef.dat*'))]
        seeds=[p.split('/')[-1].split('.')[2] for f in foldp \
                 for p in glob(os.path.join(f,'*tef.dat*'))]
        for path in glob(os.path.join(my_dir,'%s/*/*tef.dat*' % st)):
            folder = path.split('/')[-2]
            seed = path.split('/')[-1].split('.')[2]
            self.wrk[st][seed]={}
        for path in glob(os.path.join(my_dir,'%s/*/*tef.dat*' % st)):
            folder = path.split('/')[-2]
            seed = path.split('/')[-1].split('.')[2]
            sample_i = np.loadtxt(path)
            # errcheck
            if len(sample_i)!= 2001:
                os.remove(path)
                del self.wrk[st][seed]
                seeds.remove(seed)
            else:
                acc.append(sample_i)
                data_1=np.array(sample_i)
                tew,wf=calc_work(data_1,st,self.w_c) #sample_i/data => tew
                self.wrk[st][seed]=folder,tew,wf
                os.remove(path)
        data=np.array(acc)
        JA=calc_pmf(data,st,self.w_c)               # get JA
        wf_sd=dict([(self.wrk[st][s][2],s) for s in seeds])
        sel_seed=wf_sd.get(JA, wf_sd[min(wf_sd.keys(), key=lambda k: abs(k-JA))])
        work_ss=self.wrk[st][sel_seed][2]
        self.d_cp[st][sel_seed]=self.wrk[st][sel_seed][0]
        print wf_sd
        print 'JA:',JA
        print 'selected seed',sel_seed
        print 'work for selected seed:',work_ss
        self.w_c[int(st)]=work_ss
        return seeds

def cp_file(f_dir,f,d_dir,d):
    shutil.copy(os.path.join(f_dir,f),os.path.join(d_dir,d))


def main_call(st,w_c,d_cp):
    wrk_pkl={}
    c1 = asmd_calcs(wrk_pkl,w_c,d_cp)
    seed_l = c1.create_pkl(st)
    nextnum=str(int(st)+1).zfill(2)
    seed_folder = d_cp[st]
    if st == num:
        for s,f in seed_folder.items():
            if not os.path.exists(os.path.join(my_dir,nextnum)):
                os.makedirs(os.path.join(my_dir,nextnum))
            cp_file(os.path.join(my_dir,st,f),'daOut.coor.%s' % s, \
                    os.path.join(my_dir,nextnum),'00.coor')
            cp_file(os.path.join(my_dir,st,f),'daOut.vel.%s' % s, \
                    os.path.join(my_dir,nextnum),'00.vel')
            if os.path.isfile(os.path.join(my_dir,nextnum,'00.vel'))==True:
                print 'COOR VEL in place'
                for path in glob(os.path.join(my_dir,'%s/*/*daOut*' % st)):
                    any_seed = path.split('/')[-1].split('.')[-1]
                    print any_seed
                    if any_seed == s:
                        pass
                    else:
                        os.remove(path)
    pickle.dump(w_c,open('%s-wc.pkl' % st,'w'))
    pickle.dump(wrk_pkl,open('%s-sfwf.pkl' % st,'w'))
    print w_c

def load_work(st,w_c):
    if st == '00':
	w_c[0]=0
    else:
        w_c = pickle.load(open('%s-wc.pkl' % st, 'rb'))
    return w_c

# main call
w_c    = {}
d_cp   = {}

w_c = load_work(str(int(num)-1).zfill(2),w_c)

main_call(num,w_c,d_cp)

'''
dirs = []
for i in range(1,int(num)+1):
    dirs.append(str(i).zfill(2))
[main_call(st,w_c,d_cp) for st in sorted(dirs)]

# matplotlib end
plt.draw()
plotname = '%s-eval' % num
plt.savefig('%s.png' % plotname)
plt.savefig('%s.eps' % plotname)
'''
